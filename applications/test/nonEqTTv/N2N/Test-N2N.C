// heat bath 0D di N2 + N (fig. 5, 7 e 8 del paper): 5e22 particelle/m3 per
// specie, rilassamento V-T e, se richiesto, la dissociazione
// N2 + N2 -> 2N + N2 con le costanti di Park (meccanismo N2_Park, tabella 2)
//
// uso: Test-N2N <T_tr> <T_ve> <t_fine> <meccanismo> <esponente_Park> <file_csv> [tau] [C-V]
//   fig 5: Test-N2N 30000  1000 1e-5 none    0.7 output/fig5.csv
//   fig 7: Test-N2N 30000  1000 1e-3 N2_Park 0.7 output/fig7.csv
//   fig 8: Test-N2N 30000 30000 1e-4 N2_Park 0.7 output/fig8.csv
// l'esponente e' quello della temperatura di Park T^a Tv^(1-a) (eq. 29):
// 0.7 come il paper, 0.5 e' quello fisso di Mutation++
// argomenti facoltativi, per confronto con le scelte del paper (default):
//   tau: paper (eq. 9-17, default) oppure mutation (MillikanWhite di Mutation++)
//   C-V: preferential (eq. 32, alpha = 0.3, default) oppure nonPreferential (eq. 31)

#include "mutation++.h"
#include "mutationSources.H"
#include "heatBath.H"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char *argv[])
{
    if (argc < 7 || argc > 9)
    {
        std::cerr << "uso: Test-N2N <T_tr> <T_ve> <t_fine> <meccanismo> <esponente_Park> <file_csv>"
                  << " [paper|mutation] [preferential|nonPreferential]" << std::endl;
        return 1;
    }
    const double T_tr0 = std::atof(argv[1]);
    const double T_ve0 = std::atof(argv[2]);
    const double t_end = std::atof(argv[3]);
    const std::string mechanism = argv[4];
    const double parkExponent = std::atof(argv[5]);
    const char* csvName = argv[6];
    const bool chemistry = (mechanism != "none");
    // modelli del paper di default: tau delle eq. 9-17 e Q_C-V preferenziale
    const bool paperTau = (argc < 8 || std::string(argv[7]) != "mutation");
    const bool preferential = (argc < 9 || std::string(argv[8]) != "nonPreferential");
    const double alpha = 0.3;

    Mutation::MixtureOptions opts("air_5");
    opts.setStateModel("ChemNonEqTTv");
    opts.setThermodynamicDatabase("RRHO");
    opts.setMechanism(mechanism);
    Mutation::Mixture mix(opts);

    const int ns = mix.nSpecies();
    const int iN2 = mix.speciesIndex("N2");
    const int iN = mix.speciesIndex("N");

    // ---- stato iniziale: stessa densita' numerica n per N2 e N
    const double n = 5.0e22;
    std::vector<double> rho_s(ns, 0.0);
    rho_s[iN2] = n * mix.speciesMw(iN2) / Mutation::NA;
    rho_s[iN] = n * mix.speciesMw(iN) / Mutation::NA;
    const double n0 = 2.0 * n;

    const double temps0[2] = {T_tr0, T_ve0};
    mix.setState(rho_s.data(), temps0, 1);

    // energie conservate per unita' di volume (eq. 23): E resta costante,
    // anche con la chimica perche' comprende l'energia di formazione
    const double E = totalEnergy(mix, rho_s);
    double Eve = veEnergy(mix, rho_s);
    // ---- fine stato iniziale

    std::ofstream csv(csvName);
    if (!chemistry)
    {
        // senza chimica la composizione e' fissa: T finale dalla conservazione
        // dell'energia (il paper la legge dalla figura)
        const double T_eq = equilibriumTemperature(mix, rho_s, E);
        mix.setState(rho_s.data(), temps0, 1);
        csv << "# T_eq = " << T_eq << " K (conservazione dell'energia)\n";
    }
    csv << "t,Ttr,Tv,N2,N\n";

    // densita' numeriche normalizzate n_s/n0, come nelle fig. 7b e 8b
    auto nN2 = [&]() { return rho_s[iN2] / mix.speciesMw(iN2) * Mutation::NA / n0; };
    auto nN = [&]() { return rho_s[iN] / mix.speciesMw(iN) * Mutation::NA / n0; };
    csv << "0," << mix.T() << "," << mix.Tv() << "," << nN2() << "," << nN() << "\n";

    // tempi di rilassamento V-T delle molecole (qui solo N2)
    const std::vector<Vibrator> vibrators = makeVibrators(mix);

    // ---- ciclo nel tempo: passo di 1 ns come il paper
    const double dt = 1.0e-9;
    const int nSteps = int(t_end / dt + 0.5);
    // righe del csv: ogni passo all'inizio, poi a passo logaritmico
    OutputSchedule output;
    double t = 0.0;
    std::vector<double> wdot(ns, 0.0);

    for (int step = 1; step <= nSteps; step++)
    {
        // scambio V-T (eq. 8)
        double Q = sourceVT(mix, rho_s, vibrators, paperTau);

        if (chemistry)
        {
            // produzione delle specie alla temperatura di Park (eq. 27-29)
            // e energia vibro-elettronica che se ne va con esse (eq. 30)
            productionRates(mix, rho_s, parkExponent, wdot);
            Q += sourceCV(mix, wdot, vibrators, preferential, alpha);
            for (int s = 0; s < ns; s++)
            {
                rho_s[s] += wdot[s] * dt;
            }
        }

        // eq. 22: cambiano E_ve e le densita', E si conserva
        Eve += Q * dt;
        t += dt;

        // nuove temperature dalle energie
        const double energies[2] = {E, Eve};
        mix.setState(rho_s.data(), energies, 0);

        if (output.write(step, nSteps))
        {
            csv << t << "," << mix.T() << "," << mix.Tv() << "," << nN2() << "," << nN() << "\n";
        }
    }
    // ---- fine ciclo nel tempo

    std::cout << csvName << ": T_tr = " << mix.T() << " K, T_ve = " << mix.Tv()
              << " K, n_N2/n0 = " << nN2() << ", n_N/n0 = " << nN()
              << " a t = " << t << " s" << std::endl;

    return 0;
}

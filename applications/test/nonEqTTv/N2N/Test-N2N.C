// heat bath 0D di N2 + N con Mutation++ (fig. 5, 7 e 8 del paper):
// 5e22 particelle/m3 per specie, rilassamento V-T e, se richiesto, la
// dissociazione N2 + N2 -> 2N + N2 con le costanti di Park (meccanismo N2_Park)
//
// uso: Test-N2N <T_tr> <T_ve> <t_fine> <meccanismo> <file_csv>
//   fig 5: Test-N2N 30000  1000 1e-5 none    output/fig5.csv
//   fig 7: Test-N2N 30000  1000 1e-3 N2_Park output/fig7.csv
//   fig 8: Test-N2N 30000 30000 1e-4 N2_Park output/fig8.csv

#include "mutation++.h"
#include "heatBath.H"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char *argv[])
{
    if (argc != 6)
    {
        std::cerr << "uso: Test-N2N <T_tr> <T_ve> <t_fine> <meccanismo> <file_csv>" << std::endl;
        return 1;
    }
    const double T_tr0 = std::atof(argv[1]);
    const double T_ve0 = std::atof(argv[2]);
    const double t_end = std::atof(argv[3]);
    const std::string mechanism = argv[4];
    const char* csvName = argv[5];
    const bool chemistry = (mechanism != "none");

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

    // ---- ciclo nel tempo: passo di 1 ns come il paper
    const double dt = 1.0e-9;
    double t = 0.0;
    int step = 0;
    std::vector<double> Q(mix.nEnergyEqns());
    std::vector<double> wdot(ns, 0.0);

    while (t < t_end)
    {
        // Q[0]: scambio V-T (eq. 8) piu', con la chimica, l'energia
        // vibro-elettronica tolta dalle reazioni (eq. 30, modello non
        // preferenziale, l'unico di Mutation++)
        mix.energyTransferSource(Q.data());
        Eve += Q[0] * dt;

        if (chemistry)
        {
            // produzione delle specie (eq. 27), kg/m3/s
            mix.netProductionRates(wdot.data());
            for (int s = 0; s < ns; s++)
            {
                rho_s[s] += wdot[s] * dt;
            }
        }
        t += dt;
        step++;

        // nuove temperature dalle energie
        const double energies[2] = {E, Eve};
        mix.setState(rho_s.data(), energies, 0);

        if (step % 10 == 0)
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

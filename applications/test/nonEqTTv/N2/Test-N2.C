// heat bath 0D di N2 puro: rilassamento V-T (fig. 3a, 3b, 4 del paper)
//
// uso: Test-N2 <T_tr> <T_ve> <t_fine> <file_csv>
//   fig 3a: Test-N2 10000  1000 3e-5 output/fig3a.csv
//   fig 3b: Test-N2  3000 10000 1e-4 output/fig3b.csv
//   fig 4 : Test-N2 30000  1000 1e-5 output/fig4-noEl.csv   (o fig4-el.csv)
//
// con o senza energia elettronica lo decide la cartella dati di Mutation++
// (variabile MPP_DATA_DIRECTORY: mutation-data oppure mutation-data-noElectronic)

#include "mutation++.h"
#include "mutationSources.H"
#include "heatBath.H"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <vector>

int main(int argc, char *argv[])
{
    if (argc != 5)
    {
        std::cerr << "uso: Test-N2 <T_tr> <T_ve> <t_fine> <file_csv>" << std::endl;
        return 1;
    }
    const double T_tr0 = std::atof(argv[1]);
    const double T_ve0 = std::atof(argv[2]);
    const double t_end = std::atof(argv[3]);
    const char* csvName = argv[4];

    // miscela di aria a 5 specie (contiene N2), modello a due temperature,
    // energie RRHO, nessuna reazione chimica
    Mutation::MixtureOptions opts("air_5");
    opts.setStateModel("ChemNonEqTTv");
    opts.setThermodynamicDatabase("RRHO");
    opts.setMechanism("none");
    Mutation::Mixture mix(opts);

    const int ns = mix.nSpecies();
    const int iN2 = mix.speciesIndex("N2");

    // ---- stato iniziale: solo N2 a 1 atm, T_tr e T_ve date
    std::vector<double> Y(ns, 0.0);
    Y[iN2] = 1.0;
    const double P_T_Tv[3] = {Mutation::ONEATM, T_tr0, T_ve0};
    mix.setState(Y.data(), P_T_Tv, 2);

    // densita' parziali: restano costanti (scatola chiusa, niente chimica)
    std::vector<double> rho_s(ns);
    mix.densities(rho_s.data());

    // energie conservate per unita' di volume (eq. 23)
    const double E = totalEnergy(mix, rho_s);
    double Eve = veEnergy(mix, rho_s);
    // ---- fine stato iniziale

    // temperatura finale attesa: il paper non la da' sempre, la ricavo
    // dalla conservazione dell'energia
    const double T_eq = equilibriumTemperature(mix, rho_s, E);
    const double temps0[2] = {T_tr0, T_ve0};
    mix.setState(rho_s.data(), temps0, 1);

    // tempi di rilassamento V-T delle molecole (qui solo N2): formula del
    // paper (eq. 9-17); per N2 puro coincide con quella di Mutation++
    const std::vector<Vibrator> vibrators = makeVibrators(mix);
    const bool paperTau = true;

    std::ofstream csv(csvName);
    csv << "# T_eq = " << T_eq << " K (conservazione dell'energia)\n";
    csv << "t,Ttr,Tv\n";
    csv << "0," << mix.T() << "," << mix.Tv() << "\n";

    // ---- ciclo nel tempo: passo di 1 ns come il paper
    const double dt = 1.0e-9;
    const int nSteps = int(t_end / dt + 0.5);
    // righe del csv: ogni passo all'inizio, poi a passo logaritmico
    OutputSchedule output;
    double t = 0.0;

    for (int step = 1; step <= nSteps; step++)
    {
        // Q_VT: energia che passa dal serbatoio traslazionale a quello
        // vibro-elettronico (Landau-Teller, eq. 8, tau da Mutation++)
        const double Q = sourceVT(mix, rho_s, vibrators, paperTau);

        // eq. 22: cambia solo E_ve, l'energia totale E si conserva
        Eve += Q * dt;
        t += dt;

        // nuove temperature dalle energie (Mutation++ inverte E ed E_ve)
        const double energies[2] = {E, Eve};
        mix.setState(rho_s.data(), energies, 0);

        if (output.write(step, nSteps))
        {
            csv << t << "," << mix.T() << "," << mix.Tv() << "\n";
        }
    }
    // ---- fine ciclo nel tempo

    std::cout << csvName << ": T_tr = " << mix.T() << " K, T_ve = " << mix.Tv()
              << " K a t = " << t << " s; T_eq attesa = " << T_eq << " K" << std::endl;

    return 0;
}

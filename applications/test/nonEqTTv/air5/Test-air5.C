// heat bath 0D di aria a 5 specie reagente (fig. 9 del paper): 0.79 N2 e
// 0.21 O2 in moli, 0.063 atm, 10000 K, tutte le temperature in equilibrio
// fra loro (modello a una temperatura, come fa il paper per questa figura)
//
// uso: Test-air5 <meccanismo> <t_fine> <file_csv> [secondo_meccanismo]
//   fig 9: Test-air5 air5_Park 1e-3 output/fig9-park.csv
//          Test-air5 air5_QK   1e-3 output/fig9-QK.csv  air5_QK_back
// il meccanismo e' un file di mutation-data/mechanisms (air5_Park e' quello
// di Mutation++ con le costanti di Park, air5_QK il set del paper); il
// secondo meccanismo, se dato, aggiunge le sue reazioni (serve per gli
// scambi inversi del set QK, che Mutation++ non accetta nello stesso file);
// la temperatura di controllo e' sempre T, l'esponente di Park non entra

#include "mutation++.h"
#include "heatBath.H"

#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char *argv[])
{
    if (argc != 4 && argc != 5)
    {
        std::cerr << "uso: Test-air5 <meccanismo> <t_fine> <file_csv> [secondo_meccanismo]"
                  << std::endl;
        return 1;
    }
    const std::string mechanism = argv[1];
    const double t_end = std::atof(argv[2]);
    const char* csvName = argv[3];
    const std::string mechanism2 = (argc == 5) ? argv[4] : "none";

    // una sola temperatura: stato model ChemNonEq1T
    Mutation::MixtureOptions opts("air_5");
    opts.setStateModel("ChemNonEq1T");
    opts.setThermodynamicDatabase("RRHO");
    opts.setMechanism(mechanism);
    Mutation::Mixture mix(opts);

    // seconda miscela, uguale alla prima ma con le reazioni del secondo meccanismo
    opts.setMechanism(mechanism2);
    Mutation::Mixture mix2(opts);

    const int ns = mix.nSpecies();
    const std::vector<std::string> names = {"N2", "O2", "NO", "N", "O"};

    // ---- stato iniziale: frazioni molari 0.79 N2 e 0.21 O2 -> frazioni in massa
    std::vector<double> X(ns, 0.0);
    X[mix.speciesIndex("N2")] = 0.79;
    X[mix.speciesIndex("O2")] = 0.21;
    std::vector<double> Y(ns, 0.0);
    double Mmix = 0.0;
    for (int s = 0; s < ns; s++)
    {
        Mmix += X[s] * mix.speciesMw(s);
    }
    for (int s = 0; s < ns; s++)
    {
        Y[s] = X[s] * mix.speciesMw(s) / Mmix;
    }

    const double P_T[2] = {0.063 * Mutation::ONEATM, 10000.0};
    mix.setState(Y.data(), P_T, 2);

    std::vector<double> rho_s(ns);
    mix.densities(rho_s.data());
    // densita' numerica totale iniziale, per normalizzare n_s/n0
    const double n0 = mix.numberDensity();

    // energia totale per unita' di volume (con la formazione): si conserva
    const double E = totalEnergy(mix, rho_s);
    // ---- fine stato iniziale

    std::ofstream csv(csvName);
    csv << "t,T";
    for (const std::string& name : names)
    {
        csv << "," << name;
    }
    csv << "\n";

    // riga del csv: t, T e le densita' numeriche normalizzate
    auto writeLine = [&](double t)
    {
        csv << t << "," << mix.T();
        for (const std::string& name : names)
        {
            const int s = mix.speciesIndex(name);
            csv << "," << rho_s[s] / mix.speciesMw(s) * Mutation::NA / n0;
        }
        csv << "\n";
    };
    writeLine(0.0);

    // ---- ciclo nel tempo: passo di 1 ns come il paper
    const double dt = 1.0e-9;
    const int nSteps = int(t_end / dt + 0.5);
    const int writeEvery = std::max(10, nSteps / 10000);
    double t = 0.0;
    std::vector<double> wdot(ns, 0.0);
    std::vector<double> wdot2(ns, 0.0);

    for (int step = 1; step <= nSteps; step++)
    {
        // produzione delle specie (eq. 27), tutte le reazioni a T; il secondo
        // meccanismo va portato nello stesso stato e le sue velocita' sommate
        mix.netProductionRates(wdot.data());
        const double T_now[1] = {mix.T()};
        mix2.setState(rho_s.data(), T_now, 1);
        mix2.netProductionRates(wdot2.data());
        for (int s = 0; s < ns; s++)
        {
            rho_s[s] += (wdot[s] + wdot2[s]) * dt;
        }
        t += dt;

        // nuova temperatura dall'energia totale (Mutation++ la inverte)
        const double energies[1] = {E};
        mix.setState(rho_s.data(), energies, 0);

        if (step % writeEvery == 0)
        {
            writeLine(t);
        }
    }
    // ---- fine ciclo nel tempo

    std::cout << csvName << ": T = " << mix.T() << " K a t = " << t << " s" << std::endl;

    return 0;
}

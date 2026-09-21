// 0D heat bath of reacting 5-species air (fig. 9 of the paper): 0.79 N2 and
// 0.21 O2 by mole, 0.063 atm, 10000 K, all temperatures in equilibrium
// with each other (one-temperature model, as the paper does for this figure)
//
// usage: Test-air5 <mechanism> <t_end> <csv_file> [second_mechanism]
//   fig 9: Test-air5 air5_Park 1e-3 output/fig9-park.csv
//          Test-air5 air5_QK   1e-3 output/fig9-QK.csv  air5_QK_back
// the mechanism is a file of mutation-data/mechanisms (air5_Park is the
// Mutation++ one with the Park constants, air5_QK the set of the paper); the
// second mechanism, if given, adds its reactions (needed for the reverse
// exchanges of the QK set, which Mutation++ does not accept in the same file);
// the controlling temperature is always T, the Park exponent does not enter

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
        std::cerr << "usage: Test-air5 <mechanism> <t_end> <csv_file> [second_mechanism]"
                  << std::endl;
        return 1;
    }
    const std::string mechanism = argv[1];
    const double t_end = std::atof(argv[2]);
    const char* csvName = argv[3];
    const std::string mechanism2 = (argc == 5) ? argv[4] : "none";

    // a single temperature: state model ChemNonEq1T
    Mutation::MixtureOptions opts("air_5");
    opts.setStateModel("ChemNonEq1T");
    opts.setThermodynamicDatabase("RRHO");
    opts.setMechanism(mechanism);
    Mutation::Mixture mix(opts);

    // second mixture, identical to the first but with the reactions of the second mechanism
    opts.setMechanism(mechanism2);
    Mutation::Mixture mix2(opts);

    const int ns = mix.nSpecies();
    const std::vector<std::string> names = {"N2", "O2", "NO", "N", "O"};

    // ---- initial state: mole fractions 0.79 N2 and 0.21 O2 -> mass fractions
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
    // initial total number density, to normalise n_s/n0
    const double n0 = mix.numberDensity();

    // total energy per unit volume (including formation): it is conserved
    const double E = totalEnergy(mix, rho_s);
    // ---- end of initial state

    std::ofstream csv(csvName);
    csv << "t,T";
    for (const std::string& name : names)
    {
        csv << "," << name;
    }
    csv << "\n";

    // csv row: t, T and the normalised number densities
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

    // ---- time loop: time step of 1 ns as in the paper
    const double dt = 1.0e-9;
    const int nSteps = int(t_end / dt + 0.5);
    // csv rows: every step at the beginning, then at logarithmic intervals
    OutputSchedule output;
    double t = 0.0;
    std::vector<double> wdot(ns, 0.0);
    std::vector<double> wdot2(ns, 0.0);

    for (int step = 1; step <= nSteps; step++)
    {
        // species production (eq. 27), all reactions at T; the second
        // mechanism must be brought to the same state and its rates added
        mix.netProductionRates(wdot.data());
        const double T_now[1] = {mix.T()};
        mix2.setState(rho_s.data(), T_now, 1);
        mix2.netProductionRates(wdot2.data());
        for (int s = 0; s < ns; s++)
        {
            rho_s[s] += (wdot[s] + wdot2[s]) * dt;
        }
        t += dt;

        // new temperature from the total energy (Mutation++ inverts it)
        const double energies[1] = {E};
        mix.setState(rho_s.data(), energies, 0);

        if (output.write(step, nSteps))
        {
            writeLine(t);
        }
    }
    // ---- end of time loop

    std::cout << csvName << ": T = " << mix.T() << " K at t = " << t << " s" << std::endl;

    return 0;
}

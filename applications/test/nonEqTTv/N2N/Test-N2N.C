// 0D heat bath of N2 + N (fig. 5, 7 and 8 of the paper): 5e22 particles/m3 per
// species, V-T relaxation and, if requested, the dissociation
// N2 + N2 -> 2N + N2 with the Park constants (mechanism N2_Park, table 2)
//
// usage: Test-N2N <T_tr> <T_ve> <t_end> <mechanism> <Park_exponent> <csv_file> [tau] [C-V]
//   fig 5: Test-N2N 30000  1000 1e-5 none    0.7 output/fig5.csv
//   fig 7: Test-N2N 30000  1000 1e-3 N2_Park 0.7 output/fig7.csv
//   fig 8: Test-N2N 30000 30000 1e-4 N2_Park 0.7 output/fig8.csv
// the exponent is that of Park's temperature T^a Tv^(1-a) (eq. 29):
// 0.7 as in the paper, 0.5 is the fixed one of Mutation++
// optional arguments, for comparison with the choices of the paper (default):
//   tau: paper (eq. 9-17, default) or mutation (MillikanWhite of Mutation++)
//   C-V: preferential (eq. 32, alpha = 0.3, default) or nonPreferential (eq. 31)

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
        std::cerr << "usage: Test-N2N <T_tr> <T_ve> <t_end> <mechanism> <Park_exponent> <csv_file>"
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
    // models of the paper by default: tau of eq. 9-17 and preferential Q_C-V
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

    // ---- initial state: same number density n for N2 and N
    const double n = 5.0e22;
    std::vector<double> rho_s(ns, 0.0);
    rho_s[iN2] = n * mix.speciesMw(iN2) / Mutation::NA;
    rho_s[iN] = n * mix.speciesMw(iN) / Mutation::NA;
    const double n0 = 2.0 * n;

    const double temps0[2] = {T_tr0, T_ve0};
    mix.setState(rho_s.data(), temps0, 1);

    // conserved energies per unit volume (eq. 23): E stays constant,
    // even with chemistry, because it includes the energy of formation
    const double E = totalEnergy(mix, rho_s);
    double Eve = veEnergy(mix, rho_s);
    // ---- end of initial state

    std::ofstream csv(csvName);
    if (!chemistry)
    {
        // without chemistry the composition is fixed: final T from energy
        // conservation (the paper reads it from the figure)
        const double T_eq = equilibriumTemperature(mix, rho_s, E);
        mix.setState(rho_s.data(), temps0, 1);
        csv << "# T_eq = " << T_eq << " K (energy conservation)\n";
    }
    csv << "t,Ttr,Tv,N2,N\n";

    // normalised number densities n_s/n0, as in fig. 7b and 8b
    auto nN2 = [&]() { return rho_s[iN2] / mix.speciesMw(iN2) * Mutation::NA / n0; };
    auto nN = [&]() { return rho_s[iN] / mix.speciesMw(iN) * Mutation::NA / n0; };
    csv << "0," << mix.T() << "," << mix.Tv() << "," << nN2() << "," << nN() << "\n";

    // V-T relaxation times of the molecules (here N2 only)
    const std::vector<Vibrator> vibrators = makeVibrators(mix);

    // ---- time loop: time step of 1 ns as in the paper
    const double dt = 1.0e-9;
    const int nSteps = int(t_end / dt + 0.5);
    // csv rows: every step at the beginning, then at logarithmic intervals
    OutputSchedule output;
    double t = 0.0;
    std::vector<double> wdot(ns, 0.0);

    for (int step = 1; step <= nSteps; step++)
    {
        // V-T exchange (eq. 8)
        double Q = sourceVT(mix, rho_s, vibrators, paperTau);

        if (chemistry)
        {
            // species production at Park's temperature (eq. 27-29)
            // and the vibro-electronic energy that leaves with them (eq. 30)
            productionRates(mix, rho_s, parkExponent, wdot);
            Q += sourceCV(mix, wdot, vibrators, preferential, alpha);
            for (int s = 0; s < ns; s++)
            {
                rho_s[s] += wdot[s] * dt;
            }
        }

        // eq. 22: E_ve and the densities change, E is conserved
        Eve += Q * dt;
        t += dt;

        // new temperatures from the energies
        const double energies[2] = {E, Eve};
        mix.setState(rho_s.data(), energies, 0);

        if (output.write(step, nSteps))
        {
            csv << t << "," << mix.T() << "," << mix.Tv() << "," << nN2() << "," << nN() << "\n";
        }
    }
    // ---- end of time loop

    std::cout << csvName << ": T_tr = " << mix.T() << " K, T_ve = " << mix.Tv()
              << " K, n_N2/n0 = " << nN2() << ", n_N/n0 = " << nN()
              << " at t = " << t << " s" << std::endl;

    return 0;
}

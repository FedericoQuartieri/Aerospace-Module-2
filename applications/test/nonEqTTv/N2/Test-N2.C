// 0D heat bath of pure N2: V-T relaxation (fig. 3a, 3b, 4 of the paper)
//
// usage: Test-N2 <T_tr> <T_ve> <t_end> <csv_file> [vibrationalOnly]
//   fig 3a: Test-N2 10000  1000 3e-5 output/fig3a.csv
//   fig 3b: Test-N2  3000 10000 1e-4 output/fig3b.csv
//   fig 4 : Test-N2 30000  1000 1e-5 output/fig4-noEl.csv   (or fig4-el.csv)
//
// whether electronic energy is included is decided by the Mutation++ data folder
// (variable MPP_DATA_DIRECTORY: mutation-data or mutation-data-noElectronic)
// optional argument vibrationalOnly: the V-T exchange is driven by the vibrational
// energy only, as in the OmegaVT of Mutation++, instead of e_ve (paper eq. 8, default);
// it is the alternative of report table 1, which changes fig. 4 with E_el only

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
    if (argc < 5 || argc > 6)
    {
        std::cerr << "usage: Test-N2 <T_tr> <T_ve> <t_end> <csv_file> [vibrationalOnly]"
                  << std::endl;
        return 1;
    }
    const double T_tr0 = std::atof(argv[1]);
    const double T_ve0 = std::atof(argv[2]);
    const double t_end = std::atof(argv[3]);
    const char* csvName = argv[4];
    // default: e_ve drives the V-T exchange (paper eq. 8)
    const bool vibrationalOnly = (argc == 6 && std::string(argv[5]) == "vibrationalOnly");

    // 5-species air mixture (contains N2), two-temperature model,
    // RRHO energies, no chemical reactions
    Mutation::MixtureOptions opts("air_5");
    opts.setStateModel("ChemNonEqTTv");
    opts.setThermodynamicDatabase("RRHO");
    opts.setMechanism("none");
    Mutation::Mixture mix(opts);

    const int ns = mix.nSpecies();
    const int iN2 = mix.speciesIndex("N2");

    // ---- initial state: N2 only at 1 atm, given T_tr and T_ve
    std::vector<double> Y(ns, 0.0);
    Y[iN2] = 1.0;
    const double P_T_Tv[3] = {Mutation::ONEATM, T_tr0, T_ve0};
    mix.setState(Y.data(), P_T_Tv, 2);

    // partial densities: they stay constant (closed box, no chemistry)
    std::vector<double> rho_s(ns);
    mix.densities(rho_s.data());

    // conserved energies per unit volume (eq. 23)
    const double E = totalEnergy(mix, rho_s);
    double Eve = veEnergy(mix, rho_s);
    // ---- end of initial state

    // expected final temperature: the paper does not always give it, so it is
    // obtained from energy conservation
    const double T_eq = equilibriumTemperature(mix, rho_s, E);
    const double temps0[2] = {T_tr0, T_ve0};
    mix.setState(rho_s.data(), temps0, 1);

    // V-T relaxation times of the molecules (here N2 only): formula of the
    // paper (eq. 9-17); for pure N2 it coincides with the Mutation++ one
    const std::vector<Vibrator> vibrators = makeVibrators(mix);
    const bool paperTau = true;

    std::ofstream csv(csvName);
    csv << "# T_eq = " << T_eq << " K (energy conservation)\n";
    csv << "t,Ttr,Tv\n";
    csv << "0," << mix.T() << "," << mix.Tv() << "\n";

    // ---- time loop: time step of 1 ns as in the paper
    const double dt = 1.0e-9;
    const int nSteps = int(t_end / dt + 0.5);
    // csv rows: every step at the beginning, then at logarithmic intervals
    OutputSchedule output;
    double t = 0.0;

    for (int step = 1; step <= nSteps; step++)
    {
        // Q_VT: energy transferred from the translational reservoir to the
        // vibro-electronic one (Landau-Teller, eq. 8, tau from Mutation++)
        const double Q = sourceVT(mix, rho_s, vibrators, paperTau, vibrationalOnly);

        // eq. 22: only E_ve changes, the total energy E is conserved
        Eve += Q * dt;
        t += dt;

        // new temperatures from the energies (Mutation++ inverts E and E_ve)
        const double energies[2] = {E, Eve};
        mix.setState(rho_s.data(), energies, 0);

        if (output.write(step, nSteps))
        {
            csv << t << "," << mix.T() << "," << mix.Tv() << "\n";
        }
    }
    // ---- end of time loop

    std::cout << csvName << ": T_tr = " << mix.T() << " K, T_ve = " << mix.Tv()
              << " K at t = " << t << " s; expected T_eq = " << T_eq << " K" << std::endl;

    return 0;
}

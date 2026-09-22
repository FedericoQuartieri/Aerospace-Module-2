// 0D heat bath of N2 + O2 (fig. 6 of the paper): vibrational cooling with
// two vibrational temperatures, V-T exchange and, optionally, V-V exchange
//
// Mutation++ has a single vibrational temperature and no V-V exchange,
// so here the loop keeps one reservoir per molecule and uses the functions of
// mutationSources.H species by species (energies, relaxation times);
// the V-V exchange is written following eq. 18 of the paper
//
// usage: Test-N2O2 <VV: on|off> <t_end> <csv_file> [tau]
//   fig 6: Test-N2O2 off 3e-6 output/fig6-noVV.csv paper
//          Test-N2O2 on  3e-6 output/fig6-VV.csv paper
// optional tau: mutation (default, as in the thermo of the solver) or paper
// (eq. 9-17, the choice of report table 1)

#include "mutation++.h"
#include "mutationSources.H"
#include "heatBath.H"

#include <cmath>
#include <cstdlib>
#include <fstream>
#include <iostream>
#include <string>
#include <vector>

// vibro-electronic temperature of species s from its energy e_ve (J/kg),
// by bisection (e_ve grows with Tv)
static double TvFromEve(Mutation::Mixture& mix, const int s, const double e_ve)
{
    double Tlow = 100.0;
    double Thigh = 100000.0;
    for (int iter = 0; iter < 60; iter++)
    {
        const double Tv = 0.5 * (Tlow + Thigh);
        if (speciesEve(mix, s, Tv) > e_ve)
        {
            Thigh = Tv;
        }
        else
        {
            Tlow = Tv;
        }
    }
    return 0.5 * (Tlow + Thigh);
}

int main(int argc, char *argv[])
{
    if (argc != 4 && argc != 5)
    {
        std::cerr << "usage: Test-N2O2 <VV: on|off> <t_end> <csv_file> [paper|mutation]" << std::endl;
        return 1;
    }
    const bool withVV = (std::string(argv[1]) == "on");
    const double t_end = std::atof(argv[2]);
    const char* csvName = argv[3];
    const bool paperTau = (argc == 5 && std::string(argv[4]) == "paper");

    Mutation::MixtureOptions opts("air_5");
    opts.setStateModel("ChemNonEqTTv");
    opts.setThermodynamicDatabase("RRHO");
    opts.setMechanism("none");
    Mutation::Mixture mix(opts);

    const int ns = mix.nSpecies();
    const int iN2 = mix.speciesIndex("N2");
    const int iO2 = mix.speciesIndex("O2");
    const double M_N2 = mix.speciesMw(iN2);
    const double M_O2 = mix.speciesMw(iO2);
    const double R_N2 = Mutation::RU / M_N2;
    const double R_O2 = Mutation::RU / M_O2;

    // ---- initial state of the paper: T_tr = 5000 K, T_v = 30000 K, 1 atm,
    // N2 and O2 in equal parts (same number of moles)
    double T = 5000.0;
    double Tv_N2 = 30000.0;
    double Tv_O2 = 30000.0;

    std::vector<double> Y(ns, 0.0);
    Y[iN2] = M_N2 / (M_N2 + M_O2);
    Y[iO2] = M_O2 / (M_N2 + M_O2);
    const double P_T_Tv[3] = {Mutation::ONEATM, T, Tv_N2};
    mix.setState(Y.data(), P_T_Tv, 2);

    std::vector<double> rho_s(ns);
    mix.densities(rho_s.data());
    const double rho_N2 = rho_s[iN2];
    const double rho_O2 = rho_s[iO2];

    // energies per unit volume (eq. 23): the total energy E is conserved,
    // the vibro-electronic energies of the two reservoirs change in time
    const double E = totalEnergy(mix, rho_s);
    double Eve_N2 = rho_N2 * speciesEve(mix, iN2, Tv_N2);
    double Eve_O2 = rho_O2 * speciesEve(mix, iO2, Tv_O2);

    // initial translational-rotational energy: it grows with T as
    // 2.5 R T per molecule (eq. 3 and 4), and is the way to obtain T afterwards
    const double T0 = T;
    const double Etr0 = E - Eve_N2 - Eve_O2;
    const double rhoCvTr = 2.5 * (rho_N2 * R_N2 + rho_O2 * R_O2);
    // ---- end of initial state

    // final T from energy conservation (the paper reads it from the figure)
    const double T_eq = equilibriumTemperature(mix, rho_s, E);

    // ---- V-T relaxation times (eq. 9-17), one per molecule
    const std::vector<Vibrator> vibrators = makeVibrators(mix);
    const Vibrator* vibN2 = nullptr;
    const Vibrator* vibO2 = nullptr;
    for (const Vibrator& v : vibrators)
    {
        if (v.species == iN2) vibN2 = &v;
        if (v.species == iO2) vibO2 = &v;
    }
    // ---- end of relaxation times

    // ---- constants of the V-V exchange (eq. 18)
    // exchange probability recommended by the paper
    const double P_N2O2 = 0.01;
    // cross section of the N2-O2 pair: the paper does not give the value, this is
    // the one used by hy2Foam, the code of the paper (file thermo2TModel,
    // KnabCoefficients N2_O2, github.com/vincentcasseau/hyStrath, commit 984e3000a5f8)
    const double sigma_N2O2 = 2.667e-19;
    // reduced molar mass of the pair (eq. 14), kg/mol
    const double M_N2O2 = M_N2 * M_O2 / (M_N2 + M_O2);
    // ---- end of V-V constants

    std::ofstream csv(csvName);
    csv << "# T_eq = " << T_eq << " K (energy conservation)\n";
    csv << "t,Ttr,Tv_N2,Tv_O2\n";
    csv << "0," << T << "," << Tv_N2 << "," << Tv_O2 << "\n";

    // ---- time loop: time step of 1 ns as in the paper
    const double dt = 1.0e-9;
    const int nSteps = int(t_end / dt + 0.5);
    OutputSchedule output;
    double t = 0.0;

    for (int step = 1; step <= nSteps; step++)
    {
        // Mutation++ state at T_tr: needed for the relaxation times
        const double temps[2] = {T, Tv_N2};
        mix.setState(rho_s.data(), temps, 1);

        // V-T exchange of each molecule, Landau-Teller (eq. 8), with the driving
        // term e_ve(T_tr) - e_ve(T_v) of the molecule itself
        const double Q_N2_VT = rho_N2 * (speciesEve(mix, iN2, T) - speciesEve(mix, iN2, Tv_N2))
                             / relaxationTime(mix, *vibN2, paperTau);
        const double Q_O2_VT = rho_O2 * (speciesEve(mix, iO2, T) - speciesEve(mix, iO2, Tv_O2))
                             / relaxationTime(mix, *vibO2, paperTau);

        // V-V exchange (eq. 18, with the vibrational energies only), written for
        // each molecule m with its partner l as in the paper and in hy2Foam (KnabVV):
        //   Q_m = NA sigma P sqrt(8 R T / (pi M_ml)) (rho_l/M_l) rho_m
        //         * (e_v,m(T) e_v,l(Tv_l) / e_v,l(T) - e_v,m(Tv_m))
        // the two terms are not opposite (their ratio is -E_v,O2(T)/E_v,N2(T)
        // per mole): the difference goes to the translational mode, because T_tr
        // is obtained from the total energy E, which remains conserved
        double Q_N2_VV = 0.0;
        double Q_O2_VV = 0.0;
        if (withVV)
        {
            const double ev_N2_T = speciesEv(mix, iN2, T);
            const double ev_O2_T = speciesEv(mix, iO2, T);
            const double ev_N2 = speciesEv(mix, iN2, Tv_N2);
            const double ev_O2 = speciesEv(mix, iO2, Tv_O2);
            const double K = Mutation::NA * sigma_N2O2 * P_N2O2
                           * std::sqrt(8.0 * Mutation::RU * T / (Mutation::PI * M_N2O2));
            Q_N2_VV = K * (rho_O2 / M_O2) * rho_N2 * (ev_N2_T * ev_O2 / ev_O2_T - ev_N2);
            Q_O2_VV = K * (rho_N2 / M_N2) * rho_O2 * (ev_O2_T * ev_N2 / ev_N2_T - ev_O2);
        }

        // eq. 22: the two vibro-electronic reservoirs change, E is conserved
        Eve_N2 += (Q_N2_VT + Q_N2_VV) * dt;
        Eve_O2 += (Q_O2_VT + Q_O2_VV) * dt;
        t += dt;

        // new temperatures: T_v from the reservoirs, T_tr from the remaining energy
        Tv_N2 = TvFromEve(mix, iN2, Eve_N2 / rho_N2);
        Tv_O2 = TvFromEve(mix, iO2, Eve_O2 / rho_O2);
        T = T0 + (E - Eve_N2 - Eve_O2 - Etr0) / rhoCvTr;

        if (output.write(step, nSteps))
        {
            csv << t << "," << T << "," << Tv_N2 << "," << Tv_O2 << "\n";
        }
    }
    // ---- end of time loop

    std::cout << csvName << ": T_tr = " << T << " K, T_v,N2 = " << Tv_N2
              << " K, T_v,O2 = " << Tv_O2 << " K at t = " << t
              << " s; expected T_eq = " << T_eq << " K" << std::endl;

    return 0;
}

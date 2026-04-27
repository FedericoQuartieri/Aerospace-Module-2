#include "mutation++.h"

#include "argList.H"
#include "Time.H"
#include "IOstreams.H"
#include "mathematicalConstants.H"
#include "OFstream.H"

#ifdef Log
#undef Log
#endif

using namespace Foam;

int main(int argc, char *argv[])
{
    // WARNING: Set MPP_DATA_DIRECTORY env to ./local_data
    // local_data/thermo/species.xml defines custom energy levels
    // for the species, as opposed to Mutation default ones
    //
    // NOTE: The electronic data (energy levels) shown in the appendix
    // of the paper are not able to reproduce exactly their results
    // for the electronic case (Fig. 4, with or without E_el), while
    // Mutation default ones can.
    //
    // For the non-electronic case (Fig. 3), the closest we can get to their
    // result is to remove all energy levels and set the vibrational
    // temperature to 3371.0, like shown in their appendix. Lowering the
    // vibrational temperature can lead to the same exact result, but
    // there's definitely something weird happening with these constants.

    // NOTE: Electron energy (related to the number of free
    // electrons in the mixture) is not modeled, nor present
    // for N, N2-N and N2-O2

    // Initializing a 5-species air mixture
    Mutation::MixtureOptions opts("air_5");
    // We're using the RRHO two-temperature model
    opts.setStateModel("ChemNonEqTTv");
    opts.setThermodynamicDatabase("RRHO");
    opts.setMechanism("none"); // N2 already has only VT exchange
    Mutation::Mixture mix(opts);

    const int N2_idx = mix.speciesIndex("N2");

    // The mixture is composed of only N2
    std::vector<double> Y_per_specie(mix.nSpecies(), 0.0);
    Y_per_specie[N2_idx] = 1.0;
    // Pressure of the mixture in Pascal
    const double P = Mutation::ONEATM;
     // Trans-rotational temperature
    double T_tr = 30000.0; // 10000.0 for non-electronic case
    // Vibro-electronic temperature
    double T_ve = 1000.0;

    const std::vector<double> P_Ttr_Tve = { Mutation::ONEATM, T_tr, T_ve };
    mix.setState(Y_per_specie.data(), P_Ttr_Tve.data(), 2);

    std::vector<double> rho_per_specie(mix.nSpecies());
    mix.densities(rho_per_specie.data());

    const double dt = 1.0e-9;
    const double end_time = 2.0e-5;
    double t = 0.0;
    int step = 0;

    // For each step:
    // Compute Q_N2,VT                  [eq. 8 from paper]
    //   compute e_ve,N2(T_tr)          [eq. 5 or 5+7?]
    //   compute e_ve,N2(T_ve,N2)       [eq. 5 or 5+7?]
    //   compute tau_N2,VT              [eq. 9]
    //
    // Compute cv_tr and cv_ve
    // Update T_tr
    // Update T_ve
    //
    // Mutationmpp can compute Q_VT and Cvs for us :)

    OFstream out("output/results-N2.csv");
    out << "t,T_tr,T_ve" << endl;
    out << "0," << mix.T() << "," << mix.Tv() << endl;

    while (t < end_time) {
        // Computing the source term of the energy equation
        std::vector<double> Q_sources(mix.nEnergyEqns());
        mix.energyTransferSource(Q_sources.data());
        // Q_sources[0] holds Q_ve = sum_m(Q_m,VT + Q_m,VV + Q_m,CV + ..)

        // The mixture is N2, hence Q_ve = Q_N2,VT since there
        // are no other vibrationally excited molecules

        // Computing specific heat capacities
        std::vector<double> cv_per_specie(mix.nSpecies() * mix.nEnergyEqns());
        mix.getCvsMass(cv_per_specie.data());

        // Updating temperatures as dT/dt = Q_VT / (rho cv)
        // Trans-rotational source term is -Q_VT since total Q = 0
        T_tr += -Q_sources[0] * dt / (cv_per_specie[N2_idx] *
                                      rho_per_specie[N2_idx]);

        T_ve += Q_sources[0] * dt / (rho_per_specie[N2_idx] *
                                     cv_per_specie[mix.nSpecies() + N2_idx]);

        const std::vector<double> temps = { T_tr, T_ve };
        mix.setState(rho_per_specie.data(), temps.data(), 1);

        t += dt;
        step++;

        out << t << "," << mix.T() << "," << mix.Tv() << endl;
    }

    Info << "Output saved to output/results-N2.csv" << endl;

    return 0;
}

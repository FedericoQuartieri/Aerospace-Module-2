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
    // local_data/thermo/species.xml defines only ground-state
    // energy levels, as opposed to Mutation default ones

    // TODO: Test convergence temps with updated species.xml

    // Initializing a 5-species air mixture
    Mutation::MixtureOptions opts("air_5");
    // We're using the RRHO two-temperature model
    opts.setStateModel("ChemNonEqTTv");
    opts.setThermodynamicDatabase("RRHO");
    opts.setMechanism("none"); // N-N2 already has only VT exchange
    Mutation::Mixture mix(opts);

    const int N2_idx = mix.speciesIndex("N2");
    const int N_idx = mix.speciesIndex("N");

    // Number density for both N and N2 [1/m^3]
    const double n = 5.0e22;
    // Specifying density for each specie to get a N-N2 mixture
    std::vector<double> rho_per_specie(mix.nSpecies(), 0.0);
    // Converting number density [1/m^3] to mass density [kg/m^3]
    rho_per_specie[N2_idx] = n * mix.speciesMw(N2_idx) / Mutation::NA;
    rho_per_specie[N_idx] = n * mix.speciesMw(N_idx) / Mutation::NA;

    // Trans-rotational temperature is the same for the whole mixture
    double T_tr = 30000.0;
    // Each vibrationally excited molecule has its own vibro-electronic
    // temperature, here we have only one for N2 (N alone doesn't vibrate)
    double T_ve = 1000.0;

    const std::vector<double> init_temps = { T_tr, T_ve };
    mix.setState(rho_per_specie.data(), init_temps.data(), 1);

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

    OFstream out("output/results-N2N.csv");
    out << "t,T_tr,T_ve" << endl;
    out << "0," << mix.T() << "," << mix.Tv() << endl;

    while (t < end_time) {
        // Computing the source term of the energy equation
        std::vector<double> Q_sources(mix.nEnergyEqns());
        mix.energyTransferSource(Q_sources.data());
        // Q_sources[0] holds Q_ve = sum_m(Q_m,VT + Q_m,VV + Q_m,CV + ..)

        // The mixture is N-N2, hence Q_ve = Q_N2,VT since there
        // are no other vibrationally excited molecules

        // Computing specific heat capacities
        std::vector<double> cv_per_specie(mix.nSpecies() * mix.nEnergyEqns());
        mix.getCvsMass(cv_per_specie.data());

        const double cv_tr_vol =
            cv_per_specie[N2_idx] * rho_per_specie[N2_idx] +
            cv_per_specie[N_idx] * rho_per_specie[N_idx];

        const double cv_ve = cv_per_specie[mix.nSpecies() + N2_idx];
        // Updating temperatures as dT/dt = Q_VT / (rho cv)
        // Trans-rotational source term is -Q_VT since total Q = 0
        T_tr += -Q_sources[0] * dt / cv_tr_vol;
        T_ve += Q_sources[0] * dt / (rho_per_specie[N2_idx] * cv_ve);

        const std::vector<double> temps = { T_tr, T_ve };
        mix.setState(rho_per_specie.data(), temps.data(), 1);

        t += dt;
        step++;

        if (step % 10 == 0) {
            out << t << "," << mix.T() << "," << mix.Tv() << endl;
       }
    }

    Info << "Output saved to output/results-N2N.csv" << endl;

    return 0;
}

#include "argList.H"
#include "Time.H"
#include "IOstreams.H"
#include "mathematicalConstants.H"
#include "OFstream.H"

#ifdef Log
#undef Log
#endif

#include "mutation++.h"

using namespace Foam;

// Constants relevant to the N2-N case
namespace PhysConsts
{
    const double Ru = 8.314462618;
    const double theta_v_N2 = 3371.0;
    const double sigma_v_N2 = 3.0e-21;
    const double k_B = 1.380649e-23;
    const double Na = 6.02214076e23;
}

// Calculates vibrational relaxation time (Millikan-White + Park correction)
scalar calcTauVT
(
    const Mutation::Mixture& mix,
    const label specIdx,
    const scalar P_pa,
    const scalar T,
    const scalar n_tot,
    const scalar theta_v
)
{
    // Pressure in atm for Millikan-White correlation
    const scalar p_atm = P_pa / 101325.0; 
    
    // Molar mass of species i [kg/kmol] -> [kg/mol]
    const scalar Mw_i = mix.speciesMw(specIdx); 
    
    scalar tau_inv_sum = 0.0;

    for (int s = 0; s < mix.nSpecies(); ++s)
    {
        const scalar X_s = mix.X()[s];
        if (X_s < SMALL) continue;

        const scalar Mw_s = mix.speciesMw(s);

        // Reduced mass [kg/mol]
        const scalar mu = (Mw_i * Mw_s) / (Mw_i + Mw_s); 

        // 1. Millikan-White (MW)
        const scalar A = 1.16e-3 * std::sqrt(mu) * std::pow(theta_v, 4.0/3.0);
        const scalar B = 0.015 * std::pow(mu, 0.25);
        
        scalar tau_MW = 1.0;
        if (T > SMALL)
        {
            const scalar arg = A * (std::pow(T, -1.0/3.0) - B) - 18.42;
            // Standard MW is usually scaled by p_atm
            tau_MW = std::exp(arg) / (p_atm > SMALL ? p_atm : 1.0);
        }

        // 2. Park High-T Correction
        const scalar R_spec = PhysConsts::Ru / Mw_i; 
        const scalar c_bar = std::sqrt(8.0 * R_spec * T / constant::mathematical::pi);
        
        // Species specific collision cross-section
        scalar sigma = (mix.speciesName(specIdx) == "N2") ? PhysConsts::sigma_v_N2 : 1.0e-21;
        sigma *= std::pow(50000.0 / std::max(T, 1.0), 2.0);

        const scalar n_s = X_s * n_tot;
        const scalar tau_Park = 1.0 / (c_bar * sigma * std::max(n_s, 1.0e-30));

        // Combined relaxation time
        const scalar tau_combo = tau_MW + tau_Park;
        tau_inv_sum += X_s / tau_combo;
    }

    return (tau_inv_sum > SMALL) ? (1.0 / tau_inv_sum) : 1.0e30;
}


int main(int argc, char *argv[])
{
    argList::addNote("hy2Foam 0D Verification Case");
    argList args(argc, argv);
    Time runTime(Time::controlDictName, args);

    Info<< "\n--- Initializing 0D N2/N Relaxation ---" << endl;

    // --- Mutation++ Initialization
    Mutation::MixtureOptions opts("air_5");
    opts.setStateModel("ChemNonEqTTv");
    opts.setThermodynamicDatabase("RRHO");
    Mutation::Mixture mix(opts);

    Info<< "Mixture: " << mix.nSpecies() << " species loaded." << endl;

    // --- Initial Conditions
    const scalar T_tr_init = 10000.0;
    const scalar T_vib_init = 1000.0;
    const scalar n_init = 5.0e22; // [1/m3] for both species
    
    // Molar masses [kg/kmol] -> we use directly for Mutation++
    const label iN2 = mix.speciesIndex("N2");
    const label iN  = mix.speciesIndex("N");
    
    if (iN2 == -1 || iN == -1) {
        FatalErrorInFunction << "Species N2 or N not found." << exit(FatalError);
    }

    // Calculate Mass Fractions (Y) and Densities (rho)
    std::vector<double> rho_i(mix.nSpecies(), 0.0);
    
    const scalar Na_kmol = PhysConsts::Na;

    rho_i[iN]  = n_init * mix.speciesMw(iN) / Na_kmol;
    rho_i[iN2] = n_init * mix.speciesMw(iN2) / Na_kmol;
    
    scalar rho_total = rho_i[iN] + rho_i[iN2];

    // --- Energy Initialization
    std::vector<double> h_dimless(mix.nSpecies());
    std::vector<double> hv_dimless(mix.nSpecies());

    // Get enthalpies at initial T
    mix.speciesHOverRT(T_tr_init, T_vib_init, T_tr_init, T_vib_init, T_vib_init,
                       h_dimless.data(), nullptr, nullptr, hv_dimless.data(), nullptr, nullptr);

    scalar rhoE = 0.0;
    scalar rhoEv = 0.0;

    for (int i = 0; i < mix.nSpecies(); ++i)
    {
        if (rho_i[i] <= SMALL) continue;
        
        scalar R_spec = PhysConsts::Ru / mix.speciesMw(i);
        scalar h_spec = h_dimless[i] * R_spec * T_tr_init;
        scalar u_spec = h_spec - R_spec * T_tr_init;
        scalar ev_spec = hv_dimless[i] * R_spec * T_tr_init;

        rhoE  += rho_i[i] * u_spec;
        rhoEv += rho_i[i] * ev_spec;
    }

    // Set Initial State
    std::vector<double> energies = {rhoE, rhoEv};
    mix.setState(rho_i.data(), energies.data());

    Info<< "Init State: P=" << mix.P() << " Pa, T=" << mix.T() << " K, Tv=" << mix.Tv() << " K" << endl;

    // --- Output Setup
    OFstream file("results.csv");
    file << "time,T_tr,T_vib,rho_N,rho_N2" << endl;

    // --- Simulation Loop
    const scalar dt = 1e-9;
    const scalar endTime = 2.0e-5;
    scalar t = 0.0;
    label step = 0;

    // Pre-allocation
    std::vector<double> wdot_molar(mix.nSpecies());
    std::vector<double> hv_curr(mix.nSpecies());
    std::vector<double> hv_eq(mix.nSpecies());

    Info<< "\nStarting Time Loop..." << endl;

    while (t < endTime)
    {
        // Get State
        scalar T = mix.T();
        scalar Tv = mix.Tv();
        scalar P = mix.P();
        scalar nTot = mix.numberDensity();

        // Chemistry Source
        mix.netProductionRates(wdot_molar.data());
        
        // Energy Source Terms
        mix.speciesHOverRT(T, Tv, T, Tv, Tv, nullptr, nullptr, nullptr, hv_curr.data(), nullptr, nullptr);
        mix.speciesHOverRT(T, T, T, T, T, nullptr, nullptr, nullptr, hv_eq.data(), nullptr, nullptr);

        scalar Q_chem_v = 0.0;
        scalar Q_VT = 0.0;

        for (int i = 0; i < mix.nSpecies(); ++i)
        {
            scalar Mw = mix.speciesMw(i); // kg/kmol
            scalar wdot_mass = wdot_molar[i] * Mw; // kg/m3/s

            // Vibration-Chemistry coupling
            scalar ev_curr = hv_curr[i] * (PhysConsts::Ru / Mw) * T;
            if (mag(wdot_mass) > SMALL) {
                Q_chem_v += wdot_mass * ev_curr;
            }

            // VT Relaxation
            if (mag(hv_curr[i]) > 1e-10) 
            {
                scalar ev_eq_val = hv_eq[i] * (PhysConsts::Ru / Mw) * T;
                scalar tau = calcTauVT(mix, i, P, T, nTot, PhysConsts::theta_v_N2);
                
                Q_VT += rho_i[i] * (ev_eq_val - ev_curr) / tau;
            }
            
            // Update Density (Explicit Euler)
            rho_i[i] += wdot_mass * dt;
            rho_i[i] = std::max(0.0, rho_i[i]);
        }

        // Update Vibrational Energy
        rhoEv += (Q_VT + Q_chem_v) * dt;

        // Update Mutation++ State
        energies[0] = rhoE; 
        energies[1] = rhoEv;
        mix.setState(rho_i.data(), energies.data());

        // Output
        if (step % 100 == 0)
        {
            file << t << "," << mix.T() << "," << mix.Tv() << "," 
                 << rho_i[iN] << "," << rho_i[iN2] << endl;
        }

        t += dt;
        step++;
    }

    Info<< "\nDone. ExecutionTime = " << runTime.elapsedCpuTime() << " s" << endl;
    Info<< "\nResults saved to results.csv" << endl;
    return 0;
}

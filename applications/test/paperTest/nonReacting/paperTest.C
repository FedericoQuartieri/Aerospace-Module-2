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

        // Reduced mass [g/mol] - Millikan-White formula uses g/mol
        const scalar mu = (Mw_i * Mw_s) / (Mw_i + Mw_s) * 1000.0;  // kg/kmol -> g/mol 

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
    const scalar T_tr_init = 30000.0;  // Changed to 30000 K to match Figure 5
    const scalar T_vib_init = 1000.0;
    const scalar n_init = 5.0e22; // [1/m3] for both species
    
    // Vibrational characteristic temperature for N2 [K]
    const scalar theta_v_N2 = 3390.0;
    
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

    Info<< "hv_dimless[N2]=" << hv_dimless[iN2] << ", hv_dimless[N]=" << hv_dimless[iN] << endl;

    scalar rhoE = 0.0;
    scalar rhoEv = 0.0;

    for (int i = 0; i < mix.nSpecies(); ++i)
    {
        if (rho_i[i] <= SMALL) continue;
        
        scalar R_spec = PhysConsts::Ru / mix.speciesMw(i);
        scalar h_spec = h_dimless[i] * R_spec * T_tr_init;
        
        // Compute vibrational energy directly for molecules
        scalar ev_spec = 0.0;
        if (i == iN2)
        {
            // Quantum harmonic oscillator: ev = R*theta_v / (exp(theta_v/Tv) - 1)
            scalar x = theta_v_N2 / T_vib_init;
            ev_spec = R_spec * theta_v_N2 / (std::exp(x) - 1.0);
            Info<< "N2 vibrational energy: ev_spec=" << ev_spec << " J/kg at Tv=" << T_vib_init << " K" << endl;
        }
        
        scalar u_spec = h_spec - R_spec * T_tr_init - ev_spec;  // Subtract vibrational from total u!

        rhoE  += rho_i[i] * u_spec;   // Now u_spec is only trans-rot-electronic
        
        // Only molecules have vibrational energy (N2), not atoms (N)
        if (i == iN2)
        {
            rhoEv += rho_i[i] * ev_spec;  // This is vibrational only
        }
    }

    Info<< "Before setState: rhoE=" << rhoE << " J/m3, rhoEv=" << rhoEv << " J/m3, Total=" << (rhoE+rhoEv) << " J/m3" << endl;

    // Set Initial State using densities and temperatures (vars=1)
    std::vector<double> temps = {T_tr_init, T_vib_init};
    mix.setState(rho_i.data(), temps.data(), 1);  // vars=1: (rho_i, {T, Tv})
    
    // Get the energies as computed by Mutation++ for consistency
    std::vector<double> energies(2);
    mix.getEnergiesMass(energies.data());
    
    // Convert to volumetric energies (J/m³)
    scalar rho_total_val = rho_i[iN] + rho_i[iN2];
    rhoE = energies[0] * rho_total_val;   // Total internal energy
    rhoEv = energies[1] * rho_total_val;  // Vibrational energy
    
    Info<< "After setState: rhoE=" << rhoE << " J/m3, rhoEv=" << rhoEv << " J/m3, Total=" << (rhoE+rhoEv) << " J/m3" << endl;
    Info<< "Init State: P=" << mix.P() << " Pa, T=" << mix.T() << " K, Tv=" << mix.Tv() << " K" << endl;

    // --- Output Setup
    OFstream file("results.csv");
    file << "time,T_tr,T_vib,rho_N,rho_N2" << endl;

    // --- Simulation Loop
    const scalar dt = 1e-9;
    const scalar endTime = 2.0e-5;
    scalar t = 0.0;
    label step = 0;
    
    // Current temperatures (updated directly)
    scalar T = mix.T();
    scalar Tv = mix.Tv();

    Info<< "\nStarting Time Loop..." << endl;

    while (t < endTime)
    {
        // Get State
        scalar P = mix.P();
        scalar nTot = mix.numberDensity();

        // VT Relaxation - only for N2 (molecules with vibration)
        scalar Q_VT = 0.0;
        scalar cv_vib = 0.0;

        // Only N2 has vibrational energy
        {
            scalar Mw = mix.speciesMw(iN2); // kg/kmol
            scalar R_spec = PhysConsts::Ru / Mw;  // J/(kg·K)
            
            // Quantum harmonic oscillator energy: ev = R*theta_v / (exp(theta_v/T) - 1)
            scalar x_v = theta_v_N2 / std::max(Tv, 1.0);
            scalar x_eq = theta_v_N2 / std::max(T, 1.0);
            
            scalar ev_curr = R_spec * theta_v_N2 / (std::exp(x_v) - 1.0);   // at Tv
            scalar ev_eq = R_spec * theta_v_N2 / (std::exp(x_eq) - 1.0);    // at T (equilibrium)
            
            // Vibrational specific heat: cv_v = R * (theta/T)^2 * exp(theta/T) / (exp(theta/T)-1)^2
            scalar exp_x = std::exp(x_v);
            cv_vib = R_spec * x_v * x_v * exp_x / ((exp_x - 1.0) * (exp_x - 1.0));
            
            scalar tau = calcTauVT(mix, iN2, P, T, nTot, theta_v_N2);
            
            Q_VT = rho_i[iN2] * (ev_eq - ev_curr) / tau;  // [W/m³]
        }

        // Total cv: N (3/2 R) + N2 (5/2 R for trans-rot)
        scalar R_N = PhysConsts::Ru / mix.speciesMw(iN);
        scalar R_N2 = PhysConsts::Ru / mix.speciesMw(iN2);
        scalar rho_cv_tr = rho_i[iN] * 1.5 * R_N + rho_i[iN2] * 2.5 * R_N2;  // trans-rot only
        scalar rho_cv_vib = rho_i[iN2] * cv_vib;  // vibrational only

        // Update temperatures directly (energy conservation)
        // dT/dt = -Q_VT / (rho*cv_tr)  (energy leaves trans-rot)
        // dTv/dt = Q_VT / (rho_N2*cv_vib)  (energy enters vibration)
        T  -= Q_VT * dt / rho_cv_tr;
        Tv += Q_VT * dt / rho_cv_vib;

        // Update Mutation++ State with new temperatures
        std::vector<double> temps_new = {T, Tv};
        mix.setState(rho_i.data(), temps_new.data(), 1);

        // Output - More frequent at the beginning (log-spaced-like)
        bool writeOutput = false;
        if (t < 1e-7) {
            writeOutput = (step % 10 == 0);  // Every 10 steps until 1e-7 s
        } else if (t < 1e-6) {
            writeOutput = (step % 50 == 0);  // Every 50 steps until 1e-6 s
        } else {
            writeOutput = (step % 100 == 0); // Every 100 steps after
        }
        
        if (writeOutput)
        {
            file << t << "," << T << "," << Tv << "," 
                 << rho_i[iN] << "," << rho_i[iN2] << endl;
        }

        t += dt;
        step++;
    }

    Info<< "\nDone. ExecutionTime = " << runTime.elapsedCpuTime() << " s" << endl;
    Info<< "Final State: T=" << T << " K, Tv=" << Tv << " K" << endl;
    Info<< "\nResults saved to results.csv" << endl;
    return 0;
}

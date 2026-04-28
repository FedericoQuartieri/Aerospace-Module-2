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

// Converts dimensionless enthalpy to energy (h_v, h_el only)
double hToE(double h, double T, double M)
{
    // speciesHOverRT gives h/RT, moreover h_v = e_v, h_el = e_el
    // while h_tr = e_rt - RT
    return h * (Mutation::RU * T / M); // Rm here
}

// Calculates vibrational relaxation time of the i-th specie
// of a mixture (Millikan-White + Park correction, eqs. 9-17 from the paper)
double computeTauVT(const Mutation::Mixture& mix,
                    int i, // Species index
                    double p_pa,
                    double T_tr,
                    double n_tot, // Number density of the mixture
                    double theta_v)
{
    const double p_atm = p_pa / Mutation::ONEATM;
    const double M_i = mix.speciesMw(i);

    double inv_tau_vt = 0.0;
    for (int s = 0; s < mix.nSpecies(); ++s)
    {
        const double M_s = mix.speciesMw(s);
        // Reduced molar mass in g/mol
        const double M_is = (M_i * M_s) / (M_i + M_s) * 1000.0;

        // Millikan-White contribution
        double A = 1.16e-3 * std::sqrt(M_is) * std::pow(theta_v, 4.0 / 3.0);
        double B = 0.015 * std::pow(M_is, 0.25);
        double tau_MW =
            std::exp(A * (std::pow(T_tr, -1.0 / 3.0) - B) - 18.42) / p_atm;

        const double X_s = mix.X()[s];
        // Park contribution
        double R_i = Mutation::RU / M_i;
        double c_bar = std::sqrt(8.0 * R_i * T_tr / Mutation::PI);
        // NOTE: The paper used 3.0e-21 for the N2-O2 case
        double sigma = 3.0e-21 * std::pow(50000.0 / T_tr, 2.0);
        // WARNING: We assume n_i,s to be the number density of s
        double tau_Park = 1.0 / (c_bar * sigma * X_s * n_tot);

        inv_tau_vt += X_s / (tau_MW + tau_Park);
    }
    // sum_s(X_s) = 1.0
    return 1.0 / inv_tau_vt;
}

int main(int argc, char *argv[])
{
    // Initializing a 5-species air mixture
    Mutation::MixtureOptions opts("air_5");
    // We're using the RRHO two-temperature model
    opts.setStateModel("ChemNonEqTTv");
    opts.setThermodynamicDatabase("RRHO");
     // Disable chemical reactions, consider only V-T and V-V exchanges
    opts.setMechanism("none");
    Mutation::Mixture mix(opts);

    const int ns = mix.nSpecies();

    const int N2_idx = mix.speciesIndex("N2");
    const int O2_idx = mix.speciesIndex("O2");

    const int iN2 = mix.speciesIndex("N2");
    const int iO2 = mix.speciesIndex("O2");
    // Vibrational temperatures in K (from paper appendix)
    // TODO: Can we read these from the species.xml?
    const double theta_v_N2 = 3371.0;
    const double theta_v_O2 = 2256.0;

    // Molar masses per specie
    const double M_N2 = mix.speciesMw(N2_idx);
    const double M_O2 = mix.speciesMw(O2_idx);
    // The mixture is composed of N2 and O2 in equal proportions
    // (equal number of moles). Mutation needs the mass fraction Y.
    std::vector<double> Y_per_specie(ns, 0.0);
    Y_per_specie[N2_idx] = M_N2 / (M_N2 + M_O2);
    Y_per_specie[O2_idx] = M_O2 / (M_N2 + M_O2);

    // Trans-rotational temperature is the same for the whole mixture
    double T_tr = 5000.0;
    // Each vibrationally excited molecule has its own vibro-electronic
    // temperature
    double T_ve_N2 = 30000.0;
    double T_ve_O2 = 30000.0;

    // Pressure is 1 atm
    std::vector<double> state_vars = { Mutation::ONEATM, T_tr, T_ve_N2 };
    mix.setState(Y_per_specie.data(), state_vars.data(), 2);

    std::vector<double> rho_per_specie(ns);
    mix.densities(rho_per_specie.data());

    const double dt = 1.0e-9;
    const double end_time = 1.0e-6;
    double t = 0.0;
    int step = 0;

    OFstream out("output/results-N2O2-VT.csv");
    out << "t,T_tr,T_ve_N2,T_ve_O2" << endl;
    out << "0," << mix.T() << "," << mix.Tv() << "," << mix.Tv() << endl;

    // NOTE: mix.energyTransferSource() returns a unique
    // Q_ve = sum_m(Q_m,VT + ..). Since with N2-O2 we have two
    // vibrational temperatures, we need to manually compute
    // the source terms for each specie by implementing the
    // formulas from the paper.

    while (t < end_time) {
        std::vector<double> h_v(mix.nSpecies());
        std::vector<double> h_el(mix.nSpecies());

        // Compute target energies e_m,ve(T_tr) (at equilibrium)
        mix.speciesHOverRT(T_tr, T_tr, T_tr, T_tr, T_tr,
                           NULL, NULL, NULL, h_v.data(), h_el.data(), NULL);
        double e_N2_eq = hToE(h_v[N2_idx] + h_el[N2_idx], T_tr, M_N2);
        double e_O2_eq = hToE(h_v[O2_idx] + h_el[O2_idx], T_tr, M_O2);

        // Compute current energies e_m,ve(T_m,ve)
        mix.speciesHOverRT(T_tr, T_tr, T_tr, T_ve_N2, T_ve_N2,
                           NULL, NULL, NULL, h_v.data(), h_el.data(), NULL);
        double e_N2 = hToE(h_v[N2_idx] + h_el[N2_idx], T_tr, M_N2);
        mix.speciesHOverRT(T_tr, T_tr, T_tr, T_ve_O2, T_ve_O2,
                           NULL, NULL, NULL, h_v.data(), h_el.data(), NULL);
        double e_O2 = hToE(h_v[O2_idx] + h_el[O2_idx], T_tr, M_O2);

        // NOTE: For VV, e_v is used, so perhaps you should
        // compute e_v and e_el separately

        // Compute relaxation times
        double p = mix.P();
        double n_tot = mix.numberDensity();
        double tau_N2 = computeTauVT(mix, N2_idx, p, T_tr, n_tot, theta_v_N2);
        double tau_O2 = computeTauVT(mix, O2_idx, p, T_tr, n_tot, theta_v_O2);
        // Compute source terms Q_m,VT
        double Q_N2_VT = rho_per_specie[N2_idx] * (e_N2_eq - e_N2) / tau_N2;
        double Q_O2_VT = rho_per_specie[O2_idx] * (e_O2_eq - e_O2) / tau_O2;

        /*
        // VV Source Term
        // note that e_v is needed here, not e_ve
        double M_red = (mix.speciesMw(iN2) * mix.speciesMw(iO2)) /
                       (mix.speciesMw(iN2) + mix.speciesMw(iO2));
        double K_vv = Mutation::NA * 4.0e-19 * 0.01 * 
            std::sqrt(8.0 * Mutation::RU * T_tr / (Mutation::PI * M_red));

        // Detailed balance term: e_eq_N2 * e_curr_O2 / e_eq_O2
        double Q_N2_VV = K_vv * (rho_per_specie[iO2]/mix.speciesMw(iO2)) *
                         rho_per_specie[iN2] *
                         (e_eq_N2 * e_curr_O2 / std::max(e_eq_O2, 1e-10) - e_curr_N2);
        double Q_O2_VV = -Q_N2_VV;
        */

        // Currently ignoring VV exchange
        double Q_N2_VV = 0;
        double Q_O2_VV = 0;

        // Compute cv_m,ve and cv_tr from dimensionless cvs
        std::vector<double> cv_dim_N2(ns), cv_dim_O2(ns);
        mix.speciesCvOverR(T_tr, T_ve_N2, T_tr, T_ve_N2, T_ve_N2,
                           NULL, NULL, NULL, cv_dim_N2.data(), NULL);
        mix.speciesCvOverR(T_tr, T_ve_O2, T_tr, T_ve_O2, T_ve_O2,
                           NULL, NULL, NULL, cv_dim_O2.data(), NULL);

        double cv_ve_N2 = cv_dim_N2[N2_idx] * (Mutation::RU / M_N2);
        double cv_ve_O2 = cv_dim_O2[O2_idx] * (Mutation::RU / M_O2);

        // TODO: What is going on here?
        double rho_cv_tr = rho_per_specie[N2_idx] *
                           (2.5 * Mutation::RU / M_N2) +
                           rho_per_specie[O2_idx] *
                           (2.5 * Mutation::RU / M_O2);

        T_tr -= (Q_N2_VT + Q_O2_VT) * dt / rho_cv_tr;
        T_ve_N2 += (Q_N2_VT + Q_N2_VV) * dt /
                   (rho_per_specie[N2_idx] * cv_ve_N2);
        T_ve_O2 += (Q_O2_VT + Q_O2_VV) * dt / 
                   (rho_per_specie[O2_idx] * cv_ve_O2);

        // Refresh mix state for next p and n_tot lookups
        // T_ve shouldn't influence them, so any value should be fine
        std::vector<double> temps = { T_tr, T_ve_N2 };
        mix.setState(rho_per_specie.data(), temps.data(), 1);

        if (step % 10 == 0) {
            out << t << "," << T_tr << "," << T_ve_N2 << "," << T_ve_O2 << endl;
        }

        t += dt;
    }

    Info << "Output saved to output/results-N2O2-VT.csv" << endl;

    return 0;
}

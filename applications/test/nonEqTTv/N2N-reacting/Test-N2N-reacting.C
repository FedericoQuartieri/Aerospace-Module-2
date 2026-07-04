#include "mutation++.h"

#include "argList.H"
#include "Time.H"
#include "IOstreams.H"
#include "mathematicalConstants.H"
#include "OFstream.H"

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#ifdef Log
#undef Log
#endif

using namespace Foam;

// Reacting N2-N adiabatic heat bath (paper Sec. 3.4, Figs. 7 and 8).
//
// Chemistry: irreversible N2 + N2 -> 2N + N2 (Park 1993 rates), coupled
// with V-T relaxation (Millikan-White + Park) and non-preferential
// chemistry-vibration coupling (Candler): Q_CV = sum_i e_v,i * wdot_i.
//
// Integration uses the CONSERVATIVE energy formulation of the paper
// (eq. 22-26): the state vector is (rho_s, rhoEve, rhoEtot).  For an
// adiabatic constant-volume heat bath rhoEtot is constant, rho_s evolves
// with wdot_s and rhoEve with Q_ve = Q_VT + Q_CV.  Temperatures are then
// decoded with mix.setState(rho, {rhoEtot, rhoEve}, 0) (solveEnergies).
// This also pre-validates the decode path the CFD solver will use.
//
// Two modes (same integrator, different source providers):
//
//   mpp     mechanism "N2_diss_park" loaded in Mutation++;
//           wdot = mix.netProductionRates() and Q_ve =
//           mix.energyTransferSource() = OmegaVT + OmegaCV.
//           The vendored Mutation++ evaluates dissociation forward rates
//           at Tf = sqrt(T*Tv), i.e. Park exponent q = 0.5 (hard-coded
//           in kinetics/RateManager.cpp).  This is what the CFD solver
//           will use in milestone 3.
//
//   manual  mechanism "none": Mutation++ provides pure Q_VT (OmegaCV = 0
//           without reactions); wdot is computed here with the law of
//           mass action at Tf = T^q * Tv^(1-q) (q from argv, default 0.7
//           as in the paper) and Q_CV is added following the same
//           Candler formula used by Mutation++'s OmegaCV.
//
// Running manual with q = 0.5 must reproduce mpp: that cross-check
// validates both source paths against each other.
//
// Usage: Test-N2N-reacting <T_tr> <T_ve> <mpp|manual> [q]
//
// NOTE: run with MPP_DATA_DIRECTORY=mutation-data-noel for the
// non-electronic configuration (paper Sec. 3.4 compares against
// dsmcFoam, which has no electronic energy mode).

namespace
{

// Park 1993 rate for N2+N2 -> 2N+N2 (paper Table 2): A in cm^3/(mol s)
constexpr double parkA = 7.0e21;
constexpr double parkBeta = -1.6;
constexpr double parkTa = 113200.0;

// Law of mass action for the single irreversible reaction, evaluated at
// the Park controlling temperature Tc = T^q * Tv^(1-q).
// Returns wdot for every species in kg/(m^3 s) (only N2 and N nonzero).
void manualProductionRates
(
    const Mutation::Mixture& mix,
    const std::vector<double>& rho,
    int iN2,
    int iN,
    double Ttr,
    double Tve,
    double q,
    std::vector<double>& wdot
)
{
    const double Tc = std::pow(Ttr, q)*std::pow(Tve, 1.0 - q);

    // cm^3/(mol s) -> m^3/(mol s)
    const double kf = 1.0e-6*parkA*std::pow(Tc, parkBeta)*std::exp(-parkTa/Tc);

    const double concN2 = rho[iN2]/mix.speciesMw(iN2);  // mol/m^3
    const double rate = kf*concN2*concN2;               // mol/(m^3 s)

    std::fill(wdot.begin(), wdot.end(), 0.0);
    wdot[iN2] = -mix.speciesMw(iN2)*rate;
    wdot[iN] = 2.0*mix.speciesMw(iN)*rate;
}

// Non-preferential (Candler) chemistry-vibration coupling, mirroring
// Mutation++'s OmegaCV::compute_source_Candler() exactly:
//   Q_CV = sum_i (h_v,i/RT) * RU * T / Mw_i * wdot_i   [W/m^3]
double manualOmegaCV
(
    Mutation::Mixture& mix,
    const std::vector<double>& wdot
)
{
    const int ns = mix.nSpecies();
    std::vector<double> hvOverRT(ns);
    mix.speciesHOverRT(NULL, NULL, NULL, hvOverRT.data(), NULL, NULL);

    double sum = 0.0;
    for (int i = 0; i < ns; ++i)
    {
        sum += hvOverRT[i]*wdot[i]/mix.speciesMw(i);
    }

    return sum*mix.T()*Mutation::RU;
}

} // namespace


int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        Info<< "Usage: Test-N2N-reacting <T_tr> <T_ve> <mpp|manual> [q]"
            << endl;
        return 1;
    }

    const double TtrInit = std::atof(argv[1]);
    const double TveInit = std::atof(argv[2]);
    const std::string mode = argv[3];
    const double q = (argc > 4 ? std::atof(argv[4]) : 0.7);

    const bool useMpp = (mode == "mpp");
    if (!useMpp && mode != "manual")
    {
        Info<< "Unknown mode '" << mode.c_str() << "'" << endl;
        return 1;
    }

    Mutation::MixtureOptions opts("air_5");
    opts.setStateModel("ChemNonEqTTv");
    opts.setThermodynamicDatabase("RRHO");
    opts.setMechanism(useMpp ? "N2_diss_park" : "none");
    Mutation::Mixture mix(opts);

    const int ns = mix.nSpecies();
    const int iN2 = mix.speciesIndex("N2");
    const int iN = mix.speciesIndex("N");

    // Equal number densities of 5.0e22 m^-3 for N2 and N (paper Sec. 3.4)
    const double nInit = 5.0e22;
    std::vector<double> rho(ns, 0.0);
    rho[iN2] = nInit*mix.speciesMw(iN2)/Mutation::NA;
    rho[iN] = nInit*mix.speciesMw(iN)/Mutation::NA;

    const double nTotInit = 2.0*nInit;  // normalisation of Fig. 7b/8b

    const std::vector<double> TInit = {TtrInit, TveInit};
    mix.setState(rho.data(), TInit.data(), 1);

    // Conserved quantities from the initial state: rhoEtot stays constant
    // (adiabatic constant-volume bath), rhoEve is integrated in time.
    std::vector<double> eSpecies(2*ns);
    mix.getEnergiesMass(eSpecies.data());

    double rhoEtot = 0.0;
    double rhoEve = 0.0;
    for (int i = 0; i < ns; ++i)
    {
        rhoEtot += rho[i]*eSpecies[i];
        rhoEve += rho[i]*eSpecies[ns + i];
    }

    const double rhoTotInit = rho[iN2] + rho[iN];

    Info<< "Reacting N2-N heat bath, mode=" << mode.c_str();
    if (!useMpp) Info<< " (q=" << q << ")";
    Info<< nl
        << "  T_tr(0)=" << TtrInit << " K, T_ve(0)=" << TveInit << " K" << nl
        << "  p(0)=" << mix.P() << " Pa"
        << ", rhoEtot=" << rhoEtot << " J/m^3" << endl;

    // Cross-check at t=0: Mutation++ wdot vs law of mass action at the
    // exponent Mutation++ uses internally (q=0.5).  Only meaningful in
    // mpp mode; a large mismatch means broken units or wrong mechanism.
    std::vector<double> wdot(ns, 0.0);
    if (useMpp)
    {
        mix.netProductionRates(wdot.data());

        std::vector<double> wdotRef(ns, 0.0);
        manualProductionRates
        (
            mix, rho, iN2, iN, TtrInit, TveInit, 0.5, wdotRef
        );

        Info<< "  step-0 check, wdot_N2: mutation=" << wdot[iN2]
            << " manual(q=0.5)=" << wdotRef[iN2]
            << " kg/(m^3 s)" << endl;
    }

    const double dt = 1.0e-9;
    const double endTime = 1.0e-3;
    double t = 0.0;

    // File name encodes the configuration, e.g.
    //   results-N2N-reacting-30000-1000-mpp-park05.csv
    //   results-N2N-reacting-30000-1000-manual-park07.csv
    char tag[64];
    if (useMpp)
    {
        std::snprintf(tag, sizeof(tag), "mpp-park05");
    }
    else
    {
        std::snprintf(tag, sizeof(tag), "manual-park%02d",
                      static_cast<int>(std::lround(q*10.0)));
    }

    char fname[256];
    std::snprintf
    (
        fname, sizeof(fname),
        "output/results-N2N-reacting-%d-%d-%s.csv",
        static_cast<int>(std::lround(TtrInit)),
        static_cast<int>(std::lround(TveInit)),
        tag
    );

    OFstream out(fname);
    out << "t,T_tr,T_ve,nN2_over_n0,nN_over_n0" << endl;

    auto writeRow = [&]()
    {
        const double nN2 = rho[iN2]/mix.speciesMw(iN2)*Mutation::NA;
        const double nN = rho[iN]/mix.speciesMw(iN)*Mutation::NA;
        out << t << "," << mix.T() << "," << mix.Tv() << ","
            << nN2/nTotInit << "," << nN/nTotInit << endl;
    };

    writeRow();

    // Log-spaced output: ~1200 rows over six decades
    double nextWrite = dt;
    const double writeFactor = 1.012;

    std::vector<double> Q(mix.nEnergyEqns(), 0.0);

    while (t < endTime)
    {
        // --- source terms at the current state ---
        if (useMpp)
        {
            // wdot from Mutation++ (Tf = sqrt(T*Tv));
            // energyTransferSource = OmegaVT + OmegaCV
            mix.netProductionRates(wdot.data());
            mix.energyTransferSource(Q.data());
        }
        else
        {
            // pure V-T from Mutation++ (no mechanism -> OmegaCV = 0),
            // chemistry and CV coupling computed manually at exponent q
            mix.energyTransferSource(Q.data());
            manualProductionRates
            (
                mix, rho, iN2, iN, mix.T(), mix.Tv(), q, wdot
            );
            Q[0] += manualOmegaCV(mix, wdot);
        }

        // --- explicit Euler update of the conserved variables ---
        for (int i = 0; i < ns; ++i)
        {
            rho[i] = std::max(rho[i] + wdot[i]*dt, 0.0);
        }
        rhoEve += Q[0]*dt;

        const std::vector<double> energies = {rhoEtot, rhoEve};
        mix.setState(rho.data(), energies.data(), 0);

        t += dt;

        if (t >= nextWrite || t >= endTime)
        {
            writeRow();
            nextWrite = max(nextWrite*writeFactor, nextWrite + dt);
        }
    }

    // --- conservation / consistency report ---
    const double rhoTotFinal = rho[iN2] + rho[iN];

    mix.getEnergiesMass(eSpecies.data());
    double rhoEtotDecoded = 0.0;
    for (int i = 0; i < ns; ++i)
    {
        rhoEtotDecoded += rho[i]*eSpecies[i];
    }

    const double nN2 = rho[iN2]/mix.speciesMw(iN2)*Mutation::NA;
    const double nN = rho[iN]/mix.speciesMw(iN)*Mutation::NA;

    Info<< "Final state at t=" << t << " s:" << nl
        << "  T_tr=" << mix.T() << " K, T_ve=" << mix.Tv() << " K" << nl
        << "  n_N2/n0=" << nN2/nTotInit << ", n_N/n0=" << nN/nTotInit << nl
        << "  mass conservation error="
        << mag(rhoTotFinal - rhoTotInit)/rhoTotInit << nl
        << "  energy decode error="
        << mag(rhoEtotDecoded - rhoEtot)/mag(rhoEtot) << nl
        << "Output saved to " << fname << endl;

    return 0;
}

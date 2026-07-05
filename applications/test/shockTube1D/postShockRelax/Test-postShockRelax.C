#include "mutation++.h"

#include "argList.H"
#include "IOstreams.H"
#include "OFstream.H"

#include <cmath>
#include <cstdlib>
#include <string>
#include <vector>

#ifdef Log
#undef Log
#endif

using namespace Foam;

// Steady post-shock relaxation reference for the 1D shock tube validation
// (milestone 4).
//
// Given the pre-shock state (p1, T1, pure N2, Tve1 = T1) and the shock
// speed us measured from the CFD run, this tool:
//
//  1. computes the FROZEN Rankine-Hugoniot jump (eve and composition are
//     frozen across the thin shock front) using the same Mutation++
//     thermodynamics the solver bridge uses;
//
//  2. integrates the steady 1D Euler equations with relaxation source
//     terms behind the shock (shock-attached frame, x = 0 at the front):
//
//         m   = rho*u                      (constant)
//         P   = p + m*u                    (constant)
//         H0  = h + u^2/2                  (constant)
//         d(eve)/dx = Q_ve / m             (V-T + chemistry-vibration)
//         d(Y_i)/dx = wdot_i / m           (Mutation++ kinetics)
//
//     At each x the algebraic state (rho, u, p, T) is recovered from
//     (m, P, H0; eve, Y) by a bisection on u (subsonic branch).
//
// Sources are evaluated with mix.setState(rho_i, {T, Tve}, 1) followed by
// energyTransferSource (full Q_ve: OmegaVT + OmegaCV + OmegaCElec) and
// netProductionRates — exactly the physics of the CFD solver bridge.
//
// Usage: Test-postShockRelax <p1 Pa> <T1 K> <us m/s> [xEnd m] [dx m]
// Output: output/postShockRelax-<p1>-<T1>-<us>.csv
//         (x, T, Tve, rho, u_shockframe, p, Y_N2, Y_N)

namespace
{

struct MixFns
{
    Mutation::Mixture& mix;
    const int ns;
    const int iN2;
    const int iN;
    std::vector<double> w1, w2, w3, w4;

    MixFns(Mutation::Mixture& m)
    :
        mix(m),
        ns(m.nSpecies()),
        iN2(m.speciesIndex("N2")),
        iN(m.speciesIndex("N")),
        w1(ns), w2(ns), w3(ns), w4(ns)
    {}

    // Mixture gas constant [J/kg/K] for mass fractions y
    double R(const std::vector<double>& y) const
    {
        double r = 0.0;
        for (int i = 0; i < ns; ++i)
        {
            r += y[i]*Mutation::RU/mix.speciesMw(i);
        }
        return r;
    }

    // Trans-rotational internal energy [J/kg] at T (pure function)
    double etr(const std::vector<double>& y, double T)
    {
        mix.speciesHOverRT
        (
            T, T, T, T, T,
            nullptr, w1.data(), w2.data(), nullptr, nullptr, nullptr
        );
        double e = 0.0;
        for (int i = 0; i < ns; ++i)
        {
            e += y[i]*(w1[i] + w2[i] - 1.0)*Mutation::RU*T/mix.speciesMw(i);
        }
        return e;
    }

    // Vibro-electronic energy [J/kg] at (T, Tv) (pure function)
    double eve(const std::vector<double>& y, double T, double Tv)
    {
        mix.speciesHOverRT
        (
            T, Tv, T, Tv, Tv,
            nullptr, nullptr, nullptr, w3.data(), w4.data(), nullptr
        );
        double e = 0.0;
        for (int i = 0; i < ns; ++i)
        {
            e += y[i]*(w3[i] + w4[i])*Mutation::RU*T/mix.speciesMw(i);
        }
        return e;
    }

    // Formation energy [J/kg] for mass fractions y
    double ef(const std::vector<double>& y)
    {
        const double Tref = 298.15;
        mix.speciesHOverRT
        (
            Tref, Tref, Tref, Tref, Tref,
            nullptr, nullptr, nullptr, nullptr, nullptr, w1.data()
        );
        double e = 0.0;
        for (int i = 0; i < ns; ++i)
        {
            e += y[i]*w1[i]*Mutation::RU*Tref/mix.speciesMw(i);
        }
        return e;
    }

    // Trans-rotational cv [J/kg/K] (constant for RRHO)
    double cvTr(const std::vector<double>& y, double T)
    {
        mix.speciesCpOverR
        (
            T, T, T, T, T,
            nullptr, w1.data(), w2.data(), nullptr, nullptr
        );
        double cv = 0.0;
        for (int i = 0; i < ns; ++i)
        {
            cv += y[i]*(w1[i] + w2[i] - 1.0)*Mutation::RU/mix.speciesMw(i);
        }
        return cv;
    }

    // Absolute specific enthalpy [J/kg]: h = etr + eve + ef + p/rho,
    // with p/rho = R*T
    double h(const std::vector<double>& y, double T, double eveVal)
    {
        return etr(y, T) + eveVal + ef(y) + R(y)*T;
    }

    // Invert T from trans-rotational energy target
    double TFromEtr(const std::vector<double>& y, double target, double Tg)
    {
        double T = max(Tg, 50.0);
        for (int it = 0; it < 100; ++it)
        {
            const double f = etr(y, T) - target;
            const double cv = max(cvTr(y, T), 1e-8);
            double dT = -f/cv;
            dT = min(max(dT, -0.5*T), 0.5*T);
            T += dT;
            if (mag(dT) < 1e-8*T) break;
        }
        return T;
    }
};

// Given conserved (m, P, H0) and relaxing (eve, Y), recover the subsonic
// (post-shock) state by bisection on u:
//   rho = m/u,  p = P - m*u,  T = p*u/(m*R),
//   residual(u) = h(T, eve) + u^2/2 - H0
// residual is monotonic in u on the subsonic branch.
bool solveState
(
    MixFns& fns,
    const std::vector<double>& y,
    double m,
    double P,
    double H0,
    double eveVal,
    double& rho,
    double& u,
    double& p,
    double& T
)
{
    const double R = fns.R(y);

    auto residual = [&](double uu)
    {
        const double pp = P - m*uu;
        const double TT = pp*uu/(m*R);
        return fns.h(y, TT, eveVal) + 0.5*uu*uu - H0;
    };

    // The residual has two roots: the trivial pre-shock one (u = us) and
    // the subsonic post-shock one (u2 < us). Bracket the subsonic root:
    // residual < 0 for u -> 0 and > 0 between the two roots. Sample
    // upwards from small u until the sign flips, then bisect.
    const double uMax = 0.999*P/m;  // p > 0 requires u < P/m

    double uLo = 1e-3;
    double fLo = residual(uLo);

    double uHi = -1;
    double fHi = 0;
    for (int k = 1; k <= 200; ++k)
    {
        const double uu = uMax*k/201.0;
        const double ff = residual(uu);
        if (ff*fLo < 0)
        {
            uHi = uu;
            fHi = ff;
            break;
        }
        uLo = uu;
        fLo = ff;
    }

    if (uHi < 0)
    {
        return false;
    }

    for (int it = 0; it < 200; ++it)
    {
        const double uMid = 0.5*(uLo + uHi);
        const double fMid = residual(uMid);

        if (fMid*fLo <= 0)
        {
            uHi = uMid;
            fHi = fMid;
        }
        else
        {
            uLo = uMid;
            fLo = fMid;
        }

        if ((uHi - uLo) < 1e-12*uHi)
        {
            break;
        }
    }

    u = 0.5*(uLo + uHi);
    rho = m/u;
    p = P - m*u;
    T = p*u/(m*R);

    return true;
}

} // namespace


int main(int argc, char *argv[])
{
    if (argc < 4)
    {
        Info<< "Usage: Test-postShockRelax <p1 Pa> <T1 K> <us m/s>"
            << " [xEnd m] [dx m]" << endl;
        return 1;
    }

    const double p1 = std::atof(argv[1]);
    const double T1 = std::atof(argv[2]);
    const double us = std::atof(argv[3]);
    const double xEnd = (argc > 4 ? std::atof(argv[4]) : 0.2);
    const double dx = (argc > 5 ? std::atof(argv[5]) : 1e-5);

    Mutation::MixtureOptions opts("air_5");
    opts.setStateModel("ChemNonEqTTv");
    opts.setThermodynamicDatabase("RRHO");
    opts.setMechanism("N2_diss_park");
    Mutation::Mixture mix(opts);

    MixFns fns(mix);
    const int ns = fns.ns;

    // --- pre-shock state: pure N2 in thermal equilibrium ---
    std::vector<double> y(ns, 0.0);
    y[fns.iN2] = 1.0;

    const double R1 = fns.R(y);
    const double rho1 = p1/(R1*T1);
    const double eve1 = fns.eve(y, T1, T1);

    // Conserved fluxes in the shock-attached frame (pre-shock speed = us)
    const double m = rho1*us;
    const double P = p1 + m*us;
    const double H0 = fns.h(y, T1, eve1) + 0.5*us*us;

    // --- frozen post-shock state (eve and Y frozen across the front) ---
    double rho, u, p, T;
    if (!solveState(fns, y, m, P, H0, eve1, rho, u, p, T))
    {
        Info<< "Frozen RH jump: bracketing failed (shock too weak?)" << endl;
        return 1;
    }

    Info<< "Pre-shock : p=" << p1 << " Pa, T=" << T1 << " K, rho=" << rho1
        << " kg/m^3, us=" << us << " m/s (M=" << us/std::sqrt(1.4*R1*T1)
        << ")" << nl
        << "Frozen jump: p=" << p << " Pa, T=" << T << " K, rho=" << rho
        << " kg/m^3, u2(shock frame)=" << u << " m/s" << nl
        << "  frozen ratios: p2/p1=" << p/p1 << " T2/T1=" << T/T1
        << " rho2/rho1=" << rho/rho1 << endl;

    // --- integrate the relaxation zone ---
    char fname[256];
    std::snprintf
    (
        fname, sizeof(fname), "output/postShockRelax-%d-%d-%d.csv",
        static_cast<int>(std::lround(p1)),
        static_cast<int>(std::lround(T1)),
        static_cast<int>(std::lround(us))
    );

    OFstream out(fname);
    out << "x,T,Tve,rho,u_sf,p,Y_N2,Y_N" << endl;

    double eveNow = eve1;
    double Tve = T1;
    std::vector<double> rhoi(ns), wdot(ns), Q(mix.nEnergyEqns());

    double x = 0.0;
    double nextWrite = 0.0;
    const double writeDx = xEnd/2000;

    while (x <= xEnd)
    {
        // sources at the current state
        for (int i = 0; i < ns; ++i)
        {
            rhoi[i] = rho*y[i];
        }
        const std::vector<double> TT = {T, Tve};
        mix.setState(rhoi.data(), TT.data(), 1);

        mix.energyTransferSource(Q.data());   // full Q_ve (VT + CV) [W/m^3]
        mix.netProductionRates(wdot.data());  // [kg/m^3/s]

        if (x >= nextWrite)
        {
            out << x << "," << T << "," << Tve << "," << rho << ","
                << u << "," << p << "," << y[fns.iN2] << "," << y[fns.iN]
                << endl;
            nextWrite += writeDx;
        }

        // explicit Euler step in x (dx small vs relaxation length)
        eveNow += Q[0]/m*dx;
        for (int i = 0; i < ns; ++i)
        {
            y[i] = max(y[i] + wdot[i]/m*dx, 0.0);
        }

        // renormalise Y against roundoff drift
        double ySum = 0.0;
        for (int i = 0; i < ns; ++i)
        {
            ySum += y[i];
        }
        for (int i = 0; i < ns; ++i)
        {
            y[i] /= ySum;
        }

        // recover the algebraic state and decode the temperatures
        if (!solveState(fns, y, m, P, H0, eveNow, rho, u, p, T))
        {
            Info<< "State recovery failed at x=" << x << endl;
            return 1;
        }

        // Tve from eve (Newton, same functional as the solver bridge)
        for (int it = 0; it < 100; ++it)
        {
            const double f = fns.eve(y, T, Tve) - eveNow;

            // dEve/dTv via speciesCpOverR
            mix.speciesCpOverR
            (
                T, Tve, T, Tve, Tve,
                nullptr, nullptr, nullptr, fns.w3.data(), fns.w4.data()
            );
            double cvve = 0.0;
            for (int i = 0; i < ns; ++i)
            {
                cvve +=
                    y[i]*(fns.w3[i] + fns.w4[i])
                   *Mutation::RU/mix.speciesMw(i);
            }
            cvve = max(cvve, 1e-8);

            double dTv = -f/cvve;
            dTv = min(max(dTv, -0.5*Tve), 0.5*Tve);
            Tve += dTv;
            if (mag(dTv) < 1e-8*Tve) break;
        }

        x += dx;
    }

    Info<< "Post-relaxation (x=" << xEnd << " m): T=" << T
        << " K, Tve=" << Tve << " K, Y_N=" << y[fns.iN] << nl
        << "Output saved to " << fname << endl;

    return 0;
}

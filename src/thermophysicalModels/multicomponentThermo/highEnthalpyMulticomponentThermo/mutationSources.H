// source terms of the two-temperature model of the paper (Casseau et al. 2016),
// computed with Mutation++: they are used both by the 0D programs (applications/test/
// nonEqTTv) and by the solver thermo, so the two paths do the same physics
//
// Mutation++ functions only, no OpenFOAM; the mixture state (rho_s,
// T, Tv) must be set beforehand with mix.setState(...)

#ifndef mutationSources_H
#define mutationSources_H

#include "mutation++.h"
#include "HarmonicOscillator.h"
#include "MillikanWhite.h"

#include <cmath>
#include <string>
#include <vector>

// a vibrating molecule: species index, Millikan-White data of
// Mutation++ (constants A, B for each partner and Park's sigma from the file
// transfer/VT.xml) and dissociation potential D (J/kg, table A1 of the paper)
struct Vibrator
{
    int species;
    Mutation::Transfer::MillikanWhiteModel tau;
    double D;
};

// dissociation potential D_m from table A1 of the paper (J/kg); it is needed by the
// preferential Q_C-V model (eq. 32). Same values as hy2Foam (thermoDEM)
inline double dissociationPotential(const std::string& name)
{
    if (name == "N2") return 3.36e7;
    if (name == "O2") return 1.54e7;
    if (name == "NO") return 2.09e7;
    return 0.0;
}

// list of the vibrators of the mixture: all the molecules, to be created only once
inline std::vector<Vibrator> makeVibrators(Mutation::Mixture& mix)
{
    std::vector<Vibrator> vibrators;
    Mutation::Thermodynamics::HarmonicOscillatorDB oscillatorDB;
    Mutation::Transfer::MillikanWhiteModelDB relaxationDB(mix);

    for (int s = 0; s < mix.nSpecies(); s++)
    {
        if (mix.species(s).type() == Mutation::Thermodynamics::MOLECULE)
        {
            const std::string name = mix.speciesName(s);
            // theta_v of the molecule, read from thermo/species.xml
            const double thetaV =
                oscillatorDB.create(name).characteristicTemperatures()[0];
            vibrators.push_back
            (
                {s, relaxationDB.create(name, thetaV), dissociationPotential(name)}
            );
        }
    }
    return vibrators;
}

// vibrational energy e_v (J/kg) of species s at temperature Tv (eq. 5)
inline double speciesEv(Mutation::Mixture& mix, const int s, const double Tv)
{
    std::vector<double> h_v(mix.nSpecies());
    // h_v/(R Tv) depends on Tv only: all temperatures equal to Tv
    mix.speciesHOverRT(Tv, Tv, Tv, Tv, Tv, nullptr, nullptr, nullptr,
                       h_v.data(), nullptr, nullptr);
    return h_v[s] * Mutation::RU * Tv / mix.speciesMw(s);
}

// electronic energy e_el (J/kg) of species s at temperature Tv (eq. 7)
inline double speciesEel(Mutation::Mixture& mix, const int s, const double Tv)
{
    std::vector<double> h_el(mix.nSpecies());
    mix.speciesHOverRT(Tv, Tv, Tv, Tv, Tv, nullptr, nullptr, nullptr,
                       nullptr, h_el.data(), nullptr);
    return h_el[s] * Mutation::RU * Tv / mix.speciesMw(s);
}

// vibro-electronic energy e_ve = e_v + e_el (J/kg) of species s at
// temperature Tv: this is the definition of the paper (text below eq. 8,
// e_ve,m = e_v,m + e_el,m), the same quantity as the Mutation++ reservoir
inline double speciesEve(Mutation::Mixture& mix, const int s, const double Tv)
{
    const int ns = mix.nSpecies();
    std::vector<double> h_v(ns), h_el(ns);
    mix.speciesHOverRT(Tv, Tv, Tv, Tv, Tv, nullptr, nullptr, nullptr,
                       h_v.data(), h_el.data(), nullptr);
    return (h_v[s] + h_el[s]) * Mutation::RU * Tv / mix.speciesMw(s);
}

// ---- V-T relaxation time of the molecule as in the paper (eq. 9-17) and
// in the default model of hy2Foam (MillikanWhitePark): for each partner s
//   tau(m,s) = tau_MW(m,s) + tau_P(m,s)                                (eq. 10)
//   tau_MW   = exp[A (T^-1/3 - B) - 18.42] / p,  p in atm              (eq. 11)
//   tau_P    = 1 / (c_m sigma_v n_ms),  c_m = sqrt(8 R_m T / pi),
//              sigma_v = sigma (50000/T)^2                        (eq. 15-17)
// with n_ms the number density of the pair, n_m + n_s (only n_m if s = m), as
// in hy2Foam; then a harmonic average weighted by the mole fractions (eq. 9).
// Constants A, B, sigma: those of Park 1993 in VT.xml.
// MillikanWhiteModel::relaxationTime of Mutation++ instead averages tau_MW
// arithmetically and adds a single Park correction with the density of the
// molecule alone: for pure N2 the two forms coincide, for mixtures they do not
inline double relaxationTimePaper(Mutation::Mixture& mix, const Vibrator& v)
{
    const Mutation::Transfer::MillikanWhiteModelData& data = v.tau.data();
    const int m = v.species;
    const double T = mix.T();
    const double p_atm = mix.P() / Mutation::ONEATM;
    const double n = mix.numberDensity();
    const double* X = mix.X();
    const int offset = (mix.hasElectrons() ? 1 : 0);

    const double c_m = std::sqrt(8.0 * Mutation::RU * T / (Mutation::PI * data.molecularWeight()));
    const double sigma_v = data.limitingCrossSection(T);
    const double T_fac = std::pow(T, -1.0/3.0);

    double sumX = 0.0;
    double sumXoverTau = 0.0;
    for (int i = 0; i < mix.nHeavy(); i++)
    {
        const int s = i + offset;
        if (X[s] <= 0.0)
        {
            continue;
        }
        const double tau_MW = std::exp(data.a()[i] * (T_fac - data.b()[i]) - 18.42) / p_atm;
        const double n_ms = (s == m) ? n * X[m] : n * (X[m] + X[s]);
        const double tau_P = 1.0 / (c_m * sigma_v * n_ms);
        sumX += X[s];
        sumXoverTau += X[s] / (tau_MW + tau_P);
    }
    return sumX / sumXoverTau;
}
// ---- end of relaxation time

// relaxation time of the molecule: that of Mutation++ (Millikan-White with
// Park's correction, report eq. 7; the default of the thermo and of the 0D
// programs) or that of the paper (eq. 9-17), which the verification cases
// select (report table 1)
inline double relaxationTime
(
    Mutation::Mixture& mix,
    const Vibrator& v,
    const bool paperTau
)
{
    return paperTau ? relaxationTimePaper(mix, v) : v.tau.relaxationTime(mix);
}

// ---- V-T exchange of all the molecules, Landau-Teller (eq. 8), in W/m3, at the
// current state of the mixture (T = mix.T(), Tv = mix.Tv())
// the driving term is e_ve(T) - e_ve(Tv) of the relaxing molecule, as
// written in the paper (eq. 8) and as in hy2Foam (LandauTellerVT); with
// vibrationalOnly it is e_v(T) - e_v(Tv), the vibrational energy only, as in the
// OmegaVT of Mutation++: the alternative of report table 1, run by Test-N2 on
// fig. 4 with electronic levels (without them the two forms coincide)
inline double sourceVT
(
    Mutation::Mixture& mix,
    const std::vector<double>& rho_s,
    const std::vector<Vibrator>& vibrators,
    const bool paperTau,
    const bool vibrationalOnly = false
)
{
    double Q = 0.0;
    for (const Vibrator& v : vibrators)
    {
        // molecule absent: nothing to relax (and tau would be infinite)
        if (rho_s[v.species] <= 0.0)
        {
            continue;
        }
        const double tau = relaxationTime(mix, v, paperTau);
        const double deltaE = vibrationalOnly
            ? speciesEv(mix, v.species, mix.T()) - speciesEv(mix, v.species, mix.Tv())
            : speciesEve(mix, v.species, mix.T()) - speciesEve(mix, v.species, mix.Tv());
        Q += rho_s[v.species] * deltaE / tau;
    }
    return Q;
}
// ---- end of V-T exchange

// ---- species production for the dissociation reactions (eq. 27-28),
// kg/m3/s, evaluated at Park's temperature T_P = T^a Tv^(1-a) (eq. 29, with
// a = 0.7 in the paper and in hy2Foam)
// Mutation++ evaluates them at sqrt(T Tv), which cannot be changed: it is given
// the state {T_P, T_P}, so that sqrt(T_P T_P) = T_P, and then {T, Tv} is restored
// valid only for mechanisms made of dissociations alone (here N2_Park)
inline void productionRates
(
    Mutation::Mixture& mix,
    const std::vector<double>& rho_s,
    const double parkExponent,
    std::vector<double>& wdot
)
{
    const double T = mix.T();
    const double Tv = mix.Tv();
    const double Tp = std::pow(T, parkExponent) * std::pow(Tv, 1.0 - parkExponent);

    const double TpTp[2] = {Tp, Tp};
    mix.setState(rho_s.data(), TpTp, 1);
    mix.netProductionRates(wdot.data());

    const double TTv[2] = {T, Tv};
    mix.setState(rho_s.data(), TTv, 1);
}
// ---- end of species production

// ---- vibro-electronic energy removed (or added) by the reactions, W/m3
//   - non-preferential (default, report eq. 8): every species created or
//     destroyed carries its vibro-electronic energy at Tv,
//     Q_C-V = sum over all the species s of e_ve,s(Tv) omega_s
//   - preferential (eq. 30 and 32 of the paper): the vibrational part of a
//     molecule is replaced by a fraction alpha of its dissociation energy,
//     Q_C-V = sum over the molecules m of omega_m (alpha D_m + e_el,m(Tv))
// The Park configuration of the paper (fig. 7 and 8) is the preferential one with
// alpha = 0.3 (report table 1): it is the setting of the hy2Foam heatBath
// tutorial (same case as fig. 7), and the verification cases select it
inline double sourceCV
(
    Mutation::Mixture& mix,
    const std::vector<double>& wdot,
    const std::vector<Vibrator>& vibrators,
    const bool preferential,
    const double alpha
)
{
    const double Tv = mix.Tv();
    double Q = 0.0;
    if (!preferential)
    {
        for (int s = 0; s < mix.nSpecies(); s++)
        {
            if (wdot[s] != 0.0)
            {
                Q += speciesEve(mix, s, Tv) * wdot[s];
            }
        }
        return Q;
    }
    for (const Vibrator& v : vibrators)
    {
        const int s = v.species;
        if (wdot[s] == 0.0)
        {
            continue;
        }
        Q += wdot[s] * (alpha * v.D + speciesEel(mix, s, Tv));
    }
    return Q;
}
// ---- end of chemistry-vibration coupling

#endif

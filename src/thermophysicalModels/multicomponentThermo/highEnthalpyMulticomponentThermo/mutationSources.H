// sorgenti del modello a due temperature del paper (Casseau et al. 2016),
// calcolate con Mutation++: le usano sia i programmi 0D (applications/test/
// nonEqTTv) sia il thermo del solver, cosi' i due percorsi fanno la stessa fisica
//
// solo funzioni di Mutation++, niente OpenFOAM; lo stato della miscela (rho_s,
// T, Tv) va impostato prima con mix.setState(...)

#ifndef mutationSources_H
#define mutationSources_H

#include "mutation++.h"
#include "HarmonicOscillator.h"
#include "MillikanWhite.h"

#include <cmath>
#include <string>
#include <vector>

// una molecola che vibra: indice della specie, dati di Millikan-White di
// Mutation++ (costanti A, B per ogni partner e sigma di Park dal file
// transfer/VT.xml) e potenziale di dissociazione D (J/kg, tabella A1 del paper)
struct Vibrator
{
    int species;
    Mutation::Transfer::MillikanWhiteModel tau;
    double D;
};

// potenziale di dissociazione D_m della tabella A1 del paper (J/kg); serve al
// modello preferenziale di Q_C-V (eq. 32). Stessi valori di hy2Foam (thermoDEM)
inline double dissociationPotential(const std::string& name)
{
    if (name == "N2") return 3.36e7;
    if (name == "O2") return 1.54e7;
    if (name == "NO") return 2.09e7;
    return 0.0;
}

// lista dei vibratori della miscela: tutte le molecole, da creare una volta sola
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
            // theta_v della molecola, letta da thermo/species.xml
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

// energia vibrazionale e_v (J/kg) della specie s alla temperatura Tv (eq. 5)
inline double speciesEv(Mutation::Mixture& mix, const int s, const double Tv)
{
    std::vector<double> h_v(mix.nSpecies());
    // h_v/(R Tv) dipende solo da Tv: tutte le temperature uguali a Tv
    mix.speciesHOverRT(Tv, Tv, Tv, Tv, Tv, nullptr, nullptr, nullptr,
                       h_v.data(), nullptr, nullptr);
    return h_v[s] * Mutation::RU * Tv / mix.speciesMw(s);
}

// energia elettronica e_el (J/kg) della specie s alla temperatura Tv (eq. 7)
inline double speciesEel(Mutation::Mixture& mix, const int s, const double Tv)
{
    std::vector<double> h_el(mix.nSpecies());
    mix.speciesHOverRT(Tv, Tv, Tv, Tv, Tv, nullptr, nullptr, nullptr,
                       nullptr, h_el.data(), nullptr);
    return h_el[s] * Mutation::RU * Tv / mix.speciesMw(s);
}

// energia vibro-elettronica e_ve = e_v + e_el (J/kg) della specie s alla
// temperatura Tv: e' la definizione del paper (testo sotto l'eq. 8,
// e_ve,m = e_v,m + e_el,m), la stessa grandezza del serbatoio di Mutation++
inline double speciesEve(Mutation::Mixture& mix, const int s, const double Tv)
{
    const int ns = mix.nSpecies();
    std::vector<double> h_v(ns), h_el(ns);
    mix.speciesHOverRT(Tv, Tv, Tv, Tv, Tv, nullptr, nullptr, nullptr,
                       h_v.data(), h_el.data(), nullptr);
    return (h_v[s] + h_el[s]) * Mutation::RU * Tv / mix.speciesMw(s);
}

// ---- tempo di rilassamento V-T della molecola come nel paper (eq. 9-17) e
// nel modello di default di hy2Foam (MillikanWhitePark): per ogni partner s
//   tau(m,s) = tau_MW(m,s) + tau_P(m,s)                                (eq. 10)
//   tau_MW   = exp[A (T^-1/3 - B) - 18.42] / p,  p in atm              (eq. 11)
//   tau_P    = 1 / (c_m sigma_v n_ms),  c_m = sqrt(8 R_m T / pi),
//              sigma_v = sigma (50000/T)^2                        (eq. 15-17)
// con n_ms densita' numerica della coppia, n_m + n_s (solo n_m se s = m), come
// in hy2Foam; poi media armonica pesata con le frazioni molari (eq. 9).
// Costanti A, B, sigma: quelle di Park 1993 nel VT.xml.
// MillikanWhiteModel::relaxationTime di Mutation++ invece media tau_MW in modo
// aritmetico e aggiunge una sola correzione di Park con la densita' della sola
// molecola: per N2 puro le due forme coincidono, per le miscele no
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
// ---- fine tempo di rilassamento

// tempo di rilassamento della molecola: quello del paper (default) oppure
// quello di Mutation++, tenuto per confronto
inline double relaxationTime
(
    Mutation::Mixture& mix,
    const Vibrator& v,
    const bool paperTau
)
{
    return paperTau ? relaxationTimePaper(mix, v) : v.tau.relaxationTime(mix);
}

// ---- scambio V-T di tutte le molecole, Landau-Teller (eq. 8), in W/m3, allo
// stato corrente della miscela (T = mix.T(), Tv = mix.Tv())
// la forza motrice e' e_ve(T) - e_ve(Tv) della molecola che rilassa, come
// scritto nel paper (eq. 8) e come in hy2Foam (LandauTellerVT); l'OmegaVT di
// Mutation++ userebbe solo quella vibrazionale (senza livelli elettronici
// le due forme coincidono)
inline double sourceVT
(
    Mutation::Mixture& mix,
    const std::vector<double>& rho_s,
    const std::vector<Vibrator>& vibrators,
    const bool paperTau
)
{
    double Q = 0.0;
    for (const Vibrator& v : vibrators)
    {
        // molecola assente: niente da rilassare (e tau sarebbe infinito)
        if (rho_s[v.species] <= 0.0)
        {
            continue;
        }
        const double tau = relaxationTime(mix, v, paperTau);
        Q += rho_s[v.species]
           * (speciesEve(mix, v.species, mix.T()) - speciesEve(mix, v.species, mix.Tv()))
           / tau;
    }
    return Q;
}
// ---- fine scambio V-T

// ---- produzione delle specie per le reazioni di dissociazione (eq. 27-28),
// kg/m3/s, valutate alla temperatura di Park T_P = T^a Tv^(1-a) (eq. 29, con
// a = 0.7 nel paper e in hy2Foam)
// Mutation++ le valuta a sqrt(T Tv), che non si puo' cambiare: gli si passa
// lo stato {T_P, T_P}, cosi' sqrt(T_P T_P) = T_P, e poi si ripristina {T, Tv}
// vale solo per meccanismi di sola dissociazione (qui N2_Park)
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
// ---- fine produzione delle specie

// ---- energia vibro-elettronica portata via (o aggiunta) dalle reazioni, W/m3
// eq. 30: Q_C-V = somma sulle molecole m di omega_m (D'_m + e_el,m(Tv)), con
//   - preferenziale (eq. 32):     D'_m = alpha D_m, alpha costante (0.3)
//   - non preferenziale (eq. 31): D'_m = e_v,m(Tv)
// solo le molecole, come nel paper e in hy2Foam (chemistry2Model, ParkTTv).
// Il caso Park del paper (fig. 7 e 8) e' il preferenziale con alpha = 0.3: e'
// l'impostazione del tutorial heatBath di hy2Foam (stesso caso della fig. 7)
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
    for (const Vibrator& v : vibrators)
    {
        const int s = v.species;
        if (wdot[s] == 0.0)
        {
            continue;
        }
        const double Dprime = preferential ? alpha * v.D : speciesEv(mix, s, Tv);
        Q += wdot[s] * (Dprime + speciesEel(mix, s, Tv));
    }
    return Q;
}
// ---- fine accoppiamento chimica-vibrazione

#endif

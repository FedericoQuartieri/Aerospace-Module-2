// sorgenti del modello a due temperature del paper (Casseau et al. 2016),
// calcolate con Mutation++: le usano sia i programmi 0D (applications/test/
// nonEqTTv) sia il thermo del solver, cosi' i due percorsi fanno la stessa fisica
//
// solo funzioni di Mutation++, niente OpenFOAM; lo stato della miscela (rho_s,
// T, Tv) va impostato prima con mix.setState(...)
//
// thread-safety: queste funzioni non hanno stato proprio (i buffer interni sono
// thread_local), ma Mutation::Mixture si', perche' setState() la muta. Chiamarle
// da piu' thread e' sicuro solo passando a ogni thread la sua mixture e i suoi
// vibratori: e' quello che fa il thermo con i MutationWorkspace.

#ifndef mutationSources_H
#define mutationSources_H

#include "mutation++.h"
#include "HarmonicOscillator.h"
#include "MillikanWhite.h"

#include <cmath>
#include <string>
#include <vector>

// una molecola che vibra: indice della specie e modello di Mutation++ per il
// suo tempo di rilassamento V-T (Millikan-White + correzione di Park, eq. 9-17,
// con le costanti del file transfer/VT.xml)
struct Vibrator
{
    int species;
    Mutation::Transfer::MillikanWhiteModel tau;
};

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
            vibrators.push_back({s, relaxationDB.create(name, thetaV)});
        }
    }
    return vibrators;
}

// energia vibrazionale e_v (J/kg) della specie s alla temperatura Tv (eq. 5)
inline double speciesEv(Mutation::Mixture& mix, const int s, const double Tv)
{
    // buffer per thread: questa funzione e' chiamata per ogni specie di ogni
    // cella, e allocarlo sullo heap a ogni chiamata metterebbe i thread OpenMP
    // in coda sull'allocatore. thread_local lo rende anche thread-safe
    static thread_local std::vector<double> h_v;
    h_v.resize(mix.nSpecies());
    // h_v/(R Tv) dipende solo da Tv: tutte le temperature uguali a Tv
    mix.speciesHOverRT(Tv, Tv, Tv, Tv, Tv, nullptr, nullptr, nullptr,
                       h_v.data(), nullptr, nullptr);
    return h_v[s] * Mutation::RU * Tv / mix.speciesMw(s);
}

// energia vibro-elettronica e_ve = e_v + e_el (J/kg) della specie s alla
// temperatura Tv (eq. 5 e 7): e' quella che Mutation++ mette nel serbatoio
// vibro-elettronico, per gli atomi e' la sola parte elettronica
inline double speciesEve(Mutation::Mixture& mix, const int s, const double Tv)
{
    const int ns = mix.nSpecies();
    // come in speciesEv: buffer per thread, niente allocazioni nel ciclo
    static thread_local std::vector<double> h_v, h_el;
    h_v.resize(ns);
    h_el.resize(ns);
    mix.speciesHOverRT(Tv, Tv, Tv, Tv, Tv, nullptr, nullptr, nullptr,
                       h_v.data(), h_el.data(), nullptr);
    return (h_v[s] + h_el[s]) * Mutation::RU * Tv / mix.speciesMw(s);
}

// ---- scambio V-T di tutte le molecole, Landau-Teller (eq. 8), in W/m3, allo
// stato corrente della miscela (T = mix.T(), Tv = mix.Tv())
// la forza motrice e' e_ve(T) - e_ve(Tv) con l'energia elettronica della
// molecola che rilassa, come nel paper; l'OmegaVT di Mutation++ userebbe solo
// quella vibrazionale (con i dati senza livelli elettronici i due coincidono)
inline double sourceVT
(
    Mutation::Mixture& mix,
    const std::vector<double>& rho_s,
    const std::vector<Vibrator>& vibrators
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
        const double tau = v.tau.relaxationTime(mix);
        Q += rho_s[v.species]
           * (speciesEve(mix, v.species, mix.T()) - speciesEve(mix, v.species, mix.Tv()))
           / tau;
    }
    return Q;
}
// ---- fine scambio V-T

// ---- produzione delle specie per le reazioni di dissociazione (eq. 27-28),
// kg/m3/s, valutate alla temperatura di Park T_P = T^a Tv^(1-a) (eq. 29)
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

// energia vibro-elettronica portata via (o aggiunta) dalle reazioni, W/m3:
// somma_s e_ve,s(Tv) omega_s, cioe' il modello non preferenziale (eq. 30-31)
// con D' = e_v, lo stesso che Mutation++ implementa (OmegaCV + OmegaCElec)
inline double sourceCV(Mutation::Mixture& mix, const std::vector<double>& wdot)
{
    double Q = 0.0;
    for (int s = 0; s < mix.nSpecies(); s++)
    {
        if (wdot[s] != 0.0)
        {
            Q += speciesEve(mix, s, mix.Tv()) * wdot[s];
        }
    }
    return Q;
}

#endif

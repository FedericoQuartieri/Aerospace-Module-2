// funzioni comuni ai test 0D (heat bath) che usano Mutation++
//
// lo stato del gas e' descritto dalle variabili conservate del paper (eq. 23):
// densita' parziali rho_s, energia totale E ed energia vibro-elettronica E_ve,
// tutte per unita' di volume; le temperature si ricavano da queste con
// mix.setState(rho_s, {E, E_ve}, 0)

#ifndef heatBath_H
#define heatBath_H

#include "mutation++.h"

#include <algorithm>
#include <cmath>
#include <vector>

// ---- istanti in cui scrivere il csv: ogni passo all'inizio, poi a passo
// logaritmico (200 punti per decade) fino all'ultimo passo. Le figure del
// paper hanno il tempo in scala logaritmica da 1 ns: con una scrittura ogni
// 10 o 100 passi i primi nanosecondi venivano interpolati fra punti lontani
// (fra t = 0 e t = 10 ns) e gli scarti in quel tratto erano solo interpolazione
struct OutputSchedule
{
    double next = 1.0;

    // vero se al passo step (1..nSteps) va scritta una riga
    bool write(const int step, const int nSteps)
    {
        if (step < next && step != nSteps)
        {
            return false;
        }
        next = std::max(next + 1.0, next * std::pow(10.0, 1.0 / 200.0));
        return true;
    }
};
// ---- fine istanti di scrittura

// energia totale per unita' di volume, E = somma_s rho_s e_s (eq. 1 senza cinetica)
inline double totalEnergy(Mutation::Mixture& mix, const std::vector<double>& rho_s)
{
    const int ns = mix.nSpecies();
    // e_s[s] energia totale, e_s[ns + s] energia vibro-elettronica (J/kg)
    std::vector<double> e_s(2 * ns);
    mix.getEnergiesMass(e_s.data());

    double E = 0.0;
    for (int s = 0; s < ns; s++)
    {
        E += rho_s[s] * e_s[s];
    }
    return E;
}

// energia vibro-elettronica per unita' di volume, E_ve = somma_s rho_s e_ve,s (eq. 2)
inline double veEnergy(Mutation::Mixture& mix, const std::vector<double>& rho_s)
{
    const int ns = mix.nSpecies();
    std::vector<double> e_s(2 * ns);
    mix.getEnergiesMass(e_s.data());

    double Eve = 0.0;
    for (int s = 0; s < ns; s++)
    {
        Eve += rho_s[s] * e_s[ns + s];
    }
    return Eve;
}

// ---- temperatura di equilibrio dalla conservazione dell'energia
// cerca per bisezione la T tale che, con T_tr = T_ve = T e composizione fissa,
// l'energia totale sia uguale a quella iniziale E
inline double equilibriumTemperature
(
    Mutation::Mixture& mix,
    const std::vector<double>& rho_s,
    double E
)
{
    double Tlow = 100.0;
    double Thigh = 100000.0;

    for (int iter = 0; iter < 100; iter++)
    {
        const double T = 0.5 * (Tlow + Thigh);
        const double temps[2] = {T, T};
        mix.setState(rho_s.data(), temps, 1);

        // l'energia cresce con T: se e' troppa, T e' troppo alta
        if (totalEnergy(mix, rho_s) > E)
        {
            Thigh = T;
        }
        else
        {
            Tlow = T;
        }
    }
    return 0.5 * (Tlow + Thigh);
}
// ---- fine temperatura di equilibrio

#endif

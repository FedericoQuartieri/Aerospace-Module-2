// functions shared by the 0D tests (heat bath) that use Mutation++
//
// the gas state is described by the conserved variables of the paper (eq. 23):
// partial densities rho_s, total energy E and vibro-electronic energy E_ve,
// all per unit volume; the temperatures are obtained from these with
// mix.setState(rho_s, {E, E_ve}, 0)

#ifndef heatBath_H
#define heatBath_H

#include "mutation++.h"

#include <algorithm>
#include <cmath>
#include <vector>

// ---- instants at which the csv is written: every step at the beginning, then
// at logarithmic intervals (200 points per decade) up to the last step. The
// figures of the paper have time on a logarithmic scale from 1 ns: writing every
// 10 or 100 steps, the first nanoseconds were interpolated between distant points
// (between t = 0 and t = 10 ns) and the deviations there were only interpolation
struct OutputSchedule
{
    double next = 1.0;

    // true if a row must be written at step step (1..nSteps)
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
// ---- end of write instants

// total energy per unit volume, E = sum_s rho_s e_s (eq. 1 without the kinetic term)
inline double totalEnergy(Mutation::Mixture& mix, const std::vector<double>& rho_s)
{
    const int ns = mix.nSpecies();
    // e_s[s] total energy, e_s[ns + s] vibro-electronic energy (J/kg)
    std::vector<double> e_s(2 * ns);
    mix.getEnergiesMass(e_s.data());

    double E = 0.0;
    for (int s = 0; s < ns; s++)
    {
        E += rho_s[s] * e_s[s];
    }
    return E;
}

// vibro-electronic energy per unit volume, E_ve = sum_s rho_s e_ve,s (eq. 2)
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

// ---- equilibrium temperature from energy conservation
// bisection search for the T such that, with T_tr = T_ve = T and fixed
// composition, the total energy equals the initial energy E
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

        // the energy grows with T: if it is too large, T is too high
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
// ---- end of equilibrium temperature

#endif

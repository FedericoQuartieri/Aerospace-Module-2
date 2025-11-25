/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2012-2023 OpenFOAM Foundation
     \\/     M anipulation  |
-------------------------------------------------------------------------------
License
    This file is part of OpenFOAM.

    OpenFOAM is free software: you can redistribute it and/or modify it
    under the terms of the GNU General Public License as published by
    the Free Software Foundation, either version 3 of the License, or
    (at your option) any later version.

    OpenFOAM is distributed in the hope that it will be useful, but WITHOUT
    ANY WARRANTY; without even the implied warranty of MERCHANTABILITY or
    FITNESS FOR A PARTICULAR PURPOSE.  See the GNU General Public License
    for more details.

    You should have received a copy of the GNU General Public License
    along with OpenFOAM.  If not, see <http://www.gnu.org/licenses/>.

Application
    Test-thermoMixture

Description

\*---------------------------------------------------------------------------*/
#include <iostream>
#include "mutation++.h"

using namespace Mutation;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
// Main program:

int main(int argc, char *argv[])
{
    std::cout << std::fixed;
    std::cout << std::setprecision(4);

    Thermodynamics::Thermodynamics thermo("H2 O2 C4",
                                          "RRHO",
                                          "ChemNonEqTTv");

    int n_species = thermo.nSpecies();

    double T = 5000.0;
    std::vector<double> cp(n_species);
    thermo.speciesCpOverR(T, cp.data());

    std::cout << "Per-specie Cp(T=" << T << "):" << std::endl;

    for (int i = 0; i < n_species; ++i) {
        std::cout << "  " << thermo.speciesName(i) << ": "
                  << cp[i] << std::endl;
    }

    // Non-equilibrium temp Cps here

    double Tt = 5200.0; // Translational
    double Tr = 5000.0; // Rotational
    double Tv = 4000.0; // Vibrational
    double Te = 5500.0; // Electronic

    std::vector<double> cp_t(n_species);
    std::vector<double> cp_r(n_species);
    std::vector<double> cp_v(n_species);
    std::vector<double> cp_e(n_species);

    // Note that Cp = Cp_t + Cp_r + Cp_v + Cp_e

    // We set T_electron_translational = 0
    thermo.speciesCpOverR(Tt, 0, Tr, Tv, Te,
                          cp.data(), cp_t.data(), cp_r.data(),
                          cp_v.data(), cp_e.data());

    std::cout << std::endl << "Per-specie Cp(Tt=" << Tt
              << ", Tr=" << Tr << ", Tv= " << Tv << ", Te="
              << Te << "):" << std::endl;

    for (int i = 0; i < thermo.nSpecies(); ++i) {
        std::cout << "  " << thermo.speciesName(i) << ": "
                  << cp[i] << " = " << cp_t[i] << " + " << cp_r[i]
                  << " + " << cp_v[i] << " + " << cp_e[i] << std::endl;
    }

    return 0;
}


// ************************************************************************* //

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
    Mixture mix("air_5");
    std::cout << "# of elements: " << mix.nElements() << '\n';
    std::cout << "# of species: " << mix.nSpecies() << ' ';
    std::cout << mix.nGas() << " (gas) " << mix.nCondensed() << " (condensed)\n";
    std::cout << "# of reactions: " << mix.nReactions() << '\n';
    std::cout << "# of temperatures: " << mix.nEnergyEqns() << '\n';

    return 0;
}


// ************************************************************************* //

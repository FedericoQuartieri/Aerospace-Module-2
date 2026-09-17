/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2023 OpenFOAM Foundation
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

\*---------------------------------------------------------------------------*/

#include "highEnthalpyMulticomponentThermo.H"


#include "coefficientMulticomponentMixture.H"
#include "coefficientWilkeMulticomponentMixture.H"
#include "singleComponentMixture.H"

#include "forGases.H"
#include "rrhoThermo.H"

#include "makeFluidMulticomponentThermo.H"


// ---- macro per il thermo rrho del template del corso
// rrho e' una copia di janaf, ma non e' fra i thermo delle macro forGases e
// forCoeffGases del core: lo si istanzia qui con le stesse combinazioni
// (gas perfetto, trasporto const o sutherland, energia h o e)
#define forRrhoGasEqns(Mu, He, Macro, Args...)                                 \
    forThermo(Mu, He, rrhoThermo, perfectGas, specie, Macro, Args)

#define forRrhoGasEnergies(Mu, Macro, Args...)                                 \
    forRrhoGasEqns(Mu, sensibleEnthalpy, Macro, Args);                         \
    forRrhoGasEqns(Mu, sensibleInternalEnergy, Macro, Args)

#define forRrhoGases(Macro, Args...)                                           \
    forRrhoGasEnergies(constTransport, Macro, Args);                           \
    forRrhoGasEnergies(sutherlandTransport, Macro, Args)
// ---- fine macro rrho


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
    // thermo del core: hConst, eConst, janaf
    forCoeffGases
    (
        makeFluidMulticomponentThermos,
        psiThermo,
        highEnthalpyMulticomponentThermo,
        coefficientMulticomponentMixture
    );
    forCoeffGases
    (
        makeFluidMulticomponentThermos,
        psiThermo,
        highEnthalpyMulticomponentThermo,
        coefficientWilkeMulticomponentMixture
    );
    forGases
    (
        makeFluidMulticomponentThermo,
        highEnthalpyMulticomponentThermo,
        singleComponentMixture
    );

    // thermo rrho del template, in aggiunta a quelli del core
    forRrhoGases
    (
        makeFluidMulticomponentThermos,
        psiThermo,
        highEnthalpyMulticomponentThermo,
        coefficientMulticomponentMixture
    );
    forRrhoGases
    (
        makeFluidMulticomponentThermos,
        psiThermo,
        highEnthalpyMulticomponentThermo,
        coefficientWilkeMulticomponentMixture
    );
    forRrhoGases
    (
        makeFluidMulticomponentThermo,
        highEnthalpyMulticomponentThermo,
        singleComponentMixture
    );
}

// ************************************************************************* //

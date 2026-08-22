/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2023 OpenFOAM Foundation
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

#include "shockThermo.H"
#include "fvmDdt.H"
#include "fvmSup.H"
#include "fvcDiv.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::shockThermo::thermophysicalPredictor()
{
    // The two-temperature model is active when the highEnthalpyThermo
    // registered its fields: the conserved vibro-electronic energy "eve" and
    // the Landau-Teller linearisation ("eveEq", "tauVT") computed from
    // Mutation++ in thermo.correct().
    // add support to multi-specie chemistry

    tmp<fv::convectionScheme<scalar>> mvConvection
    (
        fv::convectionScheme<scalar>::New
        (
            mesh,
            fields,
            phi,
            mesh.schemes().div("div(phi,Yi_h)")
        )
    );

    reaction->correct();

    forAll(Y, i)
    {
        volScalarField& Yi = Y_[i];

        if (thermo_.solveSpecie(i))
        {
            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              + thermophysicalTransport->divj(Yi)
             ==
                reaction->R(Yi)
              + fvModels().source(rho, Yi)
            );

            YiEqn.relax();

            fvConstraints().constrain(YiEqn);

            YiEqn.solve("Yi");

            fvConstraints().constrain(Yi);
        }
        else
        {
            Yi.correctBoundaryConditions();
        }
    }

    thermo_.normaliseY();

    if (thermo_.he().name() != "e")
    {
        FatalErrorInFunction()
            << "sensible energy e is required as primary variable"
            << exit(FatalError);
    }

    //- ------------------------------------------------------------------------

    // solve the equation of energy.

    // IMPORTANT: I cannot use the function call
    // shockFluid::thermophysicalPredictor(), as it causes inconsistencies
    // between thermo classes. The thermo.correct() function must be the one
    // defined in the derived class; shockFluid::thermophysicalPredictor() is
    // pasted here below.

    volScalarField& e = thermo_.he();

    const surfaceScalarField e_pos(interpolate(e, pos, thermo.T().name()));
    const surfaceScalarField e_neg(interpolate(e, neg, thermo.T().name()));

    surfaceScalarField phiEp
    (
        "phiEp",
        aphiv_pos()*(rho_pos()*(e_pos + 0.5*magSqr(U_pos())) + p_pos())
      + aphiv_neg()*(rho_neg()*(e_neg + 0.5*magSqr(U_neg())) + p_neg())
      + aSf()*(p_pos() - p_neg())
    );

    // Make flux for pressure-work absolute
    if (mesh.moving())
    {
        phiEp += mesh.phi()*(a_pos()*p_pos() + a_neg()*p_neg());
    }

    // Solving for the total energy e,
    // for high enthalpy flows e = e_tr + e_ve
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, e) + fvc::div(phiEp)
      + fvc::ddt(rho, K)
     ==
        fvModels().source(rho, e)
    );

    if (!inviscid)
    {
        const surfaceScalarField devTauDotU
        (
            "devTauDotU",
            devTau() & (a_pos()*U_pos() + a_neg()*U_neg())
        );

        EEqn += thermophysicalTransport->divq(e) + fvc::div(devTauDotU);
    }

    EEqn.relax();

    fvConstraints().constrain(EEqn);

    EEqn.solve();

    fvConstraints().constrain(e);

    volScalarField& eve = thermo_.eve();

    tmp<volScalarField> Q_VT = thermo_.computeSourceVT(thermo_.defaultSpecie());

    // Solve for vibrational energy e_ve
    fvScalarMatrix EveEqn
    (
        fvm::ddt(rho, eve)
     ==
        Q_VT()
    );

    EveEqn.relax();

    fvConstraints().constrain(EveEqn);

    EveEqn.solve("eve");

    fvConstraints().constrain(eve);

    // Update T_ (T_tr) and Tve_ based on solved energies
    thermo_.correct();
}


// ************************************************************************* //

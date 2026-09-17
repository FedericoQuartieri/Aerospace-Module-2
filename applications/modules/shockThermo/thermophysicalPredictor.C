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

// chiamata da foamRun a ogni passo: dopo flussi e velocita', prima della pressione
void Foam::solvers::shockThermo::thermophysicalPredictor()
{
    // ---- equazioni delle specie (macchinario OpenFOAM), sorgente chimica da Mutation++
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

    forAll(Y, i)
    {
        volScalarField& Yi = Y_[i];

        if (thermo_.solveSpecie(i))
        {
            // produzione chimica della specie calcolata da Mutation++ (eq. 27)
            tmp<volScalarField> wdot = thermo_.computeSourceY(i);

            fvScalarMatrix YiEqn
            (
                fvm::ddt(rho, Yi)
              + mvConvection->fvmDiv(phi, Yi)
              + thermophysicalTransport->divj(Yi)
             ==
                wdot()
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
    // ---- fine specie

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

    // ---- energia sensibile e (copiata da shockFluid), con il calore di reazione (eq. 22)
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

    // calore di reazione: e e' l'energia sensibile, la chimica la cambia
    tmp<volScalarField> Q_chem = thermo_.computeSourceE();

    // Solving for the sensible energy e,
    // for high enthalpy flows e = e_tr + e_ve
    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, e) + fvc::div(phiEp)
      + fvc::ddt(rho, K)
     ==
        Q_chem()
      + fvModels().source(rho, e)
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
    // ---- fine energia totale

    volScalarField& eve = thermo_.eve();

    // sorgente calcolata da Mutation++ nel bridge: V-T piu' chimica (eq. 26)
    tmp<volScalarField> Q_ve = thermo_.computeSourceVe();

    // ---- equazione di eve: d(rho*eve)/dt = Q_ve (eq. 22, forma 0D: senza trasporto)
    // Solve for vibrational energy e_ve
    fvScalarMatrix EveEqn
    (
        fvm::ddt(rho, eve)
     ==
        Q_ve()
    );

    EveEqn.relax();

    fvConstraints().constrain(EveEqn);

    EveEqn.solve("eve");

    fvConstraints().constrain(eve);
    // ---- fine eve

    // decode nel bridge: dalle energie ricava T_tr e T_ve
    // Update T_ (T_tr) and Tve_ based on solved energies
    thermo_.correct();
}


// ************************************************************************* //

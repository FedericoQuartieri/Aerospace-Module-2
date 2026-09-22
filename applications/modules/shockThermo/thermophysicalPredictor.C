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

// called by foamRun at every step: after the fluxes and the velocity, before the pressure
void Foam::solvers::shockThermo::thermophysicalPredictor()
{
    // state of each cell at the start of the corrector (rho_s, T, Tve): the
    // Mutation++ sources below (chemistry, heat of reaction, V-T and Q_C-V) are
    // each evaluated at this state, the library state being reset for every
    // species and source term
    thermo_.correctSourceState();

    // ---- species equations (OpenFOAM machinery), chemical source from Mutation++
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
            // chemical production of the species computed by Mutation++ (eq. 27)
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
    // ---- end of species

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

    // ---- sensible energy e (copied from shockFluid), with the heat of reaction (eq. 22)
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

    // heat of reaction: e is the sensible energy, and chemistry changes it
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
    // ---- end of total energy

    volScalarField& eve = thermo_.eve();

    // V-T plus chemistry (eq. 26), at the state stored by correctSourceState()
    tmp<volScalarField> Q_ve = thermo_.computeSourceVe();

    // ---- eve equation: d(rho*eve)/dt = Q_ve (eq. 22, 0D form: without transport)
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
    // ---- end of eve

    // decode in the bridge: obtains T_tr and T_ve from the energies
    // Update T_ (T_tr) and Tve_ based on solved energies
    thermo_.correct();
}


// ************************************************************************* //

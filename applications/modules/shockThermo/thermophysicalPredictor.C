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
    // registered its fields: the conserved vibro-electronic energy "eve",
    // the Landau-Teller linearisation ("eveEq", "tauVT") and the chemistry
    // sources ("mutQdot", "mutQcv", "mutR_<specie>") computed from
    // Mutation++ in thermo.correct().
    //
    // CONSERVATIVE FORMULATION: the solved sensible energy "e" lives on the
    // Mutation++ datum e = e_tr(T) + e_ve(Tve). V-T exchange redistributes
    // energy between the two pools inside "e", so no relaxation term appears
    // in the total energy equation: the effect on T enters through the
    // decode T(e - eve) done by the thermo. Chemistry sources:
    //   e equation:    Qdot = -sum_i hf_i*wdot_i  (sensible datum)
    //   eve equation:  Qcv  =  sum_i e_ve,i*wdot_i (Candler, non-pref.)
    //   Yi equations:  mutR_<specie> = wdot_i from Mutation++ kinetics
    //                  (evaluated at the Park controlling temperature)
    const bool hasEve =
        mesh.foundObject<volScalarField>("eve")
     && mesh.foundObject<volScalarField>("eveEq")
     && mesh.foundObject<volScalarField>("tauVT");

    const bool hasMutChemistry =
        mesh.foundObject<volScalarField::Internal>("mutQdot")
     && mesh.foundObject<volScalarField::Internal>("mutQcv");

    bool solveEve = false;

    const dictionary& thermoProperties = thermo_.properties();
    if (hasEve && thermoProperties.found("highEnthalpyRelaxation"))
    {
        const dictionary& relaxDict =
            thermoProperties.subDict("highEnthalpyRelaxation");

        solveEve = relaxDict.lookupOrDefault<Switch>
        (
            "solveEve",
            relaxDict.lookupOrDefault<Switch>("solveTve", false)
        );
    }

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
                fvModels().source(rho, Yi)
            );

            // Mass production from Mutation++ two-temperature kinetics when
            // available; otherwise the OpenFOAM combustion model
            const word mutRName("mutR_" + Yi.name());
            if (mesh.foundObject<volScalarField::Internal>(mutRName))
            {
                YiEqn -=
                    mesh.lookupObject<volScalarField::Internal>(mutRName);
            }
            else
            {
                YiEqn -= reaction->R(Yi);
            }

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

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, e) + fvc::div(phiEp)
      + fvc::ddt(rho, K)
     ==
        fvModels().source(rho, e)
    );

    // Chemistry heat release on the sensible-energy datum
    if (hasMutChemistry)
    {
        EEqn -= mesh.lookupObject<volScalarField::Internal>("mutQdot");
    }

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

    // Decodes (T, Tve) from the solved (e, eve) and refreshes the
    // relaxation/chemistry source fields through the highEnthalpyThermo
    thermo_.correct();

    if (hasEve && solveEve)
    {
        volScalarField& eve = mesh.lookupObjectRef<volScalarField>("eve");
        const volScalarField& eveEq =
            mesh.lookupObject<volScalarField>("eveEq");
        const volScalarField& tauVT =
            mesh.lookupObject<volScalarField>("tauVT");

        // Conservative vibro-electronic energy equation: semi-implicit
        // Landau-Teller relaxation towards eveEq = e_ve(Ttr) plus the
        // chemistry-vibration coupling
        fvScalarMatrix EveEqn
        (
            fvm::ddt(rho, eve)
          + mvConvection->fvmDiv(phi, eve)
         ==
            rho*eveEq/tauVT
          - fvm::Sp(rho/tauVT, eve)
          + fvModels().source(rho, eve)
        );

        if (hasMutChemistry)
        {
            EveEqn -= mesh.lookupObject<volScalarField::Internal>("mutQcv");
        }

        EveEqn.relax();

        fvConstraints().constrain(EveEqn);

        EveEqn.solve("eve");

        fvConstraints().constrain(eve);

        eve.max(dimensionedScalar(eve.dimensions(), Zero));
        eve.correctBoundaryConditions();

        // Derive Tve from the updated eve and refresh the properties
        thermo_.correct();
    }
}


// ************************************************************************* //

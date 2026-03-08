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
#include "fvcDiv.H"
#include "fvcDdt.H"

// * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * * //

void Foam::solvers::shockThermo::thermophysicalPredictor()
{
    const bool hasTve = mesh.foundObject<volScalarField>("Tve");

    bool solveTveRelax = false;
    bool coupleEnergyRelax = false;
    scalar tauVT = scalar(1e-4);
    scalar minTve = scalar(1);
    scalar maxTve = scalar(GREAT);
    scalar cvVeOverCvTr = scalar(1);

    const dictionary& thermoProperties = thermo_.properties();
    if (hasTve && thermoProperties.found("highEnthalpyRelaxation"))
    {
        const dictionary& relaxDict =
            thermoProperties.subDict("highEnthalpyRelaxation");

        solveTveRelax = relaxDict.lookupOrDefault<Switch>("solveTve", false);
        coupleEnergyRelax =
            relaxDict.lookupOrDefault<Switch>("coupleEnergy", false);
        tauVT = relaxDict.lookupOrDefault<scalar>("tauVT", tauVT);
        minTve = relaxDict.lookupOrDefault<scalar>("minTemperature", minTve);
        maxTve = relaxDict.lookupOrDefault<scalar>("maxTemperature", maxTve);
        cvVeOverCvTr =
            relaxDict.lookupOrDefault<scalar>("cvVeOverCvTr", cvVeOverCvTr);
    }

    const dimensionedScalar tauVTDim
    (
        "tauVT",
        dimTime,
        max(tauVT, scalar(SMALL))
    );

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


    //- for high enthalpy flows, e = e_rt + e_ve.

    volScalarField eRelaxSource
    (
        IOobject
        (
            "eRelaxSource",
            mesh.time().name(),
            mesh,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh,
        dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
    );

    if (hasTve && solveTveRelax && coupleEnergyRelax)
    {
        const volScalarField& Tve = mesh.lookupObject<volScalarField>("Tve");

        // Positive source here is removed from e-equation RHS below.
        eRelaxSource =
            rho*max(cvVeOverCvTr, scalar(0))*thermo_.Cv()*(thermo_.T() - Tve)
           /tauVTDim;
    }

    fvScalarMatrix EEqn
    (
        fvm::ddt(rho, e) + fvc::div(phiEp)
      + fvc::ddt(rho, K)
     ==
        fvModels().source(rho, e) - eRelaxSource
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

    thermo_.correct();

    if (hasTve)
    {
        if (solveTveRelax)
        {
            volScalarField& Tve = mesh.lookupObjectRef<volScalarField>("Tve");

            fvScalarMatrix TveEqn
            (
                fvm::ddt(rho, Tve)
              + mvConvection->fvmDiv(phi, Tve)
             ==
                rho*(thermo_.T() - Tve)/tauVTDim
              + fvModels().source(rho, Tve)
            );

            TveEqn.relax();

            fvConstraints().constrain(TveEqn);

            TveEqn.solve("Tve");

            fvConstraints().constrain(Tve);

            Tve.max(minTve);
            Tve.min(maxTve);
            Tve.correctBoundaryConditions();

            thermo_.correct();
        }
    }

}


// ************************************************************************* //

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
    ThermoMixture

Description
    Example of coupling between OpenFOAM-13 and mutationPP

    Federico Piscaglia
    Dept. of Aerospace Science and Technology
    Politecnico di Milano

\*---------------------------------------------------------------------------*/

#include "mutation++.h"
#include <Eigen/Dense>

#include "fv.H"
#include "dictionary.H"
#include "IFstream.H"


using namespace Foam;

using namespace Mutation;

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //


void test1()
{
    // Generate the default options for the air11 mixture
    MixtureOptions opts("air_11");
    opts.setStateModel("EquilTP");

    // Rigid-rotor harmonic oscillator (RRHO) two-temperature model
    opts.setThermodynamicDatabase("RRHO");

    // Load the mixture with the new options
    Mixture mix(opts);
    int ns = mix.nSpecies();

    Eigen::MatrixXd m_Dij(ns, ns);
    Eigen::MatrixXd m_ram_Dij(ns, ns);
    Eigen::VectorXd v_Vd_ram(ns);
    Eigen::VectorXd v_Vd_sm(ns);
    Eigen::VectorXd v_Vd(ns);
    Eigen::VectorXd v_b(ns);
    Eigen::VectorXd v_rhoi(ns);

}


//- ------------------------------------------------------------------------
// Main program:

int main(int argc, char *argv[])
{
    test1();

    return 0;
}

// ************************************************************************* //

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

Class
    Foam::blottnerTransport

Description
    Transport package using the Blottner curve fit for viscosity:

    \verbatim
        mu = 0.1*exp((A*ln(T) + B)*ln(T) + C)   [kg/m/s]
    \endverbatim

    (Blottner, Johnson & Ellis 1971; coefficients tabulated for air
    species in Gupta et al. 1990). Valid over the whole 100-30000 K
    range of hypersonic boundary layers, unlike a two-parameter
    Sutherland fit which cannot track the curve from cold-wall to
    post-shock temperatures.

    The thermal conductivity (for gases) is obtained from the modified
    Eucken correction, as in sutherlandTransport:

    \verbatim
        kappa = mu*Cv*(1.32 + 1.77*R/Cv)
    \endverbatim

Usage
    \table
        Property     | Description
        A            | Blottner coefficient of ln(T)^2 [-]
        B            | Blottner coefficient of ln(T) [-]
        C            | Blottner constant coefficient [-]
    \endtable

    Example specification of blottnerTransport for N2:
    \verbatim
    transport
    {
        A               0.0268142;
        B               0.3177838;
        C               -11.3155513;
    }
    \endverbatim

SourceFiles
    blottnerTransportI.H
    blottnerTransport.C

\*---------------------------------------------------------------------------*/

#ifndef blottnerTransport_H
#define blottnerTransport_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{

// Forward declaration of friend functions and operators

template<class Thermo> class blottnerTransport;

template<class Thermo>
inline blottnerTransport<Thermo> operator+
(
    const blottnerTransport<Thermo>&,
    const blottnerTransport<Thermo>&
);

template<class Thermo>
inline blottnerTransport<Thermo> operator*
(
    const scalar,
    const blottnerTransport<Thermo>&
);

template<class Thermo>
Ostream& operator<<
(
    Ostream&,
    const blottnerTransport<Thermo>&
);


/*---------------------------------------------------------------------------*\
                      Class blottnerTransport Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class blottnerTransport
:
    public Thermo
{
    // Private Data

        // Blottner's coefficients
        scalar A_, B_, C_;


    // Private Member Functions

        //- Read coefficient from dictionary
        scalar readCoeff(const word& coeffName, const dictionary& dict);


public:

    // Constructors

        //- Construct from components
        inline blottnerTransport
        (
            const Thermo& t,
            const scalar A,
            const scalar B,
            const scalar C
        );

        //- Construct as named copy
        inline blottnerTransport(const word&, const blottnerTransport&);

        //- Construct from name and dictionary
        blottnerTransport(const word& name, const dictionary& dict);

        //- Construct from base thermo and dictionary
        blottnerTransport(const Thermo& t, const dictionary& dict);

        //- Construct and return a clone
        inline autoPtr<blottnerTransport> clone() const;


    // Member Functions

        //- Return the instantiated type name
        static word typeName()
        {
            return "blottner<" + Thermo::typeName() + '>';
        }

        //- Dynamic viscosity [kg/m/s]
        inline scalar mu(const scalar p, const scalar T) const;

        //- Thermal conductivity [W/m/K]
        inline scalar kappa(const scalar p, const scalar T) const;

        //- Write to Ostream
        void write(Ostream& os) const;


    // Member Operators

        inline void operator+=(const blottnerTransport&);
        inline void operator*=(const scalar);


    // Friend operators

        friend blottnerTransport operator+ <Thermo>
        (
            const blottnerTransport&,
            const blottnerTransport&
        );

        friend blottnerTransport operator* <Thermo>
        (
            const scalar,
            const blottnerTransport&
        );


    // Ostream Operator

        friend Ostream& operator<< <Thermo>
        (
            Ostream&,
            const blottnerTransport&
        );
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "blottnerTransportI.H"

#ifdef NoRepository
    #include "blottnerTransport.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

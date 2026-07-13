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

#include "specie.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::blottnerTransport<Thermo>::blottnerTransport
(
    const Thermo& t,
    const scalar A,
    const scalar B,
    const scalar C
)
:
    Thermo(t),
    A_(A),
    B_(B),
    C_(C)
{}


template<class Thermo>
inline Foam::blottnerTransport<Thermo>::blottnerTransport
(
    const word& name,
    const blottnerTransport& bt
)
:
    Thermo(name, bt),
    A_(bt.A_),
    B_(bt.B_),
    C_(bt.C_)
{}


template<class Thermo>
inline Foam::autoPtr<Foam::blottnerTransport<Thermo>>
Foam::blottnerTransport<Thermo>::clone() const
{
    return autoPtr<blottnerTransport<Thermo>>
    (
        new blottnerTransport<Thermo>(*this)
    );
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::scalar Foam::blottnerTransport<Thermo>::mu
(
    const scalar p,
    const scalar T
) const
{
    const scalar lnT = log(T);
    return 0.1*exp((A_*lnT + B_)*lnT + C_);
}


template<class Thermo>
inline Foam::scalar Foam::blottnerTransport<Thermo>::kappa
(
    const scalar p,
    const scalar T
) const
{
    scalar Cv_ = this->Cv(p, T);
    return mu(p, T)*Cv_*(1.32 + 1.77*this->R()/Cv_);
}


// * * * * * * * * * * * * * * * Member Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline void Foam::blottnerTransport<Thermo>::operator+=
(
    const blottnerTransport<Thermo>& bt
)
{
    scalar Y1 = this->Y();

    Thermo::operator+=(bt);

    if (mag(this->Y()) > small)
    {
        Y1 /= this->Y();
        scalar Y2 = bt.Y()/this->Y();

        A_ = Y1*A_ + Y2*bt.A_;
        B_ = Y1*B_ + Y2*bt.B_;
        C_ = Y1*C_ + Y2*bt.C_;
    }
}


template<class Thermo>
inline void Foam::blottnerTransport<Thermo>::operator*=
(
    const scalar s
)
{
    Thermo::operator*=(s);
}


// * * * * * * * * * * * * * * * Friend Operators  * * * * * * * * * * * * * //

template<class Thermo>
inline Foam::blottnerTransport<Thermo> Foam::operator+
(
    const blottnerTransport<Thermo>& bt1,
    const blottnerTransport<Thermo>& bt2
)
{
    Thermo t
    (
        static_cast<const Thermo&>(bt1) + static_cast<const Thermo&>(bt2)
    );

    if (mag(t.Y()) < small)
    {
        return blottnerTransport<Thermo>
        (
            t,
            bt1.A_,
            bt1.B_,
            bt1.C_
        );
    }
    else
    {
        scalar Y1 = bt1.Y()/t.Y();
        scalar Y2 = bt2.Y()/t.Y();

        return blottnerTransport<Thermo>
        (
            t,
            Y1*bt1.A_ + Y2*bt2.A_,
            Y1*bt1.B_ + Y2*bt2.B_,
            Y1*bt1.C_ + Y2*bt2.C_
        );
    }
}


template<class Thermo>
inline Foam::blottnerTransport<Thermo> Foam::operator*
(
    const scalar s,
    const blottnerTransport<Thermo>& bt
)
{
    return blottnerTransport<Thermo>
    (
        s*static_cast<const Thermo&>(bt),
        bt.A_,
        bt.B_,
        bt.C_
    );
}


// ************************************************************************* //

/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2016 OpenFOAM Foundation
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
    Foam::CRECK2017DevolatilisationV2

Description

\*---------------------------------------------------------------------------*/

#ifndef REACTION_STOIV5_H
#define REACTION_STOIV5_H

class reactionStoiV5
{
    // Private data

        //- Specie name
        word name_;

        //- Molar mass g/mol
        List<scalar> solidStoiData_;


public:

    // Constructors

        //- Null constructor
        reactionStoiV5()
        :
            name_(word::null),
            solidStoiData_(0)
        {}

        //- Construct from Istream
        reactionStoiV5(const word& name, const List<scalar>& l)
        :
            name_(name),
            solidStoiData_(l)
        {}

        // //- Construct as copy
        // solidMolarMass(const solidMolarMass& vd)
        // :
        //     name_(vd.name_),
        //     W_(vd.W_)
        // {}


    //- Destructor
    ~reactionStoiV5()
    {}


    // Public Member Functions

        // Access

            //- Return const access to the name
            const word& name() const
            {
                return name_;
            }

            //- Return const access to the molar mass
            // scalar W() const
            // {
            //     return W_;
            // }

        // //- Read from Istream
        // friend Istream& operator>>(Istream& is, solidMolarMass& vd)
        // {
        //     is.readBeginList("solidMolarMass");
        //     is  >>  vd.name_ >> vd.W_;
        //     is.readEndList("solidMolarMass");
        //
        //     return is;
        // }
        //
        // //- Write to Ostream
        // friend Ostream& operator<<(Ostream& os, const solidMolarMass& vd)
        // {
        //     os  << token::BEGIN_LIST
        //         << vd.name_ << token::SPACE
        //         << vd.W_
        //         << token::END_LIST;
        //
        //     return os;
        // }
};
// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

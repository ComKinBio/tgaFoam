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
    Foam::CRECK2017DevolatilisationV3::volatileData

Description

\*---------------------------------------------------------------------------*/
#ifndef VOLATILEDATAV5_H
#define VOLATILEDATAV5_H

class volatileDataV5
{
    // Private data

        //- Specie name
        word name_;


public:

    // Constructors

        //- Null constructor
        volatileDataV5()
        :
            name_(word::null)
        {}

        //- Construct from Istream
        volatileDataV5(Istream& is)
        :
            name_(is)
        {}

        //- Construct as copy
        // volatileData(const volatileData& vd)
        // :
        //     name_(vd.name_)
        // {}


    //- Destructor
    ~volatileDataV5()
    {}


    // Public Member Functions

        // Access

            //- Return const access to the name
            const word& name() const
            {
                return name_;
            }

        //- Read from Istream
        friend Istream& operator>>(Istream& is, volatileDataV5& vd)
        {
            is.readBeginList("volatileData");
            is  >>  vd.name_;
            is.readEndList("volatileData");

            return is;
        }

        //- Write to Ostream
        friend Ostream& operator<<(Ostream& os, const volatileDataV5& vd)
        {
            os  << token::BEGIN_LIST
                << vd.name_
                << token::END_LIST;

            return os;
        }
};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

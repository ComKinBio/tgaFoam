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

#ifndef READSTOI_ANCA_H
#define READSTOI_ANCA_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
/*---------------------------------------------------------------------------*\
                     Class readStoi Declaration
\*---------------------------------------------------------------------------*/
class readStoiAnca
{
  // Private data

  //- Specie name
  word name_;

  //- stoiCoeff (-)
  scalar stoi_;


public:

  // Constructors

  //- Null constructor
  readStoiAnca()
  :name_(word::null),stoi_(0.0)
  {}

    //- Construct from Istream
    readStoiAnca(Istream& is)
    :
    name_(is),
    stoi_(is)
    {}

      //- Destructor
      ~readStoiAnca(){}

      // Public Member Functions

      // Access

      //- Return const access to the name
      const word& name() const
      {
        return name_;
      }

      //- Return const access to the stoiCoeff
      scalar stoi() const
      {
        return stoi_;
      }

      //- Read from Istream
      friend Istream& operator>>(Istream& is, readStoiAnca& rs)
      {
        is.readBeginList("stoiCoeff");
        is  >>  rs.name_ >> rs.stoi_;
        is.readEndList("stoiCoeff");

        return is;
      }

      //- Write to Ostream
      friend Ostream& operator<<(Ostream& os, const readStoiAnca& rs)
      {
        os  << token::BEGIN_LIST
        << rs.name_ << token::SPACE
        << rs.stoi_
        << token::END_LIST;

        return os;
      }

};

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

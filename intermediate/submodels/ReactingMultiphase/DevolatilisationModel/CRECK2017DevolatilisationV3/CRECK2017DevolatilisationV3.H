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
    Foam::CRECK2017DevolatilisationV3

Description
    Devolatilisation model from CRECK.

\*---------------------------------------------------------------------------*/

#ifndef CRECK2017DevolatilisationV3_H
#define CRECK2017DevolatilisationV3_H

#include "DevolatilisationModel.H"
#include <unordered_map>
#include <map>

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class CRECK2017DevolatilisationV3 Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CRECK2017DevolatilisationV3
:
    public DevolatilisationModel<CloudType>
{
    #include "readStoi.H"
    #include "solidMolarMass.H"
    #include "gasMolarMass.H"
    #include "volatileData.H"
    #include "reactionStoi.H"

    // Private data

        // Model constants

            //- List of volatile data - (name)
            List<volatileData> volatileData_;

            //- List of molar mass of solid - (name, molar mass)
            List<solidMolarMass> solidMolarMass_;

            //- Map of solid molar mass map
            std::unordered_map<std::string,scalar> solidMolarMassMap_;

            //- Map of gas molar mass map
            std::unordered_map<std::string,scalar> gasMolarMassMap_;

            //- Map of solid local id map
            std::unordered_map<std::string,label> solidLocalIdMap_;

            //- Map of gas local id map
            std::unordered_map<std::string,label> gasLocalIdMap_;

            //- List of molar mass of gas - (name, molar mass)
            List<gasMolarMass> gasMolarMass_;

            // Number of reactions
            label reactionNum;

            typedef List<Tuple2<word, scalar>> wordScalarList;
            
            typedef List<Tuple2<label, scalar>> labelScalarList;

            //- List of solid stoichiomery -(name, stoi_)
            List<wordScalarList> solidStoiDatas_;
            
            //- List of id and stoi*molarMass (solid)
            List<labelScalarList> solidStoiMolarmassDatas_;

            //- List of gas stoichiomery -(name, stoi_)
            List<wordScalarList> gasStoiDatas_;
            
            //- List of id and stoi*molarMass (gas)
            List<labelScalarList> gasStoiMolarmassDatas_;
            
            //- List of id of reactants with molar Mass (name, molarMass)
            labelScalarList reactantsIdMMList_;
            
            // List of activation energy by reactions
            List<scalar> EaList;

            // List of pre-exponential factor by reactions
            List<scalar> AList;

            // List of indicator of temperature influence by reactions
            List<scalar> bList;

            // List of reactants by reaction
            List<word> reactantsList;

            //- Mapping between local and cloud gaseous species
            List<label> volatileToGasMap_;

            //- Volatile residual coefficient (0-1)
            //  When the fraction of volatiles are depleted below this
            //  threshold, combustion can occur
            const scalar residualCoeff_;


public:

    //- Runtime type information
    TypeName("CRECK2017DevolatilisationV3");


    // Constructors

    //- Construct from dictionary
    CRECK2017DevolatilisationV3(const dictionary&, CloudType& owner);

    //- Construct copy
    CRECK2017DevolatilisationV3
    (
      const CRECK2017DevolatilisationV3<CloudType>& dm
    );

    //- Construct and return a clone
    virtual autoPtr<DevolatilisationModel<CloudType>> clone() const
    {
      return autoPtr<DevolatilisationModel<CloudType>>
      (
        new CRECK2017DevolatilisationV3<CloudType>(*this)
      );
    }


    //- Destructor
    virtual ~CRECK2017DevolatilisationV3();


    // Member Functions

        //- Update model
        virtual void calculate
        (
            const scalar dt,
            const scalar age,
            const scalar mass0,
            const scalar mass,
            const scalar T,
            const scalarField& YGasEff,
            const scalarField& YLiquidEff,
            const scalarField& YSolidEff,
            label& canCombust,
            scalarField& dMassDV,
            scalarField& dMassSOLID
        ) const;

        // Read reactions from constant folder
        inline void readReactionDict
        (
          const fvMesh& mesh
        );

        void printMap() const {
          Info << "~~~~~~~~~~~~~~"<<solidMolarMassMap_.at("CELL") << nl;
        }

};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CRECK2017DevolatilisationV3.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

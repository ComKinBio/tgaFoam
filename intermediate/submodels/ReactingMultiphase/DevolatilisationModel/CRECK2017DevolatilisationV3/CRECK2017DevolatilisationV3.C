/*---------------------------------------------------------------------------*\
=========                 |
\\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
\\    /   O peration     |
\\  /    A nd           | Copyright (C) 2011-2015 OpenFOAM Foundation
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

#include "CRECK2017DevolatilisationV3.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CRECK2017DevolatilisationV3<CloudType>::CRECK2017DevolatilisationV3
(
  const dictionary& dict,
  CloudType& owner
)
:
DevolatilisationModel<CloudType>(dict, owner, typeName),
volatileData_(this->coeffDict().lookup("volatileData")),
solidMolarMass_(this->coeffDict().lookup("solidMolarMass")),
gasMolarMass_(this->coeffDict().lookup("gasMolarMass")),
volatileToGasMap_(volatileData_.size()),
residualCoeff_(readScalar(this->coeffDict().lookup("residualCoeff")))
{

  if (volatileData_.empty())
  {
    WarningInFunction
    << "Devolatilisation model selected, but no volatiles defined"
    << nl << endl;
  }
  else
  {
    Info<< "Participating volatile species:" << endl;

    // Determine mapping between active volatiles and cloud gas components
    const label idGas = owner.composition().idGas();
    forAll(volatileData_, i)
    {
      const word& specieName = volatileData_[i].name();
      const label id = owner.composition().localId(idGas, specieName);
      volatileToGasMap_[i] = id;

      Info<< "    " << specieName  << endl;
    }
  }

  if (solidMolarMass_.empty() || gasMolarMass_.empty() )
  {
    WarningInFunction
    << "Devolatilisation model selected, but no solid/gas molar mass defined"
    << nl << endl;
  }
  else
  {
    Info<< "Participating solid species:" << endl;

  // Determine mapping between active volatiles and cloud gas components
    forAll(solidMolarMass_, i)
    {
      const word& specieName1 = solidMolarMass_[i].name();
      const scalar MM = solidMolarMass_[i].W();

      Info<< "    " << specieName1  << "   "<< MM << endl;
    }

    forAll(gasMolarMass_, i)
    {
      const word& specieName1 = gasMolarMass_[i].name();
      const scalar MM = gasMolarMass_[i].W();

      Info<< "    " << specieName1  << "   "<< MM << endl;
    }
  }

  // read data from the constant/reactionProperties
  this->readReactionDict(owner.mesh());

  //this->printMap();

}


template<class CloudType>
Foam::CRECK2017DevolatilisationV3<CloudType>::CRECK2017DevolatilisationV3
(
  const CRECK2017DevolatilisationV3<CloudType>& dm
)
:
DevolatilisationModel<CloudType>(dm),
volatileData_(dm.volatileData_),
solidMolarMass_(dm.solidMolarMass_),
gasMolarMass_(dm.gasMolarMass_),
volatileToGasMap_(dm.volatileToGasMap_),
residualCoeff_(dm.residualCoeff_)
{}

// * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CRECK2017DevolatilisationV3<CloudType>::~CRECK2017DevolatilisationV3()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Main calculatioon function

template<class CloudType>
void Foam::CRECK2017DevolatilisationV3<CloudType>::calculate
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
) const
{

  scalar dmassSolid = 0.0;
  // for competitive reactions
  // scalar dmassSolid1 = 0.0;
  // scalar dmassSolid2 = 0.0;
  // scalar dmassSolid3 = 0.0;
  //scalar dmassTot = 0.0;

  scalar kappa = 0.0;
  // for competitive reactions
  // scalar kappa1 = 0.0;
  // scalar kappa2 = 0.0;
  // scalar kappa3 = 0.0;

  scalar reactantMolarMass = 1.0;
  word reactantName;
  label productId=0;

  for (label i = 0; i < reactionNum; i++)
  {
    reactantName = reactantsList[i];

    //const label idReactants = solidLocalIdMap_.at(reactantName);

    const label idReactants = reactantsIdMMList_[i].first();
    kappa = AList[i]*pow(T,bList[i])*exp(-EaList[i]*4184./(RR*T));

    dmassSolid = dt*kappa*mass*YSolidEff[idReactants];

    //reactantMolarMass = solidMolarMassMap_.at(reactantName);
    reactantMolarMass = reactantsIdMMList_[i].second();

    // mole of reactant
    scalar moleReact = dmassSolid/reactantMolarMass;

    // Because CELL, CELLA, HCE, and LIG have multi-steps reactions
    // we have to make each of them react only once

//     if(reactantName == "CELL" && reactantsList[i+1]=="CELL")
//     {
//         kappa1 = AList[i]*pow(T,bList[i])*exp(-EaList[i]*4184./(RR*T));
//         kappa2 = AList[i+1]*pow(T,bList[i+1])*exp(-EaList[i+1]*4184./(RR*T));
//         dmassSolid1 = dt*kappa1*mass*YSolidEff[idReactants];
//         dmassSolid2 = dt*kappa2*mass*YSolidEff[idReactants];
//         dmassTot = dmassSolid1 + dmassSolid2;
//         Info << "check " << reactantName << " " << reactantsList[i+1] << nl;
//
//         if(dmassTot >= mass*YSolidEff[idReactants])
//         {
//             dmassSolid = mass*YSolidEff[idReactants]*dmassSolid1/(dmassTot+ rootVSmall);
//         }
//
//     }
//
//     if(reactantName == "CELL" && reactantsList[i-1]=="CELL")
//     {
//         kappa1 = AList[i]*pow(T,bList[i])*exp(-EaList[i]*4184./(RR*T));
//         kappa2 = AList[i-1]*pow(T,bList[i-1])*exp(-EaList[i-1]*4184./(RR*T));
//         dmassSolid1 = dt*kappa1*mass*YSolidEff[idReactants];
//         dmassSolid2 = dt*kappa2*mass*YSolidEff[idReactants];
//         dmassTot = dmassSolid1 + dmassSolid2+ rootVSmall;
//         Info << "check " << reactantName << " " << reactantsList[i-1] << nl;
//
//         if(dmassTot >= mass*YSolidEff[idReactants])
//         {
//             dmassSolid = mass*YSolidEff[idReactants]*dmassSolid1/(dmassTot+ rootVSmall);
//         }
//     }
//
//     if(reactantName == "CELLA" && reactantsList[i+1]=="CELLA")
//     {
//         kappa1 = AList[i]*pow(T,bList[i])*exp(-EaList[i]*4184./(RR*T));
//         kappa2 = AList[i+1]*pow(T,bList[i+1])*exp(-EaList[i+1]*4184./(RR*T));
//         dmassSolid1 = dt*kappa1*mass*YSolidEff[idReactants];
//         dmassSolid2 = dt*kappa2*mass*YSolidEff[idReactants];
//         dmassTot = dmassSolid1 + dmassSolid2+ rootVSmall;
//         Info << "check " << reactantName << " " << reactantsList[i+1] << nl;
//
//         if(dmassTot >= mass*YSolidEff[idReactants])
//         {
//             dmassSolid = mass*YSolidEff[idReactants]*dmassSolid1/(dmassTot+ rootVSmall);
//         }
//
//     }
//
//     if(reactantName == "CELLA" && reactantsList[i-1]=="CELLA")
//     {
//         kappa1 = AList[i]*pow(T,bList[i])*exp(-EaList[i]*4184./(RR*T));
//         kappa2 = AList[i-1]*pow(T,bList[i-1])*exp(-EaList[i-1]*4184./(RR*T));
//         dmassSolid1 = dt*kappa1*mass*YSolidEff[idReactants];
//         dmassSolid2 = dt*kappa2*mass*YSolidEff[idReactants];
//         dmassTot = dmassSolid1 + dmassSolid2 + rootVSmall;
//         Info << "check " << reactantName << " " << reactantsList[i-1] << nl;
//
//         if(dmassTot >= mass*YSolidEff[idReactants])
//         {
//             dmassSolid = mass*YSolidEff[idReactants]*dmassSolid1/(dmassTot+ rootVSmall);
//         }
//     }
//
//
//     if(reactantName == "LIG" && reactantsList[i+1]=="LIG" && reactantsList[i+2]=="LIG")
//     {
//         kappa1 = AList[i]*pow(T,bList[i])*exp(-EaList[i]*4184./(RR*T));
//         kappa2 = AList[i+1]*pow(T,bList[i+1])*exp(-EaList[i+1]*4184./(RR*T));
//         kappa3 = AList[i+2]*pow(T,bList[i+2])*exp(-EaList[i+2]*4184./(RR*T));
//         dmassSolid1 = dt*kappa1*mass*YSolidEff[idReactants];
//         dmassSolid2 = dt*kappa2*mass*YSolidEff[idReactants];
//         dmassSolid3 = dt*kappa3*mass*YSolidEff[idReactants];
//         dmassTot = dmassSolid1 + dmassSolid2+ dmassSolid3+ rootVSmall;
//         Info << "check " << reactantName << " " << reactantsList[i+1]
//         << " " << reactantsList[i+2]  << nl;
//
//         if(dmassTot >= mass*YSolidEff[idReactants])
//         {
//             dmassSolid = mass*YSolidEff[idReactants]*dmassSolid1/(dmassTot+ rootVSmall);
//         }
//
//     }
//
//     if(reactantName == "LIG" && reactantsList[i-1]=="LIG" && reactantsList[i+1]=="LIG")
//     {
//         kappa1 = AList[i]*pow(T,bList[i])*exp(-EaList[i]*4184./(RR*T));
//         kappa2 = AList[i-1]*pow(T,bList[i-1])*exp(-EaList[i-1]*4184./(RR*T));
//         kappa3 = AList[i+1]*pow(T,bList[i+1])*exp(-EaList[i+1]*4184./(RR*T));
//         dmassSolid1 = dt*kappa1*mass*YSolidEff[idReactants];
//         dmassSolid2 = dt*kappa2*mass*YSolidEff[idReactants];
//         dmassSolid3 = dt*kappa3*mass*YSolidEff[idReactants];
//         dmassTot = dmassSolid1 + dmassSolid2+ dmassSolid3+ rootVSmall;
//         Info << "check " << reactantName << " " << reactantsList[i-1]
//         << " " << reactantsList[i+1]  << nl;
//
//         if(dmassTot >= mass*YSolidEff[idReactants])
//         {
//             dmassSolid = mass*YSolidEff[idReactants]*dmassSolid1/(dmassTot+ rootVSmall);
//         }
//
//     }
//
//     if(reactantName == "LIG" && reactantsList[i-1]=="LIG" && reactantsList[i-2]=="LIG")
//     {
//         kappa1 = AList[i]*pow(T,bList[i])*exp(-EaList[i]*4184./(RR*T));
//         kappa2 = AList[i-1]*pow(T,bList[i-1])*exp(-EaList[i-1]*4184./(RR*T));
//         kappa3 = AList[i-2]*pow(T,bList[i-2])*exp(-EaList[i-2]*4184./(RR*T));
//         dmassSolid1 = dt*kappa1*mass*YSolidEff[idReactants];
//         dmassSolid2 = dt*kappa2*mass*YSolidEff[idReactants];
//         dmassSolid3 = dt*kappa3*mass*YSolidEff[idReactants];
//         dmassTot = dmassSolid1 + dmassSolid2+ dmassSolid3+ rootVSmall;
//         Info << "check " << reactantName << " " << reactantsList[i-1]
//         << " " << reactantsList[i-2]  << nl;
//
//         if(dmassTot >= mass*YSolidEff[idReactants])
//         {
//             dmassSolid = mass*YSolidEff[idReactants]*dmassSolid1/(dmassTot+ rootVSmall);
//         }
//
//     }


    // reactants reduce
    if(dMassSOLID[idReactants] != mass*YSolidEff[idReactants])
    {
        dMassSOLID[idReactants] += dmassSolid;
    }

    if(dmassSolid>mass*YSolidEff[idReactants] && dmassSolid>0.)
    {
        dMassSOLID[idReactants] =
        mass*YSolidEff[idReactants];

        moleReact =
        mass*YSolidEff[idReactants]/reactantMolarMass;
    }


    // solid products
    for(auto it = solidStoiMolarmassDatas_[i].begin();
    it != solidStoiMolarmassDatas_[i].end(); it++)
    {
        productId = it->first();
        if (productId == 999)// 999 means Null
        {
          break;
        }
        dMassSOLID[productId]-=
        moleReact*it->second();
    }

    // gas products
    for(auto it = gasStoiMolarmassDatas_[i].begin();
    it != gasStoiMolarmassDatas_[i].end(); it++)
    {
        productId = it->first();
        if (productId == 999)// 999 means Null
        {
          break;
        }
        dMassDV[volatileToGasMap_[productId]]-=
        moleReact*it->second();
    }

  }

};

// read reaction data from constant folder
template<class CloudType>
void Foam::CRECK2017DevolatilisationV3<CloudType>::readReactionDict
(
  const fvMesh& mesh
)
{

  IOdictionary reactionDict_
  (
    IOobject
    (
      "reactionProperties",
      mesh.time().constant(),
      mesh,
      IOobject::MUST_READ_IF_MODIFIED,
      IOobject::NO_WRITE
    )
  );

  reactionNum = reactionDict_.size();

  word reactionNameTest;

  forAll(reactionDict_,i)
  {

    reactionNameTest = "reaction" + std::to_string(i);

    // read solid phase stoichiomery
    List<Tuple2<word, scalar>> solidStoiData_
    (
      reactionDict_.subDict(reactionNameTest).lookup("solidProducts")
    );
    solidStoiDatas_.append(solidStoiData_);

    // read gas phase stoichiomery
    List<Tuple2<word, scalar>> gasStoiData_
    (
      reactionDict_.subDict(reactionNameTest).lookup("gasProducts")
    );
    gasStoiDatas_.append(gasStoiData_);

    // reaction rate coeff
    scalar Ea = 0.0, A = 0.0, b = 0.0;
    reactionDict_.subDict(reactionNameTest).lookup("Ea") >> Ea;
    reactionDict_.subDict(reactionNameTest).lookup("A") >> A;
    reactionDict_.subDict(reactionNameTest).lookup("b") >> b;
    EaList.append(Ea);
    AList.append(A);
    bList.append(b);

    // read reactants Names in each reactions
    word reactantsName;

    reactionDict_.subDict(reactionNameTest).lookup("reactants")
    >> reactantsName;

    reactantsList.append(reactantsName);

    // Build the solid molar mass map
    forAll(solidMolarMass_, i)
    {
      word speciesName = solidMolarMass_[i].name();

      solidMolarMassMap_[speciesName] = solidMolarMass_[i].W();
    }

    // Build the gas molar mass map
    forAll(gasMolarMass_, i)
    {
      word speciesName = gasMolarMass_[i].name();

      gasMolarMassMap_[speciesName] = gasMolarMass_[i].W();
    }

    // Build solid localId map
    const label idS = this->owner().composition().idSolid();
    List<word> solidNames =
    this->owner().composition().phaseProps()[idS].names();

    forAll(solidNames, i)
    {
      solidLocalIdMap_[solidNames[i]] =
      this->owner().composition().localId(idS, solidNames[i]);
    }



    // Build gas localId map
    const label idG = this->owner().composition().idGas();
    List<word> gasNames =
    this->owner().composition().phaseProps()[idG].names();

    forAll(gasNames, i)
    {
      gasLocalIdMap_[gasNames[i]] =
      this->owner().composition().localId(idG, gasNames[i]);
    }


  }

    //buid solidStoiMolarmassDatas_ for each reaction (Solid)
    forAll(solidStoiDatas_,i)
    {
        List<Tuple2<label,scalar>> IdAndstoiMm;
        forAll(solidStoiDatas_[i], j) // loop List of Tuple<word, scalar>
        {
            const word solidName = solidStoiDatas_[i][j].first();

            // in case no solid products
            if(solidName == "Null")
            {
                const label solidId = 999;
                const scalar stoiMM = 0.0;
                Tuple2<label,scalar> IdStoiMMTuple(solidId,stoiMM);
                break;
            }

            const label solidId =solidLocalIdMap_.at(solidName); // get id
            const scalar MM = solidMolarMassMap_.at(solidName);

            const scalar stoiVal = solidStoiDatas_[i][j].second();
            const scalar stoiMM = stoiVal * MM; // get MM * stoi

            Tuple2<label,scalar> IdStoiMMTuple(solidId,stoiMM);
            IdAndstoiMm.append(IdStoiMMTuple);
        }
        solidStoiMolarmassDatas_.append(IdAndstoiMm);
    }

    //buid gasStoiMolarmassDatas_ for each reaction (gas)
    forAll(gasStoiDatas_,i)
    {
        List<Tuple2<label,scalar>> IdAndstoiMm;
        forAll(gasStoiDatas_[i], j) // loop List of Tuple<word, scalar>
        {
            const word gasName = gasStoiDatas_[i][j].first();

            // in case no solid products
            if(gasName == "Null")
            {
                const label gasId = 999;
                const scalar stoiMM = 0.0;
                Tuple2<label,scalar> IdStoiMMTuple(gasId,stoiMM);
                break;
            }

            const label gasId =gasLocalIdMap_.at(gasName); // get id
            const scalar MM = gasMolarMassMap_.at(gasName);

            const scalar stoiVal = gasStoiDatas_[i][j].second();
            const scalar stoiMM = stoiVal * MM; // get MM * stoi

            Tuple2<label,scalar> IdStoiMMTuple(gasId,stoiMM);
            IdAndstoiMm.append(IdStoiMMTuple);
        }

        gasStoiMolarmassDatas_.append(IdAndstoiMm);
    }

    //buid reactantsIdMMList
    for (label i = 0; i < reactionNum; i++)
    {
        const word reactantName = reactantsList[i];

        const label idReactants = solidLocalIdMap_.at(reactantName);
        const scalar reactantMolarMass = solidMolarMassMap_.at(reactantName);

        Tuple2<label, scalar> IdMMTuple(idReactants, reactantMolarMass);
        reactantsIdMMList_.append(IdMMTuple);
    }


}
// ************************************************************************* //

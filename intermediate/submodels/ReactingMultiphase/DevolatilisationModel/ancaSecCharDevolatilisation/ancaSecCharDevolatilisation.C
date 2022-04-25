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

#include "ancaSecCharDevolatilisation.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::ancaSecCharDevolatilisation<CloudType>::ancaSecCharDevolatilisation
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
Foam::ancaSecCharDevolatilisation<CloudType>::ancaSecCharDevolatilisation
(
  const ancaSecCharDevolatilisation<CloudType>& dm
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
Foam::ancaSecCharDevolatilisation<CloudType>::~ancaSecCharDevolatilisation()
{}

// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

// Main calculatioon function

template<class CloudType>
void Foam::ancaSecCharDevolatilisation<CloudType>::calculate
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
  scalar kappa = 0.0;

  scalar reactantMolarMass = 1.0;
  word reactantName;
  label productId=0;
  scalar reactionWeight = 0.0;

  for (label i = 0; i < reactionNum; i++)
  {
    reactantName = reactantsList[i];

    //const label idReactants = solidLocalIdMap_.at(reactantName);

    const label idReactants = reactantsIdMMList_[i].first();
    //kJ/mol -> J/mol so times 1000000
    kappa = AList[i]*pow(T,bList[i])*exp(-EaList[i]*1000000./(RR*T));

    dmassSolid = dt*kappa*mass*YSolidEff[idReactants];
    reactionWeight += kappa*YSolidEff[idReactants];

    //reactantMolarMass = solidMolarMassMap_.at(reactantName);
    reactantMolarMass = reactantsIdMMList_[i].second();

    // mole of reactant
    scalar moleReact = dmassSolid/reactantMolarMass;

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
  std::ofstream myfile;
  myfile.open("0/weightedK_",std::ios_base::app);
  myfile << T << " " << reactionWeight << "\n";
  myfile.close();

};

// read reaction data from constant folder
template<class CloudType>
void Foam::ancaSecCharDevolatilisation<CloudType>::readReactionDict
(
  const fvMesh& mesh
)
{

  // IOdictionary reactionDict_
  // (
  //   IOobject
  //   (
  //     "reactionProperties",
  //     mesh.time().constant(),
  //     mesh,
  //     IOobject::MUST_READ_IF_MODIFIED,
  //     IOobject::NO_WRITE
  //   )
  // );
  const dictionary& reactionDict_
  (
      IFstream("constant/reactionProperties")()
  );


  // number of reaction steps
  reactionNum = reactionDict_.size();

  word reactionNameTest;

  forAll(reactionDict_,i)
  {

    label reactCur = i;

    reactionNameTest = "reaction" + std::to_string(i);

    // read solid phase stoichiomery

    // primary reactions dict
    const dictionary& primaryReactionsDict_ =
    reactionDict_.subDict(reactionNameTest).subDict("primary");

    // secondary reactions dict
    const dictionary& secondaryReactionsDict_ =
    reactionDict_.subDict(reactionNameTest).subDict("secondary");

    // x fraction coeff for "this" reaction
    scalar xCoeff = 0.0;
    reactionDict_.subDict(reactionNameTest).lookup("xCoeff")>> xCoeff;

    //********************** Process primary stoi !!! ************************//
    // read stoi info for solid phase
    List<Tuple2<word, scalar>> primarySolidStoiData_
    (
      primaryReactionsDict_.lookup("solidProducts")
    );
    // read stoi info for gas phase
    List<Tuple2<word, scalar>> primaryGasStoiData_
    (
      primaryReactionsDict_.lookup("gasProducts")
    );

    // solid (primary)
    List<Tuple2<word, scalar>> primarySolidStoiDataAndX_;

    forAll(primarySolidStoiData_, i)
    {
      // calculate the x * stoi
      scalar xTimesStoi = primarySolidStoiData_[i].second()*(1-xCoeff);
      // pack the name and the x * stoi
      Tuple2<word, scalar>
      nameAndxTimesStoi(primarySolidStoiData_[i].first(), xTimesStoi);
      // put the tuple into the list
      primarySolidStoiDataAndX_.append(nameAndxTimesStoi);
    }

    // gas (primary)
    List<Tuple2<word, scalar>> primaryGasStoiDataAndX_;

    forAll(primaryGasStoiData_, i)
    {
      // calculate the x * stoi
      scalar xTimesStoi = primaryGasStoiData_[i].second()*(1-xCoeff);
      // pack the name and the x * stoi
      Tuple2<word,scalar>
      nameAndxTimesStoi(primaryGasStoiData_[i].first(), xTimesStoi);
      // put the tuple into the list
      primaryGasStoiDataAndX_.append(nameAndxTimesStoi);
    }

    //********************* Process secondary stoi !!!************************//
    // read stoi info for solid phase
    List<Tuple2<word, scalar>> secondarySolidStoiData_
    (
      secondaryReactionsDict_.lookup("solidProducts")
    );
    // read stoi info for gas phase
    List<Tuple2<word, scalar>> secondaryGasStoiData_
    (
      secondaryReactionsDict_.lookup("gasProducts")
    );

    // solid (primary)
    List<Tuple2<word, scalar>> secondarySolidStoiDataAndX_;

    forAll(secondarySolidStoiData_, i)
    {
      // calculate the x * stoi
      scalar xTimesStoi = secondarySolidStoiData_[i].second()*xCoeff;
      // pack the name and the x * stoi
      Tuple2<word, scalar>
      nameAndxTimesStoi(secondarySolidStoiData_[i].first(), xTimesStoi);
      // put the tuple into the list
      secondarySolidStoiDataAndX_.append(nameAndxTimesStoi);
    }

    // gas (primary)
    List<Tuple2<word, scalar>> secondaryGasStoiDataAndX_;

    forAll(secondaryGasStoiData_, i)
    {
      // calculate the x * stoi
      scalar xTimesStoi = secondaryGasStoiData_[i].second()*xCoeff;
      // pack the name and the x * stoi
      Tuple2<word,scalar>
      nameAndxTimesStoi(secondaryGasStoiData_[i].first(), xTimesStoi);
      // put the tuple into the list
      secondaryGasStoiDataAndX_.append(nameAndxTimesStoi);
    }
    // *******end processing primary and secondary gas stoi data *************//

    //combine the primary and secondary stoi data (solid phase)
    solidStoiDatas_.append(primarySolidStoiDataAndX_);
    forAll(secondarySolidStoiDataAndX_,i)
    {
      forAll(solidStoiDatas_[reactCur],j)
      {
        if
        (
          solidStoiDatas_[reactCur][j].first()
          == secondarySolidStoiDataAndX_[i].first()
        )
        {
          solidStoiDatas_[reactCur][j].second()
          +=secondarySolidStoiDataAndX_[i].second();
          break;
        }
        else if(
          solidStoiDatas_[reactCur][j].first()
          != secondarySolidStoiDataAndX_[i].first()
          &&
          j == solidStoiDatas_[reactCur].size()
        )
        {
          Tuple2<word, scalar> secDiff
          (secondarySolidStoiDataAndX_[i]);
          solidStoiDatas_[reactCur].append(secondarySolidStoiDataAndX_[i]);
        }
      }
    }

    //combine the primary and secondary stoi data (gas phase)
    gasStoiDatas_.append(primaryGasStoiDataAndX_);
    forAll(secondaryGasStoiDataAndX_,i)
    {
      forAll(gasStoiDatas_[reactCur],j)
      {
        if
        (
          gasStoiDatas_[reactCur][j].first()
          == secondaryGasStoiDataAndX_[i].first()
        )
        {
          gasStoiDatas_[reactCur][j].second()
          +=secondaryGasStoiDataAndX_[i].second();
          break;
        }
        else if(
          gasStoiDatas_[reactCur][j].first()
          != secondaryGasStoiDataAndX_[i].first()
          &&
          j == gasStoiDatas_[reactCur].size()
        )
        {
          Tuple2<word, scalar> secDiff
          (secondaryGasStoiDataAndX_[i]);
          gasStoiDatas_[reactCur].append(secondaryGasStoiDataAndX_[i]);
        }
      }
    }

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
    const label idS = this->owner().composition().idSolid(); // phase Id
    List<word> solidNames =
    this->owner().composition().phaseProps()[idS].names();

    forAll(solidNames, i)
    {
      solidLocalIdMap_[solidNames[i]] =
      this->owner().composition().localId(idS, solidNames[i]); // species Id
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
    // loop List of Tuple<word, scalar>
    // solidStoiDatas_[i] has type List<Tuple<word, scalar>>
    forAll(solidStoiDatas_[i], j)
    {
      // Now both primary and secondary reactions should be included
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

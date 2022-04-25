// Info<< "before the List is called " << nl;
//
// List<phaseProperties> ppst=  coalParcels.composition().coeffDict().lookup("phases");
//
// //Info<< coalParcels.composition().coeffDict().lookup("phases") << nl;
//
// Info << "print the list ***********" << nl;
//
// forAll(ppst,i)
// {
//   Info << ppst[i]<<nl;
//   Info << "the next " << nl;
// }
//
// Istream& iss = coalParcels.composition().coeffDict().lookup("phases");
//
// label parentdic = coalParcels.composition().coeffDict().startLineNumber();
//
// Info << "the topDict is " << parentdic << nl;
//
// token firstToken(iss);
//
// if (firstToken.isLabel())
// {
//   Info << "isLabel" << nl;
// }
// else if (firstToken.isCompound())
// {
//   Info << "isCompound" << nl;
// }
// else if (firstToken.isPunctuation())
// {
//   Info << "isPunctuation" << nl;
// }
// else
// {
//   FatalIOErrorInFunction(iss)
//       << "incorrect first token, expected <int> or '(', found "
//       << firstToken.info()
//       << exit(FatalIOError);
// }
//
// Info << "the line number is " << firstToken.lineNumber()<< nl;
//
// Info << "the wordToken is " << firstToken.wordToken() << nl;
//
// Info << "the labelToken is " << firstToken.labelToken() << nl;
//
// Info << "the punctuationToken is " << firstToken.pToken() << nl;
//
// Info << "the type is " << firstToken.type() << nl;


// test phaseProperties

phaseProperties myPhaseProperties = coalParcels.composition().phaseProps().props()[0];

Info << "the phaseTypeName is "<< myPhaseProperties.phaseTypeName() << nl;

List<phaseProperties> myPhasePropertiesList(coalParcels.composition().phaseProps().props());

forAll(myPhasePropertiesList, i)
{
  Info << "the phaseTypeName is "<< myPhasePropertiesList[i].phaseTypeName() << nl;
  auto listOfSpicies = myPhasePropertiesList[i].names();
  forAll(listOfSpicies,j)
  {
    Info << "the spicie " << j<< " is " << listOfSpicies[j] << nl;
  }
}

// spicies fractio
auto myY = myPhaseProperties.Y();
Info << "myY is " << myY << nl;

Info << "YMixture0() is" << coalParcels.composition().YMixture0() << nl;

// thermo
Info << "thermo carrier " << coalParcels.composition().carrier().species() << nl;
Info << "thermo liquids components " << coalParcels.composition().liquids().components() << nl;
Info << "thermo solids components" << coalParcels.composition().solids().components() << nl;
Info << "thermo solids components" << slgThermo.solids().components() << nl;

Info << "Thermo::fvMeshConstructorTablePtr_ " <<psiReactionThermo::fvMeshConstructorTablePtr_ << nl;

for (auto it = psiReactionThermo::fvMeshConstructorTablePtr_->begin() ;
          it != psiReactionThermo::fvMeshConstructorTablePtr_->end(); it++) {
  Info << it.key() << nl;
}

for (auto it = functionObject::dictionaryConstructorTablePtr_->begin() ;
          it != functionObject::dictionaryConstructorTablePtr_->end(); it++) {
  Info << it.key() << nl;
}

Info << *functionObject::dictionaryConstructorTablePtr_ << nl;

basicSpecieMixture& mcThermoTest =
    dynamic_cast<basicSpecieMixture&>(thermo);
Info << "test my dynamic_cast " << mcThermoTest.species() << nl;
Info << "what is pThermo" << pThermo << nl;
Info << "what is thermo.subDict(liquids) " << thermo.subDict("liquids") << nl;

List<word> componentsBoyao = thermo.subDict("liquids").toc();
Info << "componentsBoyao[0]" << componentsBoyao[0] << nl;
    forAll(componentsBoyao, i)
    {
        if (thermo.subDict("liquids").isDict(componentsBoyao[i]))
        {
            Info << thermo.subDict("liquids").isDict(componentsBoyao[i]) << nl;
        }
        else
        {
           Info << "it is a joke" << nl;
        }
    }

    for (auto it = liquidProperties::dictionaryConstructorTablePtr_->begin() ;
              it != liquidProperties::dictionaryConstructorTablePtr_->end(); it++) {
      Info << it.key() << nl;
    }

    Info << "~~~~~~~~~~~~~~~~~~~~" << nl;
    for (auto it = solidProperties::dictionaryConstructorTablePtr_->begin() ;
              it != solidProperties::dictionaryConstructorTablePtr_->end(); it++) {
      Info << it.key() << nl;
    }


    Info << "the all components are " << slgThermo.solids().components() << nl;
    forAll(slgThermo.solids().components(),i)
    {
      Info << "component " << i  << " of the solids is "
      << slgThermo.solids().components()[i] << nl;
    }

    Info << "typeId " << thermo.type() << nl;

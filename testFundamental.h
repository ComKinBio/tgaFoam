List<label> boyaoWordList;

for (size_t i = 0; i < 5; i++) {
  boyaoWordList.append(i);
}

forAll(boyaoWordList,i)
{
  Info << boyaoWordList[i] << nl;
}

//Istream& iss = coalParcels.composition().coeffDict().lookup("phases");
List<word> listOfPhases(thermo.subDict("solids").toc());

forAll(listOfPhases, i)
{
  Info << "list of phase " << i  << " is "<< listOfPhases[i] << nl;
}

Info << "build list of phase properties !!! " << nl; 

List<phaseProperties> listOfPhaseProperties
(
  coalParcels.composition().coeffDict().lookup("phases")
);

forAll(listOfPhaseProperties,i)
{
  Info << "list of phaseProperties " << i << listOfPhaseProperties[i]
       <<  nl;
}

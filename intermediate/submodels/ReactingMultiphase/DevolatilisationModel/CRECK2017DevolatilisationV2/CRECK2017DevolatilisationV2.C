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

#include "CRECK2017DevolatilisationV2.H"

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

template<class CloudType>
Foam::CRECK2017DevolatilisationV2<CloudType>::CRECK2017DevolatilisationV2
(
  const dictionary& dict,
  CloudType& owner
)
:
DevolatilisationModel<CloudType>(dict, owner, typeName),
volatileData_(this->coeffDict().lookup("volatileData")),
solidMolarMass_(this->coeffDict().lookup("solidMolarMass")),
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

  if (solidMolarMass_.empty())
  {
    WarningInFunction
    << "Devolatilisation model selected, but no solid molar mass defined"
    << nl << endl;
  }
  else
  {
    Info<< "Participating solid species:" << endl;

    // Determine mapping between active volatiles and cloud gas components
  //  const label idGas = owner.composition().idGas();
    forAll(solidMolarMass_, i)
    {
      const word& specieName1 = solidMolarMass_[i].name();
      const scalar MM = solidMolarMass_[i].W();

      Info<< "    " << specieName1  << "   "<< MM << endl;
    }
  }

}


template<class CloudType>
Foam::CRECK2017DevolatilisationV2<CloudType>::CRECK2017DevolatilisationV2
(
  const CRECK2017DevolatilisationV2<CloudType>& dm
)
:
DevolatilisationModel<CloudType>(dm),
volatileData_(dm.volatileData_),
solidMolarMass_(dm.solidMolarMass_),
volatileToGasMap_(dm.volatileToGasMap_),
residualCoeff_(dm.residualCoeff_)
{}


  // * * * * * * * * * * * * * * * * Destructor  * * * * * * * * * * * * * * * //

  template<class CloudType>
  Foam::CRECK2017DevolatilisationV2<CloudType>::~CRECK2017DevolatilisationV2()
  {}


    // * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

    template<class CloudType>
    void Foam::CRECK2017DevolatilisationV2<CloudType>::calculate
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
      //Info <<"start CRECK2017DevolatilisationV2"<<endl;
      //bool done = true;

      //- Get id for wood char tar
      //const label idGas = this->owner().composition().idGas();
      const label idSolid = this->owner().composition().idSolid();
      //     const label idLiquid = this->owner().composition().idLiquid();
      //const scalar YSolidTot = this->owner().composition().YMixture0()[idSolid];
      //const scalarField& YSolid = this->owner().composition().Y0(idSolid);

      //number of solid reactions
      const label nSolid = 27;

      // Model coefficients


      //2015
      //     const scalarList A={4.00000E+013,2.00000E+006,4.00000E+000,6.50000E+007,1.00000E+010,1.20000E+009,
        //                         0.15000E+000,3.00000E+000,5.00000E+009,1.33000E+015,6.70000E+012,3.30000E+008,
        //                         1.67000E+006,1.10000E+008,4.00000E+000,4.00000E+008,8.30000E-002,7.00000E+012,
        //                         5.00000E+001,1.50000E-002,1.00000E+006,5.00000E+012,5.00000E+011,5.00000E+011,
        //                         5.00000E+012,2.00000E+012,5.00000E+012}; //pre-exponential factor s-1
        //     const scalarList b={0.,0.,1.,0.,0.,0.,1.,1.,0.,0.,0.,0.,0.,0.,1.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}; //factor, -
        //     const scalarList Ea={45000.00,19100.00,10000.00,31000.00,31000.00,30000.00, 8000.00,11000.00,
          //                          33000.00,48500.00,37500.00,25500.00,31500.00,30000.00,12000.00,30000.00,
          //                           8000.00,45700.00,11000.00, 6100.00,24000.00,50000.00,71000.00,75000.00,
          //                          71667.00,50000.00,71667.00}; //activation energy cal/mol

          //2017
          const scalarList A={1.50000E+014,2.0000E+006,4.00000E+000,6.50000E+007,1.00000E+010,1.20000E+009,
            0.15000E+000,3.00000E+000,5.00000E+009,1.33000E+015,6.70000E+012,3.30000E+008,
            1.67000E+006,1.10000E+008,4.00000E+000,4.00000E+008,8.30000E-002,7.00000E+012,
            5.00000E+001,1.50000E-002,1.00000E+006,5.00000E+012,5.00000E+011,5.00000E+011,
            5.00000E+012,2.00000E+012,5.00000E+012}; //pre-exponential factor s-1
            const scalarList b={0.,0.,1.,0.,0.,0.,1.,1.,0.,0.,0.,0.,0.,0.,1.,0.,1.,0.,0.,0.,0.,0.,0.,0.,0.,0.,0.}; //factor, -
            const scalarList Ea={47000.00,19100.00,10000.00,31000.00,31000.00,30000.00, 8000.00,11000.00,
              33000.00,48500.00,37500.00,25500.00,31500.00,30000.00,12000.00,30000.00,
              8000.00,45700.00,11000.00, 6100.00,24000.00,50000.00,71000.00,75000.00,
              71667.00,50000.00,71667.00}; //activation energy cal/mol

              scalarList kappa(nSolid);
              scalarList dmassSolid(nSolid);
              scalar dmassTemp;


              const label id_CELL = this->owner().composition().localId(idSolid, "CELL");
              const label id_CELLA = this->owner().composition().localId(idSolid, "CELLA");
              const label id_HCE = this->owner().composition().localId(idSolid, "HCE");
              const label id_HCE1 = this->owner().composition().localId(idSolid, "HCE1");
              const label id_HCE2 = this->owner().composition().localId(idSolid, "HCE2");
              const label id_LIGC = this->owner().composition().localId(idSolid, "LIG-C");
              const label id_LIGH = this->owner().composition().localId(idSolid, "LIG-H");
              const label id_LIGO = this->owner().composition().localId(idSolid, "LIG-O");
              const label id_LIGCC = this->owner().composition().localId(idSolid, "LIG-CC");
              const label id_LIGOH = this->owner().composition().localId(idSolid, "LIG-OH");
              const label id_LIG = this->owner().composition().localId(idSolid, "LIG");
              const label id_TGL = this->owner().composition().localId(idSolid, "TGL");
              const label id_CTANN = this->owner().composition().localId(idSolid, "CTANN");
              const label id_ITANN = this->owner().composition().localId(idSolid, "ITANN");
              const label id_GCO2 = this->owner().composition().localId(idSolid, "GCO2");
              const label id_GCO = this->owner().composition().localId(idSolid, "GCO");
              const label id_GCOH2 = this->owner().composition().localId(idSolid, "GCOH2");
              const label id_GH2 = this->owner().composition().localId(idSolid, "GH2");
              const label id_GCH4 = this->owner().composition().localId(idSolid, "GCH4");
              const label id_GCH3OH = this->owner().composition().localId(idSolid, "GCH3OH");
              const label id_GC2H4 = this->owner().composition().localId(idSolid, "GC2H4");
              const label id_C = this->owner().composition().localId(idSolid, "C");

              // Kinetic rate
              forAll(kappa, i)
              {
                kappa[i] = A[i]*pow(T,b[i])*exp(-Ea[i]*4184./(RR*T));
              }
              // 0	H2	2	                (CELL            162.0)
              // 1	H2O	18	                (CELLA           162.0)
              // 2	CO	28	                (HCE             132.0)
              // 3	CO2	44	                (HCE1            132.0)
              // 4	CH2O	30	                (HCE2            132.0)
              // 5	HCOOH	46	                (LIGC            258.0)
              // 6	CH4	16	                (LIGH            436.0)
              // 7	CH3OH	32	                (LIGO            422.0)
              // 8	C2H2O2	58	                (LIGCC           258.0)
              // 9	C2H4	28	                (LIGOH           378.0)
              // 10	CH3CHO	44	                (LIG             208.0)
              // 11	C2H4O2	60	                (TGL           896.0)//C57H100O7
              // 12	C2H5OH	46	                (CTANN           304.0)//C15H12O7
              // 13	C2H3CHO	56	                (ITANN           210.0)//C9H6O6
              // 14	C2H5CHO	58	                (GCO2            44.0)
              // 15	C5H8O4	132	                (GCO             28.0)
              // 16	C6H5OH	94	                (GCOH2           30.0)
              // 17	C6H6O3	126	                (GH2             2.0)
              // 18	C6H10O5	162	                (GCH4              16.0)
              // 19	C6H5OCH3	108	                (GCH3OH             32.0)
              // 20	C9H10O2	150	                (GC2H4             28.0)
              // 21	C11H12O4	208	                (C               12.0)
              // delete   22	LINOACID	280	delete                (ash             60.0)
              // 22  U2ME12 (C13H22O2, double_unsaturated_methyl_ester)  210
              // 23  MLINO (C19H34O2, methyl_linoleate)  294
              scalar massKg = 1000 * mass;
              //CELL=>CELLA
              dmassSolid[0] = dt*kappa[0]*massKg*YSolidEff[id_CELL];
              //CELL=>5H2O+6C(S)
              dmassSolid[3] = dt*kappa[3]*massKg*YSolidEff[id_CELL];
              if (dmassSolid[0]+dmassSolid[3]>massKg*YSolidEff[id_CELL] && dmassSolid[0]+dmassSolid[3]>0.)
              {
                dmassTemp = dmassSolid[0]+dmassSolid[3];
                dmassSolid[0] = massKg*YSolidEff[id_CELL]*dmassSolid[0]/dmassTemp;
                dmassSolid[3] = massKg*YSolidEff[id_CELL]*dmassSolid[3]/dmassTemp;
              }
              dMassSOLID[id_CELL] = dMassSOLID[id_CELL] + dmassSolid[0]+dmassSolid[3]; //CELL
              dMassSOLID[id_CELLA] = dMassSOLID[id_CELLA] - dmassSolid[0]; //CELLA
              dMassDV[volatileToGasMap_[1]] = dMassDV[volatileToGasMap_[1]] + (dmassSolid[3]/solidMolarMass_[0].W())*5.*18.; //H2O
              dMassSOLID[id_C] = dMassSOLID[id_C] -(dmassSolid[3]/solidMolarMass_[0].W())*6.*12.; //C

              //CELLA=>  0.45 C2H4O2 + 0.2 C2H2O2 + 0.1 CH3CHO + 0.25 HMFU + 0.3 C2H5CHO + 0.15 CH3OH + 0.4 CH2O + 0.31 CO
              //		 + 0.41 CO2 + 0.05 H2 + 0.83 H2O + 0.02 HCOOH + 0.2 G{CH4} + 0.05 G{H2} + 0.61 CHAR
              dmassSolid[1] = dt*kappa[1]*massKg*YSolidEff[id_CELLA];
              //CELLA=>C6H10O5
              dmassSolid[2] = dt*kappa[2]*massKg*YSolidEff[id_CELLA];
              if (dmassSolid[1]+dmassSolid[2]>massKg*YSolidEff[id_CELLA] && dmassSolid[1]+dmassSolid[2]>0.)
              {
                dmassTemp = dmassSolid[1]+dmassSolid[2];
                dmassSolid[1] = massKg*YSolidEff[id_CELLA]*dmassSolid[1]/dmassTemp;
                dmassSolid[2] = massKg*YSolidEff[id_CELLA]*dmassSolid[2]/dmassTemp;
              }
              dMassSOLID[id_CELLA] = dMassSOLID[id_CELLA] + dmassSolid[1]+dmassSolid[2]; //CELLA
              dMassDV[volatileToGasMap_[11]] = dMassDV[volatileToGasMap_[11]] + (dmassSolid[1]/solidMolarMass_[1].W())*0.45*60.; //C2H4O2
              dMassDV[volatileToGasMap_[8]] = dMassDV[volatileToGasMap_[8]] + (dmassSolid[1]/solidMolarMass_[1].W())*0.2*58.; // C2H2O2
              dMassDV[volatileToGasMap_[10]] = dMassDV[volatileToGasMap_[10]] + (dmassSolid[1]/solidMolarMass_[1].W())*0.1*44.; // CH3CHO
              dMassDV[volatileToGasMap_[17]] = dMassDV[volatileToGasMap_[17]] +(dmassSolid[1]/solidMolarMass_[1].W())*0.25*126.; // C6H6O3 HMFU
              dMassDV[volatileToGasMap_[14]] = dMassDV[volatileToGasMap_[14]] +(dmassSolid[1]/solidMolarMass_[1].W())*0.23*58.; // C2H5CHO
              dMassDV[volatileToGasMap_[7]] = dMassDV[volatileToGasMap_[7]] + (dmassSolid[1]/solidMolarMass_[1].W())*0.15*32.; // CH3OH
              dMassDV[volatileToGasMap_[4]] = dMassDV[volatileToGasMap_[4]] + (dmassSolid[1]/solidMolarMass_[1].W())*0.4*30.; // CH2O
              dMassDV[volatileToGasMap_[2]] = dMassDV[volatileToGasMap_[2]] + (dmassSolid[1]/solidMolarMass_[1].W())*0.31*28.; // CO
              dMassDV[volatileToGasMap_[3]] = dMassDV[volatileToGasMap_[3]] + (dmassSolid[1]/solidMolarMass_[1].W())*0.41*44.; // CO2
              dMassDV[volatileToGasMap_[0]] = dMassDV[volatileToGasMap_[0]] + (dmassSolid[1]/solidMolarMass_[1].W())*0.05*2.; // H2
              dMassDV[volatileToGasMap_[1]] = dMassDV[volatileToGasMap_[1]] + (dmassSolid[1]/solidMolarMass_[1].W())*0.83*18.; // H2O
              dMassDV[volatileToGasMap_[5]] = dMassDV[volatileToGasMap_[5]] + (dmassSolid[1]/solidMolarMass_[1].W())*0.02*46.; // HCOOH
              dMassSOLID[id_GCH4] = dMassSOLID[id_GCH4] -(dmassSolid[1]/solidMolarMass_[1].W())*0.2*16.; //GCH4
              dMassSOLID[id_GH2] = dMassSOLID[id_GH2] -(dmassSolid[1]/solidMolarMass_[1].W())*0.05*2.; //GH2
              dMassSOLID[id_C] = dMassSOLID[id_C] -(dmassSolid[1]/solidMolarMass_[1].W())*0.61*12.; //C
              dMassDV[volatileToGasMap_[18]] = dMassDV[volatileToGasMap_[18]] + dmassSolid[2]; // C6H10O5
              //HCE=>0.58HCE1+0.42HCE2
              dmassSolid[4] = min(dt*kappa[4], 1.)*massKg*YSolidEff[id_HCE];
              dMassSOLID[id_HCE] = dMassSOLID[id_HCE] + dmassSolid[4]; //HCE

              //2017
              dMassSOLID[id_HCE1] = dMassSOLID[id_HCE1] - 0.58*dmassSolid[4]; //HCE1
              dMassSOLID[id_HCE2] = dMassSOLID[id_HCE2] - 0.42*dmassSolid[4]; //HCE2

              //2015
              //     dMassSOLID[id_HCE1] = dMassSOLID[id_HCE1] - 0.5*dmassSolid[4]; //HCE1
              //     dMassSOLID[id_HCE2] = dMassSOLID[id_HCE2] - 0.5*dmassSolid[4]; //HCE2

              //HCE1=>0.025 H2O + 0.5 CO2 + 0.025 HCOOH + 0.5 CO + 0.8 CH2O + 0.125 C2H5OH + 0.1 CH3OH + 0.25 C2H4 + 0.125 G{H2} + 0.275 G{CO2} + 0.4 G{COH2}
              //      + 0.45 G{CH3OH} + 0.325 G{CH4} + 0.875 CHAR
              dmassSolid[5] = dt*kappa[5]*massKg*YSolidEff[id_HCE1];
              //HCE1=>0.25 H2O + 0.8 CO2 + 0.05 HCOOH + 0.1 CO + 0.15 G{CO} + 0.15 G{CO2} + 0.2 G{H2} + 0.3 CH2O + 1.2 G{COH2} + 0.625 G{CH4} + 0.375 G{C2H4} + 0.875 CHAR
              dmassSolid[6] = dt*kappa[6]*massKg*YSolidEff[id_HCE1];
              //HCE1=>C5H8O4（XYLAN）
              dmassSolid[7] = dt*kappa[7]*massKg*YSolidEff[id_HCE1];
              if (dmassSolid[5]+dmassSolid[6]+dmassSolid[7]>massKg*YSolidEff[id_HCE1] && dmassSolid[5]+dmassSolid[6]+dmassSolid[7]>0.)
              {
                dmassTemp = dmassSolid[5]+dmassSolid[6]+dmassSolid[7];
                dmassSolid[5] = massKg*YSolidEff[id_HCE1]*dmassSolid[5]/dmassTemp;
                dmassSolid[6] = massKg*YSolidEff[id_HCE1]*dmassSolid[6]/dmassTemp;
                dmassSolid[7] = massKg*YSolidEff[id_HCE1]*dmassSolid[7]/dmassTemp;
              }
              dMassSOLID[id_HCE1] = dMassSOLID[id_HCE1]+dmassSolid[5]+dmassSolid[6]+dmassSolid[7]; //HCE1
              dMassDV[volatileToGasMap_[1]] = dMassDV[volatileToGasMap_[1]] + (dmassSolid[5]/solidMolarMass_[3].W())*0.025*18.; // H2O
              dMassDV[volatileToGasMap_[3]] = dMassDV[volatileToGasMap_[3]] + (dmassSolid[5]/solidMolarMass_[3].W())*0.5*44.; // CO2
              dMassDV[volatileToGasMap_[5]] = dMassDV[volatileToGasMap_[5]] + (dmassSolid[5]/solidMolarMass_[3].W())*0.025*46.; // HCOOH
              dMassDV[volatileToGasMap_[2]] = dMassDV[volatileToGasMap_[2]] + (dmassSolid[5]/solidMolarMass_[3].W())*0.5*28.; // CO
              dMassDV[volatileToGasMap_[4]] = dMassDV[volatileToGasMap_[4]] + (dmassSolid[5]/solidMolarMass_[3].W())*0.8*30.; // CH2O
              dMassDV[volatileToGasMap_[12]] = dMassDV[volatileToGasMap_[12]] + (dmassSolid[5]/solidMolarMass_[3].W())*0.125*46.; // C2H5OH
              dMassDV[volatileToGasMap_[7]] = dMassDV[volatileToGasMap_[7]] + (dmassSolid[5]/solidMolarMass_[3].W())*0.1*32.; // CH3OH
              dMassDV[volatileToGasMap_[9]] = dMassDV[volatileToGasMap_[9]] + (dmassSolid[5]/solidMolarMass_[3].W())*0.25*28.; // C2H4
              dMassSOLID[id_GH2] = dMassSOLID[id_GH2] - (dmassSolid[5]/solidMolarMass_[3].W())*0.125*2.; //G(H2)
              dMassSOLID[id_GCO2] = dMassSOLID[id_GCO2] - (dmassSolid[5]/solidMolarMass_[3].W())*0.275*44.; //G(CO2)
              dMassSOLID[id_GCOH2] = dMassSOLID[id_GCOH2] - (dmassSolid[5]/solidMolarMass_[3].W())*0.4*30.; //G(COH2)
              dMassSOLID[id_GCH3OH] = dMassSOLID[id_GCH3OH] - (dmassSolid[5]/solidMolarMass_[3].W())*0.45*32.; //G(CH3OH)
              dMassSOLID[id_GCH4] = dMassSOLID[id_GCH4] - (dmassSolid[5]/solidMolarMass_[3].W())*0.325*16.; //G(CH4)
              dMassSOLID[id_C] = dMassSOLID[id_C] - (dmassSolid[5]/solidMolarMass_[3].W())*0.875*12.; //C
              dMassDV[volatileToGasMap_[1]] = dMassDV[volatileToGasMap_[1]] + (dmassSolid[6]/solidMolarMass_[3].W())*0.25*18.; // H2O
              dMassDV[volatileToGasMap_[3]] = dMassDV[volatileToGasMap_[3]] + (dmassSolid[6]/solidMolarMass_[3].W())*0.8*44.; // CO2
              dMassDV[volatileToGasMap_[5]] = dMassDV[volatileToGasMap_[5]] + (dmassSolid[6]/solidMolarMass_[3].W())*0.05*46.; // HCOOH
              dMassDV[volatileToGasMap_[2]] = dMassDV[volatileToGasMap_[2]] + (dmassSolid[6]/solidMolarMass_[3].W())*0.1*28.; // CO
              dMassSOLID[id_GCO] = dMassSOLID[id_GCO] - (dmassSolid[6]/solidMolarMass_[3].W())*0.15*28.; //G(CO)
              dMassSOLID[id_GCO2] = dMassSOLID[id_GCO2] - (dmassSolid[6]/solidMolarMass_[3].W())*0.15*44.; //G(CO2)
              dMassSOLID[id_GH2] = dMassSOLID[id_GH2] - (dmassSolid[6]/solidMolarMass_[3].W())*0.2*2.; //G(H2)
              dMassDV[volatileToGasMap_[4]] = dMassDV[volatileToGasMap_[4]] + (dmassSolid[6]/solidMolarMass_[3].W())*0.3*30.; // CH2O
              dMassSOLID[id_GCOH2] = dMassSOLID[id_GCOH2] - (dmassSolid[6]/solidMolarMass_[3].W())*1.2*30.; //G(COH2)
              dMassSOLID[id_GCH4] = dMassSOLID[id_GCH4] - (dmassSolid[6]/solidMolarMass_[3].W())*0.625*16.; //G(CH4)
              dMassSOLID[id_GC2H4] = dMassSOLID[id_GC2H4] - (dmassSolid[6]/solidMolarMass_[3].W())*0.375*28.; //G(C2H4)
              dMassSOLID[id_C] = dMassSOLID[id_C] - (dmassSolid[6]/solidMolarMass_[3].W())*0.875*12.; //C
              dMassDV[volatileToGasMap_[15]] = dMassDV[volatileToGasMap_[15]] + dmassSolid[6]; // C5H8O4

              //HCE2=>0.2 H2O + 0.175 CO + 0.275 CO2 + 0.5 CH2O + 0.1 C2H5OH + 0.2 C2H4O2 + 0.025 HCOOH + 0.25 G{CH4} + 0.3 G{CH3OH} + 0.275 G{C2H4} + 0.4 G{CO2} + 0.925 G{COH2} + CHAR
              dmassSolid[8] = min(dt*kappa[8], 1.)*massKg*YSolidEff[id_HCE2];
              dMassSOLID[id_HCE2] = dMassSOLID[id_HCE2] + dmassSolid[8]; //HCE2
              dMassDV[volatileToGasMap_[1]] = dMassDV[volatileToGasMap_[1]] + (dmassSolid[8]/solidMolarMass_[4].W())*0.2*18.; // H2O
              dMassDV[volatileToGasMap_[2]] = dMassDV[volatileToGasMap_[2]] + (dmassSolid[8]/solidMolarMass_[4].W())*0.175*28.; // CO
              dMassDV[volatileToGasMap_[3]] = dMassDV[volatileToGasMap_[3]] + (dmassSolid[8]/solidMolarMass_[4].W())*0.275*44.; // CO2
              dMassDV[volatileToGasMap_[4]] = dMassDV[volatileToGasMap_[4]] + (dmassSolid[8]/solidMolarMass_[4].W())*0.5*30.; // CH2O
              dMassDV[volatileToGasMap_[12]] = dMassDV[volatileToGasMap_[12]] + (dmassSolid[8]/solidMolarMass_[4].W())*0.1*46.; // C2H5OH
              dMassDV[volatileToGasMap_[11]] = dMassDV[volatileToGasMap_[11]] + (dmassSolid[8]/solidMolarMass_[4].W())*0.2*60.; // C2H4O2
              dMassDV[volatileToGasMap_[5]] = dMassDV[volatileToGasMap_[5]] + (dmassSolid[8]/solidMolarMass_[4].W())*0.025*46.; // HCOOH
              dMassSOLID[id_GCH4] = dMassSOLID[id_GCH4] - (dmassSolid[8]/solidMolarMass_[4].W())*0.25*16.; //G(CH4)
              dMassSOLID[id_GCH3OH] = dMassSOLID[id_GCH3OH] - (dmassSolid[8]/solidMolarMass_[4].W())*0.3*32.; //G(CH3OH)
              dMassSOLID[id_GC2H4] = dMassSOLID[id_GC2H4] - (dmassSolid[8]/solidMolarMass_[4].W())*0.275*28.; //G(C2H4)
              dMassSOLID[id_GCO2] = dMassSOLID[id_GCO2] - (dmassSolid[8]/solidMolarMass_[4].W())*0.4*44.; //G(CO2)
              dMassSOLID[id_GCOH2] = dMassSOLID[id_GCOH2] - (dmassSolid[8]/solidMolarMass_[4].W())*0.925*30.; //G(COH2)
              dMassSOLID[id_C] = dMassSOLID[id_C] - (dmassSolid[8]/solidMolarMass_[4].W())*1.*12.; //C

              //LIGC=>0.35 LIGCC + 0.1 COUMARYL + 0.08 PHENOL + 0.41 C2H4 + 1.0H2O + 0.7 G{COH2} + + 0.3 CH2O + 0.32 CO + 0.495 G{CH4} + 5.735 CHAR
              dmassSolid[9] = min(dt*kappa[9], 1.)*massKg*YSolidEff[id_LIGC];
              dMassSOLID[id_LIGC] = dMassSOLID[id_LIGC] + dmassSolid[9]; //LIGC
              dMassSOLID[id_LIGCC] = dMassSOLID[id_LIGCC] - (dmassSolid[9]/solidMolarMass_[5].W())*0.35*258.; //LIGCC
              dMassDV[volatileToGasMap_[20]] = dMassDV[volatileToGasMap_[20]] + (dmassSolid[9]/solidMolarMass_[5].W())*0.1*150.; // C9H10O2 COUMARYL
              dMassDV[volatileToGasMap_[18]] = dMassDV[volatileToGasMap_[18]] + (dmassSolid[9]/solidMolarMass_[5].W())*0.08*94.; // C6H5OH PHENOL
              dMassDV[volatileToGasMap_[9]] = dMassDV[volatileToGasMap_[9]] + (dmassSolid[9]/solidMolarMass_[5].W())*0.41*28.; // C2H4
              dMassDV[volatileToGasMap_[1]] = dMassDV[volatileToGasMap_[1]] + (dmassSolid[9]/solidMolarMass_[5].W())*1.0*18.; // H2O
              dMassSOLID[id_GCOH2] = dMassSOLID[id_GCOH2] - (dmassSolid[9]/solidMolarMass_[5].W())*0.7*30.; //G(COH2)
              dMassDV[volatileToGasMap_[4]] = dMassDV[volatileToGasMap_[4]] + (dmassSolid[9]/solidMolarMass_[5].W())*0.3*30.; // CH2O
              dMassDV[volatileToGasMap_[2]] = dMassDV[volatileToGasMap_[2]] + (dmassSolid[9]/solidMolarMass_[5].W())*0.32*28.; // CO
              dMassSOLID[id_GCH4] = dMassSOLID[id_GCH4] - (dmassSolid[9]/solidMolarMass_[5].W())*0.495*16.; //G(CH4)
              dMassSOLID[id_C] = dMassSOLID[id_C] - (dmassSolid[9]/solidMolarMass_[5].W())*5.735*12.; //C

              //LIGH=>LIGOH + 0.5 C2H5CHO + 0.5 C2H4 + 0.25 C2H4O2
              dmassSolid[10] = min(dt*kappa[10], 1.)*massKg*YSolidEff[id_LIGH];
              dMassSOLID[id_LIGH] = dMassSOLID[id_LIGH] + dmassSolid[10]; //LIGH
              dMassSOLID[id_LIGOH] = dMassSOLID[id_LIGOH]-(dmassSolid[10]/solidMolarMass_[6].W())*1.*378.; //LIGOH
              dMassDV[volatileToGasMap_[14]] = dMassDV[volatileToGasMap_[14]] + (dmassSolid[10]/solidMolarMass_[6].W())*0.5*58.; // C2H5CHO
              dMassDV[volatileToGasMap_[9]] = dMassDV[volatileToGasMap_[9]] + (dmassSolid[10]/solidMolarMass_[6].W())*0.5*28.; // C2H4
              dMassDV[volatileToGasMap_[11]] = dMassDV[volatileToGasMap_[11]] + (dmassSolid[10]/solidMolarMass_[6].W())*0.25*60.; // C2H4O2

              //LIGO=>LIGOH+CO2                                                                                                   1.00000E+009     .000  25500.00
              dmassSolid[11] = min(dt*kappa[11], 1.)*massKg*YSolidEff[id_LIGO];
              dMassSOLID[id_LIGO] = dMassSOLID[id_LIGO] + dmassSolid[11]; //LIGO
              dMassSOLID[id_LIGOH] = dMassSOLID[id_LIGOH] - (dmassSolid[11]/solidMolarMass_[7].W())*1.*378.; //LIGOH
              dMassDV[volatileToGasMap_[3]] = dMassDV[volatileToGasMap_[3]] + (dmassSolid[11]/solidMolarMass_[7].W())*1.*44.; // CO2

              //LIGCC=>0.3 COUMARYL + 0.2 PHENOL + 0.35 C2H4O2 + 0.7 H2O + 0.65 G{CH4} + 0.6 G{C2H4} + G{COH2} + 0.4 CO + 0.4 G{CO} + 6.75 CHAR
              dmassSolid[12] = min(dt*kappa[12], 1.)*massKg*YSolidEff[id_LIGCC];
              dMassSOLID[id_LIGCC] = dMassSOLID[id_LIGCC] + dmassSolid[12]; //LIGCC
              dMassDV[volatileToGasMap_[20]] = dMassDV[volatileToGasMap_[20]] + (dmassSolid[12]/solidMolarMass_[8].W())*0.3*150.; // C9H10O2 COUMARYL
              dMassDV[volatileToGasMap_[18]] = dMassDV[volatileToGasMap_[18]] + (dmassSolid[12]/solidMolarMass_[8].W())*0.2*94.; // C6H5OH PHENOL
              dMassDV[volatileToGasMap_[11]] = dMassDV[volatileToGasMap_[11]] + (dmassSolid[12]/solidMolarMass_[8].W())*0.35*60.; // C2H4O2
              dMassDV[volatileToGasMap_[1]] = dMassDV[volatileToGasMap_[1]] + (dmassSolid[12]/solidMolarMass_[8].W())*0.7*18.; // H2O
              dMassSOLID[id_GCH4] = dMassSOLID[id_GCH4] - (dmassSolid[12]/solidMolarMass_[8].W())*0.65*16.; //G(CH4)
              dMassSOLID[id_GC2H4] = dMassSOLID[id_GC2H4] - (dmassSolid[12]/solidMolarMass_[8].W())*0.6*28.; //G(C2H4)
              dMassSOLID[id_GCOH2] = dMassSOLID[id_GCOH2] - (dmassSolid[12]/solidMolarMass_[8].W())*1.0*30.; //G(COH2)
              dMassDV[volatileToGasMap_[2]] = dMassDV[volatileToGasMap_[2]] + (dmassSolid[12]/solidMolarMass_[8].W())*0.4*28.; // CO
              dMassSOLID[id_GCO] = dMassSOLID[id_GCO] - (dmassSolid[12]/solidMolarMass_[8].W())*0.4*28.; //G(CO)
              dMassSOLID[id_C] = dMassSOLID[id_C] - (dmassSolid[12]/solidMolarMass_[8].W())*6.75*12.; //C

              //LIGOH=>LIG + 0.9 H2O + 0.1 CH4 + 0.6 CH3OH + 0.1 G{H2} + 0.3 G{CH3OH} + 0.05 CO2 + 0.55 CO + 0.6 G{CO} + 0.05 HCOOH + 0.85 G{COH2} + 0.35 G{CH4} + 0.2 G{C2H4} + 4.15 CHAR
              dmassSolid[13] = min(dt*kappa[13], 1.)*massKg*YSolidEff[id_LIGOH];
              dMassSOLID[id_LIGOH] = dMassSOLID[id_LIGOH] + dmassSolid[13]; //LIGOH
              dMassSOLID[id_LIG] = dMassSOLID[id_LIG] - (dmassSolid[13]/solidMolarMass_[9].W())*1.*208.; //LIG
              dMassDV[volatileToGasMap_[1]] = dMassDV[volatileToGasMap_[1]] + (dmassSolid[13]/solidMolarMass_[9].W())*0.9*18.; // H2O
              dMassDV[volatileToGasMap_[6]] = dMassDV[volatileToGasMap_[6]] + (dmassSolid[13]/solidMolarMass_[9].W())*0.1*16.; // CH4
              dMassDV[volatileToGasMap_[7]] = dMassDV[volatileToGasMap_[7]] + (dmassSolid[13]/solidMolarMass_[9].W())*0.6*32.; // CH3OH
              dMassSOLID[id_GH2] = dMassSOLID[id_GH2] - (dmassSolid[13]/solidMolarMass_[9].W())*0.1*2.; //G(H2)
              dMassSOLID[id_GCH3OH] = dMassSOLID[id_GCH3OH] - (dmassSolid[13]/solidMolarMass_[9].W())*0.3*32.; //G(CH3OH)
              dMassDV[volatileToGasMap_[3]] = dMassDV[volatileToGasMap_[3]] + (dmassSolid[13]/solidMolarMass_[9].W())*0.05*44.; // CO2
              dMassDV[volatileToGasMap_[2]] = dMassDV[volatileToGasMap_[2]] + (dmassSolid[13]/solidMolarMass_[9].W())*0.55*28.; // CO
              dMassSOLID[id_GCO] = dMassSOLID[id_GCO] - (dmassSolid[13]/solidMolarMass_[9].W())*0.6*28.; //G(CO)
              dMassDV[volatileToGasMap_[5]] = dMassDV[volatileToGasMap_[5]] + (dmassSolid[13]/solidMolarMass_[9].W())*0.05*46.; // HCOOH
              dMassSOLID[id_GCOH2] = dMassSOLID[id_GCOH2] - (dmassSolid[13]/solidMolarMass_[9].W())*0.85*30.; //G(COH2)
              dMassSOLID[id_GCH4] = dMassSOLID[id_GCH4] - (dmassSolid[13]/solidMolarMass_[9].W())*0.35*16.; //G(CH4)
              dMassSOLID[id_GC2H4] = dMassSOLID[id_GC2H4] - (dmassSolid[13]/solidMolarMass_[9].W())*0.2*28.; //G(C2H4)
              dMassSOLID[id_C] = dMassSOLID[id_C] - (dmassSolid[13]/solidMolarMass_[9].W())*4.15*12.; //C


              //LIG=>0.7 C11H12O4 + 0.3 ANISOLE + 0.3 CO + 0.3 G{CO} + 0.3 CH3CHO
              dmassSolid[14] = dt*kappa[14]*massKg*YSolidEff[id_LIG];
              //LIG=>0.95 H2O + 0.2 CH2O + 0.4 CH3OH + CO + 0.2 CH4 + 0.05 HCOOH + 0.45 G{CO} + 0.5 G{COH2} + 0.4 G{CH4} + 0.65 G{C2H4} + 0.2 CH3CHO + 0.2 C2H5CHO + 5.5 CHAR
              dmassSolid[15] = dt*kappa[15]*massKg*YSolidEff[id_LIG];
              //LIG=>0.6 H2O + 0.4 CO + 0.2 CH4 + 0.4 CH2O + 0.2 G{CO} + 0.4 G{CH4} + 0.5 G{C2H4} + 0.4 G{CH3OH} + 2 G{COH2} + 6 CHAR
              dmassSolid[16] = dt*kappa[16]*massKg*YSolidEff[id_LIG];
              if (dmassSolid[14]+dmassSolid[15]+dmassSolid[16]>massKg*YSolidEff[id_LIG] && dmassSolid[14]+dmassSolid[15]+dmassSolid[16]>0.)
              {
                dmassTemp = dmassSolid[14]+dmassSolid[15]+dmassSolid[16];
                dmassSolid[14] = massKg*YSolidEff[id_LIG]*dmassSolid[14]/dmassTemp;
                dmassSolid[15] = massKg*YSolidEff[id_LIG]*dmassSolid[15]/dmassTemp;
                dmassSolid[16] = massKg*YSolidEff[id_LIG]*dmassSolid[16]/dmassTemp;
              }
              dMassSOLID[id_LIG] = dMassSOLID[id_LIG] + dmassSolid[14] + dmassSolid[15] + dmassSolid[16]; //LIG
              dMassDV[volatileToGasMap_[21]] = dMassDV[volatileToGasMap_[21]] + (dmassSolid[14]/solidMolarMass_[10].W())*0.7*208.; // C11H12O4
              dMassDV[volatileToGasMap_[19]] = dMassDV[volatileToGasMap_[19]] + (dmassSolid[14]/solidMolarMass_[10].W())*0.3*108.; // ANISOLE C6H5OCH3
              dMassDV[volatileToGasMap_[2]] = dMassDV[volatileToGasMap_[2]] + (dmassSolid[14]/solidMolarMass_[10].W())*0.3*28.; // CO
              dMassSOLID[id_GCO] = dMassSOLID[id_GCO] - (dmassSolid[14]/solidMolarMass_[10].W())*0.3*28.; //G(CO)
              dMassDV[volatileToGasMap_[10]] = dMassDV[volatileToGasMap_[10]] + (dmassSolid[14]/solidMolarMass_[10].W())*0.3*44.; // CH3CHO

              dMassDV[volatileToGasMap_[1]] = dMassDV[volatileToGasMap_[1]] + (dmassSolid[15]/solidMolarMass_[10].W())*0.95*18.; // H2O
              dMassDV[volatileToGasMap_[4]] = dMassDV[volatileToGasMap_[4]] + (dmassSolid[15]/solidMolarMass_[10].W())*0.2*30.; // CH2O
              dMassDV[volatileToGasMap_[7]] = dMassDV[volatileToGasMap_[7]] + (dmassSolid[15]/solidMolarMass_[10].W())*0.4*32.; // CH3OH
              dMassDV[volatileToGasMap_[2]] = dMassDV[volatileToGasMap_[2]] + (dmassSolid[15]/solidMolarMass_[10].W())*1.0*28.; // CO
              dMassDV[volatileToGasMap_[6]] = dMassDV[volatileToGasMap_[6]] + (dmassSolid[15]/solidMolarMass_[10].W())*0.2*16.; // CH4
              dMassDV[volatileToGasMap_[5]] = dMassDV[volatileToGasMap_[5]] + (dmassSolid[15]/solidMolarMass_[10].W())*0.05*46.; // HCOOH
              dMassSOLID[id_GCO] = dMassSOLID[id_GCO] - (dmassSolid[15]/solidMolarMass_[10].W())*0.45*28.; //G(CO)
              dMassSOLID[id_GCOH2] = dMassSOLID[id_GCOH2] - (dmassSolid[15]/solidMolarMass_[10].W())*0.5*30.; //G(COH2)
              dMassSOLID[id_GCH4] = dMassSOLID[id_GCH4] - (dmassSolid[15]/solidMolarMass_[10].W())*0.4*16.; //G(CH4)
              dMassSOLID[id_GC2H4] = dMassSOLID[id_GC2H4] - (dmassSolid[15]/solidMolarMass_[10].W())*0.65*28.; //G(C2H4)
              dMassDV[volatileToGasMap_[10]] = dMassDV[volatileToGasMap_[10]] + (dmassSolid[15]/solidMolarMass_[10].W())*0.2*44.; // CH3CHO
              dMassDV[volatileToGasMap_[14]] = dMassDV[volatileToGasMap_[14]] + (dmassSolid[15]/solidMolarMass_[10].W())*0.2*58.; // C2H5CHO
              dMassSOLID[id_C] = dMassSOLID[id_C] - (dmassSolid[15]/solidMolarMass_[10].W())*5.5*12.; //C

              dMassDV[volatileToGasMap_[1]] = dMassDV[volatileToGasMap_[1]] + (dmassSolid[16]/solidMolarMass_[10].W())*0.6*18.; // H2O
              dMassDV[volatileToGasMap_[2]] = dMassDV[volatileToGasMap_[2]] + (dmassSolid[16]/solidMolarMass_[10].W())*0.4*28.; // CO
              dMassDV[volatileToGasMap_[6]] = dMassDV[volatileToGasMap_[6]] + (dmassSolid[16]/solidMolarMass_[10].W())*0.2*16.; // CH4
              dMassDV[volatileToGasMap_[4]] = dMassDV[volatileToGasMap_[4]] + (dmassSolid[16]/solidMolarMass_[10].W())*0.4*30.; // CH2O
              dMassSOLID[id_GCO] = dMassSOLID[id_GCO] - (dmassSolid[16]/solidMolarMass_[10].W())*0.2*28.; //G(CO)
              dMassSOLID[id_GCH4] = dMassSOLID[id_GCH4] - (dmassSolid[16]/solidMolarMass_[10].W())*0.4*16.; //G(CH4)
              dMassSOLID[id_GC2H4] = dMassSOLID[id_GC2H4] - (dmassSolid[16]/solidMolarMass_[10].W())*0.5*28.; //G(C2H4)
              dMassSOLID[id_GCH3OH] = dMassSOLID[id_GCH3OH] - (dmassSolid[16]/solidMolarMass_[10].W())*0.4*32.; //G(CH3OH)
              dMassSOLID[id_GCOH2] = dMassSOLID[id_GCOH2] - (dmassSolid[16]/solidMolarMass_[10].W())*2.0*30.; //G(COH2)
              dMassSOLID[id_C] = dMassSOLID[id_C] - (dmassSolid[16]/solidMolarMass_[10].W())*6.0*12.; //C

              //TGL=>ACROLEIN + 0.5U2ME12+2.5MLINO !!!! New !!!!
              dmassSolid[17] = min(dt*kappa[17], 1.)*massKg*YSolidEff[id_TGL];
              dMassSOLID[id_TGL] = dMassSOLID[id_TGL] + dmassSolid[17]; //TGL
              dMassDV[volatileToGasMap_[13]] = dMassDV[volatileToGasMap_[13]] + (dmassSolid[17]/solidMolarMass_[11].W())*1.0*56.; // ACROLEIN // C2H3CHO
              dMassDV[volatileToGasMap_[22]] = dMassDV[volatileToGasMap_[22]] + (dmassSolid[17]/solidMolarMass_[11].W())*0.5*210.; // U2ME12 //C13H22O2
              dMassDV[volatileToGasMap_[23]] = dMassDV[volatileToGasMap_[23]] + (dmassSolid[17]/solidMolarMass_[11].W())*2.5*294.; // MLINO //C19H34O2
              //NOTE Changed from C19H34O2, Linoleic to 0.5U2ME12+2.5MLINO

              //CTANN=>PHENOL + ITANN
              dmassSolid[18] = min(dt*kappa[18], 1.)*massKg*YSolidEff[id_CTANN];
              dMassSOLID[id_CTANN] = dMassSOLID[id_CTANN] + dmassSolid[18]; //CTANN
              dMassDV[volatileToGasMap_[18]] = dMassDV[volatileToGasMap_[18]] + (dmassSolid[18]/solidMolarMass_[12].W())*1.0*94.; // C6H5OH PHENOL
              dMassSOLID[id_ITANN] = dMassSOLID[id_ITANN] - (dmassSolid[18]/solidMolarMass_[12].W())*1.0*210.; //ITANN

              //ITANN=>6 CHAR + 3 CO + 3 H2O
              dmassSolid[19] = min(dt*kappa[19], 1.)*massKg*YSolidEff[id_ITANN];
              dMassSOLID[id_ITANN] = dMassSOLID[id_ITANN] + dmassSolid[19]; //ITANN
              dMassSOLID[id_C] = dMassSOLID[id_C] - (dmassSolid[19]/solidMolarMass_[13].W())*6.0*12.; //C
              dMassDV[volatileToGasMap_[2]] = dMassDV[volatileToGasMap_[2]] + (dmassSolid[19]/solidMolarMass_[13].W())*3.0*28.; // CO
              dMassDV[volatileToGasMap_[1]] = dMassDV[volatileToGasMap_[1]] + (dmassSolid[19]/solidMolarMass_[13].W())*3.0*18.; // H2O

              //G(CO2)=>CO2
              dmassSolid[20] = min(dt*kappa[20], 1.)*massKg*YSolidEff[id_GCO2];
              dMassSOLID[id_GCO2] = dMassSOLID[id_GCO2] + dmassSolid[20]; //G(CO2)
              dMassDV[volatileToGasMap_[3]] = dMassDV[volatileToGasMap_[3]] + dmassSolid[20]; // CO2

              //G(CO)=>CO
              dmassSolid[21] = min(dt*kappa[21], 1.)*massKg*YSolidEff[id_GCO];
              dMassSOLID[id_GCO] = dMassSOLID[id_GCO] + dmassSolid[21]; //G(CO)
              dMassDV[volatileToGasMap_[2]] = dMassDV[volatileToGasMap_[2]] + dmassSolid[21]; // CO

              //G(COH2)=>CO+H2
              dmassSolid[22] = min(dt*kappa[22], 1.)*massKg*YSolidEff[id_GCOH2];
              dMassSOLID[id_GCOH2] = dMassSOLID[id_GCOH2] + dmassSolid[22]; //G(COH2)
              dMassDV[volatileToGasMap_[2]] = dMassDV[volatileToGasMap_[2]] + (dmassSolid[22]/solidMolarMass_[16].W())*1.*28.; // CO
              dMassDV[volatileToGasMap_[0]] = dMassDV[volatileToGasMap_[0]]+ (dmassSolid[22]/solidMolarMass_[16].W())*1.*2; // H2

              //G(H2)=>H2
              dmassSolid[23] = min(dt*kappa[23], 1.)*massKg*YSolidEff[id_GH2];
              dMassSOLID[id_GH2] = dMassSOLID[id_GH2] + dmassSolid[23]; //G(H2)
              dMassDV[volatileToGasMap_[0]] = dMassDV[volatileToGasMap_[0]] + dmassSolid[23]; // H2

              //G(CH4)=>CH4
              dmassSolid[24] = min(dt*kappa[24], 1.)*massKg*YSolidEff[id_GCH4];
              dMassSOLID[id_GCH4] = dMassSOLID[id_GCH4] + dmassSolid[24]; //G(CH4)
              dMassDV[volatileToGasMap_[6]] = dMassDV[volatileToGasMap_[6]] + dmassSolid[24]; // CH4

              //G(CH3OH)=>CH3OH
              dmassSolid[25] = min(dt*kappa[25], 1.)*massKg*YSolidEff[id_GCH3OH];
              dMassSOLID[id_GCH3OH] = dMassSOLID[id_GCH3OH] + dmassSolid[25]; //G(CH3OH)
              dMassDV[volatileToGasMap_[7]] = dMassDV[volatileToGasMap_[7]] + dmassSolid[25]; // CH3OH

              //G(C2H4)=>C2H4
              dmassSolid[26] = min(dt*kappa[26], 1.)*massKg*YSolidEff[id_GC2H4];
              dMassSOLID[id_GC2H4] = dMassSOLID[id_GC2H4] + dmassSolid[26]; //G(C2H4)
              dMassDV[volatileToGasMap_[9]] = dMassDV[volatileToGasMap_[9]] + dmassSolid[26]; // C2H4

              forAll(dMassSOLID,i)
              {
                dMassSOLID[i] = dMassSOLID[i]/1000;
              }

              forAll(dMassDV,i)
              {
                dMassDV[i] = dMassDV[i]/1000;
              }

              //     Info<<"2 YSolidEff "<<YSolidEff<<endl;
              //     Info<<"2 dMassSOLID: "<<dMassSOLID<<endl;
              //     Info<<"2 dMassDV: "<<dMassDV<<endl;

              //     forAll(volatileData_, i)
              //     {
                //         const label id = volatileToGasMap_[i];
                //         const scalar massVolatile0 = mass0*YVolatile0_[i];
                //         const scalar massVolatile = mass*YGasEff[id];
                //
                //         // Combustion allowed once all volatile components evolved
                //         done = done && (massVolatile <= residualCoeff_*massVolatile0);
                //
                //         // Model coefficients
                //         const scalar A1 = volatileData_[i].A1();
                //         const scalar E = volatileData_[i].E();
                //
                //         // Kinetic rate
                //         const scalar kappa = A1*exp(-E/(RR*T));
                //
                //         // Mass transferred from particle to carrier gas phase
                //         dMassDV[id] = min(dt*kappa*massVolatile, massVolatile);
                //     }
                //
                //     if (done && canCombust != -1)
                //     {
                  //         canCombust = 1;
                  //     }
                }


                // ************************************************************************* //

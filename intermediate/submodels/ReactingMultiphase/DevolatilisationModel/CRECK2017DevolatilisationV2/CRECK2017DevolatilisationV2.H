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
    Devolatilisation model from CRECK.

    Originally in
    Extractives Extend the Applicability of Multistep Kinetic Scheme of Biomass Pyrolysis
    Paulo Eduardo Amaral Debiagi, Chiara Pecchi, Giancarlo Gentile, Alessio Frassoldati, Alberto Cuoci, Tiziano Faravelli, and Eliseo Ranzi
    Energy & Fuels 2015 29 (10), 6544-6555
    DOI: 10.1021/acs.energyfuels.5b01753

    Later version
    Paulo Eduardo Amaral Debiagi, Giancarlo Gentile, Matteo Pelucchi, Alessio Frassoldati, Alberto Cuoci, Tiziano Faravelli, Eliseo Ranzi,
    Detailed kinetic mechanism of gas-phase reactions of volatiles released from biomass pyrolysis,
    Biomass and Bioenergy,
    Volume 93,
    2016,
    Pages 60-71,
    ISSN 0961-9534,
    https://doi.org/10.1016/j.biombioe.2016.06.015.
    (http://www.sciencedirect.com/science/article/pii/S0961953416302173)

    2015 and 2016 are slightly different!!!!

    Pyrolysis reactions	Kinetic Parameters
(A (s-1), Eact (kcal/kmol))	delt_Hr (700 K) (kcal/kmol)
            n 10^n !!!!
Cellulose
CELL	->	CELLA	1.5 1014·exp(-47,000/RT)	-1300
CELLA	->	0.45 C2H4O2 + 0.2 C2H2O2 + 0.1 CH3CHO + 0.25 HMFU + 0.3 C2H5CHO + 0.15 CH3OH + 0.4 CH2O + 0.31 CO + 0.41 CO2 + 0.05 H2 + 0.83 H2O + 0.02 HCOOH + 0.2 G{CH4} + 0.05 G{H2} + 0.61 CHAR	2 106·exp(-19,100/RT)	27,100
CELLA	->	LVG	4·T·exp(-10,000/RT)	23,200
CELL	->	5H2O + 6CHAR	6.5 107·exp(-31,000/RT)	-62,700
Hemicellulose
HECELL	->	0.58 HCE1 + 0.42 HCE2	1 1010·exp(-31,000/RT)	-5000
HCE1	->	0.025 H2O + 0.5 CO2 + 0.025 HCOOH + 0.5 CO + 0.8 CH2O + 0.125 C2H5OH + 0.1 CH3OH + 0.25 C2H4 + 0.125 G{H2} + 0.275 G{CO2} + 0.4 G{COH2} + 0.45 G{CH3OH} + 0.325 G{CH4} + 0.875 CHAR	1.2 109·exp(-30,000/RT)	-500
HCE1	->	0.25 H2O + 0.8 CO2 + 0.05 HCOOH + 0.1 CO + 0.15 G{CO} + 0.15 G{CO2} + 0.2 G{H2} + 0.3 CH2O + 1.2 G{COH2} + 0.625 G{CH4} + 0.375 G{C2H4} + 0.875 CHAR	1.5 10-1·T·exp(-8000/RT)	-42,400
HCE1	->	XYLAN	3·T·exp(-11,000/RT)	17,900
HCE2	->	0.2 H2O + 0.175 CO + 0.275 CO2 + 0.5 CH2O + 0.1 C2H5OH + 0.2 C2H4O2 + 0.025 HCOOH + 0.25 G{CH4} + 0.3 G{CH3OH} + 0.275 G{C2H4} + 0.4 G{CO2} + 0.925 G{COH2} + CHAR	5 109·exp(-33,000/RT)	12,000
Lignins
LIGC	->	0.35 LIGCC + 0.1 COUMARYL + 0.08 PHENOL + 0.41 C2H4 + 1.0H2O + 0.7 G{COH2} +  + 0.3 CH2O + 0.32 CO + 0.495 G{CH4} + 5.735 CHAR	1.33 1015 exp(-48,500/RT)	-10,300
LIGH	->	LIGOH + 0.5 C2H5CHO + 0.5 C2H4 + 0.25C2H4O2	6.7 1012·exp(-37,500/RT)	30,700
LIGO	->	LIGOH + CO2	3.3 108·exp(-25,500/RT)	26,000
LIGCC	->	0.3 COUMARYL + 0.2 PHENOL + 0.35 C2H4O2 + 0.7 H2O + 0.65 G{CH4} + 0.6 G{C2H4} + G{COH2} + 0.4 CO + 0.4 G{CO} + 6.75 CHAR	1.67 106·exp(-31,500/RT)	-31,100
LIGOH	->	LIG + 0.9 H2O + 0.1 CH4 + 0.6 CH3OH + 0.1 G{H2} + 0.3 G{CH3OH} + 0.05 CO2 + 0.55 CO + 0.6 G{CO} + 0.05 HCOOH + 0.85 G{COH2} + 0.35 G{CH4} + 0.2 G{C2H4} + 4.15 CHAR	1 108·exp(-30,000/RT)	-26,100
LIG	->	0.7 C11H12O4 + 0.3 ANISOLE + 0.3 CO + 0.3 G{CO} + 0.3 CH3CHO	4·T·exp(-12,000/RT)	46,200
LIG	->	0.95 H2O + 0.2 CH2O + 0.4 CH3OH + CO + 0.2 CH4 + 0.05 HCOOH + 0.45 G{CO} + 0.5 G{COH2} + 0.4 G{CH4} + 0.65 G{C2H4} + 0.2 CH3CHO + 0.2 C2H5CHO + 5.5 CHAR	4 108·exp(-30,000/RT)	-21,100
LIG	->	0.6 H2O + 0.4 CO + 0.2 CH4 + 0.4 CH2O + 0.2 G{CO} + 0.4 G{CH4} + 0.5 G{C2H4} + 0.4 G{CH3OH} + 2 G{COH2} + 6 CHAR	8.3 10-2·T·exp(-8,000/RT)	-83,600
Extractives
!!!!!!!!!TGL	->	ACROLEIN + 3 LINOLEIC ACID	7 1012·exp(-45,700/RT)	1300// change to  0.5U2ME12+2.5MLINO
A predictive model of biochar formation and characterization P. Debiagi, G. Gentile, A. Cuoci, A. Frassoldati, E. Ranzi, T. Faravelli
CTANN	->	PHENOL + ITANN	5 101·exp(-11,000/RT)	1300
ITANN	->	6 CHAR + 3 CO + 3 H2O	1.5 10-2·exp(-6100/RT)	10,100
Metaplastic
G{CO2}	->	CO2	1 106·exp(-24,000/RT)	-29,100
G{CO}	->	CO	5 1012·exp(-50,000/RT)	-13,400
G{COH2}	->	CO+H2	5 1011·exp(-71,000/RT)	48,600
G{H2}	->	H2	5 1011·exp(-75,000/RT)	0
G{CH4}	->	CH4	5 1012·exp(-71,667/RT)	0
G{CH3OH}	->	CH3OH	2 1012·exp(-50,000/RT)	0
G{C2H4}	->	C2H4	5 1012·exp(-71,667/RT)	0




\*---------------------------------------------------------------------------*/

#ifndef CRECK2017DevolatilisationV2_H
#define CRECK2017DevolatilisationV2_H

#include "DevolatilisationModel.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //
#define declareSpecieId(speciesName) \
    const label id_##speciesName = this->owner().composition().localId(idSolid, #speciesName);

namespace Foam
{
/*---------------------------------------------------------------------------*\
                     Class CRECK2017DevolatilisationV2 Declaration
\*---------------------------------------------------------------------------*/

template<class CloudType>
class CRECK2017DevolatilisationV2
:
    public DevolatilisationModel<CloudType>
{
    class volatileData
    {
        // Private data

            //- Specie name
            word name_;


    public:

        // Constructors

            //- Null constructor
            volatileData()
            :
                name_(word::null)
            {}

            //- Construct from Istream
            volatileData(Istream& is)
            :
                name_(is)
            {}

            //- Construct as copy
            // volatileData(const volatileData& vd)
            // :
            //     name_(vd.name_)
            // {}


        //- Destructor
        ~volatileData()
        {}


        // Public Member Functions

            // Access

                //- Return const access to the name
                const word& name() const
                {
                    return name_;
                }

            //- Read from Istream
            friend Istream& operator>>(Istream& is, volatileData& vd)
            {
                is.readBeginList("volatileData");
                is  >>  vd.name_;
                is.readEndList("volatileData");

                return is;
            }

            //- Write to Ostream
            friend Ostream& operator<<(Ostream& os, const volatileData& vd)
            {
                os  << token::BEGIN_LIST
                    << vd.name_
                    << token::END_LIST;

                return os;
            }
    };

    class solidMolarMass
    {
        // Private data

            //- Specie name
            word name_;

            //- Molar mass g/mol
            scalar W_;


    public:

        // Constructors

            //- Null constructor
            solidMolarMass()
            :
                name_(word::null),
                W_(0.0)
            {}

            //- Construct from Istream
            solidMolarMass(Istream& is)
            :
                name_(is),
                W_(is)
            {}

            // //- Construct as copy
            // solidMolarMass(const solidMolarMass& vd)
            // :
            //     name_(vd.name_),
            //     W_(vd.W_)
            // {}


        //- Destructor
        ~solidMolarMass()
        {}


        // Public Member Functions

            // Access

                //- Return const access to the name
                const word& name() const
                {
                    return name_;
                }

                //- Return const access to the molar mass
                scalar W() const
                {
                    return W_;
                }

            //- Read from Istream
            friend Istream& operator>>(Istream& is, solidMolarMass& vd)
            {
                is.readBeginList("solidMolarMass");
                is  >>  vd.name_ >> vd.W_;
                is.readEndList("solidMolarMass");

                return is;
            }

            //- Write to Ostream
            friend Ostream& operator<<(Ostream& os, const solidMolarMass& vd)
            {
                os  << token::BEGIN_LIST
                    << vd.name_ << token::SPACE
                    << vd.W_
                    << token::END_LIST;

                return os;
            }
    };
    // Private data

        // Model constants

            //- List of volatile data - (name)
            List<volatileData> volatileData_;

            //- List of molar mass of solid - (name, molar mass)
            List<solidMolarMass> solidMolarMass_;


            //- Mapping between local and cloud gaseous species
            List<label> volatileToGasMap_;

            //- Volatile residual coefficient (0-1)
            //  When the fraction of volatiles are depleted below this
            //  threshold, combustion can occur
            const scalar residualCoeff_;


public:

    //- Runtime type information
    TypeName("CRECK2017DevolatilisationV2");


    // Constructors

        //- Construct from dictionary
        CRECK2017DevolatilisationV2(const dictionary&, CloudType& owner);

        //- Construct copy
        CRECK2017DevolatilisationV2(const CRECK2017DevolatilisationV2<CloudType>& dm);

        //- Construct and return a clone
        virtual autoPtr<DevolatilisationModel<CloudType>> clone() const
        {
            return autoPtr<DevolatilisationModel<CloudType>>
            (
                new CRECK2017DevolatilisationV2<CloudType>(*this)
            );
        }


    //- Destructor
    virtual ~CRECK2017DevolatilisationV2();


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
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#ifdef NoRepository
    #include "CRECK2017DevolatilisationV2.C"
#endif

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

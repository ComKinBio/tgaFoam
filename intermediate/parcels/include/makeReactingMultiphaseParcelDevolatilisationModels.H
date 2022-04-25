/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2011-2018 OpenFOAM Foundation
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

#ifndef makeReactingMultiphaseParcelDevolatilisationModels_H
#define makeReactingMultiphaseParcelDevolatilisationModels_H

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#include "ConstantRateDevolatilisation.H"
#include "NoDevolatilisation.H"
#include "SingleKineticRateDevolatilisation.H"
#include "DiBlasiDevolatilisation.H"
#include "CRECK2017DevolatilisationV2.H"
#include "CRECK2017DevolatilisationV3.H"
#include "CRECKDevolatilisationV4.H"
#include "CRECKDevolatilisationV5.H"
#include "ancaSecCharDevolatilisation.H"


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#define makeReactingMultiphaseParcelDevolatilisationModels(CloudType)          \
                                                                               \
    makeDevolatilisationModel(CloudType);                                      \
    makeDevolatilisationModelType(ConstantRateDevolatilisation, CloudType);    \
    makeDevolatilisationModelType(NoDevolatilisation, CloudType);              \
    makeDevolatilisationModelType                                              \
    (                                                                          \
        SingleKineticRateDevolatilisation,                                     \
        CloudType                                                              \
    );                                                                         \
    makeDevolatilisationModelType(DiBlasiDevolatilisation, CloudType);         \
    makeDevolatilisationModelType(CRECK2017DevolatilisationV2, CloudType);     \
    makeDevolatilisationModelType(CRECK2017DevolatilisationV3, CloudType);     \
    makeDevolatilisationModelType(CRECKDevolatilisationV4, CloudType);         \
    makeDevolatilisationModelType(CRECKDevolatilisationV5, CloudType);         \
    makeDevolatilisationModelType(ancaSecCharDevolatilisation, CloudType);

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

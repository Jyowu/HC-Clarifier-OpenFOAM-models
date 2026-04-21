/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     | Website:  https://openfoam.org
    \\  /    A nd           | Copyright (C) 2016-2023 OpenFOAM Foundation
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
    Foam::RASModels::kOmegaSST_BCS2

Description
    k-omega SST turbulence model using the shared
    BuoyantCurvatureSwirlTools2 helper for curvature, swirl and optional
    buoyancy modifications.

SourceFiles
    kOmegaSST_BCS2.C

\*---------------------------------------------------------------------------*/

#ifndef kOmegaSST_BCS2_H
#define kOmegaSST_BCS2_H

#include "kOmegaSST.H"
#include "BuoyantCurvatureSwirlTools2.H"

namespace Foam
{
namespace RASModels
{

template<class BasicMomentumTransportModel>
class kOmegaSST_BCS2
:
    public kOmegaSST<BasicMomentumTransportModel>
{
protected:

    BuoyantCurvatureSwirlTools2 curvatureSwirlTools_;

public:

    typedef typename BasicMomentumTransportModel::alphaField alphaField;
    typedef typename BasicMomentumTransportModel::rhoField rhoField;

    TypeName("kOmegaSST_BCS2");

    kOmegaSST_BCS2
    (
        const alphaField& alpha,
        const rhoField& rho,
        const volVectorField& U,
        const surfaceScalarField& alphaRhoPhi,
        const surfaceScalarField& phi,
        const viscosity& viscosity,
        const word& type = typeName
    );

    virtual ~kOmegaSST_BCS2()
    {}

    virtual bool read();

    virtual void correct();
};


} // End namespace RASModels
} // End namespace Foam

#ifdef NoRepository
    #include "RAS/kOmegaSST_BCS2/kOmegaSST_BCS2.C"
#endif

#endif

// ************************************************************************* //

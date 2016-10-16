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

Class
    Foam::incompressible::RASModels::VLESSpalartAllmaras

Description
    Spalart-Allmaras VLES model

    The default model coefficients correspond to the following:
    \verbatim
        VLESSpalartAllmarasCoeffs
        {
            curvatureCorrection false;
            Cb1                 0.1355;
            Cb2                 0.622;
            Cw2                 0.3;
            Cw3                 2.0;
            Cv1                 7.1;
            Cr1                 1.0;
            Cr2                 12.0;
            Cr3                 1.0;
            sigmaNut            0.66666;
            kappa               0.41;
        }
    \endverbatim

Original author(s)
    Aleksandar Jemcov

Reprinted by Xiaoping Qiu

SourceFiles
    VLESSpalartAllmaras.C

\*---------------------------------------------------------------------------*/

#ifndef VLESSpalartAllmaras_H
#define VLESSpalartAllmaras_H

#include "turbulentTransportModel.H"
#include "eddyViscosity.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

/*---------------------------------------------------------------------------*\
                           Class VLESSpalartAllmaras Declaration
\*---------------------------------------------------------------------------*/

class VLESSpalartAllmaras
:
    public eddyViscosity<incompressible::RASModel>
{

    protected:


        Switch curvatureCorrection_;

        Switch delayed_;

        dimensionedScalar sigmaNut_;
        dimensionedScalar kappa_;

        dimensionedScalar Cb1_;
        dimensionedScalar Cb2_;
        dimensionedScalar Cw1_;
        dimensionedScalar Cw2_;
        dimensionedScalar Cw3_;
        dimensionedScalar Cv1_;
        dimensionedScalar Cr1_;
        dimensionedScalar Cr2_;
        dimensionedScalar Cr3_;
        dimensionedScalar Cx_;

        volScalarField nuTilda_;
        volScalarField fr1_;
        volScalarField Fr_;

	//- Wall distance
	//  Note: different to wall distance in parent RASModel
	//  which is for near-wall cells only
	const volScalarField& y_;


    // Protected Member Functions

        virtual void correctNut();

        tmp<volScalarField> chi() const;

        tmp<volScalarField> fv1(const volScalarField& chi) const;

        tmp<volScalarField> fv2
        (
            const volScalarField& chi,
            const volScalarField& fv1
        ) const;

        tmp<volScalarField> fw(const volScalarField& Stilda) const;

public:

    //- Runtime type information
    TypeName("VLESSpalartAllmaras");


    // Constructors

        //- Construct from components
        VLESSpalartAllmaras
        (
            const geometricOneField& alpha,
            const geometricOneField& rho,
            const volVectorField& U,
            const surfaceScalarField& alphaRhoPhi,
            const surfaceScalarField& phi,
            const transportModel& transport,
            const word& propertiesName = turbulenceModel::propertiesName,
            const word& type = typeName
        );


    //- Destructor
    virtual ~VLESSpalartAllmaras()
    {}


    // Member Functions

        tmp<volScalarField> F1() const;

        virtual tmp<volScalarField> nut() const
        {
            return nut_;
        }

        //- Return the effective diffusivity for nuTilda
        tmp<volScalarField> DnuTildaEff() const;

        //- Return the turbulence kinetic energy
        virtual tmp<volScalarField> k() const;

        //- Return the turbulence kinetic energy dissipation rate
        virtual tmp<volScalarField> epsilon() const;

        //- Return the effective diffusivity for k

        //- Validate the turbulence fields after construction
        //  Update turbulence viscosity and other derived fields as requires
        virtual void validate();

        //- Solve the turbulence equations and correct the turbulence viscosity
        //- Read RASProperties dictionary
        virtual bool read();

        virtual void correct();
};


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

#endif

// ************************************************************************* //

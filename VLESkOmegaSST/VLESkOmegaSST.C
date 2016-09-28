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

\*---------------------------------------------------------------------------*/

#include "VLESkOmegaSST.H"
#include "bound.H"
#include "wallDist.H"
#include "addToRunTimeSelectionTable.H"

// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

namespace Foam
{
namespace incompressible
{
namespace RASModels
{

// * * * * * * * * * * * * * * Static Data Members * * * * * * * * * * * * * //

defineTypeNameAndDebug(VLESkOmegaSST, 0);
addToRunTimeSelectionTable(RASModel, VLESkOmegaSST, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> VLESkOmegaSST::F1(const volScalarField& CDkOmega) const
{
    if (!delayed_)
    {
        return tmp<volScalarField>
        (
            new volScalarField   
            (
                IOobject
                (
                    "zero",
                    mesh_.time().timeName(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                dimensionedScalar("zero", dimless, 0),
                calculatedFvPatchScalarField::typeName
            )
        );
    }    

    else
    {
	tmp<volScalarField> CDkOmegaPlus = max
	(
	    CDkOmega,
	    dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
	);

	tmp<volScalarField> arg1 = min
	(
	    max
	    (
	        (scalar(1)/betaStar_)*sqrt(k_)/(omega_*y_),
	        scalar(500)*nu()/(sqr(y_)*omega_)
	    ),
	    (4*alphaOmega2_)*k_/(CDkOmegaPlus*sqr(y_))
	);

	return tanh(pow4(arg1));
    }
}


tmp<volScalarField> VLESkOmegaSST::F2() const
{
    tmp<volScalarField> arg2 = max
    (
        (scalar(2)/betaStar_)*sqrt(k_)/(omega_*y_),
        scalar(500)*nu()/(sqr(y_)*omega_)
    );

    return tanh(sqr(arg2));
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void VLESkOmegaSST::correctNut()
{
    // Currently this function is not implemented due to the complexity of
    // evaluating nut.  Better calculate nut at the end of correct()
    NotImplemented;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

VLESkOmegaSST::VLESkOmegaSST
(
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
)
:
    eddyViscosity<incompressible::RASModel>
    (
        type,
        alpha,
        rho,
        U,
        alphaRhoPhi,
        phi,
        transport,
        propertiesName
    ),

    delayed_
    (
        coeffDict_.lookupOrDefault<Switch>
        (
            "delayed", 
            true
        )
    ),
    alphaK1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK1",
            coeffDict_,
            0.85
        )
    ),
    alphaK2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaK2",
            coeffDict_,
            1.0
        )
    ),
    alphaOmega1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega1",
            coeffDict_,
            0.5
        )
    ),
    alphaOmega2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "alphaOmega2",
            coeffDict_,
            0.856
        )
    ),
    gamma1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma1",
            coeffDict_,
            0.5555
        )
    ),
    gamma2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "gamma2",
            coeffDict_,
            0.44
        )
    ),
    beta1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta1",
            coeffDict_,
            0.075
        )
    ),
    beta2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "beta2",
            coeffDict_,
            0.0828
        )
    ),
    betaStar_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "betaStar",
            coeffDict_,
            0.09
        )
    ),
    a1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "a1",
            coeffDict_,
            0.31
        )
    ),
    c1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "c1",
            coeffDict_,
            10.0
        )
    ),
    Cx_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cx",
            coeffDict_,
            0.61
        )
    ),

    k_
    (
        IOobject
        (
            "k",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    omega_
    (
        IOobject
        (
            "omega",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    nut_
    (
        IOobject
        (
            "nut",
            runTime_.timeName(),
            mesh_,
            IOobject::MUST_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_
    ),
    Fr_
    (
        IOobject
        (
            "Fr",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("fr", dimless, 1),
        calculatedFvPatchScalarField::typeName
    ),
    y_(wallDist::New(mesh_).y())
{
    bound(k_, kMin_);
    bound(omega_, omegaMin_);

    if (type == typeName)
    {
        nut_ =
        (
            a1_*k_
          / max
            (
                a1_*omega_,
                F2()*sqrt(2.0)*mag(symm(fvc::grad(U_)))
            )
        );

        nut_.correctBoundaryConditions();

        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool VLESkOmegaSST::read()
{
    if (eddyViscosity<incompressible::RASModel>::read())
    {
        delayed_.readIfPresent("delayed", coeffDict());
        alphaK1_.readIfPresent(coeffDict());
        alphaK2_.readIfPresent(coeffDict());
        alphaOmega1_.readIfPresent(coeffDict());
        alphaOmega2_.readIfPresent(coeffDict());
        gamma1_.readIfPresent(coeffDict());
        gamma2_.readIfPresent(coeffDict());
        beta1_.readIfPresent(coeffDict());
        beta2_.readIfPresent(coeffDict());
        betaStar_.readIfPresent(coeffDict());
        a1_.readIfPresent(coeffDict());
        c1_.readIfPresent(coeffDict());
        Cx_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void VLESkOmegaSST::validate()
{}


void VLESkOmegaSST::correct()
{
    eddyViscosity<incompressible::RASModel>::correct();

    if (!turbulence_)
    {
        return;
    }

    volScalarField Lc
    (
        IOobject
        (
            "Lc",
            mesh_.time().timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("Lc", dimLength, SMALL),
        calculatedFvPatchScalarField::typeName
    );

    label nD = mesh_.nGeometricD();

    if (nD == 3)
    {
        Lc.ref() = Cx_*pow(mesh_.V(), 1.0/3.0);
    }
    else if (nD == 2)
    {
        Vector<label> const& directions = mesh_.geometricD();

        scalar thickness = 0.0;
        for (direction dir=0; dir<directions.nComponents; dir++)
        {
            if (directions[dir] == -1)
            {
                thickness = mesh_.bounds().span()[dir];
                break;
            }
        }

        Lc.ref() = Cx_*sqrt(mesh_.V()/dimensionedScalar("thickness", dimLength, thickness));
    }
    else
    {
        FatalErrorIn("VLESKOmegaSST.C")
            << "Case is not 3D or 2D, VLES is not applicable"
            << exit(FatalError);
    }  

    tmp<volScalarField> const Li = pow(k_,3.0/2.0)/(betaStar_*k_*omega_);
    tmp<volScalarField> const Lk = 
        pow(nu(),3.0/4.0)/pow(betaStar_*k_*omega_,1.0/4.0);

    volScalarField const S2(2*magSqr(symm(fvc::grad(U_))));
    volScalarField const Omega(2*magSqr(skew(fvc::grad(U_))));
    volScalarField G(GName(), nut_*S2);

    // Update omega and G at the wall
    omega_.boundaryFieldRef().updateCoeffs();

    volScalarField const CDkOmega
    (
        (2*alphaOmega2_)*(fvc::grad(k_) & fvc::grad(omega_))/omega_
    );

    volScalarField const F1(this->F1(CDkOmega));

    dimensionedScalar nutMin
    (
        "nutMin", 
        kMin_.dimensions()/omegaMin_.dimensions(), 
        SMALL
     );

    // Turbulent frequency equation
    tmp<fvScalarMatrix> omegaEqn
    (
        fvm::ddt(omega_)
      + fvm::div(phi_, omega_)
      - fvm::laplacian(DomegaEff(F1), omega_)
     ==
        gamma(F1)* min(G, c1_*betaStar_*k_*omega_)/(nut_ + nutMin)
      - fvm::Sp(beta(F1)*omega_, omega_)
      + 2.0*(1.0-F1)*alphaOmega2_*(fvc::grad(k_)&fvc::grad(omega_))/omega_
    );

    omegaEqn.ref().relax();

    omegaEqn.ref().boundaryManipulate(omega_.boundaryFieldRef());
    
    //mesh_.updateFvMatrix(omegaEqn());
    solve(omegaEqn);
    bound(omega_, omegaMin_);

    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::laplacian(DkEff(F1), k_)
     ==
        min(G, c1_*betaStar_*k_*omega_)
      - fvm::Sp(betaStar_*omega_, k_)
    );

    kEqn.ref().relax();
    //mesh_.updateFvMatrix(kEqn());
    solve(kEqn);
    bound(k_, kMin_);

    // Re-calculate viscosity

    if (delayed_)
    {
        Fr_ = max
        (
            min
            (
                scalar(1.0),
                pow
                (
                    (scalar(1.0)-(1-F1)*exp(-0.002*Lc/Lk()))
                    /
                    (scalar(1.0)-(1-F1)*exp(-0.002*Li()/Lk()) + SMALL),
                    2.0
                )
            ),
            scalar(0.0)
        );
    }
    else
    {
        Fr_ = max
        (
            min
            (
                scalar(1.0),
                pow
                (
                    (scalar(1.0)-exp(-0.002*Lc/Lk()))
                    /
                    (scalar(1.0)-exp(-0.002*Li()/Lk()) + SMALL),
                    2.0
                )
            ),
            scalar(0.0)
        );
    }

    Fr_.correctBoundaryConditions();

    nut_ = Fr_*a1_*k_/max(a1_*omega_, F2()*sqrt(S2));
    nut_.correctBoundaryConditions();
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //

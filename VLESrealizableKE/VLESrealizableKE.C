/*---------------------------------------------------------------------------*\
  =========                 |
  \\      /  F ield         | OpenFOAM: The Open Source CFD Toolbox
   \\    /   O peration     |
    \\  /    A nd           | Copyright (C) 2011-2013 OpenFOAM Foundation
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

#include "VLESrealizableKE.H"
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

defineTypeNameAndDebug(VLESrealizableKE, 0);
addToRunTimeSelectionTable(RASModel, VLESrealizableKE, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> VLESrealizableKE::F1() const
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

        const volScalarField CDkOmega
        (
            (2*0.856)*(fvc::grad(k_) 
            & 
            fvc::grad(epsilon_/k_/0.09))/(epsilon_/k_/0.09)
        );

        tmp<volScalarField> CDkOmegaPlus = max
        (
            CDkOmega,
            dimensionedScalar("1.0e-10", dimless/sqr(dimTime), 1.0e-10)
        );

        tmp<volScalarField> arg1 = min
        (
            max
            (
                (scalar(1)/0.09)*sqrt(k_)/(epsilon_/k_*y_/0.09),
                scalar(500)*nu()/(sqr(y_)*epsilon_/k_/0.09)
            ),
            (4*0.856)*k_/(CDkOmegaPlus*sqr(y_))
        );

        return tanh(pow4(arg1));
    }
}

tmp<volScalarField> VLESrealizableKE::rCmu
(
    const volTensorField& gradU,
    const volScalarField& S2,
    const volScalarField& magS
)
{
    tmp<volSymmTensorField> tS = dev(symm(gradU));
    const volSymmTensorField& S = tS();

    tmp<volScalarField> W
    (
        (2*sqrt(2.0))*((S&S)&&S)
       /(
            magS*S2
          + dimensionedScalar("small", dimensionSet(0, 0, -3, 0, 0), SMALL)
        )
    );

    tS.clear();

    tmp<volScalarField> phis
    (
        (1.0/3.0)*acos(min(max(sqrt(6.0)*W, -scalar(1)), scalar(1)))
    );
    tmp<volScalarField> As = sqrt(6.0)*cos(phis);
    tmp<volScalarField> Us = sqrt(S2/2.0 + magSqr(skew(gradU)));

    return 1.0/(A0_ + As*Us*k_/epsilon_);
}


tmp<volScalarField> VLESrealizableKE::rCmu
(
    const volTensorField& gradU
)
{
    const volScalarField S2(2*magSqr(dev(symm(gradU))));
    tmp<volScalarField> magS = sqrt(S2);
    return rCmu(gradU, S2, magS);
}

void VLESrealizableKE::correctNut()
{
    // Currently this function is not implemented due to the complexity of
    // evaluating nut.  Better calculate nut at the end of correct()
    NotImplemented;
}

// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

VLESrealizableKE::VLESrealizableKE
(
    const geometricOneField& alpha,
    const geometricOneField& rho,
    const volVectorField& U,
    const surfaceScalarField& alphaRhoPhi,
    const surfaceScalarField& phi,
    const transportModel& transport,
    const word& propertiesName,
    const word& type
):
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
    Cmu_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cmu",
            coeffDict_,
            0.09
        )
    ),
    A0_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "A0",
            coeffDict_,
            4.0
        )
    ),
    C2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "C2",
            coeffDict_,
            1.9
        )
    ),
    sigmak_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmak",
            coeffDict_,
            1.0
        )
    ),
    sigmaEps_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaEps",
            coeffDict_,
            1.2
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
    epsilon_
    (
        IOobject
        (
            "epsilon",
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
    bound(epsilon_, epsilonMin_);

    if (type == typeName)
    {
	nut_ = rCmu(fvc::grad(U_))*sqr(k_)/epsilon_;
	nut_.correctBoundaryConditions();

        printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

bool VLESrealizableKE::read()
{
    if (RASModel::read())
    {
        delayed_.readIfPresent("delayed", coeffDict());
        Cmu_.readIfPresent(coeffDict());
        A0_.readIfPresent(coeffDict());
        C2_.readIfPresent(coeffDict());
        sigmak_.readIfPresent(coeffDict());
        sigmaEps_.readIfPresent(coeffDict());
        Cx_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void VLESrealizableKE::validate()
{}

void VLESrealizableKE::correct()
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
        FatalErrorIn("realizableVLESKE.C")
            << "Case is not 3D or 2D, LES is not applicable"
            << exit(FatalError);
    }

    tmp<volScalarField> const Li = pow(k_,3.0/2.0)/epsilon_;
    tmp<volScalarField> const Lk = pow(nu(),3.0/4.0)/pow(epsilon_,1.0/4.0);

    volTensorField const gradU(fvc::grad(U_));
    volScalarField const S2(2*magSqr(dev(symm(gradU))));
    volScalarField const magS(sqrt(S2));

    volScalarField const eta(magS*k_/epsilon_);
    tmp<volScalarField> C1 = max(eta/(scalar(5) + eta), scalar(0.43));

    volScalarField G(GName(), nut_*S2);

    // Update epsilon and G at the wall
    epsilon_.boundaryFieldRef().updateCoeffs();
    

    if (delayed_)
    {
        Fr_ = max
        (
            min
            (
                scalar(1.0),
                pow
                (
                    (scalar(1.0)-(1-this->F1())*exp(-0.002*Lc/Lk()))
                    /
                    (scalar(1.0)-(1-this->F1())*exp(-0.002*Li()/Lk()) + SMALL),
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

    // Dissipation equation
    tmp<fvScalarMatrix> epsEqn
    (
        fvm::ddt(epsilon_)
      + fvm::div(phi_, epsilon_)
      - fvm::Sp(fvc::div(phi_), epsilon_)
      - fvm::laplacian(DepsilonEff(), epsilon_)
     ==
        C1*magS*epsilon_
      - fvm::Sp
        (
            C2_*epsilon_/(k_ + sqrt(nu()*epsilon_)),
            epsilon_
        )
    );

    epsEqn.ref().relax();

    epsEqn.ref().boundaryManipulate(epsilon_.boundaryFieldRef());

   // mesh_.updateFvMatrix(epsEqn());
    solve(epsEqn);
    bound(epsilon_, epsilonMin_);


    // Turbulent kinetic energy equation
    tmp<fvScalarMatrix> kEqn
    (
        fvm::ddt(k_)
      + fvm::div(phi_, k_)
      - fvm::Sp(fvc::div(phi_), k_)
      - fvm::laplacian(DkEff(), k_)
     ==
        G - fvm::Sp(epsilon_/k_, k_)
    );

    kEqn.ref().relax();
    //mesh_.updateFvMatrix(kEqn());
    solve(kEqn);
    bound(k_, kMin_);

    nut_ = Fr_*rCmu(gradU, S2, magS)*sqr(k_)/epsilon_;
    nut_.correctBoundaryConditions();

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //

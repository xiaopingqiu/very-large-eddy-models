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

#include "VLESSpalartAllmaras.H"
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

defineTypeNameAndDebug(VLESSpalartAllmaras, 0);
addToRunTimeSelectionTable(RASModel, VLESSpalartAllmaras, dictionary);

// * * * * * * * * * * * * Private Member Functions  * * * * * * * * * * * * //

tmp<volScalarField> VLESSpalartAllmaras::F1() const
{

    volScalarField const CDkOmega
    (
        (2*0.856)*(fvc::grad(k()) 
        & 
        fvc::grad(epsilon()/k()/0.09))/(epsilon()/k()/0.09)
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
            (scalar(1)/0.09)*sqrt(k())/(epsilon()/k()*y_/0.09),
            scalar(500)*nu()/(sqr(y_)*epsilon()/k()/0.09)
        ),
        (4*0.856)*k()/(CDkOmegaPlus*sqr(y_))
    );

    return tanh(pow4(arg1));
}

tmp<volScalarField> VLESSpalartAllmaras::chi() const
{
    return nuTilda_/nu();
}


tmp<volScalarField> VLESSpalartAllmaras::fv1(const volScalarField& chi) const
{
    const volScalarField chi3(pow3(chi));
    return chi3/(chi3 + pow3(Cv1_));
}


tmp<volScalarField> VLESSpalartAllmaras::fv2
(
    const volScalarField& chi,
    const volScalarField& fv1
) const
{
    return (scalar(1.0) - chi/(scalar(1.0) + chi*fv1));
}

tmp<volScalarField> VLESSpalartAllmaras::fw(volScalarField const& Stilda) const
{
    volScalarField r
    (
        min
        (
            nuTilda_
           /(
               max
               (
                   Stilda,
                   dimensionedScalar("SMALL", Stilda.dimensions(), SMALL)
               )
              *sqr(kappa_*y_)
            ),
            scalar(10.0)
        )
    );
    r.boundaryFieldRef() == 0.0;

    volScalarField const g(r + Cw2_*(pow6(r) - r));

    return g*pow((1.0 + pow6(Cw3_))/(pow6(g) + pow6(Cw3_)), 1.0/6.0);
}

// * * * * * * * * * * * * Protected Member Functions  * * * * * * * * * * * //

void VLESSpalartAllmaras::correctNut()
{
    // Currently this function is not implemented due to the complexity of
    // evaluating nut.  Better calculate nut at the end of correct()
    NotImplemented;
}


// * * * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * //

VLESSpalartAllmaras::VLESSpalartAllmaras
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
    curvatureCorrection_
    (
        coeffDict_.lookupOrDefault<Switch>
        (
            "curvatureCorrection", 
            false
        )
    ),

    delayed_(coeffDict_.lookupOrDefault<Switch>("delayed", false)),
    sigmaNut_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "sigmaNut",
            coeffDict_,
            0.66666
        )
    ),
    kappa_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "kappa",
            coeffDict_,
            0.41
        )
    ),

    Cb1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb1",
            coeffDict_,
            0.1355
        )
    ),
    Cb2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cb2",
            coeffDict_,
            0.622
        )
    ),
    Cw1_(Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_),
    Cw2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw2",
            coeffDict_,
            0.3
        )
    ),
    Cw3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cw3",
            coeffDict_,
            2.0
        )
    ),
    Cv1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cv1",
            coeffDict_,
            7.1
        )
    ),

    Cr1_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cr1",
            coeffDict_,
            1.0
        )
    ),
    Cr2_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cr2",
            coeffDict_,
            12.0
        )
    ),
    Cr3_
    (
        dimensioned<scalar>::lookupOrAddToDict
        (
            "Cr3",
            coeffDict_,
            1.0
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
    nuTilda_
    (
        IOobject
        (
            "nuTilda",
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

    fr1_
    (
        IOobject
        (
            "fr1",
            runTime_.timeName(),
            mesh_,
            IOobject::NO_READ,
            IOobject::AUTO_WRITE
        ),
        mesh_,
        dimensionedScalar("one", dimless, 1)
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
        dimensionedScalar("one", dimless, 1),
        calculatedFvPatchScalarField::typeName
    ),
    y_(wallDist::New(mesh_).y())
{
    if (type == typeName)
    {
	if (curvatureCorrection_)
	{
	    Info<<" Curvature correction modification enabled" << endl;
	}

	printCoeffs(type);
    }
}


// * * * * * * * * * * * * * * * Member Functions  * * * * * * * * * * * * * //

tmp<volScalarField> VLESSpalartAllmaras::DnuTildaEff() const
{
    return tmp<volScalarField>
    (
        new volScalarField("DnuTildaEff", (nuTilda_ + nu())/sigmaNut_)
    );
}


tmp<volScalarField> VLESSpalartAllmaras::k() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "k",
                runTime_.timeName(),
                mesh_
            ),
            (
                sqrt
                (
                    epsilon()*nut_/0.09
                    +
                    dimensionedScalar
                    (
                        "sk",
                        dimensionSet(0,4,-4,0,0,0,0),
                        scalar(VSMALL)
                    )
                )
            )
        )
    );
}


tmp<volScalarField> VLESSpalartAllmaras::epsilon() const
{
    return tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "epsilon",
                runTime_.timeName(),
                mesh_
            ),
            (
                scalar(2.0)*nuEff()*magSqr(symm(fvc::grad(U())))
                +
                dimensionedScalar
                (
                    "se",
                    dimensionSet(0,2,-3,0,0,0,0),
                    scalar(VSMALL)
                )
            )
        )
    );
}

bool VLESSpalartAllmaras::read()
{
    if (eddyViscosity<incompressible::RASModel>::read())
    {
        curvatureCorrection_.readIfPresent("curvatureCorrection", coeffDict());
        delayed_.readIfPresent("delayed", coeffDict());
        sigmaNut_.readIfPresent(coeffDict());
        kappa_.readIfPresent(coeffDict());
        Cb1_.readIfPresent(coeffDict());
        Cb2_.readIfPresent(coeffDict());
        Cw1_ = Cb1_/sqr(kappa_) + (1.0 + Cb2_)/sigmaNut_;
        Cw2_.readIfPresent(coeffDict());
        Cw3_.readIfPresent(coeffDict());
        Cv1_.readIfPresent(coeffDict());
        Cr1_.readIfPresent(coeffDict());
        Cr2_.readIfPresent(coeffDict());
        Cr3_.readIfPresent(coeffDict());
        Cx_.readIfPresent(coeffDict());

        return true;
    }
    else
    {
        return false;
    }
}


void VLESSpalartAllmaras::validate()
{}


void VLESSpalartAllmaras::correct()
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

        Lc.ref() = Cx_*sqrt(mesh_.V()/thickness);
    }
    else
    {
        FatalErrorIn("VLESSpalartAllmaras.C")
            << "Case is not 3D or 2D, VLES is not applicable"
            << exit(FatalError);
    }  


    volScalarField const chi(this->chi());
    volScalarField const fv1(this->fv1(chi));
    volTensorField const Omegaij(skew(fvc::grad(this->U_)));
    volScalarField const sqrOmega(2*magSqr(Omegaij));

    nut_ = fv1*nuTilda_;
    nut_.correctBoundaryConditions();

    tmp<volScalarField> const Li = pow(k(),3.0/2.0)/epsilon();

    tmp<volScalarField> const Lk = pow(nu(),3.0/4.0)/pow(epsilon(),1.0/4.0);

    if (delayed_)
    {
        Fr_ = min
        (
            scalar(1.0),
            pow
            (
                (scalar(1.0)-(1-F1())*exp(-0.002*Lc/Lk()))
                /
                (scalar(1.0)-(1-F1())*exp(-0.002*Li()/Lk())),
                2.0
            )
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
                    (scalar(1.0)-exp(-0.002*Li()/Lk())),
                    2.0
                )
            ),
            scalar(0.0)
        );
    }

    Fr_.correctBoundaryConditions();

    // Curvature correction terms
    if (curvatureCorrection_)
    {
       volSymmTensorField const Sij(symm(fvc::grad(this->U_)));    
       volScalarField const sqrS(2*magSqr(Sij));
       dimensionedScalar const smallOmega
       (
           "smallOmega", 
           sqrOmega.dimensions(), 
           SMALL
       );
       volScalarField const sqrD(0.5*(sqrS + sqrOmega + smallOmega));
       volScalarField const rStar(sqrt(sqrS)/sqrt(sqrOmega+smallOmega));
       volSymmTensorField const DSijDt(fvc::DDt(this->phi_,Sij));
       volScalarField const rTilda
       (  
           (scalar(2.0)/sqr(sqrD))*(Omegaij && (Sij & DSijDt))
       );
       fr1_ = 
       (scalar(1.0) + Cr1_)*scalar(2.0)*rStar/(scalar(1.0) + rStar)
            *(scalar(1.0)-Cr3_*atan(Cr2_*rTilda)) - Cr1_;
    }

    volScalarField const Stilda
    (
        sqrt(sqrOmega) + fv2(chi, fv1)*nuTilda_/sqr(kappa_*y_)
    );

    tmp<fvScalarMatrix> nuTildaEqn
    (
        fvm::ddt(nuTilda_)
      + fvm::div(phi_, nuTilda_)
      - fvm::laplacian(DnuTildaEff(), nuTilda_)
      - Cb2_/sigmaNut_*magSqr(fvc::grad(nuTilda_))
     ==
        fr1_*Cb1_*Stilda*nuTilda_
      - fvm::Sp(Cw1_*fw(Stilda)*nuTilda_/sqr(y_), nuTilda_)
    );

    nuTildaEqn.ref().relax();

    //mesh_.updateFvMatrix(nuTildaEqn());
    solve(nuTildaEqn);
    bound(nuTilda_, dimensionedScalar("0", nuTilda_.dimensions(), 0.0));
    nuTilda_.correctBoundaryConditions();

    // Re-calculate viscosity
    nut_.ref() = Fr_.internalField()*fv1*nuTilda_.internalField();
    nut_.correctBoundaryConditions();

}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * //

} // End namespace RASModels
} // End namespace incompressible
} // End namespace Foam

// ************************************************************************* //

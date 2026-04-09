#include "BuoyantCurvatureSwirlTools.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMesh.H"
#include "uniformDimensionedFields.H"
#include "fvcGrad.H"
#include "mathematicalConstants.H"

namespace Foam
{

defineTypeNameAndDebug(BuoyantCurvatureSwirlTools, 0);


// * * * * * * * * * * * * * * Private Helper Functions  * * * * * * * * * //

vector BuoyantCurvatureSwirlTools::normalisedAxis(const vector& axis)
{
    scalar magAxis = mag(axis);

    if (magAxis < SMALL)
    {
        FatalErrorInFunction
            << "axisDirection has near-zero magnitude: " << axis
            << exit(FatalError);
    }

    return axis/magAxis;
}


vectorField BuoyantCurvatureSwirlTools::cellRelPosition() const
{
    const vectorField& C = mesh_.C();


    vectorField rel(C.size(), Zero);

    forAll(C, i)
    {
        rel[i] = C[i] - axisOrigin_;
    }

    return rel;
}


tmp<volScalarField> BuoyantCurvatureSwirlTools::axialVelocity
(
    const volVectorField& U
) const
{
    const dimensionedVector eAxis
    (
        "eAxis",
        dimless,
        axisDirection_
    );

    return U & eAxis;
}


void BuoyantCurvatureSwirlTools::radialDistanceAndTangentialVelocity
(
    const volVectorField& U,
    tmp<volScalarField>& tr,
    tmp<volScalarField>& tuTheta
) const
{
    const vectorField relPos(cellRelPosition());

    volScalarField* rPtr = new volScalarField
    (
        IOobject
        (
            "r",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimLength, 0.0)
    );

    volScalarField* uThetaPtr = new volScalarField
    (
        IOobject
        (
            "uTheta",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimVelocity, 0.0)
    );

    volScalarField& uTheta = *uThetaPtr;
    volScalarField& r = *rPtr;

    forAll(U, i)
    {
        const vector& rel = relPos[i];
        const scalar axialProj = (rel & axisDirection_);
        const vector vr = rel - axialProj*axisDirection_;
        const scalar rMag = mag(vr);
        const vector er = vr/(rMag + SMALL);
        const vector eTheta = axisDirection_ ^ er;

        r[i] = rMag;
        uTheta[i] = (U[i] & eTheta);
    }

    r.correctBoundaryConditions();
    uTheta.correctBoundaryConditions();

    tr = tmp<volScalarField>(rPtr);
    tuTheta = tmp<volScalarField>(uThetaPtr);
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

BuoyantCurvatureSwirlTools::BuoyantCurvatureSwirlTools
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
    dict_(dict),
    curvatureCorrection_
    (
        dict.lookupOrDefault<Switch>("curvatureCorrection", true)
    ),
    swirlCorrection_
    (
        dict.lookupOrDefault<Switch>("swirlCorrection", true)
    ),
    buoyancyCorrection_
    (
        dict.lookupOrDefault<Switch>("buoyancyCorrection", false)
    ),
    axisOrigin_(dict.lookupOrDefault<vector>("axisOrigin", vector::zero)),
    axisDirection_
    (
        normalisedAxis
        (
            dict.lookupOrDefault<vector>("axisDirection", vector(0, 0, 1))
        )
    ),
    S0_
    (
        "S0",
        dimless,
        dict.lookupOrDefault<scalar>("S0", 0.2)
    ),
    S1_
    (
        "S1",
        dimless,
        dict.lookupOrDefault<scalar>("S1", 0.8)
    ),
    frMax_
    (
        "frMax",
        dimless,
        dict.lookupOrDefault<scalar>("frMax", 1.25)
    ),
    swirlEps_
    (
        "swirlEps",
        dimVelocity,
        dict.lookupOrDefault<scalar>("swirlEps", 1e-6)
    ),
    Cg_
    (
        "Cg",
        dimless,
        dict.lookupOrDefault<scalar>("Cg", 1.0)
    ),
    writeFields_(dict.lookupOrDefault<Switch>("writeFields", false))
{
    if (mag(axisDirection_) < SMALL)
    {
        FatalErrorInFunction
            << "Invalid axisDirection after normalization."
            << exit(FatalError);
    }

    if (S1_.value() <= S0_.value())
    {
        FatalErrorInFunction
            << "Require S1 > S0, but got S0=" << S0_.value()
            << " and S1=" << S1_.value()
            << exit(FatalError);
    }
}


bool BuoyantCurvatureSwirlTools::read(const dictionary& dict)
{
    curvatureCorrection_ =
        dict.lookupOrDefault<Switch>("curvatureCorrection", curvatureCorrection_);
    swirlCorrection_ =
        dict.lookupOrDefault<Switch>("swirlCorrection", swirlCorrection_);
    buoyancyCorrection_ =
        dict.lookupOrDefault<Switch>("buoyancyCorrection", buoyancyCorrection_);
    axisOrigin_ = dict.lookupOrDefault<vector>("axisOrigin", axisOrigin_);
    axisDirection_ = normalisedAxis
    (
        dict.lookupOrDefault<vector>("axisDirection", axisDirection_)
    );

    S0_.readIfPresent(dict);
    S1_.readIfPresent(dict);
    frMax_.readIfPresent(dict);
    swirlEps_.readIfPresent(dict);
    Cg_.readIfPresent(dict);
    writeFields_ = dict.lookupOrDefault<Switch>("writeFields", writeFields_);

    if (mag(axisDirection_) < SMALL)
    {
        FatalErrorInFunction
            << "Invalid axisDirection after normalization."
            << exit(FatalError);
    }

    if (S1_.value() <= S0_.value())
    {
        FatalErrorInFunction
            << "Require S1 > S0, but got S0=" << S0_.value()
            << " and S1=" << S1_.value()
            << exit(FatalError);
    }

    return true;
}


// * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * * *

bool BuoyantCurvatureSwirlTools::active() const
{
    return curvatureCorrection_ || swirlCorrection_ || buoyancyCorrection_;
}


curvatureSwirlData BuoyantCurvatureSwirlTools::evaluate
(
    const volVectorField& U,
    const volTensorField& gradU,
    const volScalarField& S2,
    const volScalarField& magS,
    const volScalarField& rho,
    const volScalarField& k,
    const volScalarField& nut
) const
{
    curvatureSwirlData data;

    const dimensionedScalar zero(dimless, 0.0);
    const dimensionedScalar one(dimless, 1.0);
    const dimensionedScalar zeroVelocity(dimVelocity, 0.0);
    const dimensionedScalar zeroGb(dimless/dimTime, 0.0);
    const bool needSwirlFields =
        swirlCorrection_ || writeFields_;
    const bool needCurvatureFields =
        curvatureCorrection_ || writeFields_;

    if (needSwirlFields)
    {
        radialDistanceAndTangentialVelocity(U, data.r, data.uTheta);
        data.Uaxial = axialVelocity(U);

        data.Slocal = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Slocal",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mag(data.uTheta())/(mag(data.Uaxial()) + swirlEps_)
            )
        );

        data.Fswirl = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Fswirl",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                min
                (
                    max
                    (
                        (data.Slocal() - S0_)/(S1_ - S0_),
                        zero
                    ),
                    one
                )
            )
        );

        if (!swirlCorrection_)
        {
            data.Slocal.ref() = zeroVelocity;
            data.Fswirl.ref() = zero;
        }
    }
    else
    {
        data.Fswirl = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Fswirl",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                zero
            )
        );
    }

    if (needCurvatureFields && curvatureCorrection_)
    {
        const volTensorField Omega(skew(gradU));

        const volScalarField::Internal magOmega
        (
            typedName("magOmega"),
            sqrt(2.0*magSqr(Omega)())
        );

        const volScalarField::Internal magScc
        (
            typedName("magScc"),
            magS()
        );

        const volScalarField::Internal frInternal
        (
            typedName("fr"),
            min
            (
                max
                (
                    (
                        magScc
                       /max
                        (
                            magOmega,
                            dimensionedScalar
                            (
                                "smallMagOmega",
                                magOmega.dimensions(),
                                SMALL
                            )
                        )
                    ),
                    scalar(0)
                ),
                frMax_.value()
            )
        );

        data.fr = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "fr",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                frInternal,
                data.Fswirl().boundaryField()
            )
        );
    }
    else
    {
        data.fr = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "fr",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                one
            )
        );
    }

    if (curvatureCorrection_)
    {
        if (swirlCorrection_)
        {
            data.frEff = tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "frEff",
                        mesh_.time().name(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    one + data.Fswirl()*(data.fr() - one)
                )
            );
        }
        else
        {
            data.frEff = tmp<volScalarField>
            (
                new volScalarField
                (
                    IOobject
                    (
                        "frEff",
                        mesh_.time().name(),
                        mesh_,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    data.fr()
                )
            );
        }
    }
    else
    {
        data.frEff = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "frEff",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                mesh_,
                one
            )
        );
    }

    data.GbCoef = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "GbCoef",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            zeroGb
        )
    );

    data.buoyLimiter = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "buoyLimiter",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            zero
        )
    );

    if (buoyancyCorrection_)
    {
        const uniformDimensionedVectorField& g =
            mesh_.objectRegistry::lookupObject<uniformDimensionedVectorField>("g");

        if (mag(g.value()) > small)
        {
            const vector gHat(g.value()/mag(g.value()));
            const volVectorField gradRho(fvc::grad(rho));
            const dimensionedScalar kSmall("kSmall", k.dimensions(), SMALL);
            const dimensionedScalar uSmall("uSmall", dimVelocity, SMALL);

            data.GbCoef.ref() =
                Cg_*nut*(g & gradRho)/(k + kSmall);

            const volScalarField v(typedName("v"), gHat & U);
            const volScalarField u
            (
                typedName("u"),
                mag(U - gHat*v) + uSmall
            );

            data.buoyLimiter.ref() = tanh(mag(v)/u);
        }
    }

    if (writeFields_)
    {
        if (data.r.valid()) data.r().write();
        if (data.uTheta.valid()) data.uTheta().write();
        if (data.Uaxial.valid()) data.Uaxial().write();
        if (data.Slocal.valid()) data.Slocal().write();
        if (data.Fswirl.valid()) data.Fswirl().write();
        if (data.fr.valid()) data.fr().write();
        data.frEff().write();
        data.GbCoef().write();
        data.buoyLimiter().write();
    }

    return data;
}


} // End namespace Foam

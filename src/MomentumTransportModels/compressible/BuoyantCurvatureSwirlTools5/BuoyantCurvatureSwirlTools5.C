#include "BuoyantCurvatureSwirlTools5.H"
#include "addToRunTimeSelectionTable.H"
#include "volFields.H"
#include "surfaceFields.H"
#include "fvMesh.H"
#include "uniformDimensionedFields.H"
#include "fvcGrad.H"
#include "mathematicalConstants.H"

namespace Foam
{

defineTypeNameAndDebug(BuoyantCurvatureSwirlTools5, 0);


// * * * * * * * * * * * * * * Private Helper Functions  * * * * * * * * * //

vector BuoyantCurvatureSwirlTools5::normalisedAxis(const vector& axis)
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


vectorField BuoyantCurvatureSwirlTools5::cellRelPosition() const
{
    const vectorField& C = mesh_.C();

    vectorField rel(C.size(), Zero);

    forAll(C, i)
    {
        rel[i] = C[i] - axisOrigin_;
    }

    return rel;
}


tmp<volScalarField> BuoyantCurvatureSwirlTools5::axialVelocity
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

    tmp<volScalarField> tUaxial
    (
        new volScalarField
        (
            IOobject
            (
                "Uaxial",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedScalar("zero", dimVelocity, 0.0)
        )
    );

    tUaxial.ref() = (U & eAxis);

    return tUaxial;
}


tmp<volScalarField> BuoyantCurvatureSwirlTools5::tangentialVelocity
(
    const volVectorField& U
) const
{
    const vectorField relPos(cellRelPosition());

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

    forAll(U, i)
    {
        const vector& rel = relPos[i];
        const scalar axialProj = (rel & axisDirection_);
        const vector vr = rel - axialProj*axisDirection_;
        const scalar rMag = mag(vr);
        const vector er = vr/(rMag + SMALL);
        const vector eTheta = axisDirection_ ^ er;

        uTheta[i] = (U[i] & eTheta);
    }

    uTheta.correctBoundaryConditions();

    return tmp<volScalarField>(uThetaPtr);
}


tmp<volScalarField> BuoyantCurvatureSwirlTools5::radialDistance() const
{
    const vectorField relPos(cellRelPosition());

    volScalarField* rPtr = new volScalarField
    (
        IOobject
        (
            "rRadial",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedScalar("zero", dimLength, 0.0)
    );

    volScalarField& r = *rPtr;

    forAll(r, i)
    {
        const vector& rel = relPos[i];
        const scalar axialProj = (rel & axisDirection_);
        const vector vr = rel - axialProj*axisDirection_;

        r[i] = mag(vr);
    }

    r.correctBoundaryConditions();

    return tmp<volScalarField>(rPtr);
}


tmp<volVectorField> BuoyantCurvatureSwirlTools5::radialUnitVector() const
{
    const vectorField relPos(cellRelPosition());

    volVectorField* erPtr = new volVectorField
    (
        IOobject
        (
            "eRadial",
            mesh_.time().name(),
            mesh_,
            IOobject::NO_READ,
            IOobject::NO_WRITE
        ),
        mesh_,
        dimensionedVector("zero", dimless, vector::zero)
    );

    volVectorField& er = *erPtr;

    forAll(er, i)
    {
        const vector& rel = relPos[i];
        const scalar axialProj = (rel & axisDirection_);
        const vector vr = rel - axialProj*axisDirection_;
        const scalar rMag = mag(vr);

        er[i] = vr/(rMag + SMALL);
    }

    er.correctBoundaryConditions();

    return tmp<volVectorField>(erPtr);
}


// * * * * * * * * * * * * * * Constructors  * * * * * * * * * * * * * * * //

BuoyantCurvatureSwirlTools5::BuoyantCurvatureSwirlTools5
(
    const fvMesh& mesh,
    const dictionary& dict
)
:
    mesh_(mesh),
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
    frMax_
    (
        "frMax",
        dimless,
        dict.lookupOrDefault<scalar>("frMax", 1.25)
    ),
    cCurv_
    (
        "cCurv",
        dimless,
        dict.lookupOrDefault<scalar>("cCurv", 1.0)
    ),
    cr1_
    (
        "cr1",
        dimless,
        dict.lookupOrDefault<scalar>("cr1", 1.0)
    ),
    cr2_
    (
        "cr2",
        dimless,
        dict.lookupOrDefault<scalar>("cr2", 2.0)
    ),
    cr3_
    (
        "cr3",
        dimless,
        dict.lookupOrDefault<scalar>("cr3", 1.0)
    ),
    rMin_
    (
        "rMin",
        dimLength,
        dict.lookupOrDefault<scalar>("rMin", 1e-6)
    ),
    cPhi_
    (
        "cPhi",
        dimless,
        dict.lookupOrDefault<scalar>("cPhi", 1.0)
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

    if (cCurv_.value() <= 0)
    {
        FatalErrorInFunction
            << "Require cCurv > 0, but got cCurv=" << cCurv_.value()
            << exit(FatalError);
    }

    if (rMin_.value() <= 0)
    {
        FatalErrorInFunction
            << "Require rMin > 0, but got rMin=" << rMin_.value()
            << exit(FatalError);
    }

    if (cPhi_.value() <= 0)
    {
        FatalErrorInFunction
            << "Require cPhi > 0, but got cPhi=" << cPhi_.value()
            << exit(FatalError);
    }

    Info<< "BuoyantCurvatureSwirlTools5 (BCST5): "
        << "L = r*uTheta,  Phi = (1/r^3) d(L^2)/dr,  "
        << "Fgate = clamp(-Phi/(cPhi*D2), 0, 1)  [destabilizing-only gate]"
        << nl;
}


bool BuoyantCurvatureSwirlTools5::read(const dictionary& dict)
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

    frMax_.readIfPresent(dict);
    cCurv_.readIfPresent(dict);
    cr1_.readIfPresent(dict);
    cr2_.readIfPresent(dict);
    cr3_.readIfPresent(dict);
    rMin_.readIfPresent(dict);
    cPhi_.readIfPresent(dict);
    Cg_.readIfPresent(dict);
    writeFields_ = dict.lookupOrDefault<Switch>("writeFields", writeFields_);

    if (mag(axisDirection_) < SMALL)
    {
        FatalErrorInFunction
            << "Invalid axisDirection after normalization."
            << exit(FatalError);
    }

    if (cCurv_.value() <= 0)
    {
        FatalErrorInFunction
            << "Require cCurv > 0, but got cCurv=" << cCurv_.value()
            << exit(FatalError);
    }

    if (rMin_.value() <= 0)
    {
        FatalErrorInFunction
            << "Require rMin > 0, but got rMin=" << rMin_.value()
            << exit(FatalError);
    }

    if (cPhi_.value() <= 0)
    {
        FatalErrorInFunction
            << "Require cPhi > 0, but got cPhi=" << cPhi_.value()
            << exit(FatalError);
    }

    return true;
}


bool BuoyantCurvatureSwirlTools5::active() const
{
    return curvatureCorrection_ || swirlCorrection_ || buoyancyCorrection_;
}


curvatureSwirlData5 BuoyantCurvatureSwirlTools5::evaluate
(
    const volVectorField& U,
    const volTensorField& gradU,
    const volScalarField& rho,
    const volScalarField& k,
    const volScalarField& omegaLike,
    const volScalarField& nut
) const
{
    curvatureSwirlData5 data;

    const dimensionedScalar zero(dimless, 0.0);
    const dimensionedScalar one(dimless, 1.0);
    const dimensionedScalar zeroGb("zeroGb", k.dimensions()/dimTime, 0.0);
    const bool needSwirlFields =
        swirlCorrection_ || writeFields_;
    const bool needCurvatureFields =
        curvatureCorrection_ || writeFields_;

    if (needSwirlFields)
    {
        data.uTheta = tangentialVelocity(U);
        data.Uaxial = axialVelocity(U);

        const volScalarField r(radialDistance());
        const volVectorField er(radialUnitVector());

        // L = r * uTheta  (specific angular momentum)
        data.L = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "L",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                r*data.uTheta()
            )
        );

        // Phi = (1/r^3) d(L^2)/dr, evaluated locally as
        // grad(L^2) . er  /  max(r, rMin)^3
        const volScalarField L2(typedName("L2"), sqr(data.L()));
        const volVectorField gradL2(typedName("gradL2"), fvc::grad(L2));
        const volScalarField dL2dr(typedName("dL2dr"), gradL2 & er);
        const volScalarField r3
        (
            typedName("r3"),
            pow3(max(r, rMin_))
        );

        data.Phi = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Phi",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                dL2dr/r3
            )
        );

        // Local mean-flow deformation scale used to normalize Phi into a
        // dimensionless gate: D2 = max(S^2, 0.09*omegaLike^2).
        const volSymmTensorField Sgate(symm(gradU));
        const volTensorField Omegagate(skew(gradU));

        const volScalarField strainMagGate
        (
            typedName("strainMagGate"),
            sqrt(2.0*magSqr(Sgate))
        );

        const volScalarField rotationMagGate
        (
            typedName("rotationMagGate"),
            sqrt(2.0*magSqr(Omegagate))
        );

        const volScalarField D2gate
        (
            typedName("D2gate"),
            max(sqr(strainMagGate), 0.09*sqr(omegaLike))
        );

        const dimensionedScalar phiSmall
        (
            "phiSmall",
            data.Phi().dimensions(),
            SMALL
        );

        // Fgate = clamp( -Phi / (cPhi*D2), 0, 1 )  -- destabilizing-only:
        // fires where Phi < 0 (Rayleigh-unstable), fully off where Phi >= 0.
        data.Fgate = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Fgate",
                    mesh_.time().name(),
                    mesh_,
                    IOobject::NO_READ,
                    IOobject::NO_WRITE
                ),
                min
                (
                    max
                    (
                        -data.Phi()/(cPhi_*D2gate + phiSmall),
                        zero
                    ),
                    one
                )
            )
        );

        if (!swirlCorrection_)
        {
            data.Fgate.ref() = zero;
        }
    }
    else
    {
        data.Fgate = tmp<volScalarField>
        (
            new volScalarField
            (
                IOobject
                (
                    "Fgate",
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
        const volSymmTensorField S(symm(gradU));
        const volTensorField Omega(skew(gradU));

        const volScalarField strainMag
        (
            typedName("strainMag"),
            sqrt(2.0*magSqr(S))
        );

        const volScalarField rotationMag
        (
            typedName("rotationMag"),
            sqrt(2.0*magSqr(Omega))
        );

        const scalar deltaTValue = max(mesh_.time().deltaTValue(), SMALL);
        const dimensionedScalar deltaT("deltaT", dimTime, deltaTValue);
        const dimensionedScalar omegaSmall
        (
            "omegaSmall",
            omegaLike.dimensions(),
            SMALL
        );
        const dimensionedScalar d2Small(sqr(omegaSmall));
        const dimensionedScalar d4Small(sqr(d2Small));

        const volScalarField D2
        (
            typedName("D2"),
            max(sqr(strainMag), 0.09*sqr(omegaLike))
        );

        const volScalarField D
        (
            typedName("D"),
            sqrt(D2)
        );

        const volSymmTensorField Sold(symm(fvc::grad(U.oldTime())));

        const volScalarField DSxxDt
        (
            typedName("DSxxDt"),
            (S.component(symmTensor::XX) - Sold.component(symmTensor::XX))/deltaT
          + (U & fvc::grad(S.component(symmTensor::XX)))
        );

        const volScalarField DSxyDt
        (
            typedName("DSxyDt"),
            (S.component(symmTensor::XY) - Sold.component(symmTensor::XY))/deltaT
          + (U & fvc::grad(S.component(symmTensor::XY)))
        );

        const volScalarField DSxzDt
        (
            typedName("DSxzDt"),
            (S.component(symmTensor::XZ) - Sold.component(symmTensor::XZ))/deltaT
          + (U & fvc::grad(S.component(symmTensor::XZ)))
        );

        const volScalarField DSyyDt
        (
            typedName("DSyyDt"),
            (S.component(symmTensor::YY) - Sold.component(symmTensor::YY))/deltaT
          + (U & fvc::grad(S.component(symmTensor::YY)))
        );

        const volScalarField DSyzDt
        (
            typedName("DSyzDt"),
            (S.component(symmTensor::YZ) - Sold.component(symmTensor::YZ))/deltaT
          + (U & fvc::grad(S.component(symmTensor::YZ)))
        );

        const volScalarField DSzzDt
        (
            typedName("DSzzDt"),
            (S.component(symmTensor::ZZ) - Sold.component(symmTensor::ZZ))/deltaT
          + (U & fvc::grad(S.component(symmTensor::ZZ)))
        );

        volSymmTensorField DSDt
        (
            IOobject
            (
                "DSDt",
                mesh_.time().name(),
                mesh_,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh_,
            dimensionedSymmTensor("zero", DSxxDt.dimensions(), symmTensor::zero)
        );

        DSDt.replace(symmTensor::XX, DSxxDt);
        DSDt.replace(symmTensor::XY, DSxyDt);
        DSDt.replace(symmTensor::XZ, DSxzDt);
        DSDt.replace(symmTensor::YY, DSyyDt);
        DSDt.replace(symmTensor::YZ, DSyzDt);
        DSDt.replace(symmTensor::ZZ, DSzzDt);

        const volScalarField rStar
        (
            typedName("rStar"),
            strainMag/max(rotationMag, omegaSmall)
        );

        const volScalarField rTilde
        (
            typedName("rTilde"),
            (twoSymm(Omega & S) && DSDt)/max(rotationMag*D2*D, d4Small)
        );

        const volScalarField fRot
        (
            typedName("fRot"),
            ((one + cr1_)*((2.0*rStar)/(one + rStar))*(one - cr3_*atan(cr2_*rTilde)))
          - cr1_
        );

        const volScalarField frTilda
        (
            typedName("frTilda"),
            max(min(fRot, frMax_), zero)
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
                mesh_,
                one
            )
        );

        data.fr.ref() = max(zero, one + cCurv_*(frTilda - one));
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
                    one + data.Fgate()*(data.fr() - one)
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

    data.Gb = tmp<volScalarField>
    (
        new volScalarField
        (
            IOobject
            (
                "Gb",
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
            const dimensionedScalar rhoSmall("rhoSmall", rho.dimensions(), SMALL);
            const dimensionedScalar uSmall("uSmall", dimVelocity, SMALL);

            data.Gb.ref() =
                Cg_*nut*(g & gradRho)/max(rho, rhoSmall);

            const volScalarField v(typedName("v"), gHat & U);
            const volScalarField u
            (
                typedName("u"),
                mag(U - gHat*v) + uSmall
            );

            data.buoyLimiter.ref() = tanh(mag(v)/u);
        }
    }

    if (writeFields_ && mesh_.time().writeTime())
    {
        if (data.uTheta.valid()) data.uTheta().write();
        if (data.Uaxial.valid()) data.Uaxial().write();
        if (data.L.valid()) data.L().write();
        if (data.Phi.valid()) data.Phi().write();
        if (data.Fgate.valid()) data.Fgate().write();
        if (data.fr.valid()) data.fr().write();
        data.frEff().write();
        data.Gb().write();
        data.buoyLimiter().write();
    }

    return data;
}


} // End namespace Foam

# realizableKE_BCS and BuoyantCurvatureSwirlTools: Step-by-Step Explanation

## 1. Include Block

```cpp
#include "realizableKE_BCS.H"
#include "fvModels.H"
#include "fvConstraints.H"
#include "bound.H"
#include "BuoyantCurvatureSwirlTools.H"
```

This block pulls in the model class itself, OpenFOAM finite-volume source and constraint infrastructure, the field bounding utility, and the custom helper that supplies curvature, swirl, and buoyancy corrections.

There is no direct equation here, but it tells us the model is built from:

realizableKE_BCS = base realizable RAS model + FV source framework + FV constraint framework + custom correction tool

## 2. `boundEpsilon()`

```cpp
epsilon_ = max(epsilon_, 0.09*sqr(k_)/(this->nutMaxCoeff_*this->nu()));
```

This enforces a lower bound on the dissipation rate so the turbulent viscosity does not blow up.

epsilon <- max(epsilon, 0.09 k^2 / (nutMaxCoeff nu))

Since the model later uses:

nu_t ~ k^2 / epsilon

this bound prevents unphysically small `epsilon` from causing excessively large turbulent viscosity.

## 3. `rCmu(gradU, S2, magS)`

```cpp
tmp<volSymmTensorField> tS = dev(symm(gradU));
const volSymmTensorField& S = tS();

const volScalarField W
(
    (2*sqrt(2.0))*((S&S)&&S)
   /(
        magS*S2
      + dimensionedScalar(dimensionSet(0, 0, -3, 0, 0), small)
    )
);

tS.clear();

const volScalarField phis
(
    (1.0/3.0)*acos(min(max(sqrt(6.0)*W, -scalar(1)), scalar(1)))
);
const volScalarField As(sqrt(6.0)*cos(phis));
const volScalarField Us(sqrt(S2/2.0 + magSqr(skew(gradU))));

return 1.0/(A0_ + As*Us*k_/epsilon_);
```

This is the realizable-model replacement for a constant `Cmu`.

S = dev(symm(grad U))

W = 2 sqrt(2) (S_ij S_jk S_ki) / (|S| S^2 + epsilon_small)

phi_s = (1/3) acos(clip(sqrt(6) W, -1, 1))

A_s = sqrt(6) cos(phi_s)

U_s = sqrt(S^2/2 + |skew(grad U)|^2)

rCmu = 1 / (A0 + A_s U_s k / epsilon)

Then the turbulent viscosity is computed later as:

nu_t = rCmu k^2 / epsilon

## 4. `correctNut()`

```cpp
boundEpsilon();
this->nut_ = rCmu(gradU, S2, magS)*sqr(k_)/epsilon_;
this->nut_.correctBoundaryConditions();
fvConstraints::New(this->mesh_).constrain(this->nut_);
```

This computes the eddy viscosity using the realizable coefficient and then applies boundary conditions and FV constraints.

nu_t = rCmu k^2 / epsilon

The overload that takes no arguments simply rebuilds:

grad U, S^2 = 2 |dev(symm(grad U))|^2, |S| = sqrt(S^2)

and forwards them to the main version.

## 5. `kSource()` and `epsilonSource()`

These functions currently return empty FV matrices, so at present:

S_k,extra = 0

S_epsilon,extra = 0

They exist as extension hooks for future model-specific sources.

## 6. Constructor

```cpp
A0_       = 4.0
C2_       = 1.9
sigmak_   = 1.0
sigmaEps_ = 1.2
curvatureSwirlTools_(this->mesh_, this->coeffDict_)
k_ and epsilon_ are MUST_READ/AUTO_WRITE fields
```

The constructor initializes realizable-model coefficients, constructs the helper tool, reads `k` and `epsilon` from the case, and applies the initial bounds:

k <- max(k, k_min)

epsilon <- max(epsilon, 0.09 k^2 / (nutMaxCoeff nu))

## 7. `read()`

```cpp
A0_.readIfPresent(this->coeffDict());
C2_.readIfPresent(this->coeffDict());
sigmak_.readIfPresent(this->coeffDict());
sigmaEps_.readIfPresent(this->coeffDict());
curvatureSwirlTools_.read(this->coeffDict());
```

This re-reads model coefficients and all helper settings at runtime. It changes parameters inside the equations, not the form of the equations themselves.

## 8. Start of `correct()`

```cpp
volScalarField::Internal divU
(
    typedName("divU"),
    fvc::div(fvc::absolute(this->phi(), U))()
);

const volTensorField gradU(fvc::grad(U));
const volScalarField S2(typedName("S2"), 2*magSqr(dev(symm(gradU))));
const volScalarField magS(typedName("magS"), sqrt(S2));

const volScalarField::Internal eta
(
    typedName("eta"), magS()*k_()/epsilon_()
);
const volScalarField::Internal C1
(
    typedName("C1"),
    max(eta/(scalar(5) + eta), scalar(0.43))
);

const volScalarField::Internal G
(
    this->GName(),
    nut*(gradU.v() && dev(twoSymm(gradU.v())))
);
```

This block prepares the local flow measures used in both the realizable model and the helper.

div U = grad . U

grad U = nabla U

S^2 = 2 |dev(symm(grad U))|^2

|S| = sqrt(S^2)

eta = |S| k / epsilon

C1 = max(eta / (5 + eta), 0.43)

G = nu_t [grad U : dev(2 symm(grad U))]

## 9. Coupling to `BuoyantCurvatureSwirlTools`

```cpp
tmp<volScalarField> tFrEff;
const volScalarField* frEffPtr = nullptr;
tmp<volScalarField> tGbCoef;
const volScalarField* GbCoefPtr = nullptr;
tmp<volScalarField> tBuoyLimiter;
const volScalarField* buoyLimiterPtr = nullptr;

volScalarField oneFrEff(... 1.0);
volScalarField zeroGbCoef(... 0.0);
volScalarField zeroBuoyLimiter(... 0.0);

if (curvatureSwirlTools_.active())
{
    curvatureSwirlData swirlData = curvatureSwirlTools_.evaluate
    (
        U,
        gradU,
        S2,
        magS,
        rho,
        k_,
        nut
    );
    tFrEff = swirlData.frEff;
    frEffPtr = &tFrEff();
    tGbCoef = swirlData.GbCoef;
    GbCoefPtr = &tGbCoef();
    tBuoyLimiter = swirlData.buoyLimiter;
    buoyLimiterPtr = &tBuoyLimiter();
}
```

Default fields are created so the model cleanly falls back to the base realizable form when corrections are off:

f_r,eff^default = 1

G_b^coef,default = 0

f_b^default = 0

If the helper is active, it returns the computed correction/source fields.

## 10. What the Helper Computes

The important helper outputs are:

S_local = |u_theta| / (|U_axial| + epsilon_swirl)

F_swirl = clip((S_local - S0) / (S1 - S0), 0, 1)

f_r = clip(|S| / max(|Omega|, epsilon_small), 0, f_r,max)

f_r,eff = 1 + F_swirl (f_r - 1)

G_b^coef = C_g nu_t (g . grad rho) / (k + k_small)

g-hat = g / |g|

v = g-hat . U

u = |U - g-hat v| + u_small

f_b = tanh(|v| / u)

## 11. Epsilon Equation

```cpp
epsilon_.boundaryFieldRef().updateCoeffs();

tmp<fvScalarMatrix> epsEqn
(
    fvm::ddt(alpha, rho, epsilon_)
  + fvm::div(alphaRhoPhi, epsilon_)
  - fvm::laplacian(alpha*rho*DepsilonEff(), epsilon_)
 ==
    C1*alpha()*rho()*magS()*epsilon_()
  - fvm::Sp
    (
        C2_*alpha()*rho()*epsilon_()/(k_() + sqrt(this->nu()()*epsilon_())),
        epsilon_
    )
  - fvm::SuSp
    (
        alpha*C1*(buoyLimiterPtr ? *buoyLimiterPtr : zeroBuoyLimiter)
       *(GbCoefPtr ? *GbCoefPtr : zeroGbCoef),
        epsilon_
    )
  + epsilonSource()
  + fvModels.source(alpha, rho, epsilon_)
);
```

The solved dissipation equation is approximately:

d(alpha rho epsilon)/dt + div(alpha rho U epsilon) - div(alpha rho D_epsilon^eff grad epsilon)

= C1 alpha rho |S| epsilon - C2 alpha rho epsilon^2 / (k + sqrt(nu epsilon)) - alpha C1 f_b G_b^coef epsilon + S_epsilon,extra

D_epsilon^eff = nu_t / sigma_epsilon + nu

After assembly, the matrix is relaxed, constrained, solved, and `epsilon` is bounded again.

## 12. k Equation

```cpp
tmp<fvScalarMatrix> kEqn
(
    fvm::ddt(alpha, rho, k_)
  + fvm::div(alphaRhoPhi, k_)
  - fvm::laplacian(alpha*rho*DkEff(), k_)
 ==
    alpha()*rho()*G*(frEffPtr ? *frEffPtr : oneFrEff)
  - fvm::SuSp(alpha*(GbCoefPtr ? *GbCoefPtr : zeroGbCoef), k_)
  - fvm::SuSp(2.0/3.0*alpha()*rho()*divU, k_)
  - fvm::Sp(alpha()*rho()*epsilon_()/k_(), k_)
  + kSource()
  + fvModels.source(alpha, rho, k_)
);

kEqn.ref().relax();
fvConstraints.constrain(kEqn.ref());
solve(kEqn);
fvConstraints.constrain(k_);
bound(k_, this->kMin_);
```

This is the transport equation for turbulence kinetic energy.

d(alpha rho k)/dt + div(alpha rho U k) - div(alpha rho D_k^eff grad k)

= alpha rho G f_r,eff - alpha G_b^coef k - (2/3) alpha rho (div U) k - alpha rho epsilon + S_k,extra

D_k^eff = nu_t / sigma_k + nu

The helper modifies this equation in two direct ways:

G -> G f_r,eff

- alpha G_b^coef k

with

f_r,eff = 1 + F_swirl (f_r - 1)

G_b^coef = C_g nu_t (g . grad rho) / (k + k_small)

After assembly, the matrix is relaxed, constrained, solved, and `k` is bounded from below.

## 13. BuoyantCurvatureSwirlTools Includes and Registration

```cpp
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
```

This block sets up the helper on OpenFOAM finite-volume infrastructure and gravity-gradient support. The equations it eventually supports are:

f_r,eff = 1 + F_swirl (f_r - 1)

G_b^coef = C_g nu_t (g . grad rho) / (k + k_small)

## 14. `normalisedAxis(...)`

```cpp
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
```

This converts the user axis into a unit vector:

e = axis / |axis|

with validity requirement:

|axis| > SMALL

## 15. `cellRelPosition()`

```cpp
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
```

This builds the relative position field:

rel = x - x0

where `x0 = axisOrigin`.

## 16. `axialVelocity(U)`

```cpp
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
```

This gives the scalar velocity along the chosen axis:

U_axial = U . e

## 17. `radialDistanceAndTangentialVelocity(...)`

```cpp
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
```

This computes in one pass:

rel = x - x0

v_r = rel - (rel . e) e

r = |v_r|

e_r = v_r / (r + epsilon)

e_theta = e x e_r

u_theta = U . e_theta

## 18. BuoyantCurvatureSwirlTools Constructor

This reads and initializes:

curvatureCorrection

swirlCorrection

buoyancyCorrection

axisOrigin = x0

axisDirection = e

S0, S1

frMax

swirlEps

Cg

with checks:

|e| > SMALL

S1 > S0

## 19. BuoyantCurvatureSwirlTools `read()`

This re-reads the same settings at runtime and updates the parameters inside the later equations:

S_local = |u_theta| / (|U_axial| + epsilon_swirl)

F_swirl = clip((S_local - S0) / (S1 - S0), 0, 1)

f_r = clip(|S| / max(|Omega|, epsilon_small), 0, frMax)

G_b^coef = C_g nu_t (g . grad rho) / (k + k_small)

## 20. `active()`

```cpp
return curvatureCorrection_ || swirlCorrection_ || buoyancyCorrection_;
```

Logical form:

active = curvatureCorrection OR swirlCorrection OR buoyancyCorrection

If all are off, the caller can use fallback fields:

f_r,eff = 1

G_b^coef = 0

f_b = 0

## 21. Start of `evaluate(...)`

The helper receives:

U

grad U

S^2 = 2 |dev(symm(grad U))|^2

|S| = sqrt(S^2)

rho

k

nu_t

and prepares default constants:

0, 1, zero velocity, zero buoyancy coefficient

It also decides lazily whether swirl and curvature intermediate fields need to be computed.

## 22. Swirl Fields in `evaluate()`

When needed, the helper computes:

r = |v_r|

u_theta = U . e_theta

U_axial = U . e

S_local = |u_theta| / (|U_axial| + epsilon_swirl)

F_swirl = clip((S_local - S0) / (S1 - S0), 0, 1)

If swirl correction is disabled, it forces:

S_local = 0

F_swirl = 0

## 23. Curvature Factor in `evaluate()`

The raw curvature factor is:

f_r = clip(|S| / max(|Omega|, epsilon_small), 0, frMax)

The effective factor is:

if curvature off:

f_r,eff = 1

if curvature on, swirl off:

f_r,eff = f_r

if curvature on, swirl on:

f_r,eff = 1 + F_swirl (f_r - 1)

## 24. Buoyancy Block in `evaluate()`

When buoyancy correction is on and gravity is nonzero:

g-hat = g / |g|

grad rho = grad(rho)

G_b^coef = C_g nu_t (g . grad rho) / (k + k_small)

v = g-hat . U

u = |U - g-hat v| + u_small

f_b = tanh(|v| / u)

If buoyancy is off:

G_b^coef = 0

f_b = 0

## 25. Final Writing and Return

If `writeFields` is enabled, the helper writes valid diagnostic fields and always writes:

frEff

GbCoef

buoyLimiter

Then it returns the `curvatureSwirlData` struct to `realizableKE_BCS`.

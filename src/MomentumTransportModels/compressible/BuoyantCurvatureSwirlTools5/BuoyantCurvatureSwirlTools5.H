#ifndef BuoyantCurvatureSwirlTools5_H
#define BuoyantCurvatureSwirlTools5_H

#include "volFields.H"
#include "fvMesh.H"
#include "dictionary.H"
#include "Switch.H"
#include "vector.H"
#include "tmp.H"

namespace Foam
{

struct curvatureSwirlData5
{
    tmp<volScalarField> uTheta;
    tmp<volScalarField> Uaxial;
    tmp<volScalarField> L;
    tmp<volScalarField> Phi;
    tmp<volScalarField> Fgate;
    tmp<volScalarField> fr;
    tmp<volScalarField> frEff;
    tmp<volScalarField> Gb;
    tmp<volScalarField> buoyLimiter;
};


// BuoyantCurvatureSwirlTools5 (BCST5)
//
// Generation-5 gate: replaces the BCS3/BCS4 magnitude-only swirl fraction
// F_swirl = |u_theta| / sqrt(u_theta^2 + u_z^2 + eps) with a signed Rayleigh
// (angular-momentum) stability discriminant:
//
//     L          = r * u_theta                     (specific angular momentum)
//     Phi(r)     = (1/r^3) * d(L^2)/dr              (Rayleigh discriminant)
//
//     Phi > 0  -> stabilizing   (solid-body-like core)
//     Phi ~ 0  -> marginal      (free vortex)
//     Phi < 0  -> destabilizing (near concave wall)
//
// The curvature-correction gate is destabilizing-only: it fires where the
// local rotation profile is Rayleigh-unstable (Phi < 0) and is fully
// suppressed where Phi >= 0 (stabilizing or marginal), normalized by the
// local mean-flow deformation scale D2 = max(S^2, 0.09*omegaLike^2):
//
//     Fgate = clamp( -Phi / (cPhi * D2), 0, 1 )
//
//     frEff = 1 + Fgate*(fr - 1)      when swirlCorrection = true
//     frEff = fr                      when swirlCorrection = false
//
// fr itself (the Spalart-Shur-style curvature/rotation correction factor)
// and the buoyancy-production machinery (Gb, buoyLimiter) are unchanged
// from BuoyantCurvatureSwirlTools3.
//
// Not implemented in this generation (flagged as open design choices in
// the source equation set, no fixed functional form given there):
//   - Section 4: circumferential-variance jet-merge applicability gate.
//   - Section 5: near-wall F1/F2 (Menter SST) attenuation of the gate.

class BuoyantCurvatureSwirlTools5
{
    // Private data

        const fvMesh& mesh_;

        Switch curvatureCorrection_;
        Switch swirlCorrection_;
        Switch buoyancyCorrection_;
        vector axisOrigin_;
        vector axisDirection_;

        dimensionedScalar frMax_;
        dimensionedScalar cCurv_;
        dimensionedScalar cr1_;
        dimensionedScalar cr2_;
        dimensionedScalar cr3_;
        dimensionedScalar rMin_;
        dimensionedScalar cPhi_;
        dimensionedScalar Cg_;

        Switch writeFields_;


    // Private Member Functions

        static vector normalisedAxis(const vector& axis);

        vectorField cellRelPosition() const;

        tmp<volScalarField> axialVelocity
        (
            const volVectorField& U
        ) const;

        tmp<volScalarField> tangentialVelocity
        (
            const volVectorField& U
        ) const;

        tmp<volScalarField> radialDistance() const;

        tmp<volVectorField> radialUnitVector() const;


public:

    TypeName("BuoyantCurvatureSwirlTools5");


    // Constructors

        BuoyantCurvatureSwirlTools5
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        virtual ~BuoyantCurvatureSwirlTools5() = default;


    // Member Functions

        bool active() const;

        curvatureSwirlData5 evaluate
        (
            const volVectorField& U,
            const volTensorField& gradU,
            const volScalarField& rho,
            const volScalarField& k,
            const volScalarField& omegaLike,
            const volScalarField& nut
        ) const;

        bool read(const dictionary& dict);
};


} // End namespace Foam

#endif

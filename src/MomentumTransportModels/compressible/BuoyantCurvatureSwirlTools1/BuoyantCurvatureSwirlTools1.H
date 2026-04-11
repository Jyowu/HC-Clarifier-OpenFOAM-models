#ifndef BuoyantCurvatureSwirlTools1_H
#define BuoyantCurvatureSwirlTools1_H

#include "volFields.H"
#include "fvMesh.H"
#include "dictionary.H"
#include "Switch.H"
#include "vector.H"
#include "tmp.H"

namespace Foam
{

struct curvatureSwirlData1
{
    tmp<volScalarField> r;
    tmp<volScalarField> uTheta;
    tmp<volScalarField> Uaxial;
    tmp<volScalarField> Slocal;
    tmp<volScalarField> Fswirl;
    tmp<volScalarField> fr;
    tmp<volScalarField> frEff;
    tmp<volScalarField> GbCoef;
    tmp<volScalarField> buoyLimiter;
};


class BuoyantCurvatureSwirlTools1
{
    // Private data

        const fvMesh& mesh_;
        const dictionary& dict_;

        Switch curvatureCorrection_;
        Switch swirlCorrection_;
        Switch buoyancyCorrection_;
        vector axisOrigin_;
        vector axisDirection_;

        dimensionedScalar S0_;
        dimensionedScalar S1_;
        dimensionedScalar frMax_;
        dimensionedScalar cCurv_;
        dimensionedScalar cr1_;
        dimensionedScalar cr2_;
        dimensionedScalar cr3_;
        dimensionedScalar swirlEps_;
        dimensionedScalar Cg_;

        Switch writeFields_;


    // Private Member Functions

        static vector normalisedAxis(const vector& axis);

        vectorField cellRelPosition() const;

        tmp<volScalarField> axialVelocity
        (
            const volVectorField& U
        ) const;

        void radialDistanceAndTangentialVelocity
        (
            const volVectorField& U,
            tmp<volScalarField>& tr,
            tmp<volScalarField>& tuTheta
        ) const;


public:

    TypeName("BuoyantCurvatureSwirlTools1");


    // Constructors

        BuoyantCurvatureSwirlTools1
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        virtual ~BuoyantCurvatureSwirlTools1() = default;


    // Member Functions

        bool active() const;

        curvatureSwirlData1 evaluate
        (
            const volVectorField& U,
            const volTensorField& gradU,
            const volScalarField& S2,
            const volScalarField& magS,
            const volScalarField& rho,
            const volScalarField& k,
            const volScalarField& epsilon,
            const volScalarField& nut
        ) const;

        bool read(const dictionary& dict);
};


} // End namespace Foam

#endif

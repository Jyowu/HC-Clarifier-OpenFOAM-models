#ifndef BuoyantCurvatureSwirlTools_H
#define BuoyantCurvatureSwirlTools_H

#include "volFields.H"
#include "fvMesh.H"
#include "dictionary.H"
#include "Switch.H"
#include "vector.H"
#include "tmp.H"

namespace Foam
{

struct curvatureSwirlData
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


class BuoyantCurvatureSwirlTools
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

    TypeName("BuoyantCurvatureSwirlTools");


    // Constructors

        BuoyantCurvatureSwirlTools
        (
            const fvMesh& mesh,
            const dictionary& dict
        );

        virtual ~BuoyantCurvatureSwirlTools() = default;


    // Member Functions

        bool active() const;

        curvatureSwirlData evaluate
        (
            const volVectorField& U,
            const volTensorField& gradU,
            const volScalarField& S2,
            const volScalarField& magS,
            const volScalarField& rho,
            const volScalarField& k,
            const volScalarField& nut
        ) const;

        bool read(const dictionary& dict);
};


} // End namespace Foam

#endif

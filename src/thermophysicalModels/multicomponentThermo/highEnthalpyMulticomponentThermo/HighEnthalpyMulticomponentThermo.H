#ifndef HighEnthalpyMulticomponentThermo_H
#define HighEnthalpyMulticomponentThermo_H

#include "fluidMulticomponentThermo.H"
#include "multicomponentThermo.H"
#include "psiThermo.H"
#include "basicThermo.H"
#include "fvMesh.H"
#include "zeroGradientFvPatchFields.H"

#ifdef Log
#undef Log
#endif

#include "mutation++.h"

#include <algorithm>
#include <exception>
#include <memory>
#include <vector>


namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class HighEnthalpyMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

template<class Thermo>
class HighEnthalpyMulticomponentThermo
:
    public FluidMulticomponentThermo<Thermo>
{
private:

    //- Relaxation time used where V-T exchange is negligible/unavailable [s]
    static constexpr scalar tauInf_ = 1e10;

    //- Vibro-electronic energy per unit mass [J/kg]. Conserved variable
    //  solved by the flow solver; Tve_ is derived from it in correct().
    volScalarField eve_;

    //- Equilibrium vibro-electronic energy e_ve(Ttr) [J/kg], source field
    //  for the Landau-Teller relaxation solved by the flow solver
    volScalarField eveEq_;

    //- Effective V-T relaxation time [s] such that
    //  Q_VT = rho*(eveEq - eve)/tauVT matches Mutation++'s source term
    volScalarField tauVT_;

    //- Ratio between the base-thermo mixture Cv (which links the solved
    //  energy to T through the polynomial inversion) and the trans-rotational
    //  Cv from Mutation++. Scales the V-T sink in the total energy equation
    //  so that dT/dt = -Q_VT/(rho*cv_tr) as in the two-temperature model.
    volScalarField cvTrRatio_;

    //- Optional Mutation++ bridge. If not available, fallback to base thermo.
    std::unique_ptr<Mutation::Mixture> mutationMixturePtr_;

    //- Mapping OpenFOAM specie index -> Mutation++ specie index
    List<label> foamToMutationSpecie_;

    word mutationMixtureName_;
    word mutationStateModel_;
    word mutationThermoDatabase_;
    word mutationMechanism_;

    scalar mutationMinT_;
    scalar mutationMaxT_;
    scalar mutationMinP_;
    Switch mutationForceThermalEquilibrium_;

    bool mutationAvailable_;

    //- Cached Mutation++ mixture data
    int nMutationSpecies_;
    int nMutationEnergyEqns_;
    std::vector<double> mutationMw_;

    //- Work buffers (avoid per-cell allocations)
    mutable std::vector<double> mutW1_;
    mutable std::vector<double> mutW2_;
    mutable std::vector<double> mutW3_;
    mutable std::vector<double> mutW4_;
    mutable std::vector<double> mutY_;
    mutable std::vector<double> mutQ_;
    mutable std::vector<double> mutCv_;

    //- Load Mutation++ options and build species mapping
    void initialiseMutation();

    //- Initialise eve from the current (T, Tve, Y) fields
    void initialiseEve();

    //- Fill mutY_ with normalised Mutation-ordered mass fractions taken
    //  from an accessor yFoam(speciei). Returns false if no mass is mapped.
    template<class YAccessor>
    inline bool mutationMassFractions(const YAccessor& yFoam) const;

    //- Mixture vibro-electronic energy [J/kg] at (Ttr, Tve) for the
    //  composition currently stored in mutY_. Pure function (no setState).
    inline scalar mixtureEve(const scalar T, const scalar Tv) const;

    //- Mixture vibro-electronic specific heat [J/kg/K] at (Ttr, Tve) for the
    //  composition currently stored in mutY_. Pure function (no setState).
    inline scalar mixtureCvVe(const scalar T, const scalar Tv) const;

    //- Invert Tve from a target eve [J/kg] with Newton iterations for the
    //  composition currently stored in mutY_
    inline scalar TveFromEve
    (
        const scalar T,
        const scalar eveTarget,
        const scalar TvGuess
    ) const;

    //- Update thermo/transport properties from Mutation++ for one state point
    inline bool updateFromMutation
    (
        const std::vector<double>& yMutation,
        const scalar p,
        scalar& T,
        scalar& Tv,
        scalar& Cp,
        scalar& Cv,
        scalar& psi,
        scalar& mu,
        scalar& kappa
    );

    //- Compute Landau-Teller linearisation (eveEq, tauVT) and the
    //  trans-rotational Cv at the state currently set in the mixture.
    //  Must be called right after a successful updateFromMutation().
    inline void relaxationAtCurrentState
    (
        const scalar T,
        const scalar Tv,
        scalar& eveNow,
        scalar& eveEqv,
        scalar& tau,
        scalar& cvTr
    ) const;

public:

    // Constructors
    HighEnthalpyMulticomponentThermo
    (
        const fvMesh& mesh,
        const word& phaseName = word::null
    )
    :
        FluidMulticomponentThermo<Thermo>(mesh, phaseName),
        eve_
        (
            IOobject
            (
                this->phasePropertyName("eve"),
                mesh.time().name(),
                mesh,
                IOobject::NO_READ,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimEnergy/dimMass, Zero),
            zeroGradientFvPatchScalarField::typeName
        ),
        eveEq_
        (
            IOobject
            (
                this->phasePropertyName("eveEq"),
                mesh.time().name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimEnergy/dimMass, Zero),
            zeroGradientFvPatchScalarField::typeName
        ),
        tauVT_
        (
            IOobject
            (
                this->phasePropertyName("tauVT"),
                mesh.time().name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimTime, scalar(tauInf_)),
            zeroGradientFvPatchScalarField::typeName
        ),
        cvTrRatio_
        (
            IOobject
            (
                this->phasePropertyName("cvTrRatio"),
                mesh.time().name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimless, 1),
            zeroGradientFvPatchScalarField::typeName
        ),
        mutationMixturePtr_(nullptr),
        foamToMutationSpecie_(),
        mutationMixtureName_("air_11"),
        mutationStateModel_("ChemNonEqTTv"),
        mutationThermoDatabase_("RRHO"),
        mutationMechanism_(word::null),
        mutationMinT_(scalar(200)),
        mutationMaxT_(scalar(50000)),
        mutationMinP_(scalar(1)),
        mutationForceThermalEquilibrium_(false),
        mutationAvailable_(false),
        nMutationSpecies_(0),
        nMutationEnergyEqns_(0)
    {
        // Fallback to translational temperature when Tve is not provided.
        scalarField& TveCells = this->Tve_.primitiveFieldRef();
        const scalarField& TCells = this->T_.primitiveField();
        forAll(TveCells, celli)
        {
            if (TveCells[celli] <= SMALL)
            {
                TveCells[celli] = TCells[celli];
            }
        }

        volScalarField::Boundary& TveBf = this->Tve_.boundaryFieldRef();
        const volScalarField::Boundary& TBf = this->T_.boundaryField();
        forAll(TveBf, patchi)
        {
            forAll(TveBf[patchi], facei)
            {
                if (TveBf[patchi][facei] <= SMALL)
                {
                    TveBf[patchi][facei] = TBf[patchi][facei];
                }
            }
        }

        initialiseMutation();
        initialiseEve();
    }

    virtual ~HighEnthalpyMulticomponentThermo() = default;

    //- Access to vibro-electronic energy [J/kg]
    const volScalarField& eve() const
    {
        return eve_;
    }

    //- Non-const access to vibro-electronic energy [J/kg]
    volScalarField& eve()
    {
        return eve_;
    }

    // Member Functions
    virtual void correct() override
    {
        // First execute the default OpenFOAM update to preserve robustness.
        FluidMulticomponentThermo<Thermo>::correct();

        if (!mutationAvailable_ || !mutationMixturePtr_)
        {
            return;
        }

        // Snapshot of the base-thermo mixture Cv before the bridge overwrites
        // it below: it links the solved energy to T through the polynomial
        // inversion and is needed for the cvTrRatio source scaling.
        const volScalarField cvPoly("cvPoly", this->Cv_);

        const PtrList<volScalarField>& Y = this->Y();

        const scalarField& pCells = this->p_.primitiveField();
        scalarField& TCells = this->T_.primitiveFieldRef();
        scalarField& TveCells = this->Tve_.primitiveFieldRef();
        const scalarField& eveCells = eve_.primitiveField();
        scalarField& eveEqCells = eveEq_.primitiveFieldRef();
        scalarField& tauCells = tauVT_.primitiveFieldRef();
        scalarField& ratioCells = cvTrRatio_.primitiveFieldRef();
        scalarField& CpCells = this->Cp_.primitiveFieldRef();
        scalarField& CvCells = this->Cv_.primitiveFieldRef();
        scalarField& psiCells = this->psi_.primitiveFieldRef();
        scalarField& muCells = this->mu_.primitiveFieldRef();
        scalarField& kappaCells = this->kappa_.primitiveFieldRef();

        bool hadFailure = false;

        forAll(TCells, celli)
        {
            const auto yFoam = [&](const label speciei)
            {
                return Y[speciei][celli];
            };

            mutationMassFractions(yFoam);

            scalar T = max(TCells[celli], mutationMinT_);

            // Derive Tve from the conserved eve field
            scalar Tv = mutationForceThermalEquilibrium_
                ? T
                : TveFromEve
                  (
                      T,
                      eveCells[celli],
                      max(TveCells[celli], mutationMinT_)
                  );

            scalar Cp = CpCells[celli];
            scalar Cv = CvCells[celli];
            scalar psi = psiCells[celli];
            scalar mu = muCells[celli];
            scalar kappa = kappaCells[celli];

            if
            (
                updateFromMutation
                (
                    mutY_,
                    max(pCells[celli], mutationMinP_),
                    T,
                    Tv,
                    Cp,
                    Cv,
                    psi,
                    mu,
                    kappa
                )
            )
            {
                if (mutationForceThermalEquilibrium_)
                {
                    Tv = T;
                }

                TCells[celli] = T;
                TveCells[celli] = Tv;
                CpCells[celli] = Cp;
                CvCells[celli] = Cv;
                psiCells[celli] = psi;
                muCells[celli] = mu;
                kappaCells[celli] = kappa;

                scalar eveNow, eveEqv, tau, cvTr;
                relaxationAtCurrentState(T, Tv, eveNow, eveEqv, tau, cvTr);

                eveEqCells[celli] = eveEqv;
                tauCells[celli] = tau;
                ratioCells[celli] =
                    cvTr > SMALL ? max(cvPoly[celli]/cvTr, scalar(1)) : 1;
            }
            else
            {
                // Disable relaxation where the bridge failed
                eveEqCells[celli] = eveCells[celli];
                tauCells[celli] = tauInf_;
                ratioCells[celli] = 1;
                hadFailure = true;
            }
        }

        const volScalarField::Boundary& pBf = this->p_.boundaryField();
        volScalarField::Boundary& TBf = this->T_.boundaryFieldRef();
        volScalarField::Boundary& TveBf = this->Tve_.boundaryFieldRef();
        const volScalarField::Boundary& eveBf = eve_.boundaryField();
        volScalarField::Boundary& eveEqBf = eveEq_.boundaryFieldRef();
        volScalarField::Boundary& tauBf = tauVT_.boundaryFieldRef();
        volScalarField::Boundary& ratioBf = cvTrRatio_.boundaryFieldRef();
        volScalarField::Boundary& CpBf = this->Cp_.boundaryFieldRef();
        volScalarField::Boundary& CvBf = this->Cv_.boundaryFieldRef();
        volScalarField::Boundary& psiBf = this->psi_.boundaryFieldRef();
        volScalarField::Boundary& muBf = this->mu_.boundaryFieldRef();
        volScalarField::Boundary& kappaBf = this->kappa_.boundaryFieldRef();

        forAll(TBf, patchi)
        {
            const fvPatchScalarField& pp = pBf[patchi];
            fvPatchScalarField& pT = TBf[patchi];
            fvPatchScalarField& pTv = TveBf[patchi];
            const fvPatchScalarField& pEve = eveBf[patchi];
            fvPatchScalarField& pEveEq = eveEqBf[patchi];
            fvPatchScalarField& pTau = tauBf[patchi];
            fvPatchScalarField& pRatio = ratioBf[patchi];
            fvPatchScalarField& pCp = CpBf[patchi];
            fvPatchScalarField& pCv = CvBf[patchi];
            fvPatchScalarField& pPsi = psiBf[patchi];
            fvPatchScalarField& pMu = muBf[patchi];
            fvPatchScalarField& pKappa = kappaBf[patchi];

            forAll(pT, facei)
            {
                const auto yFoam = [&](const label speciei)
                {
                    return Y[speciei].boundaryField()[patchi][facei];
                };

                mutationMassFractions(yFoam);

                scalar T = max(pT[facei], mutationMinT_);

                scalar Tv = mutationForceThermalEquilibrium_
                    ? T
                    : TveFromEve
                      (
                          T,
                          pEve[facei],
                          max(pTv[facei], mutationMinT_)
                      );

                scalar Cp = pCp[facei];
                scalar Cv = pCv[facei];
                scalar psi = pPsi[facei];
                scalar mu = pMu[facei];
                scalar kappa = pKappa[facei];

                if
                (
                    updateFromMutation
                    (
                        mutY_,
                        max(pp[facei], mutationMinP_),
                        T,
                        Tv,
                        Cp,
                        Cv,
                        psi,
                        mu,
                        kappa
                    )
                )
                {
                    if (mutationForceThermalEquilibrium_)
                    {
                        Tv = T;
                    }

                    pT[facei] = T;
                    pTv[facei] = Tv;
                    pCp[facei] = Cp;
                    pCv[facei] = Cv;
                    pPsi[facei] = psi;
                    pMu[facei] = mu;
                    pKappa[facei] = kappa;

                    scalar eveNow, eveEqv, tau, cvTr;
                    relaxationAtCurrentState(T, Tv, eveNow, eveEqv, tau, cvTr);

                    pEveEq[facei] = eveEqv;
                    pTau[facei] = tau;
                    pRatio[facei] =
                        cvTr > SMALL
                      ? max(cvPoly.boundaryField()[patchi][facei]/cvTr, scalar(1))
                      : 1;
                }
                else
                {
                    pEveEq[facei] = pEve[facei];
                    pTau[facei] = tauInf_;
                    pRatio[facei] = 1;
                    hadFailure = true;
                }
            }
        }

        if (hadFailure)
        {
            WarningInFunction
                << "Mutation++ update failed for at least one location. "
                << "Base thermo values were kept for failed entries."
                << endl;
        }
    }

};


template<class Thermo>
void HighEnthalpyMulticomponentThermo<Thermo>::initialiseMutation()
{
    const dictionary& propertiesDict = this->properties();

    if (propertiesDict.found("highEnthalpyMutation"))
    {
        const dictionary& mutationDict =
            propertiesDict.subDict("highEnthalpyMutation");

        mutationMixtureName_ =
            mutationDict.lookupOrDefault<word>("mixture", mutationMixtureName_);
        mutationStateModel_ =
            mutationDict.lookupOrDefault<word>("stateModel", mutationStateModel_);
        mutationThermoDatabase_ =
            mutationDict.lookupOrDefault<word>
            (
                "thermodynamicDatabase",
                mutationThermoDatabase_
            );
        mutationMechanism_ =
            mutationDict.lookupOrDefault<word>("mechanism", mutationMechanism_);
        mutationMinT_ =
            mutationDict.lookupOrDefault<scalar>("minTemperature", mutationMinT_);
        mutationMaxT_ =
            mutationDict.lookupOrDefault<scalar>("maxTemperature", mutationMaxT_);
        mutationMinP_ =
            mutationDict.lookupOrDefault<scalar>("minPressure", mutationMinP_);
        mutationForceThermalEquilibrium_ =
            mutationDict.lookupOrDefault<Switch>
            (
                "forceThermalEquilibrium",
                mutationForceThermalEquilibrium_
            );
    }

    if (mutationStateModel_.find("Equil") == 0)
    {
        WarningInFunction
            << "State model '" << mutationStateModel_
            << "' is not supported in highEnthalpyMutation bridge. "
            << "Falling back to base thermo."
            << endl;
        mutationAvailable_ = false;
        return;
    }

    try
    {
        Mutation::MixtureOptions opts(mutationMixtureName_.c_str());
        opts.setStateModel(mutationStateModel_.c_str());
        opts.setThermodynamicDatabase(mutationThermoDatabase_.c_str());

        if (!mutationMechanism_.empty())
        {
            opts.setMechanism(mutationMechanism_.c_str());
        }

        mutationMixturePtr_.reset(new Mutation::Mixture(opts));

        nMutationSpecies_ = mutationMixturePtr_->nSpecies();
        nMutationEnergyEqns_ = mutationMixturePtr_->nEnergyEqns();

        mutationMw_.resize(nMutationSpecies_);
        for (int i = 0; i < nMutationSpecies_; ++i)
        {
            mutationMw_[i] = mutationMixturePtr_->speciesMw(i);
        }

        mutW1_.resize(nMutationSpecies_);
        mutW2_.resize(nMutationSpecies_);
        mutW3_.resize(nMutationSpecies_);
        mutW4_.resize(nMutationSpecies_);
        mutY_.resize(nMutationSpecies_);
        mutQ_.resize(std::max(nMutationEnergyEqns_, 1));
        mutCv_.resize(nMutationSpecies_*nMutationEnergyEqns_);

        const speciesTable& species = this->species();
        foamToMutationSpecie_.setSize(species.size(), -1);

        forAll(species, speciei)
        {
            const label mutationSpeciei =
                mutationMixturePtr_->speciesIndex(species[speciei].c_str());

            foamToMutationSpecie_[speciei] = mutationSpeciei;

            if (mutationSpeciei < 0)
            {
                WarningInFunction
                    << "Specie '" << species[speciei]
                    << "' not present in Mutation++ mixture '"
                    << mutationMixtureName_ << "'." << endl;
            }
        }

        mutationAvailable_ = true;
    }
    catch (const std::exception& e)
    {
        WarningInFunction
            << "Unable to initialize Mutation++ mixture '"
            << mutationMixtureName_ << "': " << e.what() << nl
            << "Falling back to base thermo."
            << endl;
        mutationMixturePtr_.reset(nullptr);
        mutationAvailable_ = false;
    }
    catch (...)
    {
        WarningInFunction
            << "Unknown error during Mutation++ initialization. "
            << "Falling back to base thermo."
            << endl;
        mutationMixturePtr_.reset(nullptr);
        mutationAvailable_ = false;
    }
}


template<class Thermo>
void HighEnthalpyMulticomponentThermo<Thermo>::initialiseEve()
{
    if (!mutationAvailable_ || !mutationMixturePtr_)
    {
        return;
    }

    const PtrList<volScalarField>& Y = this->Y();

    const scalarField& TCells = this->T_.primitiveField();
    const scalarField& TveCells = this->Tve_.primitiveField();
    scalarField& eveCells = eve_.primitiveFieldRef();

    forAll(eveCells, celli)
    {
        const auto yFoam = [&](const label speciei)
        {
            return Y[speciei][celli];
        };

        mutationMassFractions(yFoam);

        eveCells[celli] = mixtureEve
        (
            max(TCells[celli], mutationMinT_),
            max(TveCells[celli], mutationMinT_)
        );
    }

    volScalarField::Boundary& eveBf = eve_.boundaryFieldRef();
    const volScalarField::Boundary& TBf = this->T_.boundaryField();
    const volScalarField::Boundary& TveBf = this->Tve_.boundaryField();

    forAll(eveBf, patchi)
    {
        forAll(eveBf[patchi], facei)
        {
            const auto yFoam = [&](const label speciei)
            {
                return Y[speciei].boundaryField()[patchi][facei];
            };

            mutationMassFractions(yFoam);

            eveBf[patchi][facei] = mixtureEve
            (
                max(TBf[patchi][facei], mutationMinT_),
                max(TveBf[patchi][facei], mutationMinT_)
            );
        }
    }

    eveEq_.primitiveFieldRef() = eve_.primitiveField();
    eveEq_.boundaryFieldRef() = eve_.boundaryField();
}


template<class Thermo>
template<class YAccessor>
inline bool HighEnthalpyMulticomponentThermo<Thermo>::mutationMassFractions
(
    const YAccessor& yFoam
) const
{
    std::fill(mutY_.begin(), mutY_.end(), 1e-8);

    const label nFoamSpecies = foamToMutationSpecie_.size();

    double ySum = 0.0;
    for (label speciei = 0; speciei < nFoamSpecies; ++speciei)
    {
        const label mutationSpeciei = foamToMutationSpecie_[speciei];
        if (mutationSpeciei >= 0)
        {
            const double y = max(yFoam(speciei), scalar(0));
            mutY_[mutationSpeciei] += y;
            ySum += y;
        }
    }

    if (ySum <= SMALL)
    {
        mutY_[0] = 1.0;
        return false;
    }

    for (int i = 0; i < nMutationSpecies_; ++i)
    {
        mutY_[i] /= ySum;
    }

    return true;
}


template<class Thermo>
inline Foam::scalar HighEnthalpyMulticomponentThermo<Thermo>::mixtureEve
(
    const scalar T,
    const scalar Tv
) const
{
    // Same temperature layout as Mutation++'s ChemNonEqTTv state model:
    // (Th, Te, Tr, Tv, Tel) = (T, Tv, T, Tv, Tv); h is normalised by Th.
    mutationMixturePtr_->speciesHOverRT
    (
        T, Tv, T, Tv, Tv,
        mutW1_.data(),
        nullptr,
        nullptr,
        mutW3_.data(),
        mutW4_.data(),
        nullptr
    );

    double e = 0.0;
    for (int i = 0; i < nMutationSpecies_; ++i)
    {
        e += mutY_[i]*(mutW3_[i] + mutW4_[i])*Mutation::RU*T/mutationMw_[i];
    }

    return e;
}


template<class Thermo>
inline Foam::scalar HighEnthalpyMulticomponentThermo<Thermo>::mixtureCvVe
(
    const scalar T,
    const scalar Tv
) const
{
    mutationMixturePtr_->speciesCpOverR
    (
        T, Tv, T, Tv, Tv,
        nullptr,
        mutW1_.data(),
        mutW2_.data(),
        mutW3_.data(),
        mutW4_.data()
    );

    double cv = 0.0;
    for (int i = 0; i < nMutationSpecies_; ++i)
    {
        cv += mutY_[i]*(mutW3_[i] + mutW4_[i])*Mutation::RU/mutationMw_[i];
    }

    return cv;
}


template<class Thermo>
inline Foam::scalar HighEnthalpyMulticomponentThermo<Thermo>::TveFromEve
(
    const scalar T,
    const scalar eveTarget,
    const scalar TvGuess
) const
{
    scalar Tv = min(max(TvGuess, mutationMinT_), mutationMaxT_);

    for (int iter = 0; iter < 50; ++iter)
    {
        const scalar eveNow = mixtureEve(T, Tv);
        const scalar cvVe = max(mixtureCvVe(T, Tv), scalar(1e-8));

        scalar dTv = (eveTarget - eveNow)/cvVe;

        // Damp large Newton steps to remain within the bracket
        dTv = min(max(dTv, -0.5*Tv), 0.5*Tv);

        Tv = min(max(Tv + dTv, mutationMinT_), mutationMaxT_);

        if (mag(dTv) < 1e-4)
        {
            break;
        }
    }

    return Tv;
}


template<class Thermo>
inline bool HighEnthalpyMulticomponentThermo<Thermo>::updateFromMutation
(
    const std::vector<double>& yMutation,
    const scalar p,
    scalar& T,
    scalar& Tv,
    scalar& Cp,
    scalar& Cv,
    scalar& psi,
    scalar& mu,
    scalar& kappa
)
{
    if (!mutationMixturePtr_)
    {
        return false;
    }

    const double state[3] =
    {
        max(p, mutationMinP_),
        max(T, mutationMinT_),
        max(Tv, mutationMinT_)
    };

    try
    {
        mutationMixturePtr_->setState(yMutation.data(), state, 2);

        T = max(mutationMixturePtr_->T(), mutationMinT_);
        Tv = max(mutationMixturePtr_->Tv(), mutationMinT_);

        Cp = max(mutationMixturePtr_->mixtureFrozenCpMass(), SMALL);
        Cv = max(mutationMixturePtr_->mixtureFrozenCvMass(), SMALL);

        mu = max(mutationMixturePtr_->viscosity(), SMALL);
        kappa = max(mutationMixturePtr_->frozenThermalConductivity(), SMALL);

        const scalar rho = max(mutationMixturePtr_->density(), SMALL);
        psi = max(rho/max(p, mutationMinP_), SMALL);
    }
    catch (...)
    {
        return false;
    }

    return true;
}


template<class Thermo>
inline void HighEnthalpyMulticomponentThermo<Thermo>::relaxationAtCurrentState
(
    const scalar T,
    const scalar Tv,
    scalar& eveNow,
    scalar& eveEqv,
    scalar& tau,
    scalar& cvTr
) const
{
    eveNow = mixtureEve(T, Tv);
    eveEqv = mixtureEve(T, T);

    // Trans-rotational Cv of the mixture [J/kg/K] (first block of getCvsMass)
    mutationMixturePtr_->getCvsMass(mutCv_.data());

    cvTr = 0.0;
    for (int i = 0; i < nMutationSpecies_; ++i)
    {
        cvTr += mutY_[i]*mutCv_[i];
    }

    // Mutation++ V-T (and other ve-pool) energy transfer source [W/m^3];
    // recast into the Landau-Teller form Q = rho*(eveEq - eve)/tau
    tau = tauInf_;

    try
    {
        mutationMixturePtr_->energyTransferSource(mutQ_.data());

        const double Q = mutQ_[0];
        const double dE = eveEqv - eveNow;
        const double rho = mutationMixturePtr_->density();

        if (mag(Q) > SMALL && dE*Q > 0)
        {
            tau = min(max(rho*dE/Q, scalar(1e-12)), scalar(tauInf_));
        }
    }
    catch (...)
    {
        tau = tauInf_;
    }
}

} // End namespace Foam

#endif

// ************************************************************************* //

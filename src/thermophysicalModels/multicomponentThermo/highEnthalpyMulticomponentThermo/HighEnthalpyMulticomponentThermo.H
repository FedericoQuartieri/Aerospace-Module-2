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

    Two-temperature thermodynamics bridged to Mutation++.

    CONSERVATIVE FORMULATION (milestone 3): the solved sensible energy "e"
    uses the Mutation++ datum e = e_tr(T) + e_ve(Tve) (formation enthalpies
    excluded, exactly as in the 0D reacting tests).  Temperatures are decoded
    from the conserved variables:
        Tve  from  eve                 (Newton on e_ve)
        T    from  e - eve = e_tr(T)   (Newton on e_tr; linear for RRHO)
    V-T exchange moves energy between the tr and ve pools inside "e", so no
    relaxation term appears in the total energy equation; chemistry adds
        Qdot = -sum_i hf_i * wdot_i        to the e equation, and
        Qcv  =  sum_i e_ve,i * wdot_i      to the eve equation (Candler).
    The V-T relaxation time is recast from Mutation++'s
    energyTransferSource() after subtracting the chemistry-vibration part.
\*---------------------------------------------------------------------------*/

template<class Thermo>
class HighEnthalpyMulticomponentThermo
:
    public FluidMulticomponentThermo<Thermo>
{
private:

    //- Relaxation time used where V-T exchange is negligible/unavailable [s]
    static constexpr scalar tauInf_ = 1e10;

    //- Vibrational-electronic temperature field [K]
    volScalarField Tve_;

    //- Vibro-electronic energy per unit mass [J/kg]. Conserved variable
    //  solved by the flow solver; Tve_ is derived from it in correct().
    volScalarField eve_;

    //- Equilibrium vibro-electronic energy e_ve(Ttr) [J/kg], source field
    //  for the Landau-Teller relaxation solved by the flow solver
    volScalarField eveEq_;

    //- Effective V-T relaxation time [s] such that
    //  Q_VT = rho*(eveEq - eve)/tauVT matches Mutation++'s pure V-T source
    volScalarField tauVT_;

    //- Chemistry heat release for the sensible energy equation [W/m^3]:
    //  Qdot = -sum_i hf_i*wdot_i (formation enthalpies excluded from "e")
    volScalarField::Internal Qdot_;

    //- Chemistry-vibration coupling for the eve equation [W/m^3]:
    //  Qcv = sum_i (e_v,i + e_el,i)(Tve)*wdot_i (non-preferential Candler,
    //  mirrors Mutation++'s OmegaCV + OmegaCElec)
    volScalarField::Internal Qcv_;

    //- Net mass production rate per FOAM specie [kg/m^3/s], registered as
    //  "mutR_<specie>" for use as source in the specie transport equations
    PtrList<volScalarField::Internal> RY_;

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

    //- Formation enthalpy per unit mass of each Mutation++ specie [J/kg]
    std::vector<double> mutationHfMass_;

    //- Work buffers (avoid per-cell allocations)
    mutable std::vector<double> mutW1_;
    mutable std::vector<double> mutW2_;
    mutable std::vector<double> mutW3_;
    mutable std::vector<double> mutW4_;
    mutable std::vector<double> mutY_;
    mutable std::vector<double> mutQ_;
    mutable std::vector<double> mutWdot_;

    //- Load Mutation++ options and build species mapping
    void initialiseMutation();

    //- Initialise eve and re-initialise the solved energy "e" on the
    //  Mutation++ datum e = e_tr(T) + e_ve(Tve) from the (T, Tve, Y) fields
    void initialiseEnergies();

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

    //- Mixture trans-rotational sensible energy [J/kg] at Ttr for the
    //  composition currently stored in mutY_. Pure function (no setState).
    inline scalar mixtureEtr(const scalar T) const;

    //- Mixture trans-rotational specific heat [J/kg/K] at Ttr for the
    //  composition currently stored in mutY_. Pure function (no setState).
    inline scalar mixtureCvTr(const scalar T) const;

    //- Invert Tve from a target eve [J/kg] with Newton iterations for the
    //  composition currently stored in mutY_
    inline scalar TveFromEve
    (
        const scalar T,
        const scalar eveTarget,
        const scalar TvGuess
    ) const;

    //- Invert T from a target trans-rotational energy [J/kg] with Newton
    //  iterations for the composition currently stored in mutY_
    inline scalar TFromEtr
    (
        const scalar etrTarget,
        const scalar TGuess
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

    //- Compute the relaxation linearisation (eveEq, tauVT) and the chemistry
    //  sources (wdot in mutWdot_, Qcv, Qdot) at the state currently set in
    //  the mixture. Must be called right after a successful
    //  updateFromMutation(). tauVT is recast from the pure V-T part of
    //  energyTransferSource (total minus Candler CV+CElec).
    inline void sourcesAtCurrentState
    (
        const scalar T,
        const scalar Tv,
        scalar& eveNow,
        scalar& eveEqv,
        scalar& tau,
        scalar& qcv,
        scalar& qdot
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
        Tve_
        (
            IOobject
            (
                this->phasePropertyName("Tve"),
                mesh.time().name(),
                mesh,
                IOobject::READ_IF_PRESENT,
                IOobject::AUTO_WRITE
            ),
            mesh,
            dimensionedScalar(dimTemperature, Zero)
        ),
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
        Qdot_
        (
            IOobject
            (
                this->phasePropertyName("mutQdot"),
                mesh.time().name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
        ),
        Qcv_
        (
            IOobject
            (
                this->phasePropertyName("mutQcv"),
                mesh.time().name(),
                mesh,
                IOobject::NO_READ,
                IOobject::NO_WRITE
            ),
            mesh,
            dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
        ),
        RY_(),
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
        scalarField& TveCells = Tve_.primitiveFieldRef();
        const scalarField& TCells = this->T_.primitiveField();
        forAll(TveCells, celli)
        {
            if (TveCells[celli] <= SMALL)
            {
                TveCells[celli] = TCells[celli];
            }
        }

        volScalarField::Boundary& TveBf = Tve_.boundaryFieldRef();
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

        // One production-rate field per FOAM specie, named "mutR_<specie>"
        const speciesTable& species = this->species();
        RY_.setSize(species.size());
        forAll(species, speciei)
        {
            RY_.set
            (
                speciei,
                new volScalarField::Internal
                (
                    IOobject
                    (
                        this->phasePropertyName("mutR_" + species[speciei]),
                        mesh.time().name(),
                        mesh,
                        IOobject::NO_READ,
                        IOobject::NO_WRITE
                    ),
                    mesh,
                    dimensionedScalar(dimMass/dimVolume/dimTime, Zero)
                )
            );
        }

        initialiseMutation();
        initialiseEnergies();
    }

    virtual ~HighEnthalpyMulticomponentThermo() = default;

    //- Access to vibrational-electronic temperature [K]
    const volScalarField& Tve() const
    {
        return Tve_;
    }

    //- Non-const access to vibrational-electronic temperature [K]
    volScalarField& Tve()
    {
        return Tve_;
    }

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

    //- Bring the remaining he() overloads (field accessors, (p,T) forms)
    //  into scope: the overrides below would otherwise hide them
    using FluidMulticomponentThermo<Thermo>::he;

    //- Enthalpy/Internal energy for patch [J/kg] on the Mutation++ datum
    //  e = e_tr(T) + eve. Energy boundary conditions (fixedEnergy,
    //  mixedEnergy, ...) evaluate the energy from the temperature BC via
    //  this call: without the override they would use the base-thermo
    //  polynomial datum (~3e5 J/kg off for N2), breaking fixed-temperature
    //  walls under the conservative decode T(e - eve).
    virtual tmp<scalarField> he
    (
        const scalarField& T,
        const label patchi
    ) const override
    {
        if (!mutationAvailable_ || !mutationMixturePtr_)
        {
            return FluidMulticomponentThermo<Thermo>::he(T, patchi);
        }

        const PtrList<volScalarField>& Y = this->Y();
        const fvPatchScalarField& pEve = eve_.boundaryField()[patchi];

        tmp<scalarField> the(new scalarField(T.size()));
        scalarField& heP = the.ref();

        forAll(heP, facei)
        {
            const auto yFoam = [&](const label speciei)
            {
                return Y[speciei].boundaryField()[patchi][facei];
            };

            mutationMassFractions(yFoam);

            heP[facei] =
                mixtureEtr(max(T[facei], mutationMinT_)) + pEve[facei];
        }

        return the;
    }

    //- Enthalpy/Internal energy for cell-set [J/kg] on the Mutation++ datum
    virtual tmp<scalarField> he
    (
        const scalarField& T,
        const labelList& cells
    ) const override
    {
        if (!mutationAvailable_ || !mutationMixturePtr_)
        {
            return FluidMulticomponentThermo<Thermo>::he(T, cells);
        }

        const PtrList<volScalarField>& Y = this->Y();

        tmp<scalarField> the(new scalarField(T.size()));
        scalarField& heC = the.ref();

        forAll(heC, i)
        {
            const label celli = cells[i];

            const auto yFoam = [&](const label speciei)
            {
                return Y[speciei][celli];
            };

            mutationMassFractions(yFoam);

            heC[i] = mixtureEtr(max(T[i], mutationMinT_)) + eve_[celli];
        }

        return the;
    }

    // Member Functions
    virtual void correct() override
    {
        if (!mutationAvailable_ || !mutationMixturePtr_)
        {
            // No bridge: plain OpenFOAM update
            FluidMulticomponentThermo<Thermo>::correct();
            return;
        }

        // NOTE: the base correct() is intentionally NOT called here: the
        // solved energy uses the Mutation++ datum (see initialiseEnergies),
        // so the base polynomial inversion of T from "e" would be
        // inconsistent. All fields the base update would touch are set
        // below from the conservative decode.

        const PtrList<volScalarField>& Y = this->Y();

        const scalarField& pCells = this->p_.primitiveField();
        const scalarField& eCells = this->he().primitiveField();
        scalarField& TCells = this->T_.primitiveFieldRef();
        scalarField& TveCells = Tve_.primitiveFieldRef();
        const scalarField& eveCells = eve_.primitiveField();
        scalarField& eveEqCells = eveEq_.primitiveFieldRef();
        scalarField& tauCells = tauVT_.primitiveFieldRef();
        // DimensionedField publicly derives from Field: direct indexing
        scalarField& QdotCells = Qdot_;
        scalarField& QcvCells = Qcv_;
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

            // --- conservative decode: Tve from eve, T from e - eve ---
            scalar Tv = mutationForceThermalEquilibrium_
                ? max(TCells[celli], mutationMinT_)
                : TveFromEve
                  (
                      max(TCells[celli], mutationMinT_),
                      eveCells[celli],
                      max(TveCells[celli], mutationMinT_)
                  );

            scalar T = TFromEtr
            (
                eCells[celli] - eveCells[celli],
                max(TCells[celli], mutationMinT_)
            );

            if (mutationForceThermalEquilibrium_)
            {
                Tv = T;
            }

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
                TCells[celli] = T;
                TveCells[celli] = Tv;
                CpCells[celli] = Cp;
                CvCells[celli] = Cv;
                psiCells[celli] = psi;
                muCells[celli] = mu;
                kappaCells[celli] = kappa;

                scalar eveNow, eveEqv, tau, qcv, qdot;
                sourcesAtCurrentState(T, Tv, eveNow, eveEqv, tau, qcv, qdot);

                eveEqCells[celli] = eveEqv;
                tauCells[celli] = tau;
                QcvCells[celli] = qcv;
                QdotCells[celli] = qdot;

                forAll(RY_, speciei)
                {
                    const label mutationSpeciei =
                        foamToMutationSpecie_[speciei];
                    RY_[speciei][celli] =
                        mutationSpeciei >= 0 ? mutWdot_[mutationSpeciei] : 0;
                }
            }
            else
            {
                // Disable relaxation and chemistry where the bridge failed
                eveEqCells[celli] = eveCells[celli];
                tauCells[celli] = tauInf_;
                QcvCells[celli] = 0;
                QdotCells[celli] = 0;

                forAll(RY_, speciei)
                {
                    RY_[speciei][celli] = 0;
                }

                hadFailure = true;
            }
        }

        const volScalarField::Boundary& pBf = this->p_.boundaryField();
        const volScalarField::Boundary& eBf = this->he().boundaryField();
        volScalarField::Boundary& TBf = this->T_.boundaryFieldRef();
        volScalarField::Boundary& TveBf = Tve_.boundaryFieldRef();
        const volScalarField::Boundary& eveBf = eve_.boundaryField();
        volScalarField::Boundary& eveEqBf = eveEq_.boundaryFieldRef();
        volScalarField::Boundary& tauBf = tauVT_.boundaryFieldRef();
        volScalarField::Boundary& CpBf = this->Cp_.boundaryFieldRef();
        volScalarField::Boundary& CvBf = this->Cv_.boundaryFieldRef();
        volScalarField::Boundary& psiBf = this->psi_.boundaryFieldRef();
        volScalarField::Boundary& muBf = this->mu_.boundaryFieldRef();
        volScalarField::Boundary& kappaBf = this->kappa_.boundaryFieldRef();

        forAll(TBf, patchi)
        {
            const fvPatchScalarField& pp = pBf[patchi];
            const fvPatchScalarField& pE = eBf[patchi];
            fvPatchScalarField& pT = TBf[patchi];
            fvPatchScalarField& pTv = TveBf[patchi];
            const fvPatchScalarField& pEve = eveBf[patchi];
            fvPatchScalarField& pEveEq = eveEqBf[patchi];
            fvPatchScalarField& pTau = tauBf[patchi];
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

                scalar Tv = mutationForceThermalEquilibrium_
                    ? max(pT[facei], mutationMinT_)
                    : TveFromEve
                      (
                          max(pT[facei], mutationMinT_),
                          pEve[facei],
                          max(pTv[facei], mutationMinT_)
                      );

                scalar T = TFromEtr
                (
                    pE[facei] - pEve[facei],
                    max(pT[facei], mutationMinT_)
                );

                if (mutationForceThermalEquilibrium_)
                {
                    Tv = T;
                }

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
                    pT[facei] = T;
                    pTv[facei] = Tv;
                    pCp[facei] = Cp;
                    pCv[facei] = Cv;
                    pPsi[facei] = psi;
                    pMu[facei] = mu;
                    pKappa[facei] = kappa;

                    scalar eveNow, eveEqv, tau, qcv, qdot;
                    sourcesAtCurrentState
                    (
                        T, Tv, eveNow, eveEqv, tau, qcv, qdot
                    );

                    pEveEq[facei] = eveEqv;
                    pTau[facei] = tau;
                }
                else
                {
                    pEveEq[facei] = pEve[facei];
                    pTau[facei] = tauInf_;
                    hadFailure = true;
                }
            }
        }

        if (hadFailure)
        {
            WarningInFunction
                << "Mutation++ update failed for at least one location. "
                << "Previous values were kept for failed entries."
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
        mutWdot_.resize(nMutationSpecies_);

        // Formation enthalpies per unit mass [J/kg]: hf/RT is returned
        // normalised by the passed temperature, so hf = (hf/RT)*RU*T is
        // temperature-independent. Evaluated once here.
        {
            const double Tref = 298.15;
            mutationMixturePtr_->speciesHOverRT
            (
                Tref, Tref, Tref, Tref, Tref,
                nullptr, nullptr, nullptr, nullptr, nullptr,
                mutW1_.data()
            );

            mutationHfMass_.resize(nMutationSpecies_);
            for (int i = 0; i < nMutationSpecies_; ++i)
            {
                mutationHfMass_[i] =
                    mutW1_[i]*Mutation::RU*Tref/mutationMw_[i];
            }
        }

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
void HighEnthalpyMulticomponentThermo<Thermo>::initialiseEnergies()
{
    if (!mutationAvailable_ || !mutationMixturePtr_)
    {
        return;
    }

    const PtrList<volScalarField>& Y = this->Y();

    const scalarField& TCells = this->T_.primitiveField();
    const scalarField& TveCells = Tve_.primitiveField();
    scalarField& eveCells = eve_.primitiveFieldRef();
    scalarField& eCells = this->he().primitiveFieldRef();

    forAll(eveCells, celli)
    {
        const auto yFoam = [&](const label speciei)
        {
            return Y[speciei][celli];
        };

        mutationMassFractions(yFoam);

        const scalar T = max(TCells[celli], mutationMinT_);
        const scalar Tv = max(TveCells[celli], mutationMinT_);

        eveCells[celli] = mixtureEve(T, Tv);

        // Solved energy on the Mutation++ datum: e = e_tr(T) + e_ve(Tve)
        eCells[celli] = mixtureEtr(T) + eveCells[celli];
    }

    volScalarField::Boundary& eveBf = eve_.boundaryFieldRef();
    volScalarField::Boundary& eBf = this->he().boundaryFieldRef();
    const volScalarField::Boundary& TBf = this->T_.boundaryField();
    const volScalarField::Boundary& TveBf = Tve_.boundaryField();

    forAll(eveBf, patchi)
    {
        forAll(eveBf[patchi], facei)
        {
            const auto yFoam = [&](const label speciei)
            {
                return Y[speciei].boundaryField()[patchi][facei];
            };

            mutationMassFractions(yFoam);

            const scalar T = max(TBf[patchi][facei], mutationMinT_);
            const scalar Tv = max(TveBf[patchi][facei], mutationMinT_);

            eveBf[patchi][facei] = mixtureEve(T, Tv);
            eBf[patchi][facei] = mixtureEtr(T) + eveBf[patchi][facei];
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
    // Floor on every Mutation++ specie: exactly-zero mole fractions make
    // Mutation++'s OmegaVT produce infinities (tau_park = 1/(n*X*...) for
    // each vibrator of the mixture). Those are benign in IEEE arithmetic
    // (zero net contribution) but trip OpenFOAM's floating point trapping.
    // A 1e-30 floor is physically indistinguishable from zero.
    constexpr double yFloor = 1e-30;

    std::fill(mutY_.begin(), mutY_.end(), yFloor);

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

    ySum += nMutationSpecies_*yFloor;

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
inline Foam::scalar HighEnthalpyMulticomponentThermo<Thermo>::mixtureEtr
(
    const scalar T
) const
{
    // ht and hr only depend on the trans-rotational temperature; the -1
    // converts enthalpy to internal energy (h = e + RT per mole)
    mutationMixturePtr_->speciesHOverRT
    (
        T, T, T, T, T,
        nullptr,
        mutW1_.data(),
        mutW2_.data(),
        nullptr,
        nullptr,
        nullptr
    );

    double e = 0.0;
    for (int i = 0; i < nMutationSpecies_; ++i)
    {
        e +=
            mutY_[i]*(mutW1_[i] + mutW2_[i] - 1.0)
           *Mutation::RU*T/mutationMw_[i];
    }

    return e;
}


template<class Thermo>
inline Foam::scalar HighEnthalpyMulticomponentThermo<Thermo>::mixtureCvTr
(
    const scalar T
) const
{
    mutationMixturePtr_->speciesCpOverR
    (
        T, T, T, T, T,
        nullptr,
        mutW1_.data(),
        mutW2_.data(),
        nullptr,
        nullptr
    );

    double cv = 0.0;
    for (int i = 0; i < nMutationSpecies_; ++i)
    {
        cv +=
            mutY_[i]*(mutW1_[i] + mutW2_[i] - 1.0)
           *Mutation::RU/mutationMw_[i];
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
inline Foam::scalar HighEnthalpyMulticomponentThermo<Thermo>::TFromEtr
(
    const scalar etrTarget,
    const scalar TGuess
) const
{
    scalar T = min(max(TGuess, mutationMinT_), mutationMaxT_);

    // e_tr is linear in T for RRHO: this converges in very few iterations
    for (int iter = 0; iter < 50; ++iter)
    {
        const scalar f = mixtureEtr(T) - etrTarget;
        const scalar cvTr = max(mixtureCvTr(T), scalar(1e-8));

        scalar dT = -f/cvTr;

        dT = min(max(dT, -0.5*T), 0.5*T);

        T = min(max(T + dT, mutationMinT_), mutationMaxT_);

        if (mag(dT) < 1e-4)
        {
            break;
        }
    }

    return T;
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
inline void HighEnthalpyMulticomponentThermo<Thermo>::sourcesAtCurrentState
(
    const scalar T,
    const scalar Tv,
    scalar& eveNow,
    scalar& eveEqv,
    scalar& tau,
    scalar& qcv,
    scalar& qdot
) const
{
    eveNow = mixtureEve(T, Tv);
    eveEqv = mixtureEve(T, T);

    tau = tauInf_;
    qcv = 0;
    qdot = 0;

    try
    {
        // Net mass production rates [kg/m^3/s] at the current state
        mutationMixturePtr_->netProductionRates(mutWdot_.data());

        // Candler chemistry-vibration (+ chemistry-electronic) coupling and
        // reaction heat release, mirroring OmegaCV/OmegaCElec exactly:
        // e_ve,i(Tv) = (hv_i + hel_i)/RT * RU*T / Mw_i  [J/kg]
        mutationMixturePtr_->speciesHOverRT
        (
            nullptr,
            nullptr,
            nullptr,
            mutW1_.data(),
            mutW2_.data(),
            nullptr
        );

        const double Tstate = mutationMixturePtr_->T();

        for (int i = 0; i < nMutationSpecies_; ++i)
        {
            const double eVei =
                (mutW1_[i] + mutW2_[i])*Mutation::RU*Tstate/mutationMw_[i];

            qcv += eVei*mutWdot_[i];
            qdot -= mutationHfMass_[i]*mutWdot_[i];
        }

        // Total ve-pool source = pure V-T + CV + CElec; recast only the pure
        // V-T part into the Landau-Teller form Q = rho*(eveEq - eve)/tau
        mutationMixturePtr_->energyTransferSource(mutQ_.data());

        const double Qvt = mutQ_[0] - qcv;
        const double dE = eveEqv - eveNow;
        const double rho = mutationMixturePtr_->density();

        if (mag(Qvt) > SMALL && dE*Qvt > 0)
        {
            tau = min(max(rho*dE/Qvt, scalar(1e-12)), scalar(tauInf_));
        }
    }
    catch (...)
    {
        tau = tauInf_;
        qcv = 0;
        qdot = 0;
        std::fill(mutWdot_.begin(), mutWdot_.end(), 0.0);
    }
}

} // End namespace Foam

#endif

// ************************************************************************* //

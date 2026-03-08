#ifndef HighEnthalpyMulticomponentThermo_H
#define HighEnthalpyMulticomponentThermo_H

#include "fluidMulticomponentThermo.H"
#include "multicomponentThermo.H"
#include "psiThermo.H"
#include "basicThermo.H"
#include "fvMesh.H"

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

    //- Vibrational-electronic temperature field [K]
    volScalarField Tve_;

    //- Optional Mutation++ bridge. If not available, fallback to base thermo.
    std::unique_ptr<Mutation::Mixture> mutationMixturePtr_;

    //- Mapping OpenFOAM specie index -> Mutation++ specie index
    List<label> foamToMutationSpecie_;

    word mutationMixtureName_;
    word mutationStateModel_;
    word mutationThermoDatabase_;

    scalar mutationMinT_;
    scalar mutationMinP_;
    Switch mutationForceThermalEquilibrium_;

    bool mutationAvailable_;

    //- Load Mutation++ options and build species mapping
    void initialiseMutation();

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
        mutationMixturePtr_(nullptr),
        foamToMutationSpecie_(),
        mutationMixtureName_("air_11"),
        mutationStateModel_("ChemNonEqTTv"),
        mutationThermoDatabase_("RRHO"),
        mutationMinT_(scalar(200)),
        mutationMinP_(scalar(1)),
        mutationForceThermalEquilibrium_(false),
        mutationAvailable_(false)
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

        initialiseMutation();
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

    // Member Functions
    virtual void correct() override
    {
        // First execute the default OpenFOAM update to preserve robustness.
        FluidMulticomponentThermo<Thermo>::correct();

        if (!mutationAvailable_ || !mutationMixturePtr_)
        {
            return;
        }

        const PtrList<volScalarField>& Y = this->Y();

        const scalarField& pCells = this->p_.primitiveField();
        scalarField& TCells = this->T_.primitiveFieldRef();
        scalarField& TveCells = Tve_.primitiveFieldRef();
        scalarField& CpCells = this->Cp_.primitiveFieldRef();
        scalarField& CvCells = this->Cv_.primitiveFieldRef();
        scalarField& psiCells = this->psi_.primitiveFieldRef();
        scalarField& muCells = this->mu_.primitiveFieldRef();
        scalarField& kappaCells = this->kappa_.primitiveFieldRef();

        std::vector<double> yMutation(mutationMixturePtr_->nSpecies(), 0.0);

        bool hadFailure = false;

        forAll(TCells, celli)
        {
            for (std::size_t i = 0; i < yMutation.size(); ++i)
            {
                yMutation[i] = 0.0;
            }

            forAll(Y, speciei)
            {
                const label mutationSpeciei = foamToMutationSpecie_[speciei];
                if (mutationSpeciei >= 0)
                {
                    yMutation[mutationSpeciei] += max(Y[speciei][celli], scalar(0));
                }
            }

            scalar T = max(TCells[celli], mutationMinT_);
            scalar Tv = mutationForceThermalEquilibrium_
                ? T
                : max(TveCells[celli], mutationMinT_);
            scalar Cp = CpCells[celli];
            scalar Cv = CvCells[celli];
            scalar psi = psiCells[celli];
            scalar mu = muCells[celli];
            scalar kappa = kappaCells[celli];

            if
            (
                updateFromMutation
                (
                    yMutation,
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
            }
            else
            {
                hadFailure = true;
            }
        }

        const volScalarField::Boundary& pBf = this->p_.boundaryField();
        volScalarField::Boundary& TBf = this->T_.boundaryFieldRef();
        volScalarField::Boundary& TveBf = Tve_.boundaryFieldRef();
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
            fvPatchScalarField& pCp = CpBf[patchi];
            fvPatchScalarField& pCv = CvBf[patchi];
            fvPatchScalarField& pPsi = psiBf[patchi];
            fvPatchScalarField& pMu = muBf[patchi];
            fvPatchScalarField& pKappa = kappaBf[patchi];

            forAll(pT, facei)
            {
                for (std::size_t i = 0; i < yMutation.size(); ++i)
                {
                    yMutation[i] = 0.0;
                }

                forAll(Y, speciei)
                {
                    const label mutationSpeciei = foamToMutationSpecie_[speciei];
                    if (mutationSpeciei >= 0)
                    {
                        yMutation[mutationSpeciei] +=
                            max(Y[speciei].boundaryField()[patchi][facei], scalar(0));
                    }
                }

                scalar T = max(pT[facei], mutationMinT_);
                scalar Tv = mutationForceThermalEquilibrium_
                    ? T
                    : max(pTv[facei], mutationMinT_);
                scalar Cp = pCp[facei];
                scalar Cv = pCv[facei];
                scalar psi = pPsi[facei];
                scalar mu = pMu[facei];
                scalar kappa = pKappa[facei];

                if
                (
                    updateFromMutation
                    (
                        yMutation,
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
                }
                else
                {
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
        mutationMinT_ =
            mutationDict.lookupOrDefault<scalar>("minTemperature", mutationMinT_);
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

        mutationMixturePtr_.reset(new Mutation::Mixture(opts));

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
bool HighEnthalpyMulticomponentThermo<Thermo>::updateFromMutation
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

    std::vector<double> y = yMutation;
    scalar ySum = 0;
    for (std::size_t i = 0; i < y.size(); ++i)
    {
        y[i] = max(y[i], 0.0);
        ySum += y[i];
    }

    if (ySum <= SMALL)
    {
        y.assign(y.size(), 0.0);
        y[0] = 1.0;
    }
    else
    {
        for (std::size_t i = 0; i < y.size(); ++i)
        {
            y[i] /= ySum;
        }
    }

    const double state[3] =
    {
        max(p, mutationMinP_),
        max(T, mutationMinT_),
        max(Tv, mutationMinT_)
    };

    try
    {
        mutationMixturePtr_->setState(y.data(), state, 2);

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

} // End namespace Foam

#endif

// ************************************************************************* //

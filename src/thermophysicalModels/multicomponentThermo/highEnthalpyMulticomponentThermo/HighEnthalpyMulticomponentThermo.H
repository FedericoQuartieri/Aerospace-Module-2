#ifndef HighEnthalpyMulticomponentThermo_H
#define HighEnthalpyMulticomponentThermo_H

#include "fluidMulticomponentThermo.H"
#include "multicomponentThermo.H"
#include "psiThermo.H"
#include "basicThermo.H"
#include "fvMesh.H"
#include "zeroGradientFvPatchFields.H"
#include "constants.H"

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

    //- Optional Mutation++ bridge.
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

    //- Load Mutation++ options and build species mapping
    void initMutation();

    inline scalar hToE(scalar h, scalar T, scalar M) const
    {
        // TODO: Use openfoam constant

        // speciesHOverRT gives h/RT, moreover h_v = e_v, h_el = e_el
        // while h_tr = e_rt - RT
        return h * (Mutation::RU * T / M); // Rm here
    }


public:

    // Constructors
    HighEnthalpyMulticomponentThermo
    (
        const fvMesh& mesh,
        const word& phaseName = word::null
    )
    :
        FluidMulticomponentThermo<Thermo>(mesh, phaseName),
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
        initMutation();
    }

    virtual ~HighEnthalpyMulticomponentThermo() = default;

    // Member Functions

    virtual tmp<volScalarField> computeSourceVT(const label speciei) const override
    {
        tmp<volScalarField> Q_VT
        (
            volScalarField::New
            (
                "Q_VT",
                this->mesh(),
                dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
            )
        );

        scalarField& QCells = Q_VT.ref().primitiveFieldRef();
        const scalarField& TCells = this->T_.primitiveField();
        const scalarField& TveCells = this->Tve_.primitiveField();
        const scalarField& PCells = this->p_.primitiveField();

        // Mixture density
        tmp<volScalarField> tRho = this->rho();
        const scalarField& RhoCells = tRho().primitiveField();
        // Mixture molecular weight field
        tmp<volScalarField> tWmix = this->W();
        const scalarField& WmixCells = tWmix().primitiveField();

        // Mass fraction field per specie
        const PtrList<volScalarField>& Y = this->Y();
        const label nSpecies = this->species().size();
        // Molecular weight of given speciei
        const scalar Wi = this->WiValue(speciei);

        // TODO: Perhaps query mutation
        constexpr scalar theta_v_N2 = 3371.0;

        std::vector<scalar> h_v(nSpecies);
        std::vector<scalar> h_el(nSpecies);

        // Loop over all cells of the mesh
        forAll(TCells, celli)
        {
            const scalar p = PCells[celli];
            const scalar Ttr = TCells[celli];

            // Mixture molecular weight
            const scalar Wmix = WmixCells[celli];

            // Compute tau_vt
            scalar inv_tau_vt = 0;
            // Loop over the species
            for (label s = 0; s < nSpecies; s++) {
                // TODO: Precompute cell independent values
                const scalar Ws = this->WiValue(s);
                const scalar Xs = Y[s][celli] * (Wmix / Ws);

                if (Xs > 1e-15)  {
                    const scalar Wis = (Wi * Ws) / (Wi + Ws);

                    // Millikan-white contribution
                    const scalar A = 1.16e-3 * sqrt(Wis)
                                             * pow(theta_v_N2, 4.0 / 3);
                    const scalar B = 1.5e-2 * pow(Wis, 0.25);

                    const scalar tauMW = exp(A * (
                        pow(Ttr, -1.0 / 3) - B) - 18.42) /
                        (p / 101325); // convert to atm

                    // TODO: Get raw pointer to YCells?
                    // Mole fraction of specie s in celli
                    // Park contribution
                    const scalar Ri = constant::thermodynamic::RR / Wi;
                    const scalar cbar = sqrt(8.0 * Ri * Ttr /
                                             constant::mathematical::pi);
                    const scalar sigma = 3e-21 * pow(50000.0 / Ttr, 2.0);
                    const scalar nTot = p / (constant::physicoChemical::k.value() * Ttr);
                    const double tauPark = 1.0 / (cbar * sigma * Xs * nTot);

                    inv_tau_vt += Xs / (tauMW + tauPark);
                }
            }

            // Partial density of speciei
            const scalar rho = RhoCells[celli] * Y[speciei][celli];
            const scalar Tve = TveCells[celli];

            const label mIdx = foamToMutationSpecie_[speciei];

            // Compute target energy e_i,ve(T_tr) (at equilibrium)
            mutationMixturePtr_->speciesHOverRT(
                Ttr, Ttr, Ttr, Ttr, Ttr,
                nullptr, nullptr, nullptr, h_v.data(), h_el.data(), nullptr);
            scalar eve_eq = hToE(h_v[mIdx] + h_el[mIdx], Ttr, Wi);

            // Compute current energy e_i,ve(T_i,ve)
            mutationMixturePtr_->speciesHOverRT(
                Ttr, Ttr, Ttr, Tve, Tve,
                nullptr, nullptr, nullptr, h_v.data(), h_el.data(), nullptr);
            double eve = hToE(h_v[mIdx] + h_el[mIdx], Ttr, Wi);

            // rho (e_ve(T_tr) - e_ve(T_ve)) / tau_vt
            QCells[celli] = rho * (eve_eq - eve) * inv_tau_vt;
        }

        return Q_VT;
    }

    virtual void correct() override
    {
        // First execute the default OpenFOAM update to preserve robustness.
        FluidMulticomponentThermo<Thermo>::correct();

        // WARNING: Currently deactivated.
        return;

        /*
        const PtrList<volScalarField>& Y = this->Y();

        const scalarField& pCells = this->p_.primitiveField();
        scalarField& TCells = this->T_.primitiveFieldRef();
        scalarField& TveCells = this->Tve_.primitiveFieldRef();
        const scalarField& eveCells = this->eve_.primitiveField();
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

        if (hadFailure)
        {
            WarningInFunction
                << "sparati"
                << endl;
        }
    */
    }

};

template<class Thermo>
void HighEnthalpyMulticomponentThermo<Thermo>::initMutation()
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

    Mutation::MixtureOptions opts(mutationMixtureName_.c_str());
    opts.setStateModel(mutationStateModel_.c_str());
    opts.setThermodynamicDatabase(mutationThermoDatabase_.c_str());

    if (!mutationMechanism_.empty())
    {
        opts.setMechanism(mutationMechanism_.c_str());
    }

    mutationMixturePtr_.reset(new Mutation::Mixture(opts));

    nMutationSpecies_ = mutationMixturePtr_->nSpecies();

    mutationMw_.resize(nMutationSpecies_);
    for (int i = 0; i < nMutationSpecies_; ++i)
    {
        mutationMw_[i] = mutationMixturePtr_->speciesMw(i);
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

} // End namespace Foam

#endif

// ************************************************************************* //

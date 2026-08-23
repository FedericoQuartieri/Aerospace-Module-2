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

    //- Mutation++ bridge
    std::unique_ptr<Mutation::Mixture> mutationMixturePtr_;

    //- Mapping OpenFOAM specie index -> Mutation++ specie index
    List<label> foamToMutationSpecie_;

    word mutationMixtureName_;
    word mutationStateModel_;
    word mutationThermoDatabase_;
    word mutationMechanism_;

    int nMutationSpecies_;
    int nMutationEnergyEqns_;

    void initMutation();

    inline scalar hToE(scalar h, scalar T, scalar M) const
    {
        // speciesHOverRT gives h/RT, moreover h_v = e_v, h_el = e_el
        // while h_tr = e_rt - RT
        return h * (constant::physicoChemical::RR.value() * T / M); // Rm here
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
        mutationMechanism_("none"),
        nMutationSpecies_(0),
        nMutationEnergyEqns_(0)
    {
        initMutation();

        // Initialize eve_ and he_ from Tve_ and T_
        const scalarField& TCells = this->T_.primitiveField();
        const scalarField& TveCells = this->Tve_.primitiveField();
        scalarField& eveCells = this->eve_.primitiveFieldRef();
        scalarField& heCells = this->he().primitiveFieldRef();

        const label defaultIdx = this->defaultSpecie();
        const scalar Wi = this->WiValue(defaultIdx);
        const scalar Ru = constant::physicoChemical::RR.value();
        const scalar Rs = Ru / Wi;
        const scalar Cv_tr = 2.5 * Rs;

        const label mIdx = foamToMutationSpecie_[defaultIdx];
        std::vector<double> h_v(nMutationSpecies_, 0.0);
        std::vector<double> h_el(nMutationSpecies_, 0.0);

        forAll(TCells, celli)
        {
            const scalar Ttr = TCells[celli];
            const scalar Tve = TveCells[celli];

            mutationMixturePtr_->speciesHOverRT(
                Ttr, Ttr, Ttr, Tve, Tve,
                nullptr, nullptr, nullptr, h_v.data(), h_el.data(), nullptr);

            eveCells[celli] = hToE(h_v[mIdx] + h_el[mIdx], Ttr, Wi);
            heCells[celli]  = Cv_tr * Ttr + eveCells[celli];
        }
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
        const scalarField& rhoCells = tRho().primitiveField();
        // Mixture molecular weight field
        tmp<volScalarField> tWmix = this->W();
        const scalarField& WmixCells = tWmix().primitiveField();

        // Mass fraction field per specie
        const PtrList<volScalarField>& Y = this->Y();
        const label nSpecies = this->species().size();
        // Molecular weight of given speciei
        const scalar Wi = this->WiValue(speciei);

        constexpr scalar theta_v_N2 = 3371.0;
        const scalar Ru = constant::physicoChemical::RR.value();
        const scalar kb = constant::physicoChemical::k.value();
        const scalar PI = constant::mathematical::pi;

        const label mIdx = foamToMutationSpecie_[speciei];
        std::vector<double> h_v(nMutationSpecies_, 0.0);
        std::vector<double> h_el(nMutationSpecies_, 0.0);

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
                // Mole fraction of specie s in celli
                const scalar Xs = Y[s][celli] * (Wmix / Ws);

                if (Xs > 1e-15)  {
                    const scalar Wis = (Wi * Ws) / (Wi + Ws);

                    // Millikan-white contribution
                    const scalar A = 1.16e-3 * sqrt(Wis)
                                             * pow(theta_v_N2, 4.0 / 3.0);
                    const scalar B = 1.5e-2 * pow(Wis, 0.25);

                    const scalar tauMW = exp(A * (
                        pow(Ttr, -1.0 / 3.0) - B) - 18.42) /
                        (p / 101325.0); // convert to atm

                    // TODO: Get raw pointer to YCells?
                    // Park contribution
                    const scalar Ri = Ru / Wi;
                    const scalar cbar = sqrt(8.0 * Ri * Ttr / PI);
                    const scalar sigma = 3e-21 * pow(50000.0 / Ttr, 2.0);
                    const scalar nTot = p / (kb * Ttr);
                    const double tauPark = 1.0 / (cbar * sigma * Xs * nTot);

                    inv_tau_vt += Xs / (tauMW + tauPark);
                }
            }

            // Partial density of speciei
            const scalar rho = rhoCells[celli] * Y[speciei][celli];
            const scalar Tve = TveCells[celli];

            // Compute target equilibrium energy and current energy from Mutation++
            mutationMixturePtr_->speciesHOverRT(
                Ttr, Ttr, Ttr, Ttr, Ttr,
                nullptr, nullptr, nullptr, h_v.data(), h_el.data(), nullptr);
            const scalar eve_eq = hToE(h_v[mIdx] + h_el[mIdx], Ttr, Wi);

            mutationMixturePtr_->speciesHOverRT(
                Ttr, Ttr, Ttr, Tve, Tve,
                nullptr, nullptr, nullptr, h_v.data(), h_el.data(), nullptr);
            const scalar eve = hToE(h_v[mIdx] + h_el[mIdx], Ttr, Wi);

            // rho (e_ve(T_tr) - e_ve(T_ve)) / tau_vt
            QCells[celli] = rho * (eve_eq - eve) * inv_tau_vt;
        }

        return Q_VT;
    }

    virtual void correct() override
    {
        // TODO: Should this be called?
        //FluidMulticomponentThermo<Thermo>::correct();

        const scalarField& eTotCells = this->he().primitiveField();
        const scalarField& eveCells  = this->eve_.primitiveField();

        scalarField& TCells   = this->T_.primitiveFieldRef();
        scalarField& TveCells = this->Tve_.primitiveFieldRef();
        scalarField& psiCells = this->psi_.primitiveFieldRef();

        const label defaultIdx = this->defaultSpecie();
        const scalar Wi = this->WiValue(defaultIdx);
        const scalar Ru = constant::physicoChemical::RR.value();
        const scalar Rs = Ru / Wi;
        const scalar Cv_tr = 2.5 * Rs;
        constexpr scalar theta_v_N2 = 3371.0;

        const label mIdx = foamToMutationSpecie_[defaultIdx];

        std::vector<double> h_v(nMutationSpecies_, 0.0);
        std::vector<double> h_el(nMutationSpecies_, 0.0);
        std::vector<double> cv_dim(nMutationSpecies_, 0.0);

        forAll(TCells, celli)
        {
            const scalar eTot = eTotCells[celli];
            const scalar eve  = eveCells[celli];

            const scalar e_tr = eTot - eve;
            const scalar Ttr  = e_tr / Cv_tr;

            // Harmonic oscillator initial guess for Tve
            scalar Tv = theta_v_N2 / Foam::log(1.0 + (Rs * theta_v_N2) / eve);

            // Newton refinement using Mutation++ (accounting for electronic excitation)
            for (int iter = 0; iter < 2; ++iter)
            {
                mutationMixturePtr_->speciesHOverRT(
                    Ttr, Ttr, Ttr, Tv, Tv,
                    nullptr, nullptr, nullptr, h_v.data(), h_el.data(), nullptr);
                const scalar eve_calc = (h_v[mIdx] + h_el[mIdx]) * (Rs * Ttr);

                mutationMixturePtr_->speciesCvOverR(
                    Ttr, Tv, Ttr, Tv, Tv,
                    nullptr, nullptr, nullptr, cv_dim.data(), nullptr);
                const scalar cv_ve = cv_dim[mIdx] * Rs;

                if (cv_ve > 1e-6)
                {
                    const scalar deltaT = (eve - eve_calc) / cv_ve;
                    Tv += deltaT;
                    if (Foam::mag(deltaT) < 1e-2) break;
                }
            }

            TCells[celli]   = Ttr;
            TveCells[celli] = Tv;
            psiCells[celli] = 1.0 / (Rs * Ttr);
        }
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
    }

    Mutation::MixtureOptions opts(mutationMixtureName_);
    opts.setStateModel(mutationStateModel_);
    opts.setThermodynamicDatabase(mutationThermoDatabase_);
    opts.setMechanism(mutationMechanism_);

    mutationMixturePtr_ = std::make_unique<Mutation::Mixture>(opts);
    nMutationSpecies_ = mutationMixturePtr_->nSpecies();
    nMutationEnergyEqns_ = mutationMixturePtr_->nEnergyEqns();

    const wordList& foamSpecies = this->species();
    foamToMutationSpecie_.setSize(foamSpecies.size(), -1);

    forAll(foamSpecies, i)
    {
        foamToMutationSpecie_[i] = mutationMixturePtr_->speciesIndex(foamSpecies[i].c_str());
    }
}

} // End namespace Foam

#endif

// ************************************************************************* //

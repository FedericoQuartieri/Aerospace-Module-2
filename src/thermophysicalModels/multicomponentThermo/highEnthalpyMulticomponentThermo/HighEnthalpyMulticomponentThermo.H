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
#include "mutationSources.H"

#include <algorithm>
#include <memory>
#include <vector>


namespace Foam
{

/*---------------------------------------------------------------------------*\
            Class HighEnthalpyMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

// bridge between OpenFOAM and Mutation++: in each cell the gas state is given
// by the conserved variables of the paper (eq. 23), partial densities rho_s,
// total energy and vibro-electronic energy; Mutation++ obtains from these
// the two temperatures, and the sources (V-T and chemistry) are those of
// mutationSources.H, the same as the 0D programs
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

    // exponent a of Park's temperature T^a Tv^(1-a) (eq. 29)
    scalar parkExponent_;

    // models of the paper: V-T tau of eq. 9-17 (false: the Mutation++ one)
    // and preferential Q_C-V with factor alpha (eq. 32; false: eq. 31)
    bool paperTau_;
    bool preferentialCV_;
    scalar alphaCV_;

    int nMutationSpecies_;
    int nMutationEnergyEqns_;

    // vibrating molecules, with their V-T relaxation time
    std::vector<Vibrator> vibrators_;

    // enthalpy of formation of each Mutation++ species (J/kg): the energy e
    // of the solver is sensible, the Mutation++ one includes the formation
    std::vector<double> hf_;

    // buffers filled cell by cell, in the order of the Mutation++ species
    mutable std::vector<double> yS_;    // mass fractions
    mutable std::vector<double> rhoS_;  // partial densities
    mutable std::vector<double> eS_;    // energies per species (total, then ve)
    mutable std::vector<double> wS_;    // species production

    // sources of each cell computed by correctSources()
    List<scalarField> wdotCells_;       // production of each Mutation++ species
    scalarField QchemCells_;            // heat of reaction for e
    scalarField QveCells_;              // source of eve

    void initMutation();

    void updatePsi();

    // mass fractions of the cell (or of the face) in the Mutation++ order
    template<class Field>
    void fillMassFractions(const Field& Yvalues, const label i) const
    {
        std::fill(yS_.begin(), yS_.end(), 0.0);
        forAll(foamToMutationSpecie_, s)
        {
            yS_[foamToMutationSpecie_[s]] = Yvalues[s][i];
        }
    }

    // mean energy of formation of the mixture (J/kg) from yS_
    double mixtureHf() const
    {
        double hf = 0.0;
        for (int s = 0; s < nMutationSpecies_; s++)
        {
            hf += yS_[s]*hf_[s];
        }
        return hf;
    }

    // brings Mutation++ to the state of the cell: partial densities and the two
    // temperatures (fills yS_ and rhoS_)
    void setCellState(const label celli, const scalar rho) const
    {
        fillMassFractions(this->Y(), celli);
        for (int s = 0; s < nMutationSpecies_; s++)
        {
            rhoS_[s] = rho*yS_[s];
        }
        const double TTv[2] = {this->T_[celli], this->Tve_[celli]};
        mutationMixturePtr_->setState(rhoS_.data(), TTv, 1);
    }

public:

    // ---- start-up: creates Mutation++ and turns the initial temperatures into energies
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
        mutationMixtureName_("air_5"),
        mutationStateModel_("ChemNonEqTTv"),
        mutationThermoDatabase_("RRHO"),
        mutationMechanism_("none"),
        parkExponent_(0.7),
        paperTau_(true),
        preferentialCV_(true),
        alphaCV_(0.3),
        nMutationSpecies_(0),
        nMutationEnergyEqns_(0)
    {
        // creates the Mutation++ mixture (see initMutation at the end)
        initMutation();

        Mutation::Mixture& mix = *mutationMixturePtr_;
        const int ns = nMutationSpecies_;

        // ---- initial energies of each cell from p, T and Tve read from the case
        // the base thermo (janaf) has already recomputed T from the energy, limiting it
        // to its validity range (20000 K): T is read again from the case file
        {
            volScalarField Tcase
            (
                IOobject
                (
                    this->T_.name(),
                    mesh.time().name(),
                    mesh,
                    IOobject::MUST_READ,
                    IOobject::NO_WRITE,
                    false
                ),
                mesh
            );
            this->T_ = Tcase;
        }

        const scalarField& pCells = this->p_.primitiveField();
        const scalarField& TCells = this->T_.primitiveField();
        const scalarField& TveCells = this->Tve_.primitiveField();
        scalarField& eveCells = this->eve_.primitiveFieldRef();
        scalarField& eCells = this->he().primitiveFieldRef();

        forAll(TCells, celli)
        {
            fillMassFractions(this->Y(), celli);
            const double pTTv[3] = {pCells[celli], TCells[celli], TveCells[celli]};
            mix.setState(yS_.data(), pTTv, 2);

            // eS_[s]: total energy of species s, eS_[ns + s]: the vibro-electronic one
            mix.getEnergiesMass(eS_.data());
            double e = 0.0;
            double eve = 0.0;
            for (int s = 0; s < ns; s++)
            {
                e += yS_[s]*(eS_[s] - hf_[s]);
                eve += yS_[s]*eS_[ns + s];
            }
            // e = e_tr + e_ve (eq. 1, 3-7, without formation), eve = e_ve (eq. 2)
            eCells[celli] = e;
            eveCells[celli] = eve;
        }
        this->he().correctBoundaryConditions();
        this->eve_.correctBoundaryConditions();
        // ---- end of initial energies

        // psi consistent with the T read again (the base thermo had computed it at 20000 K)
        updatePsi();
    }
    // ---- end of start-up

    virtual ~HighEnthalpyMulticomponentThermo() = default;

    // Member Functions

    // ---- sources of each cell, computed ONCE per corrector (called at the
    // start of thermophysicalPredictor) at the same state: Y, T, Tve and rho at
    // the start of the corrector. The equations of the species, of the total
    // energy and of eve then use the same production rates, so the mass produced,
    // the heat of reaction and Q_C-V are consistent with each other (previously each
    // species recomputed the chemistry after the previous species had already advanced)
    virtual void correctSources() override
    {
        const label nCells = this->T_.size();

        wdotCells_.setSize(nMutationSpecies_);
        forAll(wdotCells_, s)
        {
            wdotCells_[s] = scalarField(nCells, 0.0);
        }
        QchemCells_ = scalarField(nCells, 0.0);
        QveCells_ = scalarField(nCells, 0.0);

        // keeps the tmp of rho (= psi*p) alive as long as it is needed
        const tmp<volScalarField> tRho = this->rho();
        const scalarField& rhoCells = tRho().primitiveField();

        Mutation::Mixture& mix = *mutationMixturePtr_;
        const bool chemistry = (mutationMechanism_ != "none");

        forAll(QveCells_, celli)
        {
            setCellState(celli, rhoCells[celli]);

            // V-T exchange (eq. 8)
            double Qve = sourceVT(mix, rhoS_, vibrators_, paperTau_);

            if (chemistry)
            {
                // species production at Park's temperature (eq. 27-29)
                productionRates(mix, rhoS_, parkExponent_, wS_);

                // vibro-electronic energy removed by the reactions (eq. 30-32)
                Qve += sourceCV(mix, wS_, vibrators_, preferentialCV_, alphaCV_);

                // heat of reaction for the sensible energy e: the reactions
                // convert sensible energy into energy of formation,
                // -sum_s hf_s omega_s (eq. 1, last term)
                double Qchem = 0.0;
                for (int s = 0; s < nMutationSpecies_; s++)
                {
                    wdotCells_[s][celli] = wS_[s];
                    Qchem -= hf_[s]*wS_[s];
                }
                QchemCells_[celli] = Qchem;
            }
            QveCells_[celli] = Qve;
        }
    }
    // ---- end of sources

    // the three functions below return the values of correctSources()
    void checkSources() const
    {
        if (QveCells_.size() != this->T_.size())
        {
            FatalErrorInFunction
                << "call correctSources() before the source terms"
                << exit(FatalError);
        }
    }

    // source of the vibro-electronic energy (eq. 26): V-T plus Q_C-V [W/m3]
    virtual tmp<volScalarField> computeSourceVe() const override
    {
        checkSources();
        tmp<volScalarField> tQ
        (
            volScalarField::New
            (
                "Q_ve",
                this->mesh(),
                dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
            )
        );
        tQ.ref().primitiveFieldRef() = QveCells_;
        return tQ;
    }

    // chemical production of the OpenFOAM species speciei (eq. 27) [kg/m3/s]
    virtual tmp<volScalarField> computeSourceY(const label speciei) const override
    {
        checkSources();
        tmp<volScalarField> tw
        (
            volScalarField::New
            (
                "wdot_" + this->Y()[speciei].name(),
                this->mesh(),
                dimensionedScalar(dimDensity/dimTime, Zero)
            )
        );
        tw.ref().primitiveFieldRef() = wdotCells_[foamToMutationSpecie_[speciei]];
        return tw;
    }

    // heat of reaction for the sensible energy e, -sum hf_s omega_s [W/m3]
    virtual tmp<volScalarField> computeSourceE() const override
    {
        checkSources();
        tmp<volScalarField> tQ
        (
            volScalarField::New
            (
                "Q_chem",
                this->mesh(),
                dimensionedScalar(dimEnergy/dimVolume/dimTime, Zero)
            )
        );
        tQ.ref().primitiveFieldRef() = QchemCells_;
        return tQ;
    }

    // ---- decode: obtains T_tr, T_ve and psi from the solved energies
    virtual void correct() override
    {
        const scalarField& eCells = this->he().primitiveField();
        const scalarField& eveCells = this->eve_.primitiveField();
        const tmp<volScalarField> tRho = this->rho();
        const scalarField& rhoCells = tRho().primitiveField();

        scalarField& TCells = this->T_.primitiveFieldRef();
        scalarField& TveCells = this->Tve_.primitiveFieldRef();

        Mutation::Mixture& mix = *mutationMixturePtr_;

        forAll(TCells, celli)
        {
            const double rho = rhoCells[celli];
            fillMassFractions(this->Y(), celli);
            for (int s = 0; s < nMutationSpecies_; s++)
            {
                rhoS_[s] = rho*yS_[s];
            }

            // energies per unit volume: Mutation++ wants the total one
            // including formation, and the vibro-electronic one
            const double energies[2] =
            {
                rho*(eCells[celli] + mixtureHf()),
                rho*eveCells[celli]
            };
            // Mutation++ inverts the energies (Newton) and finds T_tr and T_ve
            mix.setState(rhoS_.data(), energies, 0);

            TCells[celli] = mix.T();
            TveCells[celli] = mix.Tv();
        }

        this->T_.correctBoundaryConditions();
        this->Tve_.correctBoundaryConditions();

        updatePsi();
    }
    // ---- end of decode

};

// ---- initMutation: reads highEnthalpyMutation from physicalProperties and creates the mixture
template<class Thermo>
void HighEnthalpyMulticomponentThermo<Thermo>::initMutation()
{
    const dictionary& propertiesDict = this->properties();

    // reads mixture, stateModel, thermodynamicDatabase, mechanism, parkExponent,
    // relaxationTime, chemistryVibration and preferentialFactor
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
        parkExponent_ =
            mutationDict.lookupOrDefault<scalar>("parkExponent", parkExponent_);

        // V-T tau: paper (eq. 9-17, default) or mutation
        const word tauModel =
            mutationDict.lookupOrDefault<word>("relaxationTime", "paper");
        // chemistry-vibration coupling: preferential (eq. 32, default)
        // or nonPreferential (eq. 31)
        const word cvModel =
            mutationDict.lookupOrDefault<word>("chemistryVibration", "preferential");
        if
        (
            (tauModel != "paper" && tauModel != "mutation")
         || (cvModel != "preferential" && cvModel != "nonPreferential")
        )
        {
            FatalErrorInFunction
                << "relaxationTime must be paper or mutation, "
                << "chemistryVibration preferential or nonPreferential"
                << exit(FatalError);
        }
        paperTau_ = (tauModel == "paper");
        preferentialCV_ = (cvModel == "preferential");
        alphaCV_ =
            mutationDict.lookupOrDefault<scalar>("preferentialFactor", alphaCV_);
    }

    Info<< "highEnthalpyMutation: tau V-T " << (paperTau_ ? "paper" : "mutation")
        << ", Q_C-V " << (preferentialCV_ ? "preferential" : "nonPreferential")
        << " (alpha " << alphaCV_ << "), Park exponent " << parkExponent_
        << endl;

    Mutation::MixtureOptions opts(mutationMixtureName_);
    opts.setStateModel(mutationStateModel_);
    opts.setThermodynamicDatabase(mutationThermoDatabase_);
    opts.setMechanism(mutationMechanism_);

    mutationMixturePtr_ = std::make_unique<Mutation::Mixture>(opts);
    nMutationSpecies_ = mutationMixturePtr_->nSpecies();
    nMutationEnergyEqns_ = mutationMixturePtr_->nEnergyEqns();

    // links each OpenFOAM species to the same Mutation++ species
    const wordList& foamSpecies = this->species();
    foamToMutationSpecie_.setSize(foamSpecies.size(), -1);

    forAll(foamSpecies, i)
    {
        foamToMutationSpecie_[i] =
            mutationMixturePtr_->speciesIndex(foamSpecies[i].c_str());

        if (foamToMutationSpecie_[i] < 0)
        {
            FatalErrorInFunction
                << "species " << foamSpecies[i]
                << " does not exist in the Mutation++ mixture " << mutationMixtureName_
                << exit(FatalError);
        }
    }

    // vibrating molecules and their relaxation times
    vibrators_ = makeVibrators(*mutationMixturePtr_);

    // work buffers
    yS_.assign(nMutationSpecies_, 0.0);
    rhoS_.assign(nMutationSpecies_, 0.0);
    eS_.assign(2*nMutationSpecies_, 0.0);
    wS_.assign(nMutationSpecies_, 0.0);

    // enthalpies of formation (J/kg): hf/RT from Mutation++, at any T
    const double T0 = 1000.0;
    std::vector<double> hfOverRT(nMutationSpecies_, 0.0);
    mutationMixturePtr_->speciesHOverRT
    (
        T0, T0, T0, T0, T0,
        nullptr, nullptr, nullptr, nullptr, nullptr, hfOverRT.data()
    );
    hf_.assign(nMutationSpecies_, 0.0);
    for (int s = 0; s < nMutationSpecies_; s++)
    {
        hf_[s] = hfOverRT[s]*Mutation::RU*T0/mutationMixturePtr_->speciesMw(s);
    }
}
// ---- end of initMutation

// ---- updatePsi: psi = rho/p in cells and boundary faces, with rho from Mutation++
// for given p, T, Tve and composition (eq. 24: p is the sum of the partial pressures)
template<class Thermo>
void HighEnthalpyMulticomponentThermo<Thermo>::updatePsi()
{
    Mutation::Mixture& mix = *mutationMixturePtr_;

    const scalarField& pCells = this->p_.primitiveField();
    const scalarField& TCells = this->T_.primitiveField();
    const scalarField& TveCells = this->Tve_.primitiveField();
    scalarField& psiCells = this->psi_.primitiveFieldRef();

    forAll(psiCells, celli)
    {
        fillMassFractions(this->Y(), celli);
        const double pTTv[3] = {pCells[celli], TCells[celli], TveCells[celli]};
        mix.setState(yS_.data(), pTTv, 2);
        psiCells[celli] = mix.density()/pCells[celli];
    }

    volScalarField::Boundary& psiBf = this->psi_.boundaryFieldRef();
    forAll(psiBf, patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        const fvPatchScalarField& pTve = this->Tve_.boundaryField()[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];

        forAll(ppsi, facei)
        {
            std::fill(yS_.begin(), yS_.end(), 0.0);
            forAll(foamToMutationSpecie_, s)
            {
                yS_[foamToMutationSpecie_[s]] =
                    this->Y()[s].boundaryField()[patchi][facei];
            }
            const double pTTv[3] = {pp[facei], pT[facei], pTve[facei]};
            mix.setState(yS_.data(), pTTv, 2);
            ppsi[facei] = mix.density()/pp[facei];
        }
    }
}
// ---- end of updatePsi

} // End namespace Foam

#endif

// ************************************************************************* //

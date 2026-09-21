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

// ponte fra OpenFOAM e Mutation++: in ogni cella lo stato del gas e' dato
// dalle variabili conservate del paper (eq. 23), densita' parziali rho_s,
// energia totale ed energia vibro-elettronica; Mutation++ ricava da queste
// le due temperature, e le sorgenti (V-T e chimica) sono quelle di
// mutationSources.H, le stesse dei programmi 0D
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

    // esponente a della temperatura di Park T^a Tv^(1-a) (eq. 29)
    scalar parkExponent_;

    // modelli del paper: tau V-T delle eq. 9-17 (falso: quello di Mutation++)
    // e Q_C-V preferenziale con fattore alpha (eq. 32; falso: eq. 31)
    bool paperTau_;
    bool preferentialCV_;
    scalar alphaCV_;

    int nMutationSpecies_;
    int nMutationEnergyEqns_;

    // molecole che vibrano, con il loro tempo di rilassamento V-T
    std::vector<Vibrator> vibrators_;

    // entalpia di formazione di ogni specie di Mutation++ (J/kg): l'energia e
    // del solver e' sensibile, quella di Mutation++ comprende la formazione
    std::vector<double> hf_;

    // buffer riempiti cella per cella, nell'ordine delle specie di Mutation++
    mutable std::vector<double> yS_;    // frazioni in massa
    mutable std::vector<double> rhoS_;  // densita' parziali
    mutable std::vector<double> eS_;    // energie per specie (totale, poi ve)
    mutable std::vector<double> wS_;    // produzione delle specie

    // sorgenti di ogni cella calcolate da correctSources()
    List<scalarField> wdotCells_;       // produzione di ogni specie di Mutation++
    scalarField QchemCells_;            // calore di reazione per e
    scalarField QveCells_;              // sorgente di eve

    void initMutation();

    void updatePsi();

    // frazioni in massa della cella (o della faccia) nell'ordine di Mutation++
    template<class Field>
    void fillMassFractions(const Field& Yvalues, const label i) const
    {
        std::fill(yS_.begin(), yS_.end(), 0.0);
        forAll(foamToMutationSpecie_, s)
        {
            yS_[foamToMutationSpecie_[s]] = Yvalues[s][i];
        }
    }

    // energia di formazione media della miscela (J/kg) da yS_
    double mixtureHf() const
    {
        double hf = 0.0;
        for (int s = 0; s < nMutationSpecies_; s++)
        {
            hf += yS_[s]*hf_[s];
        }
        return hf;
    }

    // porta Mutation++ nello stato della cella: densita' parziali e le due
    // temperature (riempie yS_ e rhoS_)
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

    // ---- avvio: crea Mutation++ e trasforma le temperature iniziali in energie
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
        // crea la miscela di Mutation++ (vedi initMutation in fondo)
        initMutation();

        Mutation::Mixture& mix = *mutationMixturePtr_;
        const int ns = nMutationSpecies_;

        // ---- energie di partenza di ogni cella da p, T e Tve letti dal caso
        // il thermo base (janaf) ha gia' ricalcolato T dall'energia limitandola
        // al suo intervallo di validita' (20000 K): rileggo T dal file del caso
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

            // eS_[s]: energia totale della specie s, eS_[ns + s]: quella vibro-elettronica
            mix.getEnergiesMass(eS_.data());
            double e = 0.0;
            double eve = 0.0;
            for (int s = 0; s < ns; s++)
            {
                e += yS_[s]*(eS_[s] - hf_[s]);
                eve += yS_[s]*eS_[ns + s];
            }
            // e = e_tr + e_ve (eq. 1, 3-7, senza formazione), eve = e_ve (eq. 2)
            eCells[celli] = e;
            eveCells[celli] = eve;
        }
        this->he().correctBoundaryConditions();
        this->eve_.correctBoundaryConditions();
        // ---- fine energie di partenza

        // psi coerente con la T riletta (il thermo base l'aveva calcolata a 20000 K)
        updatePsi();
    }
    // ---- fine avvio

    virtual ~HighEnthalpyMulticomponentThermo() = default;

    // Member Functions

    // ---- sorgenti di ogni cella, calcolate UNA volta per correttore (chiamata
    // all'inizio di thermophysicalPredictor) allo stesso stato: Y, T, Tve e rho
    // di inizio correttore. Le equazioni delle specie, dell'energia totale e di
    // eve usano poi le stesse velocita' di produzione, quindi la massa prodotta,
    // il calore di reazione e Q_C-V sono coerenti fra loro (prima ogni specie
    // ricalcolava la chimica dopo che le specie precedenti erano gia' avanzate)
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

        // tengo vivo il tmp di rho (= psi*p) finche' serve
        const tmp<volScalarField> tRho = this->rho();
        const scalarField& rhoCells = tRho().primitiveField();

        Mutation::Mixture& mix = *mutationMixturePtr_;
        const bool chemistry = (mutationMechanism_ != "none");

        forAll(QveCells_, celli)
        {
            setCellState(celli, rhoCells[celli]);

            // scambio V-T (eq. 8)
            double Qve = sourceVT(mix, rhoS_, vibrators_, paperTau_);

            if (chemistry)
            {
                // produzione delle specie alla temperatura di Park (eq. 27-29)
                productionRates(mix, rhoS_, parkExponent_, wS_);

                // energia vibro-elettronica portata via dalle reazioni (eq. 30-32)
                Qve += sourceCV(mix, wS_, vibrators_, preferentialCV_, alphaCV_);

                // calore di reazione per l'energia sensibile e: le reazioni
                // convertono energia sensibile in energia di formazione,
                // -somma_s hf_s omega_s (eq. 1, ultimo termine)
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
    // ---- fine sorgenti

    // le tre funzioni sotto restituiscono i valori di correctSources()
    void checkSources() const
    {
        if (QveCells_.size() != this->T_.size())
        {
            FatalErrorInFunction
                << "chiamare correctSources() prima delle sorgenti"
                << exit(FatalError);
        }
    }

    // sorgente dell'energia vibro-elettronica (eq. 26): V-T piu' Q_C-V [W/m3]
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

    // produzione chimica della specie di OpenFOAM speciei (eq. 27) [kg/m3/s]
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

    // calore di reazione per l'energia sensibile e, -somma hf_s omega_s [W/m3]
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

    // ---- decode: dalle energie risolte ricava T_tr, T_ve e psi
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

            // energie per unita' di volume: Mutation++ vuole quella totale
            // con la formazione, e quella vibro-elettronica
            const double energies[2] =
            {
                rho*(eCells[celli] + mixtureHf()),
                rho*eveCells[celli]
            };
            // Mutation++ inverte le energie (Newton) e trova T_tr e T_ve
            mix.setState(rhoS_.data(), energies, 0);

            TCells[celli] = mix.T();
            TveCells[celli] = mix.Tv();
        }

        this->T_.correctBoundaryConditions();
        this->Tve_.correctBoundaryConditions();

        updatePsi();
    }
    // ---- fine decode

};

// ---- initMutation: legge highEnthalpyMutation da physicalProperties e crea la miscela
template<class Thermo>
void HighEnthalpyMulticomponentThermo<Thermo>::initMutation()
{
    const dictionary& propertiesDict = this->properties();

    // legge mixture, stateModel, thermodynamicDatabase, mechanism, parkExponent,
    // relaxationTime, chemistryVibration e preferentialFactor
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

        // tau V-T: paper (eq. 9-17, default) oppure mutation
        const word tauModel =
            mutationDict.lookupOrDefault<word>("relaxationTime", "paper");
        // accoppiamento chimica-vibrazione: preferential (eq. 32, default)
        // oppure nonPreferential (eq. 31)
        const word cvModel =
            mutationDict.lookupOrDefault<word>("chemistryVibration", "preferential");
        if
        (
            (tauModel != "paper" && tauModel != "mutation")
         || (cvModel != "preferential" && cvModel != "nonPreferential")
        )
        {
            FatalErrorInFunction
                << "relaxationTime deve essere paper o mutation, "
                << "chemistryVibration preferential o nonPreferential"
                << exit(FatalError);
        }
        paperTau_ = (tauModel == "paper");
        preferentialCV_ = (cvModel == "preferential");
        alphaCV_ =
            mutationDict.lookupOrDefault<scalar>("preferentialFactor", alphaCV_);
    }

    Info<< "highEnthalpyMutation: tau V-T " << (paperTau_ ? "paper" : "mutation")
        << ", Q_C-V " << (preferentialCV_ ? "preferential" : "nonPreferential")
        << " (alpha " << alphaCV_ << "), esponente di Park " << parkExponent_
        << endl;

    Mutation::MixtureOptions opts(mutationMixtureName_);
    opts.setStateModel(mutationStateModel_);
    opts.setThermodynamicDatabase(mutationThermoDatabase_);
    opts.setMechanism(mutationMechanism_);

    mutationMixturePtr_ = std::make_unique<Mutation::Mixture>(opts);
    nMutationSpecies_ = mutationMixturePtr_->nSpecies();
    nMutationEnergyEqns_ = mutationMixturePtr_->nEnergyEqns();

    // collega ogni specie di OpenFOAM alla stessa specie di Mutation++
    const wordList& foamSpecies = this->species();
    foamToMutationSpecie_.setSize(foamSpecies.size(), -1);

    forAll(foamSpecies, i)
    {
        foamToMutationSpecie_[i] =
            mutationMixturePtr_->speciesIndex(foamSpecies[i].c_str());

        if (foamToMutationSpecie_[i] < 0)
        {
            FatalErrorInFunction
                << "la specie " << foamSpecies[i]
                << " non esiste nella miscela Mutation++ " << mutationMixtureName_
                << exit(FatalError);
        }
    }

    // molecole che vibrano e loro tempi di rilassamento
    vibrators_ = makeVibrators(*mutationMixturePtr_);

    // buffer di lavoro
    yS_.assign(nMutationSpecies_, 0.0);
    rhoS_.assign(nMutationSpecies_, 0.0);
    eS_.assign(2*nMutationSpecies_, 0.0);
    wS_.assign(nMutationSpecies_, 0.0);

    // entalpie di formazione (J/kg): hf/RT da Mutation++, a una T qualsiasi
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
// ---- fine initMutation

// ---- updatePsi: psi = rho/p in celle e facce di bordo, con rho da Mutation++
// per p, T, Tve e composizione date (eq. 24: p somma delle pressioni parziali)
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
// ---- fine updatePsi

} // End namespace Foam

#endif

// ************************************************************************* //

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
#include <exception>
#include <memory>
#include <string>
#include <vector>

#ifdef _OPENMP
#include <omp.h>
#endif


namespace Foam
{

// quanti thread OpenMP e quale sta girando; senza -fopenmp sono 1 e 0, e il
// codice resta esattamente quello di prima
inline int highEnthalpyMaxThreads()
{
#ifdef _OPENMP
    return omp_get_max_threads();
#else
    return 1;
#endif
}

inline int highEnthalpyThreadIndex()
{
#ifdef _OPENMP
    return omp_get_thread_num();
#else
    return 0;
#endif
}


/*---------------------------------------------------------------------------*\
                       Struct MutationWorkspace Declaration
\*---------------------------------------------------------------------------*/

// tutto cio' che Mutation++ modifica mentre si passa da una cella all'altra.
//
// Mutation::Mixture non e' thread-safe: setState() muta l'oggetto e tutte le
// interrogazioni successive (T(), Tv(), density(), netProductionRates(), ...)
// leggono quello stato. Due thread sulla stessa mixture si sovrascriverebbero
// a vicenda. Percio' ogni thread ha la sua mixture e i suoi buffer, creati una
// volta sola all'avvio in initMutation().
struct MutationWorkspace
{
    std::unique_ptr<Mutation::Mixture> mix;

    // i modelli di Millikan-White sono costruiti a partire dalla mixture qui
    // sopra e vanno usati solo con quella
    std::vector<Vibrator> vibrators;

    // buffer riempiti cella per cella, nell'ordine delle specie di Mutation++
    std::vector<double> yS;    // frazioni in massa
    std::vector<double> rhoS;  // densita' parziali
    std::vector<double> eS;    // energie per specie (totale, poi ve)
    std::vector<double> wS;    // produzione delle specie
};


/*---------------------------------------------------------------------------*\
            Class HighEnthalpyMulticomponentThermo Declaration
\*---------------------------------------------------------------------------*/

// ponte fra OpenFOAM e Mutation++: in ogni cella lo stato del gas e' dato
// dalle variabili conservate del paper (eq. 23), densita' parziali rho_s,
// energia totale ed energia vibro-elettronica; Mutation++ ricava da queste
// le due temperature, e le sorgenti (V-T e chimica) sono quelle di
// mutationSources.H, le stesse dei programmi 0D
//
// i cicli sulle celle sono indipendenti fra loro (ogni cella e' un problema a
// se') e sono distribuiti sui thread con OpenMP: vedi parallelFor()
template<class Thermo>
class HighEnthalpyMulticomponentThermo
:
    public FluidMulticomponentThermo<Thermo>
{
private:

    //- Mutation++ bridge: una copia per thread OpenMP
    mutable std::vector<std::unique_ptr<MutationWorkspace>> workspaces_;

    //- Mapping OpenFOAM specie index -> Mutation++ specie index
    List<label> foamToMutationSpecie_;

    word mutationMixtureName_;
    word mutationStateModel_;
    word mutationThermoDatabase_;
    word mutationMechanism_;

    // esponente a della temperatura di Park T^a Tv^(1-a) (eq. 29)
    scalar parkExponent_;

    int nMutationSpecies_;
    int nMutationEnergyEqns_;

    // entalpia di formazione di ogni specie di Mutation++ (J/kg): l'energia e
    // del solver e' sensibile, quella di Mutation++ comprende la formazione
    std::vector<double> hf_;

    // ---- sorgenti calcolate da updateSources() in un solo giro sulle celle
    // e poi lette dalle computeSource*: prima ogni computeSource* rifaceva il
    // giro per conto suo, e productionRates() veniva ricalcolata sei volte
    // (una per specie risolta, una per il calore di reazione, una per eve)
    // per ottenere numeri che escono tutti dallo stesso calcolo

    // omega_s di ogni specie, nell'ordine di Mutation++ [kg/m3/s]
    List<scalarField> wCache_;

    // sorgente dell'energia vibro-elettronica, V-T piu' chimica [W/m3]
    scalarField QveCache_;

    // calore di reazione per l'energia sensibile [W/m3]
    scalarField QchemCache_;

    // vero fra updateSources() e la correct() che segue: le computeSource*
    // senza un updateSources() davanti leggerebbero sorgenti di un altro stato
    bool sourcesValid_;
    // ---- fine sorgenti in cache

    void initMutation();

    void updatePsi();

    // le computeSource* non hanno niente da leggere se il solver non ha
    // chiamato updateSources(): meglio fermarsi che restituire numeri vecchi
    void checkSources() const
    {
        if (!sourcesValid_)
        {
            FatalErrorInFunction
                << "updateSources() non e' stata chiamata per lo stato"
                   " corrente: il solver deve chiamarla all'inizio del"
                   " predictor, prima di usare le computeSource*"
                << exit(FatalError);
        }
    }

    // lo spazio di lavoro del thread corrente
    MutationWorkspace& workspace() const
    {
        return *workspaces_[highEnthalpyThreadIndex()];
    }

    // ---- parallelFor: esegue body(i, ws) per i = 0 ... n-1
    //
    // ogni thread riceve il proprio MutationWorkspace, quindi il corpo puo'
    // chiamare setState() e le sorgenti senza sincronizzazione: le celle non
    // si parlano fra loro.
    //
    // schedule(static) e non runtime: il default di libgomp per runtime e'
    // dynamic con fette da UNA cella, che qui e' il peggiore dei mondi (molta
    // sincronizzazione e ordine variabile). Con static le celle, che costano
    // tutte uguale, sono ripartite in fette fisse: nessun overhead e stesso
    // risultato bit a bit a ogni corsa, che e' come si scoprono le corse
    // critiche (vedi applications/test/openmp).
    //
    // Su una CPU con core disuguali si guadagna un 10-15% passando a
    // schedule(dynamic, 64): su Apple Silicon quattro core sono performance e
    // quattro efficiency, e con le fette fisse tutti aspettano i thread finiti
    // sui core lenti. Misurato su 20000 celle a 8 thread: static 2.1x,
    // guided 2.2x, dynamic,64 2.4x. Si paga in riproducibilita': assegnando le
    // fette a chi si libera, l'ordine cambia a ogni corsa, e siccome il Newton
    // di Mutation++ riparte dalle temperature dell'ultima cella trattata dal
    // thread, i risultati ballano nell'ultima cifra.
    //
    // il numero di thread e' limitato dalla taglia del ciclo: aprire una
    // regione parallela costa, e i thread che restano senza lavoro poi girano
    // a vuoto in attesa (GOMP fa spin). Senza questo limite i casi heat bath a
    // una cella, e soprattutto le patch di bordo da due facce che updatePsi()
    // percorre a ogni passo, lancerebbero otto thread per niente
    template<class Body>
    void parallelFor(const label n, const Body& body) const
    {
        bool failed = false;
        std::string failMsg;

        // almeno minCellsPerThread celle a testa, altrimenti meno thread
        const label minCellsPerThread = 32;
        const int nThreads = std::min<label>
        (
            highEnthalpyMaxThreads(),
            std::max<label>(1, n/minCellsPerThread)
        );

        #pragma omp parallel num_threads(nThreads) if (nThreads > 1)
        {
            MutationWorkspace& ws = workspace();

            #pragma omp for schedule(static)
            for (label i = 0; i < n; i++)
            {
                try
                {
                    body(i, ws);
                }
                catch (const std::exception& e)
                {
                    // un'eccezione che esce da una regione OpenMP fa terminare
                    // il programma senza messaggio: la si ferma qui e si
                    // segnala dopo la regione (tipico: il Newton di
                    // setState() che non converge in una cella)
                    #pragma omp critical(highEnthalpyMutationFailure)
                    {
                        if (!failed)
                        {
                            failed = true;
                            failMsg = e.what();
                        }
                    }
                }
                catch (...)
                {
                    #pragma omp critical(highEnthalpyMutationFailure)
                    {
                        if (!failed)
                        {
                            failed = true;
                            failMsg = "eccezione non identificata";
                        }
                    }
                }
            }
        }

        if (failed)
        {
            FatalErrorInFunction
                << "Mutation++ ha fallito in una cella: " << failMsg.c_str()
                << exit(FatalError);
        }
    }
    // ---- fine parallelFor

    // frazioni in massa della cella (o della faccia) nell'ordine di Mutation++
    template<class Field>
    void fillMassFractions
    (
        const Field& Yvalues,
        const label i,
        MutationWorkspace& ws
    ) const
    {
        std::fill(ws.yS.begin(), ws.yS.end(), 0.0);
        forAll(foamToMutationSpecie_, s)
        {
            ws.yS[foamToMutationSpecie_[s]] = Yvalues[s][i];
        }
    }

    // energia di formazione media della miscela (J/kg) da ws.yS
    double mixtureHf(const MutationWorkspace& ws) const
    {
        double hf = 0.0;
        for (int s = 0; s < nMutationSpecies_; s++)
        {
            hf += ws.yS[s]*hf_[s];
        }
        return hf;
    }

    // porta la mixture del thread nello stato della cella: densita' parziali
    // e le due temperature (riempie ws.yS e ws.rhoS)
    void setCellState
    (
        const label celli,
        const scalar rho,
        MutationWorkspace& ws
    ) const
    {
        fillMassFractions(this->Y(), celli, ws);
        for (int s = 0; s < nMutationSpecies_; s++)
        {
            ws.rhoS[s] = rho*ws.yS[s];
        }
        const double TTv[2] = {this->T_[celli], this->Tve_[celli]};
        ws.mix->setState(ws.rhoS.data(), TTv, 1);
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
        workspaces_(),
        foamToMutationSpecie_(),
        mutationMixtureName_("air_5"),
        mutationStateModel_("ChemNonEqTTv"),
        mutationThermoDatabase_("RRHO"),
        mutationMechanism_("none"),
        parkExponent_(0.7),
        nMutationSpecies_(0),
        nMutationEnergyEqns_(0),
        wCache_(),
        QveCache_(),
        QchemCache_(),
        sourcesValid_(false)
    {
        // crea le miscele di Mutation++, una per thread (vedi initMutation in fondo)
        initMutation();

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

        parallelFor
        (
            TCells.size(),
            [&](const label celli, MutationWorkspace& ws)
            {
                fillMassFractions(this->Y(), celli, ws);
                const double pTTv[3] =
                    {pCells[celli], TCells[celli], TveCells[celli]};
                ws.mix->setState(ws.yS.data(), pTTv, 2);

                // ws.eS[s]: energia totale della specie s,
                // ws.eS[ns + s]: quella vibro-elettronica
                ws.mix->getEnergiesMass(ws.eS.data());
                double e = 0.0;
                double eve = 0.0;
                for (int s = 0; s < ns; s++)
                {
                    e += ws.yS[s]*(ws.eS[s] - hf_[s]);
                    eve += ws.yS[s]*ws.eS[ns + s];
                }
                // e = e_tr + e_ve (eq. 1, 3-7, senza formazione), eve = e_ve (eq. 2)
                eCells[celli] = e;
                eveCells[celli] = eve;
            }
        );
        this->he().correctBoundaryConditions();
        this->eve_.correctBoundaryConditions();
        // ---- fine energie di partenza

        // psi coerente con la T riletta (il thermo base l'aveva calcolata a 20000 K)
        updatePsi();
    }
    // ---- fine avvio

    virtual ~HighEnthalpyMulticomponentThermo() = default;

    // Member Functions

    // ---- updateSources: un solo giro sulle celle per tutte le sorgenti
    //
    // Prima ogni computeSource* faceva il suo giro, e per ogni cella
    // ricalcolava productionRates(): quattro volte per le specie risolte, una
    // per il calore di reazione, una dentro la sorgente di eve. Sei calcoli
    // identici della stessa cosa, il 65% del tempo dell'intera corsa.
    //
    // Qui si calcola una volta sola, e ne escono insieme le omega_s di tutte
    // le specie, il calore di reazione e l'accoppiamento chimica-vibrazione,
    // perche' vengono tutti dal medesimo vettore omega_s.
    //
    // Cambia anche il risultato, e in meglio: prima la sorgente della specie i
    // vedeva le specie 0...i-1 gia' avanzate dalle loro YiEqn, e il calore di
    // reazione le vedeva tutte avanzate e rinormalizzate. Le eq. 27 e 30 del
    // paper valutano tutti gli omega_s allo STESSO stato, che e' quello che si
    // fa adesso.
    virtual void updateSources() override
    {
        const label nCells = this->mesh().nCells();
        const bool chemistry = (mutationMechanism_ != "none");
        const int ns = nMutationSpecies_;

        if (QveCache_.size() != nCells)
        {
            QveCache_.setSize(nCells);
            QchemCache_.setSize(nCells);

            // le omega_s servono solo con un meccanismo di reazione
            wCache_.setSize(chemistry ? ns : 0);
            forAll(wCache_, s)
            {
                wCache_[s].setSize(nCells);
            }
        }

        // tengo vivo il tmp di rho (= psi*p) finche' serve
        const tmp<volScalarField> tRho = this->rho();
        const scalarField& rhoCells = tRho().primitiveField();

        parallelFor
        (
            nCells,
            [&](const label celli, MutationWorkspace& ws)
            {
                setCellState(celli, rhoCells[celli], ws);

                // scambio V-T di Landau-Teller (eq. 8)
                scalar Qve = sourceVT(*ws.mix, ws.rhoS, ws.vibrators);
                scalar Qchem = 0;

                if (chemistry)
                {
                    // l'unico calcolo delle velocita' di reazione
                    productionRates(*ws.mix, ws.rhoS, parkExponent_, ws.wS);

                    // energia vibro-elettronica portata via dalle reazioni
                    // (eq. 30)
                    Qve += sourceCV(*ws.mix, ws.wS);

                    for (int s = 0; s < ns; s++)
                    {
                        // produzione della specie s (eq. 27)
                        wCache_[s][celli] = ws.wS[s];
                        // calore di reazione, -somma_s hf_s omega_s
                        Qchem -= hf_[s]*ws.wS[s];
                    }
                }

                QveCache_[celli] = Qve;
                QchemCache_[celli] = Qchem;
            }
        );

        sourcesValid_ = true;
    }
    // ---- fine updateSources

    // ---- sorgente dell'energia vibro-elettronica in ogni cella (eq. 26):
    // scambio V-T (eq. 8) piu', con la chimica, l'energia vibro-elettronica
    // portata via dalle reazioni (eq. 30)
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
        tQ.ref().primitiveFieldRef() = QveCache_;

        return tQ;
    }
    // ---- fine sorgente di eve

    // ---- produzione chimica della specie di OpenFOAM speciei (kg/m3/s), eq. 27
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

        // senza meccanismo di reazione non c'e' nulla da produrre
        if (mutationMechanism_ != "none")
        {
            tw.ref().primitiveFieldRef() =
                wCache_[foamToMutationSpecie_[speciei]];
        }

        return tw;
    }
    // ---- fine produzione chimica

    // ---- calore di reazione per l'energia sensibile e del solver (W/m3):
    // le reazioni convertono energia sensibile in energia di formazione,
    // -somma_s hf_s * omega_s (eq. 1, ultimo termine)
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
        tQ.ref().primitiveFieldRef() = QchemCache_;

        return tQ;
    }
    // ---- fine calore di reazione

    // ---- decode: dalle energie risolte ricava T_tr, T_ve e psi
    virtual void correct() override
    {
        const scalarField& eCells = this->he().primitiveField();
        const scalarField& eveCells = this->eve_.primitiveField();
        const tmp<volScalarField> tRho = this->rho();
        const scalarField& rhoCells = tRho().primitiveField();

        scalarField& TCells = this->T_.primitiveFieldRef();
        scalarField& TveCells = this->Tve_.primitiveFieldRef();

        const int ns = nMutationSpecies_;

        parallelFor
        (
            TCells.size(),
            [&](const label celli, MutationWorkspace& ws)
            {
                const double rho = rhoCells[celli];
                fillMassFractions(this->Y(), celli, ws);
                for (int s = 0; s < ns; s++)
                {
                    ws.rhoS[s] = rho*ws.yS[s];
                }

                // energie per unita' di volume: Mutation++ vuole quella totale
                // con la formazione, e quella vibro-elettronica
                const double energies[2] =
                {
                    rho*(eCells[celli] + mixtureHf(ws)),
                    rho*eveCells[celli]
                };
                // Mutation++ inverte le energie (Newton) e trova T_tr e T_ve
                ws.mix->setState(ws.rhoS.data(), energies, 0);

                TCells[celli] = ws.mix->T();
                TveCells[celli] = ws.mix->Tv();
            }
        );

        this->T_.correctBoundaryConditions();
        this->Tve_.correctBoundaryConditions();

        updatePsi();

        // lo stato e' cambiato: le sorgenti in cache sono di quello vecchio
        sourcesValid_ = false;
    }
    // ---- fine decode

};

// ---- initMutation: legge highEnthalpyMutation da physicalProperties e crea
// una miscela per thread OpenMP
template<class Thermo>
void HighEnthalpyMulticomponentThermo<Thermo>::initMutation()
{
    const dictionary& propertiesDict = this->properties();

    // legge mixture, stateModel, thermodynamicDatabase, mechanism e parkExponent
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
    }

    Mutation::MixtureOptions opts(mutationMixtureName_);
    opts.setStateModel(mutationStateModel_);
    opts.setThermodynamicDatabase(mutationThermoDatabase_);
    opts.setMechanism(mutationMechanism_);

    // una miscela per thread: costruirle costa (ognuna rilegge gli xml dei dati)
    // ma si paga una volta sola all'avvio, non a ogni passo temporale
    const int nThreads = highEnthalpyMaxThreads();
    workspaces_.resize(nThreads);

    for (int t = 0; t < nThreads; t++)
    {
        workspaces_[t] = std::make_unique<MutationWorkspace>();
        workspaces_[t]->mix = std::make_unique<Mutation::Mixture>(opts);
    }

    Mutation::Mixture& mix0 = *workspaces_[0]->mix;
    nMutationSpecies_ = mix0.nSpecies();
    nMutationEnergyEqns_ = mix0.nEnergyEqns();

    // collega ogni specie di OpenFOAM alla stessa specie di Mutation++
    const wordList& foamSpecies = this->species();
    foamToMutationSpecie_.setSize(foamSpecies.size(), -1);

    forAll(foamSpecies, i)
    {
        foamToMutationSpecie_[i] = mix0.speciesIndex(foamSpecies[i].c_str());

        if (foamToMutationSpecie_[i] < 0)
        {
            FatalErrorInFunction
                << "la specie " << foamSpecies[i]
                << " non esiste nella miscela Mutation++ " << mutationMixtureName_
                << exit(FatalError);
        }
    }

    // molecole che vibrano e loro tempi di rilassamento, e buffer di lavoro:
    // legati alla miscela del proprio thread, mai condivisi
    for (int t = 0; t < nThreads; t++)
    {
        MutationWorkspace& ws = *workspaces_[t];
        ws.vibrators = makeVibrators(*ws.mix);
        ws.yS.assign(nMutationSpecies_, 0.0);
        ws.rhoS.assign(nMutationSpecies_, 0.0);
        ws.eS.assign(2*nMutationSpecies_, 0.0);
        ws.wS.assign(nMutationSpecies_, 0.0);
    }

    // entalpie di formazione (J/kg): hf/RT da Mutation++, a una T qualsiasi
    const double T0 = 1000.0;
    std::vector<double> hfOverRT(nMutationSpecies_, 0.0);
    mix0.speciesHOverRT
    (
        T0, T0, T0, T0, T0,
        nullptr, nullptr, nullptr, nullptr, nullptr, hfOverRT.data()
    );
    hf_.assign(nMutationSpecies_, 0.0);
    for (int s = 0; s < nMutationSpecies_; s++)
    {
        hf_[s] = hfOverRT[s]*Mutation::RU*T0/mix0.speciesMw(s);
    }
}
// ---- fine initMutation

// ---- updatePsi: psi = rho/p in celle e facce di bordo, con rho da Mutation++
// per p, T, Tve e composizione date (eq. 24: p somma delle pressioni parziali)
template<class Thermo>
void HighEnthalpyMulticomponentThermo<Thermo>::updatePsi()
{
    const scalarField& pCells = this->p_.primitiveField();
    const scalarField& TCells = this->T_.primitiveField();
    const scalarField& TveCells = this->Tve_.primitiveField();
    scalarField& psiCells = this->psi_.primitiveFieldRef();

    parallelFor
    (
        psiCells.size(),
        [&](const label celli, MutationWorkspace& ws)
        {
            fillMassFractions(this->Y(), celli, ws);
            const double pTTv[3] =
                {pCells[celli], TCells[celli], TveCells[celli]};
            ws.mix->setState(ws.yS.data(), pTTv, 2);
            psiCells[celli] = ws.mix->density()/pCells[celli];
        }
    );

    volScalarField::Boundary& psiBf = this->psi_.boundaryFieldRef();
    forAll(psiBf, patchi)
    {
        const fvPatchScalarField& pp = this->p_.boundaryField()[patchi];
        const fvPatchScalarField& pT = this->T_.boundaryField()[patchi];
        const fvPatchScalarField& pTve = this->Tve_.boundaryField()[patchi];
        fvPatchScalarField& ppsi = psiBf[patchi];

        parallelFor
        (
            ppsi.size(),
            [&](const label facei, MutationWorkspace& ws)
            {
                std::fill(ws.yS.begin(), ws.yS.end(), 0.0);
                forAll(foamToMutationSpecie_, s)
                {
                    ws.yS[foamToMutationSpecie_[s]] =
                        this->Y()[s].boundaryField()[patchi][facei];
                }
                const double pTTv[3] = {pp[facei], pT[facei], pTve[facei]};
                ws.mix->setState(ws.yS.data(), pTTv, 2);
                ppsi[facei] = ws.mix->density()/pp[facei];
            }
        );
    }
}
// ---- fine updatePsi

} // End namespace Foam

#endif

// ************************************************************************* //

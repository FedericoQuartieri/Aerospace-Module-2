# Milestone 5 — Blunted cone Mach 11.3: prerequisiti, parallelo e primo 2D

**Data:** 5-6 luglio 2026
**Branch:** `milestone-5-blunted-cone`
**Obiettivo:** portare il solver sul primo caso 2D del paper Part Two (blunted cone Mach 11.3, §3.1) — passando dai due prerequisiti individuati alla fine della M4: la consistenza del datum energetico a parete e il funzionamento in parallelo (per il cluster di Federico: 4 job × 28 core).

---

## 1. Prerequisito A — BC di parete a temperatura fissa (datum dell'energia)

**Il problema.** Le BC di energia di OpenFOAM (`fixedEnergy`, `mixedEnergy` — scelte automaticamente in base alla BC di T) valutano l'energia di parete chiamando `thermo.he(Tp, patchi)`. L'implementazione base usa i polinomi del thermo (datum janaf, zero a 298 K), mentre dal refactor M3 il campo `e` vive sul datum Mutation (`e_tr+e_ve`, da 0 K): per N₂ lo scarto è **~3×10⁵ J/kg**. Su una parete isoterma la decodifica `T(e−eve)` avrebbe prodotto temperature di parete completamente sbagliate — e con esse il flusso termico, che è LA grandezza da validare sul cono. Nei casi precedenti non si vedeva solo perché tutte le pareti erano `zeroGradient`.

**Il fix.** Override nel thermo di `he(T, patchi)` e `he(T, cells)` sul datum Mutation: `e = e_tr(T) + eve` (usando il campo `eve` corrente, così il decode restituisce esattamente la T imposta). Un dettaglio C++: gli overload nuovi nascondono gli altri `he()` della gerarchia → serve `using FluidMulticomponentThermo<Thermo>::he;`.

**Il test.** Nuovo caso `applications/test/wallSlab1D/`: slab 1D di N₂ a 1000 K con parete isoterma a 300 K, viscoso. Risultato: **T di parete decodificata = 300.0002 K** (errore 0.2 mK), profilo di raffreddamento monotono. PASS.

## 2. Prerequisito B — parallelo MPI

Il bridge Mutation++ è per-cella (ogni rank costruisce la propria `Mixture` e lavora sulle celle locali), quindi la decomposizione di dominio standard doveva funzionare — ma è stato **verificato**, non assunto: il caso `shockTube1D/relaxing` (2000 celle, chimica attiva) girato seriale e con `mpirun -np 4` (nuovo flag `./Allrun -parallel [n]`) produce:

| Campo | Max diff relativa seriale vs parallelo |
|---|---|
| T, Tve, p | ~1–2×10⁻⁵ (roundoff MPI atteso) |
| N | 2.7×10⁻³ relativo, ma su valori ~10⁻⁷ (assoluto trascurabile) |

Metriche di validazione identiche a tutte le cifre stampate. **Il codice è parallelo e pronto per il cluster.**

## 3. Bug sistemati lungo la strada

1. **`setRDeltaT` dimensionalmente rotto** (mai esercitato: nessun caso usava LTS prima): la formula era un ibrido tra la versione di shockFluid (`amaxSf` volumetrico) e quella di multicomponentFluid (flusso di massa, con `rho`): il `*rho()` extra rompeva le dimensioni. Rimosso.
2. **`smoluchowskiJumpT` legge `mixture/transport/Pr`** hard-coded dal layout pureMixture (stile rhoCentralFoam): con il thermo multicomponente il subdict non esiste → aggiunto al `physicalProperties` del caso un subdict `mixture { transport { Pr 0.72; } }` che il thermo multicomponente ignora (documentato nel dict).
3. **`limitTemperature` in OF-13 richiede `cellZone all;`**.
4. Avvio impulsivo a M11: con maxCo 0.4 il transitorio esplode (FPE nello smoother di U); **maxCo 0.15** con `rDeltaTDampingCoeff 0.9` + `limitTemperature` [50, 60000] K è stabile. Il restart con maxCo 0.4 su campo sviluppato è ANCH'ESSO instabile su questa mesh coarse: restare a 0.15.

## 4. Il caso `bluntedCone2D`

Geometria del paper (§3.1): naso sferico Rn=6.35 mm + cono 25°, estensione streamwise 5 cm, assialsimmetrico (wedge 4°). La mesh è **generata da script** (`makeBlockMeshDict.py coarse|fine`):
- **coarse**: (40+80)×60 = 7200 celle, prima cella ~25 µm — per la prova locale della pipeline;
- **fine**: (200+400)×200 = 120k celle, prima cella ~2 µm — grado-paper (600×200), per il cluster.

Condizioni (Tabella 1 del paper): U∞=2764.5 m/s, p∞=21.9139 Pa, T∞=Tve∞=144.4 K, N₂ puro non-reagente. Parete: `maxwellSlipU` + `smoluchowskiJumpT` (accomodazione 1.0, Tw=297.2 K), Tve a parete `zeroGradient` (coerente con l'assenza di diffusione del pool ve nella EveEqn — limite noto, irrilevante qui: la vibrazione resta quasi congelata). Regime stazionario via **LTS** (`localEuler`).

### Risultati coarse, a convergenza (vedi §6 per la storia della convergenza)

| Grandezza | Valore | Riferimento |
|---|---|---|
| Picco T sulla linea di ristagno | **25.6×T∞** | salto RH frozen teorico M11.3: **25.8** ✓ |
| Standoff dello shock | 1.01 mm (0.159 Rn) | correlazioni strong-shock 0.1–0.15 Rn ✓ (smearing coarse) |
| Tve massima | ~460 K | "vibrational mode barely excited" (paper §3.1) ✓ |
| ρ a parete | ~20×ρ∞ | strato limite freddo, fisico ✓ |
| Cp di ristagno | **1.30 (convergente, non transitorio)** | Rayleigh pitot 1.83 — **29% sotto, da capire con mesh fine (§6)** |

Nota metodologica: lo standoff va misurato sulla **temperatura** (half-rise), non sul gradiente di densità — il picco di ρ nello strato limite a parete fredda domina il gradiente e falsa la misura.

## 5. Pacchetto cluster (`bluntedCone2D/cluster/`)

- `job-cone-fine.sh`: template PBS/Torque (`qsub`) 28 core, adattato al cluster di Federico (non SLURM)
- `README.md`: setup una tantum (build Mutation++ da `thirdParty/` — **gitignorata, va copiata a mano** —, catena Allwmake), stima tempi, nota sul dynamicCode (serve g++ almeno sul nodo di lancio), cosa riportare indietro.
- I 4 job disponibili si prestano a varianti in parallelo (mesh/maxCo/accomodazione).

## 6. Run coarse portato a convergenza: risultato onesto

Il run è stato completato fino a `endTime` (pseudo-tempo 3×10⁻³, interrotto una volta da uno spegnimento imprevisto e ripreso senza perdita — l'ultimo checkpoint scritto era integro su tutti i rank). **Verificato che la soluzione è davvero stazionaria**, non solo terminata per timeout: i campi T e p differiscono meno dello 0.5% tra pseudo-tempo 10⁻³ e 3×10⁻³, coerente con i residui che oscillano in un ciclio limite (~5×10⁻⁶) senza scendere ulteriormente.

**Il Cp di ristagno converge a 1.304 — non a ~1.8.** Non è un run incompleto: è il risultato vero di questa mesh coarse, e resta **~29% sotto** la stima Rayleigh-Pitot ideale (1.833). Questo va segnalato onestamente, non minimizzato. La causa più probabile è la qualità di mesh: `checkMesh` riporta non-ortogonalità massima **34.6°** (media 16.1°), su una griglia solo 40×60 con stiramento radiale 30:1 — abbastanza da degradare sensibilmente la ricostruzione dei flussi vicino al naso curvo, e lo schema Kurganov (dissipativo) su risoluzione circonferenziale grossa può ulteriormente smerare il recupero di pressione post-shock. Il resto della fisica (T di picco 25.6×T∞ vs 25.8 teorico, standoff 0.159 Rn, Tve quasi congelata) resta in ottimo accordo — è **specificamente** la pressione di parete, la grandezza più sensibile alla risoluzione locale del naso, a essere fuori.

**Conseguenza pratica**: il run fine sul cluster (200×200 in più, non-ortogonalità attesa molto minore, prima cella 2 µm) non è solo "più preciso" — è il test che stabilisce se questo è davvero un problema di mesh (il Cp dovrebbe salire verso ~1.8 infittendo) o qualcosa di più profondo nel bridge/schema. Se il Cp fine restasse basso, andrebbe indagato il termine di pressione nel solver o lo schema di ricostruzione vicino a pareti curve non ortogonali.

## 7. Stato e prossimi passi

- Pipeline 2D **dimostrata end-to-end** in locale (mesh generata → run parallelo LTS → postProcess con metriche fisiche corrette), inclusa la resilienza a un'interruzione (checkpoint LTS, ripresa senza perdita).
- Run coarse **concluso**: fisica di shock/temperatura validata qualitativamente; Cp di parete **sotto stima del 29%**, imputato a mesh coarse/non-ortogonalità — da confermare/smentire con la mesh fine.
- **Da fare da Federico**: run fine (`./Allrun -fine 28`) sul cluster via `qsub`; controllare **per primo** se il Cp di ristagno sale verso 1.8 (test diagnostico prioritario, non solo un affinamento); poi digitalizzazione dei riferimenti del paper (Fig 2: linea di ristagno Wang&Boyd/MONACO, Cp/Cf/St con esperimenti CUBRC run 31) per il confronto quantitativo — pattern Engauge già usato per la Fig 5.
- Rimandati noti: diffusione del pool ve nella EveEqn (κ_ve separata, eq. 7-8 Part Two), cilindro Mach 20 reagente (secondo caso del paper), Fig 9/CVDV/Tve multiple.

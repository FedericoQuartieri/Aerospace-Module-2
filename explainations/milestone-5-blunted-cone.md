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

**Il Cp di ristagno converge a 1.304 — non a ~1.8.** Non è un run incompleto: è il risultato vero di questa mesh coarse, e resta **~29% sotto** la stima Rayleigh-Pitot ideale (1.833). Questo va segnalato onestamente, non minimizzato.

**Prima ipotesi (SBAGLIATA, poi corretta):** inizialmente attribuita alla non-ortogonalità di mesh (34.6° max riportata da `checkMesh`). Verifica diretta: la mesh *fine* (120k celle) ha **la stessa identica non-ortogonalità** (34.9° max) — refinement non la riduce, quindi è una proprietà geometrica del blocco, non della risoluzione. Localizzando le 342 facce oltre 30° (via `checkMesh -writeSets` + lettura di owner/faces/points), risultano **tutte nell'angolo esterno a valle del blocco cono** (x=29–47 mm, r=41–55 mm) — lontano dal naso e dal ristagno. Non può essere la causa del Cp basso: è geometricamente irrilevante per quella regione.

**Diagnosi corretta, per ispezione diretta del profilo di pressione** cella-per-cella lungo la linea di ristagno (dalla parete fino a valle dello shock): il Cp non supera mai **1.35** in nessun punto — sale con continuità dal fronte d'urto verso parete, senza il salto quasi discontinuo atteso dietro un bow shock forte seguito da un plateau in strato limite. Lo shock è **spalmato su gran parte del millimetro di stand-off**, non catturato come discontinuità netta. Causa: solo **~20 celle radiali** coprono l'intero stand-off su questa mesh coarse (60 celle radiali totali, stiramento 30:1) — pochissimo per uno schema centrale (Kurganov, intrinsecamente dissipativo) che tipicamente smeraggia un urto su 3–5 celle. Il resto della fisica (T di picco 25.6×T∞ vs 25.8 teorico, standoff 0.159 Rn, Tve quasi congelata) resta in ottimo accordo — è specificamente la **cattura dell'urto/recupero di pressione**, sensibile alla risoluzione radiale locale, a essere insufficiente.

**Conseguenza pratica, favorevole al run fine**: con 200 celle radiali (stiramento 136:1, prima cella 2.3 µm) la stessa distanza di stand-off (~1 mm) è coperta da **~100 celle** — stimato da `r=136^(1/199)≈1.025` e la somma della serie geometrica fino a 1 mm — cioè **5 volte più risoluzione** proprio dove serve. Il run fine non è solo "più preciso": è il test diretto se questo è un problema di risoluzione del blocco d'urto (il Cp dovrebbe salire sensibilmente verso ~1.8) o qualcosa di più radicato nello schema numerico. Se restasse comunque basso con 100 celle nello strato d'urto, andrebbe indagato lo schema di ricostruzione/i limitatori (vanAlbada) piuttosto che la mesh.

## 7. Il run fine sul cluster e la risoluzione del giallo del Cp

Il run fine (120k celle, 28 core, coda `cpu` del cluster PBS — dopo la trafila di ambiente documentata nel `cluster/README.md`: sourcing di `/opt/openfoam13`, bug `pop_var_context` aggirato con sourcing in sottoprocesso, log auto-gestito via `exec`) ha dato il verdetto del test diagnostico:

**Cp fine = 1.335 vs coarse 1.304** — praticamente invariato con 17× più celle e ~5× risoluzione nello strato d'urto. **Anche la seconda diagnosi (sotto-risoluzione dello shock) era quindi sbagliata**: due mesh drasticamente diverse, stesso risultato → errore sistematico, non di discretizzazione.

**La causa vera, terza e definitiva**: il bridge Mutation++ ha un clamp di sicurezza `minTemperature 200` (default nato per i casi caldi delle milestone precedenti), ma il free-stream del paper è a **T∞ = 144.4 K < 200 K**. Il clamp riscaldava silenziosamente tutto il flusso indisturbato a 200.00 K esatti (verificato sui campi convergiuti). Conseguenza: il solver simulava — correttamente! — un flusso a **Mach 9.59 invece di 11.3** (ρ∞ effettiva 3.69×10⁻⁴ invece di 5.113×10⁻⁴). La prova numerica che il solver era esatto: al Mach effettivo, il Rayleigh-Pitot dà Cp = 1.831 e il CFD misurava 1.806 — **accordo all'1.4%**. Il "29% di errore" era interamente l'artefatto del confronto tra il flusso a M9.6 (simulato) e la teoria a M11.3 (nominale).

**Fix e conferma**: `minTemperature 50` nel dict del caso, rerun coarse locale (3h su 4 core — più lungo dei run precedenti perché ora il free-stream è davvero a 144 K):

| Grandezza | Prima (clamp attivo) | Dopo il fix | Teorico |
|---|---|---|---|
| T∞ effettiva | 200.0 K | **144.4 K** | 144.4 K |
| M effettivo | 9.59 | **11.29** | 11.3 |
| Cp di ristagno | 1.304 | **1.802** | 1.833 (**err. 1.7%**) |
| Standoff | 0.159 Rn | **0.147 Rn** | 0.1–0.15 Rn |

**Guard-rail permanente**: `postProcess-cone.py` ora misura il free-stream *effettivamente simulato* (T, p, M nelle celle a monte dello shock), lo stampa accanto al nominale, avvisa esplicitamente se differiscono, e calcola ogni confronto teorico al Mach effettivo. La lezione: **mai confrontare col nominale senza verificare cosa è stato davvero simulato** — tre diagnosi sbagliate (non-ortogonalità, risoluzione shock) prima di quella giusta, tutte perché l'assunzione di base ("il free-stream è quello imposto") non era stata controllata.

## 8. Stato e prossimi passi

- Pipeline 2D **dimostrata end-to-end** locale + cluster (PBS), inclusa la resilienza a interruzioni (checkpoint LTS).
- **Cp di ristagno validato all'1.7%** dal teorico sulla mesh coarse con free-stream corretto; standoff nella banda delle correlazioni.
- **Da fare da Federico**: `git pull` sul cluster + rilancio `qsub job-cone-fine.sh` (il fix `minTemperature 50` è committato) per il run fine definitivo alle condizioni giuste; poi digitalizzazione dei riferimenti del paper (Fig 2: linea di ristagno Wang&Boyd/MONACO, Cp/Cf/St con esperimenti CUBRC run 31) per il confronto quantitativo — pattern Engauge già usato per la Fig 5.
- Rimandati noti: diffusione del pool ve nella EveEqn (κ_ve separata, eq. 7-8 Part Two), cilindro Mach 20 reagente (secondo caso del paper), Fig 9/CVDV/Tve multiple.

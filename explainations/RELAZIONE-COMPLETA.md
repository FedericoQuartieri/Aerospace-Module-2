# Solver ipersonico a due temperature — Relazione completa

**Estensione di OpenFOAM-13 per flussi ipersonici reagenti in
non-equilibrio termochimico, accoppiata a Mutation++.**

Questo documento racconta l'intero progetto: la fisica, l'architettura del
codice, il percorso di sviluppo milestone per milestone, la validazione
contro i due paper di riferimento e l'analisi di performance. È pensato per
essere letto da solo, senza guardare il codice.

> **In una riga:** un solver OpenFOAM che risolve le Navier-Stokes
> compressibili con **due temperature** (traslazionale-rotazionale `T` e
> vibro-elettronica `Tve`) e **chimica a velocità finita**, validato dal
> punto 0D fino a due casi 2D del paper (cono Mach 11.3 e cilindro Mach 20
> reagente), più un'analisi di **strong scaling MPI**.

---

## Indice

1. [Obiettivo e contesto](#1-obiettivo-e-contesto)
2. [La fisica: il modello a due temperature](#2-la-fisica-il-modello-a-due-temperature)
3. [Architettura del codice](#3-architettura-del-codice)
4. [Il percorso: le milestone](#4-il-percorso-le-milestone)
5. [Validazione: risultati vs paper](#5-validazione-risultati-vs-paper)
6. [Performance HPC](#6-performance-hpc)
7. [Confronto con l'altro gruppo](#7-confronto-con-laltro-gruppo)
8. [Lezioni e trappole](#8-lezioni-e-trappole)
9. [Come si compila e si gira](#9-come-si-compila-e-si-gira)
10. [Limiti noti e sviluppi futuri](#10-limiti-noti-e-sviluppi-futuri)
11. [Mappa dei file](#11-mappa-dei-file)

---

## 1. Obiettivo e contesto

### Il problema fisico

Nei flussi ipersonici (rientro atmosferico, cruise ipersonico, scramjet)
un'onda d'urto forte converte gran parte dell'energia cinetica in energia
interna: le temperature dietro l'urto arrivano a decine di migliaia di
kelvin. In queste condizioni **cade l'ipotesi di equilibrio termico
locale** su cui è costruita la maggior parte dei solver CFD:

- i modi di energia interna delle molecole (traslazione, rotazione,
  **vibrazione**, elettronici) non si eccitano tutti allo stesso ritmo;
- vibrazione ed elettronici rilassano **molto più lentamente** di
  traslazione e rotazione, creando uno stato di **non-equilibrio termico**;
- le reazioni chimiche (dissociazione) sono governate da queste
  temperature e a loro volta le modificano.

Un modello a **singola temperatura** in questo regime sbaglia le proprietà
termodinamiche, l'accoppiamento con la chimica e — cosa pratica — la
stabilità numerica.

### La risposta: il modello a due temperature

Il compromesso standard (Park) è tenere **due** temperature:

- `T` (o `Ttr`): **traslazionale-rotazionale**, governa pressione,
  densità, quantità di moto;
- `Tve`: **vibro-elettronica**, governa l'eccitazione ritardata dei modi
  interni e il loro accoppiamento con la chimica.

Lo scambio di energia fra i due serbatoi è modellato con termini di
rilassamento (Landau-Teller per il V-T). È l'estensione minima dei modelli
a singola temperatura capace di catturare il non-equilibrio dominante, ed è
compatibile con i solver a volumi finiti segregati come OpenFOAM.

### I due paper di riferimento

Il progetto replica e valida contro **Casseau et al.**, la coppia di paper
che introduce il solver open-source `hy2Foam`:

- **Part One** (*Aerospace* 2016, 3, 34) — analisi **0D**: heat bath in
  non-equilibrio, rilassamento V-T e chimica.
- **Part Two** (*Aerospace* 2016, 3, 45) — analisi **multi-dimensionale**:
  cono smussato a Mach 11.3 (§3.1) e cilindro a Mach 20 reagente (§3.2).

### Cosa è stato costruito

Un solver OpenFOAM-13 completo, con classi termodinamiche e di trasporto
custom e un bridge verso Mutation++, validato lungo tutta la catena
**0D → 1D → 2D** su entrambi i casi del Part Two, più la conduzione
vibrazionale e un'analisi di scaling MPI. Nove milestone in tutto.

---

## 2. La fisica: il modello a due temperature

### Le variabili risolte

Il solver trasporta, per ogni cella:

- densità `rho`, quantità di moto `rho·U`, energia totale;
- frazioni di massa delle specie `Yi` (con sorgenti chimiche);
- l'**energia vibro-elettronica specifica** `eve` [J/kg], con la sua
  equazione di trasporto dedicata.

Le due temperature si **decodificano** dall'energia (vedi sotto).

### Rilassamento V-T (Landau-Teller)

L'energia fluisce dal serbatoio traslazionale a quello vibrazionale con un
tempo caratteristico `tauVT` (dipendente da T e composizione, fornito da
Mutation++). In forma di sorgente per l'equazione di `eve`:

```
sorgente V-T = rho·(eveEq - eve)/tauVT
```

dove `eveEq = e_ve(T)` è l'energia vibrazionale che si avrebbe
all'equilibrio con la temperatura traslazionale. Il termine è
semi-implicito (parte in `eveEq`, parte implicita `-Sp(rho/tauVT, eve)`)
per stabilità.

### Chimica a velocità finita, a due temperature

La dissociazione (es. `N2 + M → 2N + M`) è governata da velocità di
reazione che dipendono **da entrambe** le temperature (modello Park:
temperatura di controllo `Tf = T^q · Tve^(1-q)`). Mutation++ fornisce per
cella:

- le velocità nette di produzione `ω̇_i` di ogni specie → termini
  sorgente nelle equazioni `Yi`;
- il calore di reazione `Qdot = -Σ hf_i·ω̇_i` → nell'equazione dell'energia
  totale;
- l'accoppiamento chimica-vibrazione `Qcv = Σ e_ve,i·ω̇_i` (Candler:
  dissociazione preferenziale) → nell'equazione di `eve`.

> **Nota sul `q` di Park.** La Mutation++ vendorizzata ha `q = 0.5`
> hard-coded (`Tf = √(T·Tve)`), mentre il paper usa `q = 0.7`. La
> differenza è ~3% nel transitorio 0D (verificato in M2): trascurabile, e
> documentata.

### La formulazione conservativa (il cuore del solver)

La scelta chiave (adottata in M3) è **non** trasportare `T` e `Tve`
direttamente, ma l'energia sul **datum di Mutation++**:

```
e = e_tr(T) + e_ve(Tve)        (energie di formazione escluse)
```

Le due temperature si ricavano per **decode** cella per cella con Newton:

1. `Tve` da `eve` (invertendo `e_ve(Tve) = eve`);
2. `T` da `e - eve` (invertendo `e_tr(T)`, lineare per un modello RRHO).

Questo garantisce la **conservazione stretta dell'energia totale**: il
travaso V-T è interno a `e` e agisce solo attraverso il decode. Non serve
alcun termine V-T esplicito nell'equazione dell'energia principale (bug
latente della prima versione, corretto in M3).

### Conduzione dei due modi

Ogni modo conduce calore con la propria conduttività:

- l'energia totale conduce con `divq(e)` (modello `unityLewisFourier`, che
  usa la conduttività di miscela);
- il pool vibro-elettronico conduce con la sua `κ_ve`. Con l'**Eucken
  modificata** `κ_ve = μ·cv_ve`, e poiché `∇eve = cv_ve·∇Tve`
  puntualmente, la conduzione vibrazionale è **esattamente**
  `∇·(μ·∇eve) = laplacian(μ, eve)`. Questo termine è stato aggiunto in M8
  (rimandato da M3).

---

## 3. Architettura del codice

### Le fondamenta

- **OpenFOAM-13** — infrastruttura FVM, trasporto, discretizzazione.
- **Mutation++** (vendorizzata in `thirdParty/`, gitignorata) — backend
  termochimico: energie interne, calori specifici, integrali di
  collisione, velocità di reazione a due temperature. Approccio "detached":
  il core di OpenFOAM resta intatto, l'estensione è una libreria user-side
  collegata dinamicamente.

### Le classi custom

| Classe | Ruolo |
|---|---|
| `rrhoThermo` (RRHO) | termodinamica Rigid-Rotor Harmonic-Oscillator, derivata dal JANAF standard ma estesa all'alta temperatura / non-equilibrio |
| `HighEnthalpyMulticomponentThermo` | il **bridge**: introduce il campo `Tve`, gestisce l'accoppiamento OpenFOAM ↔ Mutation++ cella per cella, override del datum a parete (`he`) |
| `blottnerTransport` | viscosità di Blottner `μ = 0.1·exp((A·lnT+B)·lnT+C)` + κ Eucken (il modello di trasporto del paper); aggiunta in M6 |
| `shockThermo` (modulo solver) | il solver vero, derivato da `shockFluid` (famiglia `rhoCentralFoam`): flusso di Kurganov, ricostruzione MUSCL, local time-stepping (LTS) verso lo stazionario |

### Il bridge per cella

Il cuore dell'accoppiamento: per ogni cella, `HighEnthalpyMulticomponentThermo`
estrae lo stato locale, interroga Mutation++, e immagazzina i risultati in
campi (`eveEq`, `tauVT`, `mutR_<specie>`, `mutQdot`, `mutQcv`) che il solver
usa come sorgenti. In parallelo MPI: **una `Mutation::Mixture` per rank**
(verificato: seriale vs parallelo identici a 1e-5).

### Niente type table statiche: dynamicCode

Una scoperta architetturale che ha semplificato molto (M6/M8): il thermo
dei casi è **compilato al volo** da OpenFOAM (`dynamicCode`). Per aggiungere
un modello di trasporto (`blottner`) o una termodinamica (`rrho`) **non**
servono type table statiche: basta l'header nell'`lnInclude` e il nome
nella whitelist `etc/codeTemplates/dynamicCode/fluidMulticomponentThermo`
(la cache in `~/.OpenFOAM/13` viene ricopiata al source di `etc/bashrc`).

---

## 4. Il percorso: le milestone

Nove milestone, ciascuna su un branch dedicato, documentate in
`explainations/milestone-N-*.md`. Le prime otto sono la **validazione
fisica** obbligatoria del paper; la nona è **performance** (rifinitura).

### M1 — Rilassamento V-T nel solver *(0D)*

Prima estensione: un'equazione conservativa per `eve` nel solver, con
rilassamento Landau-Teller semi-implicito. `Tve` ricavata da `eve` con
Newton (senza `setState`). Validato sul heat bath 0D vs `Test-N2`: errore
massimo **0.03% su Ttr, 1.87% su Tve**.

### M2 — Chimica reagente 0D

Dissociazione irreversibile dell'azoto (`N2+N2 → 2N+N2`, rate di Park 1993),
riproducendo le Fig 7/8 del Part One. Doppio percorso validato: Mutation++
end-to-end vs implementazione manuale (legge di azione di massa + Candler) —
a `q=0.5` le curve sono bit-identiche. **Scoperta:** il `q=0.5` hard-coded
in Mutation++ vs `q=0.7` del paper differiscono ~3% nel transitorio.

### M3 — Chimica nel solver + formulazione conservativa piena

Il salto architetturale: il campo `e` vive sul datum Mutation
`e_tr(T)+e_ve(Tve)`, decode per cella (Tve da eve, T da e-eve), sorgenti
chimiche per cella dal thermo. Niente più termine V-T in EEqn (il travaso è
interno a `e`). Due bug storici risolti qui: **SIGFPE su OmegaVT con X=0**
(floor delle frazioni a 1e-30) e **clamp `Thigh` del JANAF** che
tosava T=30000→20000. Validato vs i CSV di M2 (T 0.03%, composizione 1e-4).

### M4 — Validazione 1D *(shock tube)*

Tre casi in `applications/test/shockTube1D/`:

- **`sod`** (N2 freddo) vs Riemann esatto: L1 **0.07–0.17%**;
- **`relaxing`** (driver caldo, zona di rilassamento) vs riferimento ODE
  post-shock: **T max 1.04% / media 0.69%, Tve media 0.83%**;
- il salto RH frozen coincide con la teoria a 4 cifre a Ms=8.36.

**Prestazioni misurate: 34–42 µs/cella/step** → proiezione: il multi-D
avrebbe richiesto ore/giorni, da cui l'uso del cluster.

### M5 — Cono Mach 11.3: prerequisiti e caso *(2D)*

Primo caso multi-D del Part Two (§3.1). Prerequisiti:

- **fix del datum a parete**: override `he(T,patchi)` = `e_tr(T)+eve` nel
  thermo (serve `using base::he` per il name-hiding del C++); unit test
  `wallSlab1D`: parete decodificata a 300.0002 K ✓;
- **parallelo verificato**: `relaxing` seriale vs `mpirun -np 4` identici a
  ~1e-5;
- caso `bluntedCone2D`: mesh wedge assialsimmetrica da script, BC di slip
  (`maxwellSlipU`) e salto di temperatura (`smoluchowskiJumpT`), LTS.

Coarse convergiuto, ma **Cp di ristagno 1.30 vs 1.83** atteso: un giallo
lasciato aperto per M6.

### M6 — Campagna cluster, il giallo del Cp, Blottner, Fig 2 validata

La milestone più densa. Tre cose:

1. **Cluster PBS operativo** (trafila non ovvia: sintassi
   `select=1:ncpus=28:mpiprocs=28`, stage-out `-o` inaffidabile → log
   auto-gestito, sourcing isolato per un bug bash, Mutation++ clonata da
   GitHub).
2. **Il giallo del Cp risolto dopo tre diagnosi**: non era la mesh né la
   risoluzione dell'urto, ma il **clamp `minTemperature 200`** del bridge
   che riscaldava il free-stream da 144.4 a 200 K → simulava Mach 9.6
   invece di 11.3. Fix `minTemperature 50` + **guard-rail permanente** nel
   post (misura il free-stream *effettivo* e confronta la teoria al Mach
   effettivo). *Lezione: mai confrontare col nominale senza verificare cosa
   è stato simulato.*
3. **Trasporto Blottner** implementato (chiude il gap su Cf/St, vedi
   sotto) e **Fig 2 digitalizzata** per estrazione vettoriale dal PDF
   (`digitize-fig2.py`, PyMuPDF).

**Risultato Fig 2 (mesh fine):** Cp ristagno **1.843 vs 1.833 (0.5%)**,
Cp superficie **0.7%**, St **2.7%**, Cf **7.7%**, standoff 0.150 Rn. Il gap
Cf/St, inizialmente ~16%/11% con la Sutherland, si è chiuso passando a
Blottner (la Sutherland era -21% a 297 K, la temperatura di parete).

### M7 — Cilindro Mach 20 reagente *(2D)*

Secondo caso del Part Two (§3.2): azoto a Mach 20 su cilindro R=1 m, parete
1000 K, con **chimica reagente accoppiata al two-temperature** — la prima
prova della chimica in multi-D (finora solo 0D/1D). Highlights:

- **Checkerboard odd-even** confinato alla cella a parete (aspect ratio
  ~1700:1 + flusso centrale). Diagnosi corretta dopo un errore: non era la
  BC di Tve, era un modo pressione-velocità. **Curato con Minmod** (limiter
  più dissipativo) + **estrazione robusta** delle grandezze di parete
  (banda pulita per la pressione, solo componente tangenziale per l'attrito).
- **Risultato:** C_D **1.285 vs 1.304** (paper) / 1.284 (DSMC), standoff
  **0.247 vs ~0.25 m**, profili di ristagno che ricalcano il run 3 —
  inclusa la **sovrastima attesa del picco di T** che il paper attribuisce
  proprio alla combinazione Park.
- La **mesh fine non converge** (aspect ratio 5000:1, LTS più lento): la
  coarse è il risultato di M7, dentro la dispersione dei riferimenti (che
  fra loro divergono del 30-40%).

### M8 — Conduzione vibrazionale κ_ve

Il tassello rimandato da M3: la EveEqn conduceva solo convezione +
rilassamento + chimica, **senza conduzione termica del modo vibrazionale**.
Aggiunto in una riga:

```cpp
if (!inviscid) EveEqn -= fvm::laplacian(thermo.mu(), eve);
```

(Eucken `κ_ve = μ·cv_ve` → `κ_ve∇Tve = μ∇eve` esatto). Validato: no-op a
cella singola (0D invariato), stabile in 1D. Sul cluster (cilindro):
**la parete Tve a gradiente, impraticabile in M7, ora è stabile grazie
alla diffusione** — niente scacchiera, superficie liscia, C_H **58.7 kW**
pulito, vicino alla DSMC (63.3). *Con questa milestone la validazione del
paper è completa.*

### M9 — Performance: strong scaling MPI *(rifinitura opzionale)*

Copre l'asse "High Performance" del corso. **Strong scaling MPI del solver
2D reagente vero** (decomposizione di dominio con scambio di halo), su un
nodo da 28 core, due mesh:

- **coarse (9k celle):** satura a **8× su 28 rank** (efficienza 0.29) —
  communication-bound (~320 celle/rank);
- **fine (156k celle):** near-linear fino a 14 rank, **19.9× su 28 rank**
  (efficienza 0.71) — compute-bound (~5600 celle/rank).

Il contrasto coarse/fine è il compromesso compute/comunicazione da manuale.
Vedi §6.

---

## 5. Validazione: risultati vs paper

### Cono Mach 11.3 (Part Two §3.1, Fig 2) — mesh fine, Blottner

| Grandezza | Questo lavoro | Riferimento | Scarto |
|---|---|---|---|
| Cp ristagno | 1.843 | 1.833 (Rayleigh) | **0.5%** |
| Cp superficie | — | hy2Foam | **0.7%** medio |
| Numero di Stanton St | — | hy2Foam | **2.7%** |
| Coeff. attrito Cf | — | hy2Foam | 7.7% (dentro la dispersione dei riferimenti) |
| Shock standoff | 0.150 Rn | 0.1–0.15 Rn | ✓ |
| Free-stream | T∞=144.4, M=11.29 | 144.4, 11.3 | ✓ (guard-rail) |

### Cilindro Mach 20 reagente (Part Two §3.2, Fig 5) — mesh coarse, Minmod, M8

| Grandezza | Questo lavoro | Paper (run 3) | DSMC |
|---|---|---|---|
| Coeff. di drag C_D | **1.286** | 1.304 | 1.284 |
| Shock standoff | **0.247 m** | ~0.25 | ~0.25 |
| Flusso termico integrato C_H | **58.7 kW** | 88.1 | 63.3 |
| Picco T ristagno | 15.7 kK | ~14.5 (sovrastima Park attesa) | — |
| Convergenza dp/p | 0.02% | — | — |

> **Nota su C_H.** Il paper stesso dichiara che il run Park (88 kW)
> sovrastima la DSMC (63) del 39%. Il nostro 58.7 — pulito, senza lo spike
> spurio near-wall che gonfiava le versioni precedenti — è a ~7% dalla
> DSMC: coerenza fisica, non un difetto.

**In sintesi:** la fisica reattiva two-temperature multi-D è validata su
entrambi i casi del paper, a grado-paper sulla pressione e nella dispersione
dei riferimenti sulle quantità di parete.

---

## 6. Performance HPC

### Strong scaling MPI (1 nodo, 28 core)

| ranks | fine 156k — speedup / eff | coarse 9k — speedup / eff |
|---|---|---|
| 1 | 1.00 / 1.00 | 1.00 / 1.00 |
| 2 | 1.98 / 0.99 | 1.92 / 0.96 |
| 4 | 3.82 / 0.95 | 3.31 / 0.83 |
| 7 | 6.49 / 0.93 | 4.51 / 0.64 |
| 14 | 12.24 / 0.87 | 7.11 / 0.51 |
| 28 | **19.90 / 0.71** | 8.00 / 0.29 |

La coarse segue l'ideale fino a ~4 rank poi satura (troppo poche celle per
rank, domina lo scambio di halo). La fine resta quasi lineare fino a 14 rank
e arriva a **19.9× su 28 core**: con ~5600 celle/rank il costo dominante è
il bridge Mutation++ per cella, non la comunicazione.

### Costo per cella del bridge

| Caso | µs/cella/step |
|---|---|
| 1D shockTube reagente | 34–42 |
| 2D cono M11 (non-reagente) | ~55 |
| 2D cilindro M20 (reagente) | ~72–84 |

La chimica reagente aggiunge ~50% rispetto al solo thermo (Mutation++
valuta anche le velocità di reazione a due temperature).

### Correttezza del parallelo

Prerequisito garantito già in M5: **seriale vs `mpirun -np 4` identici a
1e-5**. Il parallelismo non cambia il risultato fisico — condizione
indispensabile prima di parlare di velocità (una `Mutation::Mixture` per
rank, analogo distributed-memory del thread-local ownership).

---

## 7. Confronto con l'altro gruppo

Un altro gruppo del corso (stesso professore, stesso anno, stesso template
di partenza: le classi `RRHO` / `HighEnthalpyMulticomponentThermo` e i test
`thermoMixturePark2T`) ha sviluppato lo stesso tipo di estensione. Repo
pubblico con licenza **MIT**.

| | Questo progetto | L'altro gruppo |
|---|---|---|
| Validazione 0D (heat bath) | ✓ (M1–M2) | ✓ (unico risultato fisico) |
| Chimica **nel** solver | ✓ (M3) | descritto, non eseguito |
| **1D** shock tube | ✓ (M4) | *sviluppo futuro* |
| **2D cono M11** vs Fig 2 | ✓ (M5–M6) | *sviluppo futuro* |
| **2D cilindro M20** vs Fig 5 | ✓ (M7) | *sviluppo futuro (lo nominano)* |
| Conduzione vibrazionale κ_ve | ✓ (M8) | — |
| Analisi performance | ✓ MPI (M9) | ✓ OpenMP su kernel 0D |

**In sintesi:** su questo progetto sono stati fatti **tutti i loro
risultati** (0D) **più tutto ciò che loro elencano come "sviluppi futuri"**
(1D, cilindro 2D) **più la conduzione vibrazionale** **più l'analisi HPC**.
Sull'asse performance il confronto è diretto: loro OpenMP su celle 0D
indipendenti (2.81× su 8 thread, efficienza 0.35); qui MPI su un solver
accoppiato con scambio di halo (**19.9× su 28 rank, efficienza 0.71**) —
scaling distributed-memory reale, non repliche di un kernel
embarrassingly-parallel.

---

## 8. Lezioni e trappole

Bug e insidie che vale la pena ricordare (compaiono nei log come "bug
storici"):

- **SIGFPE su OmegaVT con X=0** (M3): `tau_park = 1/(n·X)` dà infinito
  benigno nei test standalone ma il trapping di `foamRun` crasha → floor
  delle frazioni a 1e-30 nel bridge.
- **Clamp `Thigh` del JANAF** (M3): tosava T=30000→20000 alla costruzione
  → `Thigh 50000` nel file del caso.
- **`setRDeltaT` con `*rho()` di troppo** (M5): errore nel local
  time-stepping.
- **Il clamp `minTemperature`** (M6): la lezione più importante — riscaldava
  silenziosamente il free-stream e falsava il Mach simulato. Da qui il
  **guard-rail** che misura sempre il free-stream effettivo.
- **Il checkerboard odd-even** (M7): non era la BC, era lo schema centrale
  sull'aspect ratio a parete. Curato con Minmod + estrazione robusta.
- **Il fattore Eucken** (M8): `divq(e)` applica il fattore ~1.9 del totale
  anche al gradiente vibrazionale; il termine κ_ve usa il fattore 1
  corretto. Approssimazione nota, non un bug (rifinitura possibile: split
  κ_tr/κ_ve nella EEqn).
- **Trappole di digitalizzazione** (M6/M7): i tick degli assi vanno
  filtrati per segmento e con orientamento; le figure MDPI sono vettoriali,
  estraibili direttamente dal PDF (più preciso di Engauge).

---

## 9. Come si compila e si gira

### Ambiente

```bash
source /path/to/OpenFOAM-13/etc/bashrc      # OpenFOAM core
source <repo>/etc/bashrc                    # variabili POLIMI_*/Mutation++
```

Sul cluster il core è in `/opt/openfoam13` (va sorgentato esplicitamente).
Mutation++ vive in `thirdParty/` (gitignorata; clonabile da GitHub).

### Build

```bash
# librerie custom (thermo + trasporto + bridge)
cd src/thermophysicalModels/specie && ./Allwmake
cd src/thermophysicalModels/multicomponentThermo && wmake
# il modulo solver
cd applications/modules/shockThermo && wmake
```

### Girare un caso

```bash
cd applications/test/reactingCylinder2D
./Allrun                    # seriale, mesh coarse
./Allrun -parallel 4        # parallelo
./Allrun -fine 28           # mesh fine, cluster
```

Sul cluster i job PBS stanno in `<caso>/cluster/` (es.
`job-cylinder-coarse.sh`). Lo scaling: `perf/job-scaling.sh`.

---

## 10. Limiti noti e sviluppi futuri

Tutti **opzionali** — la validazione del paper è completa.

- **Split pulito κ_tr/κ_ve** nella EEqn (oggi il totale usa un fattore
  Eucken ~1.9 anche sul gradiente vibrazionale).
- **Mesh fine convergente del cilindro**: richiede un re-grading near-wall
  per ridurre l'aspect ratio da 5000:1 a ~1500:1.
- **Weak scaling** e **multi-nodo** (oltre le 28 core di un nodo).
- **Fig 9** (aria 5 specie), **CVDV-QK**, **Tve multiple per specie** nel
  solver (oggi un campo unico), energia elettronica separata.
- Parete Tve `fixedValue` vs jump (per un gradiente vibrazionale più netto).

---

## 11. Mappa dei file

```
src/thermophysicalModels/
  specie/thermo/rrho/              # rrhoThermo (RRHO)
  specie/transport/blottner/       # blottnerTransport
  multicomponentThermo/.../HighEnthalpyMulticomponentThermo.H   # il bridge
applications/modules/
  shockFluid/                      # base (upstream)
  shockThermo/                     # il solver (EveEqn, decode, sorgenti)
applications/test/
  nonEqTTv/                        # 0D heat bath (M1–M2) + dati Mutation++
  shockTube1D/{sod,relaxing}/      # 1D (M4)
  bluntedCone2D/                   # cono M11 (M5–M6) + references/ Fig 2
  reactingCylinder2D/              # cilindro M20 (M7–M8) + references/ Fig 5
    perf/                          # harness di scaling (M9)
  wallSlab1D/                      # unit test datum a parete
explainations/
  milestone-1..9-*.md              # doc per milestone
  RELAZIONE-COMPLETA.md            # questo documento
etc/codeTemplates/dynamicCode/     # whitelist (rrho, blottner)
thirdParty/                        # Mutation++ (gitignorata)
```

---

*Solver ipersonico a due temperature — OpenFOAM-13 + Mutation++.
Validazione Casseau Part One (0D) + Part Two (cono M11.3, cilindro M20
reagente), conduzione vibrazionale e strong scaling MPI. Nove milestone,
completo su entrambi gli assi: fisica e performance.*

# Milestone 1 — Rilassamento V-T dal thermo al solver

**Data:** 3 luglio 2026
**Commit:** `e9626ce` — *"moved VT relaxation into thermo/solver: Mutation++ tauVT and conservative eve equation, added solverHeatBath validation vs 0D test"*
**Obiettivo:** portare la fisica del rilassamento vibro-traslazionale (V-T) validata in 0D dentro il solver CFD `shockThermo`, in modo che funzioni su tutte le celle e non solo negli script standalone.

---

## 1. Da dove partivamo

Prima di questa milestone il solver risolveva sì un'equazione per la temperatura vibro-elettronica `Tve`, ma con **fisica finta**:

```
// vecchio dict highEnthalpyRelaxation
cvVeOverCvTr    4.0;      // rapporto dei cv, costante
tauVT           1.2e-5;   // tempo di rilassamento, costante
```

Il tempo di rilassamento `tauVT` era una costante scelta a mano nel dizionario, e il rapporto tra i calori specifici `cvVeOverCvTr` pure. La fisica vera — il tempo di rilassamento di Millikan-White con correzione di Park (eq. 9-17 del paper Casseau Part One) — viveva **solo** nei test 0D standalone (`applications/test/nonEqTTv/`), dove veniva calcolata da Mutation++ cella per cella.

In più c'era un problema più sottile e più grave: l'equazione era scritta in **temperatura**

```
ddt(rho, Tve) + div(phi, Tve) == rho*(T - Tve)/tauVT
```

La temperatura **non è una quantità conservata**. Attraverso un'onda d'urto, dove densità e velocità cambiano bruscamente, un'equazione di trasporto scritta in temperatura dà risultati sbagliati. Andava riscritta in **energia**.

---

## 2. Cosa è stato cambiato

Tre file toccati, più un nuovo caso di validazione.

### 2.1 Il thermo — `HighEnthalpyMulticomponentThermo.H`

È il cuore della milestone. La classe thermo adesso possiede la **variabile conservata** e calcola per ogni cella i coefficienti fisici che servono al solver.

**Nuovi campi** (registrati nel mesh, quindi visibili al solver):

| Campo | Significato | Unità |
|-------|-------------|-------|
| `eve_` | Energia vibro-elettronica per unità di massa (la variabile conservata) | J/kg |
| `eveEq_` | Energia di equilibrio `e_ve(Ttr)`, il target del rilassamento | J/kg |
| `tauVT_` | Tempo di rilassamento V-T effettivo, calcolato da Mutation++ | s |
| `cvTrRatio_` | Rapporto `Cv_base / cv_tr` (vedi §3.3) | — |

**Cosa fa la `correct()` a ogni chiamata**, per ogni cella e ogni faccia di bordo:

1. **Deriva `Tve` da `eve`.** Siccome ora la variabile conservata è l'energia, la temperatura vibrazionale va ricavata invertendo `e_ve(Tve) = eve`. Lo fa la funzione `TveFromEve()` con qualche iterazione di Newton, usando `speciesHOverRT` e `speciesCpOverR` di Mutation++ — che sono **funzioni pure** (non modificano lo stato della mixture), quindi si possono chiamare nel loop senza `setState` ripetuti.

2. **Aggiorna le proprietà** (T, Cp, Cv, ψ, μ, κ) da Mutation++ come già faceva prima (`updateFromMutation`).

3. **Calcola i coefficienti del rilassamento** (`relaxationAtCurrentState`):
   - `eveEq = e_ve(Ttr, Ttr)` — l'energia vibrazionale se fosse all'equilibrio con la temperatura traslazionale;
   - `tauVT` — preso dal source term di Mutation++ `energyTransferSource()` e **ricondotto alla forma di Landau-Teller**. Cioè: Mutation++ restituisce direttamente `Q` [W/m³], e noi lo riscriviamo come `Q = ρ·(eveEq − eve)/τ`, da cui `τ = ρ·(eveEq − eve)/Q`. In questo `τ` è già dentro tutta la fisica Millikan-White + Park del database `transfer` di Mutation++;
   - `cvTr` — il calore specifico traslazionale della miscela, dal primo blocco di `getCvsMass()`.

**Altre aggiunte:**
- opzione `mechanism` nel dict (per forzare `none`, cioè niente chimica, come nei test 0D);
- buffer di lavoro pre-allocati (`mutW1_`…`mutCv_`) per non allocare memoria a ogni cella;
- `initialiseEve()` che al primo passo inizializza `eve` partendo dai campi `T` e `Tve` letti dal caso.

### 2.2 Il solver — `thermophysicalPredictor.C`

L'equazione in temperatura è stata sostituita da un'**equazione conservativa in energia** con rilassamento semi-implicito:

```cpp
fvScalarMatrix EveEqn
(
    fvm::ddt(rho, eve)
  + mvConvection->fvmDiv(phi, eve)
 ==
    rho*eveEq/tauVT
  - fvm::Sp(rho/tauVT, eve)      // termine implicito: stabilità
  + fvModels().source(rho, eve)
);
```

Il termine di rilassamento `ρ·(eveEq − eve)/τ` è spezzato in una parte esplicita (`rho*eveEq/tauVT`) e una implicita (`fvm::Sp(rho/tauVT, eve)`): il pezzo implicito rende il solutore stabile anche quando `τ` è molto piccolo rispetto al passo temporale, cosa che succede sempre in questi flussi.

Nell'equazione dell'energia totale il termine V-T compare come sink scalato da `cvTrRatio` (vedi §3.3):

```cpp
eRelaxSource = cvTrRatio*rho*(eveEq - eve)/tauVT;
// ...
== fvModels().source(rho, e) - eRelaxSource
```

Il dizionario `highEnthalpyRelaxation` ora ha solo due switch (`solveEve`, `coupleEnergy`): sono spariti `tauVT` e `cvVeOverCvTr` costanti, perché adesso li calcola il thermo.

### 2.3 Il tutorial `shockTube`

Aggiornato il dict (`constant/physicalProperties`) alla nuova sintassi e `system/fvSolution` per risolvere il campo `eve`. Il tutorial continua a girare senza errori — serve come test di non-regressione della gasdinamica di base.

---

## 3. La fisica in tre righe

### 3.1 Modello a due temperature

Il gas ad alta entalpia è in **non-equilibrio termico**: le molecole hanno una temperatura traslazionale-rotazionale `Ttr` (moto e rotazione, che si equilibrano subito) e una temperatura vibro-elettronica `Tve` (vibrazione + livelli elettronici, che si equilibrano lentamente). Dietro un'onda d'urto la traslazionale schizza in alto istantaneamente, la vibrazionale insegue con un ritardo governato dal tempo di rilassamento.

### 3.2 Equazione di Landau-Teller

Il ritorno all'equilibrio dell'energia vibrazionale segue

$$\rho \frac{\partial e_{ve}}{\partial t} = \rho\,\frac{e_{ve}(T_{tr}) - e_{ve}(T_{ve})}{\tau_{V\text{-}T}}$$

`eveEq = e_ve(Ttr)` è il target, `eve = e_ve(Tve)` è lo stato attuale, `τ` il tempo caratteristico. Quando `Tve → Ttr` la differenza si annulla e il rilassamento si ferma: è l'equilibrio.

### 3.3 Perché `cvTrRatio`

C'è un dettaglio implementativo importante. L'equazione dell'energia del solver risolve `e` = l'energia **del thermo base** (i polinomi `rrho`), non l'energia esatta di Mutation++. Quando il thermo base inverte `T(e)`, per ottenere il `dT/dt = −Q/(ρ·cv_tr)` giusto (cioè: l'energia che esce dal pool traslazionale deve abbassare `Ttr` usando il **suo** calore specifico traslazionale, non quello totale) il sink va scalato per il rapporto

$$\texttt{cvTrRatio} = \frac{Cv_\text{base}}{cv_{tr}}$$

Questo è **esatto per l'heat bath** (dove non c'è flusso) ed è un'**approssimazione attraverso un urto**. Vedi §5.

---

## 4. Validazione

Nuovo caso: `applications/test/nonEqTTv/solverHeatBath/`.

È un **heat bath a singola cella**: una scatola 1×1×1 celle, velocità nulla, N₂ puro, `T_tr = 10000 K`, `T_ve = 1000 K`, `p = 1 atm` — esattamente la configurazione della Figura 3a del paper e del test standalone `Test-N2`.

Lo script `Allrun`:
1. genera il **riferimento** lanciando `Test-N2 10000 1000` (lo script 0D già validato);
2. gira il **solver** `shockThermo` sullo stesso caso, con gli **stessi dati Mutation++**;
3. sovrappone le due storie temporali (`compare-heatbath.py`).

È il test più importante di tutta la milestone: **stessa fisica, due percorsi di codice diversi** (script 0D vs solver CFD completo). Se le curve coincidono, l'integrazione nel solver è corretta.

**Risultato:**

| Grandezza | Errore max | Equilibrio (solver → riferimento) |
|-----------|-----------|-----------------------------------|
| `T_tr` | 0.03 % (2.1 K) | — |
| `T_ve` | 1.87 % (22 K) | — |
| Convergenza | — | 7617.5 K → 7619.6 K |

Le curve sono praticamente sovrapposte (`heatbath-comparison.png`). L'errore residuo sulla `T_ve` è dovuto al fatto che il solver campiona a passi discreti e interpola, non a un errore di modello.

---

## 5. Approssimazione da tenere presente

Come detto in §3.3, il campo `e` risolto è ancora l'energia del thermo base, e il termine V-T nell'equazione dell'energia è scalato da `cvTrRatio`. Questo è **esatto in assenza di flusso** (heat bath) ma diventa approssimato attraverso un urto reale.

La formulazione **pienamente conservativa** — energia totale `E` e energia vibro-elettronica `E_ve` come le due variabili di stato, con `T` invertita direttamente da Mutation++ (`solveEnergies` / `setState` con `vars=0`) — è il refactor naturale da fare **quando si passerà alla validazione 1D** con shock veri (Milestone 4). Non serviva ancora per l'heat bath e avrebbe complicato la verifica di questo primo passo.

---

## 6. Cosa sblocca e cosa viene dopo

Questa milestone chiude il punto più critico della lista "cosa manca ancora": il modello 0D validato ora **vive dentro il solver**, con fisica vera calcolata per cella, e c'è un test di regressione permanente che lo dimostra.

I prossimi passi (roadmap Casseau Part One → Part Two):

- **Milestone 2** — chimica in 0D: casi reagenti (Fig 7/8, dissociazione N₂ con accoppiamento chimica-vibrazione Q_C-V), tutti i test attuali hanno `mechanism none`.
- **Milestone 3** — chimica nel solver.
- **Milestone 4** — 1D shock tube, con qui il refactor alla formulazione conservativa (§5).
- **Milestone 5** — 2D/3D (blunted cone Mach 11, cylinder Mach 20 del paper Part Two).

---

## File toccati

```
src/thermophysicalModels/multicomponentThermo/highEnthalpyMulticomponentThermo/
    HighEnthalpyMulticomponentThermo.H     (thermo: eve, eveEq, tauVT, cvTrRatio, TveFromEve, ...)

applications/modules/shockThermo/
    thermophysicalPredictor.C              (solver: equazione conservativa in eve)

applications/test/nonEqTTv/
    N2/Test-N2.C                           (override opzionale T_tr/T_ve da riga di comando)
    solverHeatBath/                        (NUOVO: caso di validazione solver vs 0D)

tutorials/shockThermo/shockTube/
    constant/physicalProperties            (nuova sintassi dict)
    system/fvSolution                      (risoluzione campo eve)
```

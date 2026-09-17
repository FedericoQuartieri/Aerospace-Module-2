# Milestone 1 — Rilassamento V-T dal thermo al solver

**Data:** 3 luglio 2026 (commit `e9626ce`), **aggiornato il 17 settembre 2026** dopo la
riscrittura del ponte con Mutation++ (vedi `milestone-2-heat-bath-paper.md`).
**Obiettivo:** portare la fisica del rilassamento vibro-traslazionale (V-T) validata in
0D dentro il solver CFD `shockThermo`, in modo che funzioni su tutte le celle e non
solo negli script standalone.

> Nota: la prima versione di questa milestone teneva nel thermo i campi `eveEq_`,
> `tauVT_` e `cvTrRatio_` e calcolava il tempo di rilassamento a mano. Quei campi non
> esistono più: il ponte è stato riscritto e questo documento descrive la versione
> attuale del codice.

---

## 1. Da dove partivamo

Prima di questa milestone il solver risolveva sì un'equazione per la temperatura
vibro-elettronica `Tve`, ma con **fisica finta**: il tempo di rilassamento e il
rapporto fra i calori specifici erano costanti scelte a mano nel dizionario. La
fisica vera viveva **solo** nei test 0D standalone (`applications/test/nonEqTTv/`).

In più l'equazione era scritta in **temperatura**, che non è una quantità
conservata: attraverso un urto avrebbe dato risultati sbagliati. Andava riscritta in
**energia**.

---

## 2. Com'è fatto adesso

### 2.1 Il thermo — `HighEnthalpyMulticomponentThermo.H`

È un **ponte sottile su Mutation++**: in ogni cella lo stato del gas è dato dalle
variabili conservate del paper (eq. 23), cioè densità parziali, energia totale ed
energia vibro-elettronica; Mutation++ ricava da queste le due temperature.

**Campi** (in `highEnthalpyMulticomponentThermo::composite`, registrati nel mesh):

| Campo | Significato | Unità |
|-------|-------------|-------|
| `Tve_` | Temperatura vibro-elettronica, letta da `0/Tve` | K |
| `eve_` | Energia vibro-elettronica per unità di massa (la variabile conservata) | J/kg |

Il campo `e` (energia sensibile, `he()` del thermo base) è la seconda variabile
conservata; il thermo base (janaf) serve solo a OpenFOAM per costruire i campi, la
sua polinomiale non entra nella fisica.

**Costruttore**: rilegge `T` dal file del caso (il thermo base l'aveva ricalcolata
dall'energia limitandola a 20000 K, il suo intervallo di validità), poi per ogni
cella chiede a Mutation++ le energie a `(p, T, Tve)` e riempie `e` ed `eve`.

**Sorgenti** (usate dal solver, calcolate cella per cella con le funzioni di
`mutationSources.H`, le stesse dei programmi 0D):

- `computeSourceVe()`: scambio V-T (eq. 8) più, con la chimica, l'energia
  vibro-elettronica portata via dalle reazioni (eq. 30);
- `computeSourceY(i)`: produzione chimica della specie `i` (eq. 27);
- `computeSourceE()`: calore di reazione per l'energia sensibile, `-Σ hf_s ω_s`.

**`correct()`**: per ogni cella passa a Mutation++ `ρ(e + Σ Y_s hf_s)` e `ρ eve`
(`setState` con `vars = 0`), che inverte le energie con un Newton e restituisce
`T` e `Tve`; poi `updatePsi()` ricalcola `psi = ρ/p` con `p` somma delle pressioni
parziali (eq. 24), in celle e facce di bordo.

### 2.2 Le sorgenti — `mutationSources.H`

Funzioni di solo Mutation++ (niente OpenFOAM), condivise fra thermo e programmi 0D:

- `makeVibrators(mix)`: per ogni molecola, il modello di Millikan-White con
  correzione di Park di Mutation++ (eq. 9-17, costanti del suo `VT.xml`);
- `speciesEv`, `speciesEve`: energia vibrazionale e vibro-elettronica di una specie a
  una data `Tv` (eq. 5 e 7);
- `sourceVT`: Landau-Teller (eq. 8) con forza motrice `e_ve(T) − e_ve(Tv)` della
  molecola, come nel paper. L'`OmegaVT` di Mutation++ userebbe solo l'energia
  vibrazionale: senza livelli elettronici i due coincidono;
- `productionRates`: velocità di reazione alla temperatura di Park `T^a Tv^(1−a)`
  (eq. 29). Mutation++ le valuta a `sqrt(T Tv)` e non si può cambiare: gli si passa
  lo stato `{T_P, T_P}` e poi si ripristina `{T, Tv}`. Vale per i meccanismi di sola
  dissociazione;
- `sourceCV`: `Σ e_ve,s(Tv) ω_s`, il modello non preferenziale (eq. 30-31), lo stesso
  che Mutation++ implementa in `OmegaCV + OmegaCElec`.

### 2.3 Il solver — `thermophysicalPredictor.C`

```cpp
// specie: sorgente chimica di Mutation++
YiEqn: ddt(rho, Yi) + div(phi, Yi) + divj(Yi) == wdot_i + fvModels
// energia sensibile: calore di reazione
EEqn:  ddt(rho, e) + div(phiEp) + ddt(rho, K) == Q_chem + fvModels
// energia vibro-elettronica: V-T + chimica-vibrazione
EveEqn: ddt(rho, eve) == Q_ve
thermo_.correct();   // decode: T, Tve, psi da Mutation++
```

Le sorgenti sono esplicite. Nell'heat bath con passo di 1 ns è la stessa
integrazione dei programmi 0D, ed è per questo che i due percorsi coincidono.

---

## 3. La fisica in tre righe

### 3.1 Modello a due temperature

Il gas ad alta entalpia è in **non-equilibrio termico**: `Ttr` (traslazione e
rotazione, che si equilibrano subito) e `Tve` (vibrazione e livelli elettronici, che
si equilibrano lentamente). Dietro un urto la traslazionale schizza in alto,
la vibrazionale insegue con un ritardo governato dal tempo di rilassamento.

### 3.2 Equazione di Landau-Teller (eq. 8)

$$\rho_m \frac{\partial e_{ve,m}}{\partial t} = \rho_m\,\frac{e_{ve,m}(T_{tr}) - e_{ve,m}(T_{ve})}{\tau_{m,V\text{-}T}}$$

Quando `Tve → Ttr` la differenza si annulla: è l'equilibrio.

### 3.3 Perché tutto in energia

`e` ed `eve` sono conservate e la loro decodifica in temperature la fa Mutation++
con le stesse funzioni termodinamiche usate per costruirle: non ci sono rapporti di
calori specifici da approssimare, e il solver e i programmi 0D condividono le stesse
formule. Questo è esatto anche attraverso un urto (a meno del trasporto di `eve`,
che nell'heat bath non c'è e nel solver è ancora da aggiungere).

---

## 4. Validazione

Caso `applications/test/nonEqTTv/solverHeatBath/`, **heat bath a singola cella**
(1×1×1 celle, velocità nulla), lanciato con `./Allrun <figura>` per ognuna delle
figure del paper riprodotte nel solver (3a, 3b, 4, 5, 7). Il confronto con il
programma 0D e con le curve del paper è in `milestone-2-heat-bath-paper.md`: solver e
programma 0D coincidono entro 0.03-0.2 % su tutte le curve, chimica compresa.

---

## 5. Approssimazioni da tenere presenti

- Il solver risolve `eve` senza il termine di trasporto `div(phi, eve)`: va aggiunto
  prima di usare il modello in un caso con flusso (tubo d'urto, milestone 4).
- Il thermo base (janaf) è valido fino a 20000 K: le sue funzioni (`mu`, `Cp`, ...)
  vengono limitate sopra quella temperatura e OpenFOAM lo segnala con un avviso; le
  temperature e `psi` non ne risentono perché le calcola Mutation++.
- Le proprietà di trasporto restano quelle del thermo base: i casi validati sono
  inviscidi.

---

## 6. Cosa viene dopo

- **Milestone 3** — trasporto di `eve` e chimica nel tubo d'urto 1D.
- **Milestone 4** — 2D/3D (blunted cone Mach 11, cylinder Mach 20 del paper Part Two).

# Milestone 3 — Chimica nel solver + formulazione conservativa

**Data:** 4 luglio 2026
**Branch:** `milestone-3-chemistry-solver`
**Obiettivo:** portare la chimica two-temperature dentro il solver CFD e, insieme (opzione "B" concordata), sostituire l'approssimazione energetica della M1 con la **formulazione conservativa piena**: temperature decodificate dalle energie conservate via Mutation++, niente più `cvTrRatio`.

---

## 1. La scelta "B": perché il refactor conservativo subito

La M1 aveva un'approssimazione dichiarata: il campo `e` risolto era l'energia sensibile del thermo base (polinomi rrho), la temperatura veniva invertita dai polinomi, e il sink V-T nell'equazione dell'energia andava scalato di `cvTrRatio = Cv_base/cv_tr` per ottenere il dT/dt giusto. Funzionava per l'heat bath, ma:
- il riferimento M2 (con chimica) è **conservativo esatto** — validare un'approssimazione contro un riferimento esatto rende ambiguo ogni scostamento (bug o modello?);
- la M4 (shock veri) avrebbe comunque richiesto il refactor → farlo in M3 evita di scrivere due volte l'aggancio all'energia.

## 2. La nuova formulazione

### 2.1 Datum dell'energia

Il campo risolto `e` ora vive sul **datum Mutation++**:

$$e = e_{tr}(T) + e_{ve}(T_{ve}) \qquad \text{(entalpie di formazione escluse)}$$

Scomposizione esatta: `e_abs = e_tr + e_ve + hf/Mw` per ogni specie (verificato: `speciesHOverRT` espone `ht`, `hr`, `hv`, `hel`, `hf` separatamente). All'avvio `initialiseEnergies()` sovrascrive `e` con questo datum a partire da (T, Tve, Y) letti dal caso.

### 2.2 Decode delle temperature (il cuore del refactor)

A ogni `correct()`, per ogni cella e faccia di bordo:
1. **Tve da eve** — Newton su `e_ve` (come M1);
2. **T da (e − eve)** — Newton su `e_tr(T)`, che per RRHO è **lineare in T** → converge in 1-2 iterazioni.

I polinomi del thermo base **non entrano più nella fisica**: servono solo da bootstrap alla costruzione. La `correct()` non chiama più la base quando il bridge è attivo.

### 2.3 Conseguenza elegante: niente più termine V-T nell'equazione dell'energia

Il travaso V-T sposta energia tra i pool `e_tr` e `e_ve` **dentro** `e`: l'energia totale non cambia. L'effetto su T passa automaticamente dal decode `T(e − eve)`. Spariti `eRelaxSource`, `cvTrRatio` e lo switch `coupleEnergy`.

## 3. La chimica nel solver

Tutte le sorgenti vengono da Mutation++ (coerente con M2 — la kinetics OpenFOAM non conosce la temperatura di controllo di Park), calcolate in `sourcesAtCurrentState()` e pubblicate come campi:

| Campo | Equazione | Formula |
|---|---|---|
| `mutR_<specie>` [kg/m³/s] | `YiEqn` | ω̇ᵢ da `netProductionRates` |
| `mutQdot` [W/m³] | `EEqn` | −Σ hfᵢ·ω̇ᵢ (datum sensibile) |
| `mutQcv` [W/m³] | `EveEqn` | Σ (e_v+e_el)ᵢ(Tve)·ω̇ᵢ (Candler, replica OmegaCV+OmegaCElec) |
| `tauVT` [s] | `EveEqn` (LT semi-implicito) | da **Q_VT puro** = `energyTransferSource` − Q_cv |

L'ultima riga corregge un bug latente della M1: con la chimica accesa `energyTransferSource` restituisce VT+CV insieme, e attribuire il CV al rilassamento verso `eveEq` sarebbe sbagliato (il CV non scala con `eveEq − eve`). Lo split è possibile perché la parte CV+CElec è esattamente la formula di Candler, replicata con gli stessi identici call Mutation++ (validato in M2: percorsi bit-identici).

Nel solver: `YiEqn` usa `mutR_<specie>` quando esiste (altrimenti fallback a `reaction->R`), `EEqn` prende `mutQdot`, `EveEqn` prende `mutQcv`. Fallback per cella: sorgenti a zero dove il bridge fallisce.

## 4. I due bug trovati (e le loro lezioni)

### 4.1 SIGFPE in OmegaVT con frazioni molari esattamente zero

`tau_park = 1/(n·X·σ)` viene calcolato per **ogni** vibratore della miscela air_5 (N₂, O₂, NO): con X(O₂)=0 esatto → divisione per zero → +Inf. Nei test standalone l'Inf è **benigno** (il contributo del vibratore assente si annulla), ma `foamRun` ha il trapping FPE attivo → crash al primo step.

La bisezione ha mostrato che **anche il codice M1 crasha oggi**: il caso M1 non era mai stato rilanciato con `dynamicCode` pulito (l'Allclean non lo cancellava — ora lo fa). Bug pre-esistente, mai colpito per caso.

**Fix alla radice:** floor 1e-30 sulle frazioni massiche nel bridge (`mutationMassFractions`), su *tutte* le specie Mutation. Fisicamente indistinguibile da zero, mantiene il trapping attivo (che nei casi multi-D serve a beccare bug veri).

### 4.2 Clamp janaf a Thigh=20000

Il caso reagente parte da 30000 K, ma l'entry `N` del `speciesThermo.janaf` ha `Thigh 20000`: l'inversione T(e) della **costruzione base** clampa T a 20000 prima ancora che il bridge parta → tutto lo stato iniziale sbagliato (ρ compresa, via ψ). Il caso M1 a 10000 K non lo vedeva.

**Fix:** `Thigh 50000` nel janaf del caso (uniforme con `maxTemperature` del bridge). Innocuo: col decode conservativo i polinomi base sono solo bootstrap. *Lezione: leggere subito il testo dei FOAM Warning — c'erano, ripetuti, e dicevano esattamente questo.*

## 5. Validazione

### 5.1 Regressione M1 (non-reagente, N₂ 10000/1000)

| Grandezza | M1 | M3 |
|---|---|---|
| T_tr max err | 0.03% | **0.02%** |
| T_ve max err | 1.87% | 1.87% (residuo da campionamento probe + semi-implicito vs Euler esplicito) |

### 5.2 Caso nuovo: `solverHeatBathReacting` (N₂/N 30000/1000, meccanismo Park)

Confronto con il CSV M2 `results-N2N-reacting-30000-1000-mpp-park05.csv` — stessa fisica, due percorsi di codice (0D esplicito vs solver CFD completo):

| Grandezza | Errore max | Finale a 1e-4 s (solver / rif.) |
|---|---|---|
| T_tr | **5.78 K (0.03%)** | 10228.7 / 10225.0 K |
| T_ve | 56 K (4.17%, solo nel transitorio ripido iniziale) | 10253.1 / 10249.3 K |
| n_N2/n₀ | **1e-4** | 0.1783 / 0.1782 |
| n_N/n₀ | **2e-4** | 1.1434 / 1.1435 |

Curve sovrapposte punto per punto, overshoot di T_ve incluso (`heatbath-reacting-comparison.png`). Run: 1e5 step, ~4.5 min.

### 5.3 Smoke test shockTube

Il tutorial gira senza errori né warning nuovi con la formulazione conservativa (frozen, T basse: qui si verifica solo che la gasdinamica di base non si sia rotta).

## 6. File toccati

```
src/.../HighEnthalpyMulticomponentThermo.H   (decode conservativo, sorgenti chimiche,
                                              floor Y, via cvTrRatio)
applications/modules/shockThermo/
    thermophysicalPredictor.C                (YiEqn+mutR, EEqn+mutQdot, EveEqn+mutQcv,
                                              via eRelaxSource/coupleEnergy)
applications/test/nonEqTTv/
    solverHeatBathReacting/                  (NUOVO caso di validazione reagente)
    solverHeatBath/Allclean                  (pulisce anche dynamicCode)
    solverHeatBath/constant/physicalProperties (via coupleEnergy)
tutorials/shockThermo/shockTube/constant/physicalProperties (via coupleEnergy)
.gitignore                                   (pattern per il nuovo caso)
```

## 7. Cosa resta / prossimi passi

- La M4 (1D shock tube con validazione fisica) eredita la formulazione già conservativa: resta "solo" il caso, la validazione contro Sod/rilassamento post-shock, e il problema prestazioni (Mutation++ per cella per step — il reagente 1-cella fa ~370 step/s).
- BC a parete con T fissa (per M5) valuteranno `he(p,T)` col datum base ≠ datum Mutation: da sistemare quando serviranno.
- CFL/chimica: le sorgenti sono esplicite in tempo (valide a dt=1e-9); per M4/M5 servirà valutare un limitatore o l'accoppiamento semi-implicito delle Y.

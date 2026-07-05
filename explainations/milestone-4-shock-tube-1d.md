# Milestone 4 — Shock tube 1D: la fisica validata in un flusso vero

**Data:** 5 luglio 2026
**Branch:** `milestone-4-shock-tube`
**Obiettivo:** portare la validazione dal 0D (heat bath) all'**1D**: verificare che gasdinamica, rilassamento V-T e chimica funzionino insieme attraverso onde d'urto reali. La formulazione conservativa era già pronta dalla M3 — qui si costruiscono i casi e i riferimenti.

---

## 1. Scoperta preliminare: il tutorial non era un tubo d'urto

Il tutorial `tutorials/shockThermo/shockTube` era in realtà **una singola cella** ("homogeneous heat-bath style test") — utile all'epoca per sviluppare il thermo, ma fuorviante. In questa milestone è stato ripristinato a vero tubo 1D (200 celle, condizioni Sod), mentre la validazione quantitativa vive in `applications/test/shockTube1D/`.

## 2. I tre pezzi della milestone

### 2.1 `postShockRelax` — il riferimento ODE stazionario

Nuovo test standalone Mutation++ che costruisce il riferimento per la zona di rilassamento dietro un urto normale:

1. **Salto di Rankine-Hugoniot frozen** (eve e composizione congelate attraverso il fronte sottile), calcolato con le stesse funzioni termodinamiche del bridge del solver;
2. **Integrazione delle equazioni 1D stazionarie** in shock frame: `m=ρu`, `P=p+mu`, `H₀=h+u²/2` costanti; `d(eve)/dx = Q_ve/m`, `dY_i/dx = ω̇_i/m` con sorgenti Mutation++ (VT+CV complete); a ogni x lo stato algebrico è recuperato con una bisezione sul ramo subsonico.

**Auto-validazione del riferimento**: per Ms=8.355 il salto frozen riproduce la teoria RH a gas ideale (γ=1.4) **a 4 cifre**: p₂/p₁=81.28, ρ₂/ρ₁=5.599, T₂/T₁=14.52. A valle: equilibrio termico T=Tve ✓. Runtime: 0.2 s.

### 2.2 Caso A — `sod/`: la gasdinamica contro la soluzione esatta

Sod classico in N₂ puro freddo (γ=1.4, vibrazione congelata): p 10⁵/10⁴ Pa, ρ ratio 8, dominio ±5 m, 1000 celle, bridge Mutation++ **attivo** (mechanism none). Confronto a t=5 ms con il **risolutore di Riemann esatto** (Toro, scritto in `exact-riemann.py`):

| Metrica | Valore |
|---|---|
| L1(ρ)/ρ_L | **0.12%** |
| L1(u)/u* | **0.17%** |
| L1(p)/p_L | **0.07%** |
| Posizione shock | 2.815 vs 2.817 m (errore 2 mm su 10 m) |

Ventaglio di espansione, contatto e shock tutti al posto giusto, discontinuità su 2-3 celle (`sod-comparison.png`). **Il bridge conservativo non rompe la gasdinamica.**

### 2.3 Caso B — `relaxing/`: la zona di rilassamento contro l'ODE

Urto forte in N₂: driver 5 MPa/3000 K → driven 1 kPa/300 K; dominio 1 m, 2000 celle (dx=0.5 mm), meccanismo `N2_diss_park` attivo. Lo shock **misurato** dal CFD viaggia a **2950 m/s (Ms=8.36)**.

Pipeline di validazione (tutta in `Allrun`):
1. misura della velocità dello shock dagli ultimi due snapshot;
2. `Test-postShockRelax` genera il riferimento ODE **alla velocità misurata**;
3. il profilo CFD dietro il fronte viene mappato in shock frame e sovrapposto.

Risultato sulla zona di rilassamento (2–34.5 mm dietro il fronte, fino all'equilibrio termico del riferimento):

| Grandezza | Errore max | Errore medio |
|---|---|---|
| T | **1.04%** | 0.69% |
| Tve | 6.27% (nel gradiente iniziale ripido) | **0.83%** |

Il quadro fisico è da manuale (`relaxation-comparison.png`): dietro il fronte T salta a ~4355 K (frozen), Tve resta ~300 K, poi la zona di Landau-Teller li porta all'equilibrio ~3740 K in ~35 mm — CFD e ODE sovrapposti. Oltre ~45 mm c'è la superficie di contatto (fisica di un'altra onda, correttamente esclusa dalle metriche).

## 3. Problemi trovati e sistemati

- **Igiene dei casi**: `setFields` modifica `0/` sul posto → i casi ora usano il pattern `0.bak/` (sorgente versionato) + copia in `Allrun`, come il tutorial. Il primo run del caso B è fallito proprio per questo (campi da 1000 celle copiati in una mesh da 2000).
- **`OFstream` non crea directory**: l'output del test ODE spariva silenziosamente se `output/` non esisteva → `mkdir -p` nell'Allrun (e `|| exit 1` sul `cp`).
- **Locale it_IT**: `printf '%.0f' 2950.0` fallisce con la locale italiana → lo script `--speed` ora stampa direttamente l'intero, coerente con il `lround` del C++.

## 4. Prestazioni (misurate)

| Caso | Celle × step | Tempo solve | µs/cella/step |
|---|---|---|---|
| Sod (mechanism none) | 1000 × 2500 | ~85 s | **34** |
| Relaxing (kinetics attiva) | 2000 × 4800 | ~410 s | **42** |

Per l'1D è comodo. Proiezione M5 (cilindro Mach 20 del paper: ~155k celle): ~6.5 s/step → migliaia di step = **ore/giorni**. Candidati per l'ottimizzazione (da fare in M5): salto delle celle in equilibrio/fredde, riuso dello stato Mutation fra celle simili, parallelizzazione (il bridge è per-cella → domain decomposition standard di OpenFOAM dovrebbe scalare).

## 5. File toccati

```
applications/test/shockTube1D/          (NUOVO gruppo di validazione M4)
    postShockRelax/                     (riferimento ODE: Make + Test-postShockRelax.C)
    sod/                                (caso A + exact-riemann.py + 0.bak pattern)
    relaxing/                           (caso B + compare-relaxation.py)
tutorials/shockThermo/shockTube/        (ripristinato a vero tubo 1D, 200 celle Sod)
    system/{blockMeshDict,setFieldsDict,controlDict}
.gitignore                              (pattern shockTube1D)
```

## 6. Prossimi passi (M5)

- Blunted cone Mach 11.3 (assialsimmetrico, non-reagente) vs dati Michigan/MONACO/CUBRC del paper Part Two — richiede: mesh, BC di parete (attenzione: `fixedEnergy` valuta `he(p,T)` col datum del thermo base ≠ datum Mutation — da sistemare **prima**), e il lavoro sulle prestazioni.
- Cilindro Mach 20 reagente vs dsmcFoam.
- Rimandati noti: Fig 9 (aria 5 specie), CVDV-QK, Tve multiple per specie.

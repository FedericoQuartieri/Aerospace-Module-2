# Milestone 2 — Chimica in 0D: heat bath reagente N₂–N (Fig 7/8)

**Data:** 4 luglio 2026
**Branch:** `milestone-2-reacting-0d`
**Obiettivo:** estendere la validazione 0D ai casi **reagenti** del paper (Casseau et al., Part One, §3.4, Figure 7 e 8): dissociazione dell'azoto accoppiata al rilassamento vibrazionale. Fino a qui tutti i test avevano `mechanism("none")` — la chimica era spenta.

---

## 1. Il caso fisico

Heat bath adiabatico a volume costante, miscela N₂–N in ugual densità numerica (5.0×10²² m⁻³ ciascuna), con **una sola reazione irreversibile**:

$$N_2 + N_2 \xrightarrow{k_f} 2N + N_2$$

con rate di Park 1993 (Tabella 2 del paper): A = 7.0×10²¹ cm³/(mol·s), β = −1.6, Ta = 113200 K.

Due condizioni iniziali:
- **Fig 7:** T_tr = 30000 K, T_ve = 1000 K (non-equilibrio termico forte)
- **Fig 8:** T_tr = T_ve = 30000 K (equilibrio termico iniziale)

Tre fenomeni si contendono l'energia:
1. **V-T** (Landau-Teller, Millikan-White+Park): travasa energia tra pool traslazionale e vibrazionale — già validato nelle milestone precedenti;
2. **Dissociazione**: consuma N₂, assorbe l'energia di legame (Ta) dal pool traslazionale;
3. **Accoppiamento chimica-vibrazione (Q_C-V)**: ogni molecola dissociata porta via anche la propria energia vibrazionale media — modello **non-preferenziale di Candler** (eq. 30-31 del paper): Q_CV = Σᵢ e_v,i·ω̇ᵢ.

## 2. Cosa sa fare Mutation++ (verificato nei sorgenti, non a memoria)

| Componente | Dove | Esito |
|---|---|---|
| Reazioni irreversibili `=>` | `kinetics/Reaction.cpp` | ✅ supportate |
| Q_CV non-preferenziale Candler | `transfer/OmegaCV.cpp` | ✅ identico all'eq. 31 del paper |
| `energyTransferSource()` | `ChemNonEqTTvStateModel.cpp` | ✅ = Ω_VT + Ω_CV (+Ω_CElec) |
| `netProductionRates()` | `kinetics/Kinetics.h` | ✅ ω̇ in kg/(m³·s) |
| decodifica T da energie | `setState(ρᵢ, {ρe, ρe_ve}, 0)` | ✅ `solveEnergies` (Newton) |
| Temperatura di controllo dissociazioni | `kinetics/RateManager.cpp:58` | ⚠️ `sqrt(T·Tv)` **hard-coded** (Park q=0.5), il paper usa q=0.7 |
| Classificazione `N2+N2=>2N+N2` | `Reaction.cpp:375` | ✅ `DISSOCIATION_M` → forward rate alla T di Park |

Il ⚠️ è la ragione del **doppio percorso** (§4).

## 3. La formulazione: conservativa in energia

A differenza dei test precedenti (update di temperatura via `dT/dt = Q/ρcv`), qui il vettore di stato è quello **conservativo del paper** (eq. 22-26):

```
U = (ρ_s, ρe_ve, ρe_tot)
```

Per un heat bath adiabatico a volume costante:
- `ρe_tot` = **costante** (include le entalpie di formazione → l'energia chimica assorbita dalla dissociazione è contabilizzata automaticamente, senza termini sorgente scritti a mano);
- `dρ_s/dt = ω̇_s`;
- `dρe_ve/dt = Q_ve = Q_VT + Q_CV`.

A ogni passo le temperature vengono **decodificate** dalle energie con `mix.setState(ρ_s, {ρe_tot, ρe_ve}, vars=0)` — lo stesso percorso (`solveEnergies`) che il solver CFD userà nella milestone 3: questo test lo pre-valida.

## 4. Doppio percorso di integrazione

Stesso integratore, due fornitori dei termini sorgente:

| Modo | ω̇ | Q_ve | Esponente Park |
|---|---|---|---|
| `mpp` | `netProductionRates()` (Mutation++) | `energyTransferSource()` = VT+CV | **0.5** (hard-coded) |
| `manual` | legge di azione di massa scritta nel test | VT da Mutation++ (mechanism `none` → CV=0) + Candler manuale | **parametrico** (default 0.7, come il paper) |

Il modo `mpp` è ciò che il solver userà in M3; il modo `manual` a q=0.7 riproduce fedelmente il paper. **Girando `manual` con q=0.5 le curve devono coincidere con `mpp`**: è la cross-validazione dei due percorsi.

## 5. Dati: root non-elettronico versionato

Il confronto del paper (§3.4) è contro dsmcFoam, che **non ha il modo elettronico** → i run vanno fatti senza energia elettronica. Finora la selezione elettronico/non-elettronico si faceva **editando a mano** `species.xml` (non riproducibile: i commenti nel file lo testimoniano). Ora esiste un secondo data root versionato:

```
mutation-data-noel/        # selezionato con MPP_DATA_DIRECTORY=mutation-data-noel
    thermo/species.xml     # tutte le specie a solo ground state; θv(N2)=3371 K (Tab. A1)
    mechanisms/, mixtures/, transfer/, transport/   # copie degli altri file
```

Il meccanismo nuovo è `mutation-data/mechanisms/N2_diss_park.xml` (presente in entrambi i root). Coerenza verificata: `air5_Park.xml` ha N2+M con A=3.0×10²² e M(N2)=0.2333 → 7.0×10²¹ ✓.

## 6. Risultati e validazione

Cinque run (`Allrun` aggiornato): Fig 7 in `mpp`/`manual 0.7`/`manual 0.5`, Fig 8 in `mpp`/`manual 0.7`.

### 6.1 Validazioni interne (le più forti)

| Check | Risultato |
|---|---|
| ω̇ Mutation++ vs legge di azione di massa a t=0 (unità!) | identici a tutte le cifre (−0.00149097 kg/m³s) |
| **Curve `mpp` vs `manual q=0.5`** (1200 punti, 6 decadi) | **differenza max = 0** su T_tr, T_ve, n/n₀ |
| Conservazione massa | ~10⁻¹⁴ (precisione macchina) |
| Errore decodifica energia a fine run | 0 ÷ 4×10⁻¹⁶ |

La seconda riga è la validazione incrociata: due implementazioni indipendenti della chimica (kinetics di Mutation++ vs legge di azione di massa + Candler scritti nel test) producono traiettorie bit-identiche.

### 6.2 Confronto col paper (config Park, Fig 7)

| Grandezza a t=10⁻³ s | Questo test | Paper (Fig 7) |
|---|---|---|
| T_tr ≈ T_ve | 8298 / 8301 K | ~8000 K |
| n_N2/n₀ | 0.137 | ~0.13 |
| n_N/n₀ | 1.225 | ~1.22 |

Comportamento qualitativo riprodotto: overshoot di T_ve a ~20 kK attorno a 10⁻⁶ s (paper: ~18-19 kK), poi raffreddamento congiunto. Per la Fig 8, riprodotto il **lag di T_ve** dietro T_tr citato esplicitamente nel paper ("the Park configuration first initiates a thermal relaxation... which results in a lag in the vibrational temperature decrease"); stato finale a 10⁻³ s: 10108 K, n_N2/n₀ = 0.047.

### 6.3 La scoperta utile: q=0.5 vs q=0.7 quasi indistinguibili

Il timore iniziale (T_controllo 5477 K vs 10820 K a t=0 → dinamiche diverse) si è rivelato **quasi irrilevante**:

| t [s] | T_tr(q=0.5) | T_tr(q=0.7) |
|---|---|---|
| 10⁻⁶ | 20620 | 19979 |
| 10⁻⁵ | 13464 | 13504 |
| 10⁻³ | 8298 | 8299 |

Differenza massima ~640 K (~3%) nel transitorio, poi convergenza. Il motivo fisico: il rilassamento V-T è più veloce della chimica, quindi quando la dissociazione diventa importante T_ve è già salita e le due temperature di controllo si avvicinano. **Conseguenza pratica per la M3:** si può usare Mutation++ end-to-end nel solver (q=0.5 hard-coded) senza patch al thirdParty, con errore trascurabile rispetto al paper.

## 7. File toccati

```
applications/test/nonEqTTv/
    mutation-data/mechanisms/N2_diss_park.xml      (NUOVO: reazione irreversibile Park)
    mutation-data-noel/                            (NUOVO: data root non-elettronico)
    N2N-reacting/                                  (NUOVO: Make + Test-N2N-reacting.C)
    plot-temps-reacting.py                         (NUOVO)
    Allwmake, Allrun                               (aggiornati)
    output/results-N2N-reacting-*.csv              (5 CSV canonici)
    output/{temp,density}-curves-N2N-reacting-*.png (4 plot)
```

## 8. Cosa resta della milestone 2 (rimandato di comune accordo)

- **Fig 9** (aria 5 specie, 19 reazioni QK, equilibrio termico): soprattutto data-entry dei rate dal paper di Scanlon; le fondamenta (formulazione conservativa, doppio percorso) sono già pronte.
- **Config CVDV-QK** (Fig 7/8, linee tratteggiate del paper): il modello CVDV di Marrone-Treanor non è in Mutation++; implementabile a mano nel percorso `manual` (eq. 33-40, oscillatore armonico) se servirà.

## 9. Prossimo passo naturale (M3)

Portare la chimica nel solver: `reaction->correct()` OpenFOAM → sorgenti da Mutation++ (`netProductionRates` per le Y_i, `energyTransferSource` completo per l'equazione E_ve), e heat bath reagente su singola cella nel solver contro questi CSV — lo stesso pattern solverHeatBath della M1.

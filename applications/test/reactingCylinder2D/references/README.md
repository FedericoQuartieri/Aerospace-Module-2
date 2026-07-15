# Riferimenti digitalizzati — Fig. 5 del paper Part Two

CSV estratti dalla **Figura 5** (cilindro Mach 20 reagente, §3.2) di:

> Casseau, V.; Espinoza, D.E.R.; Scanlon, T.J.; Brown, R.E.
> *A Two-Temperature Open-Source CFD Model for Hypersonic Reacting Flows,
> Part Two: Multi-Dimensional Analysis.* Aerospace **2016**, 3, 45.
> doi:10.3390/aerospace3040045 — licenza **CC-BY 4.0**. Il PDF non è
> duplicato: `digitize-fig5.py` legge la copia in
> `../../bluntedCone2D/references/`.

Stessa tecnica di estrazione vettoriale della Fig 2 (vedi il README del
cono per metodo e trappole generali). Novità di questa figura:

- **asse y logaritmico** nel pannello (c): calibrazione in log10 e
  de-log in uscita;
- **serie per colore** (caption): run 1 = nero, run 2 = rosso,
  **run 3 = blu (Park TTv + rates Park = il NOSTRO setup)**;
- **filtro tick con orientamento**: qui molte curve corrono appiattite
  sull'asse x (Mach post-shock ~0, Cp/Cf/q → 0 a θ alti); un tick è
  corto E perpendicolare al suo bordo, un tratto di curva lungo il bordo
  è parallelo e sopravvive (il filtro della Fig 2 lo mangiava).

## Avvertenze note

- Nei pannelli (d), (e), (f) il file `run1` contiene **anche** la curva
  non-reagente solida (stesso colore nero nel PDF): in (d) NR +
  "all reacting runs" (coincidono quasi ovunque), in (e) NR + "runs 1
  and 2", in (f) NR (quella alta, ~12 W/cm² al ristagno) + run 1
  dash-dot. La serie **run3 è pulita** ed è quella da usare per il
  confronto quantitativo.
- In (b) ogni run contiene sia T_tr che T_v (due curve per file): si
  separano per valore di y, come per la Fig 2a.
- `dsmc` = croci dsmcFoam; `qk` = cerchi QK (centri dei marker).
- Le bande grigie Kn_GLL del paper delimitano dove il confronto
  CFD-DSMC perde senso (transizione): tenerne conto nelle metriche.

## File

| File | Pannello | Contenuto | x | y |
|---|---|---|---|---|
| `fig5a-*.csv` | (a) | Mach, linea di ristagno | m (centro cilindro) | — |
| `fig5b-*.csv` | (b) | T_tr e T_v | m | kK |
| `fig5c-*.csv` | (c) | densità numerica N2, N | m | 1/m³ |
| `fig5d-*.csv` | (d) | coefficiente di pressione | θ [deg] | — |
| `fig5e-*.csv` | (e) | coefficiente d'attrito | θ [deg] | — |
| `fig5f-*.csv` | (f) | flusso termico a parete | θ [deg] | W/cm² |

QA: `fig5-overlay.png` — **da ispezionare dopo ogni rigenerazione**.
Rigenerazione: `python3 digitize-fig5.py` (richiede `pip install pymupdf`).

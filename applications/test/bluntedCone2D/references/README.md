# Riferimenti digitalizzati — Fig. 2 del paper Part Two

CSV estratti dalla **Figura 2** di:

> Casseau, V.; Espinoza, D.E.R.; Scanlon, T.J.; Brown, R.E.
> *A Two-Temperature Open-Source CFD Model for Hypersonic Reacting Flows,
> Part Two: Multi-Dimensional Analysis.* Aerospace **2016**, 3, 45.
> doi:10.3390/aerospace3040045 — licenza **CC-BY 4.0** (il PDF è incluso
> qui per riproducibilità: `aerospace-03-00045.pdf`).

## Metodo: estrazione vettoriale, non digitalizzazione manuale

Le figure MDPI sono grafica **vettoriale**: le curve esistono nel PDF come
comandi di disegno con coordinate esatte. `digitize-fig2.py` le estrae con
PyMuPDF (`get_drawings`), le calibra dai riquadri degli assi (range noti
dalle etichette) e le classifica per serie:

- **colore blu** → hy2Foam con prima cella 10 µm (pannelli e, f)
- **path neri connessi e lunghi** (contiguità >80%) → hy2Foam (linea solida)
- **path neri a tratti staccati** → CFD Michigan (dash-dot)
- **simboli pieni** (fill senza stroke) → triangoli DSMC MONACO (centroide)
- **piccoli path a 2-3 segmenti** (gambo+cap) → barre d'errore degli
  esperimenti CUBRC run 31 (punto = centro barra)

L'accuratezza è limitata solo dalla calibrazione degli assi, non da
ricalco manuale. QA: `fig2-overlay.png` mostra i punti estratti sopra il
render della pagina — **da ispezionare dopo ogni rigenerazione**.

## Avvertenze note

- Nei pannelli (b), (e) la curva Michigan è sovrapposta a hy2Foam nel
  paper stesso ("shown to be superimposed"): non emerge come serie
  separata dove coincide al pixel.
- Gli `experiments` contengono solo i punti con barra d'errore
  distinguibile nel PDF (5 in (d), 9 in (f)); i restanti simboli
  sperimentali del paper sono graficamente indistinguibili dai triangoli
  DSMC.
- In (a) le due famiglie DSMC (Tt e Tv) e le due Michigan (Ttr e Tv) sono
  in un unico CSV per serie: si separano per valore di y se serve.

## File

| File | Pannello | Contenuto | Unità x | Unità y |
|---|---|---|---|---|
| `fig2a-*.csv` | (a) | T/T∞ linea di ristagno | mm (0=parete) | — |
| `fig2b-*.csv` | (b) | ρ/ρ∞ linea di ristagno | mm | — |
| `fig2c-*.csv` | (c) | U/U∞ linea di ristagno | mm | — |
| `fig2d-*.csv` | (d) | coefficiente di pressione | cm dal ristagno | — |
| `fig2e-*.csv` | (e) | coefficiente d'attrito | cm | — |
| `fig2f-*.csv` | (f) | numero di Stanton | cm | — |

Rigenerazione: `python3 digitize-fig2.py` (richiede `pip install pymupdf`).

# Curve di riferimento estratte dal paper

Curve di hy2Foam prese dalle figure del capitolo 3 di
Casseau et al., *A Two-Temperature Open-Source CFD Model for Hypersonic
Reacting Flows, Part One: Zero-Dimensional Analysis*, Aerospace 2016, 3, 34
(il pdf e' in `../aerospace-03-00034-1.pdf`).

Le figure del pdf sono vettoriali: `extract-paper-curves.py` converte la
pagina in svg con `mutool`, legge i tratti delle curve e li riporta in
coordinate fisiche usando i tick maggiori degli assi e le etichette numeriche.
I file csv sono il risultato: due colonne, tempo in secondi e valore
(temperatura in K, oppure densita' numerica normalizzata n/n0).

| file | figura | curva del paper |
|------|--------|-----------------|
| fig3a-Ttr, fig3a-Tv | 3a | hy2Foam default (riscaldamento N2) |
| fig3b-Ttr, fig3b-Tv | 3b | hy2Foam (raffreddamento N2) |
| fig4-noEl-*, fig4-el-* | 4 | hy2Foam senza / con energia elettronica |
| fig5-* | 5 | hy2Foam default (N2 + N) |
| fig5-lemans-* | 5 | hy2Foam con la convenzione di LeMANS (curva rossa) |
| fig6-noVV-*, fig6-VV-* | 6 | hy2Foam senza / con scambio V-V (N2 + O2) |
| fig7a-*, fig7b-* | 7 | hy2Foam Park (N2 + N reagente) |
| fig8a-*, fig8b-* | 8 | hy2Foam Park (N2 + N reagente, equilibrio termico iniziale) |
| fig9a-T, fig9b-* | 9 | hy2Foam (aria a 5 specie, 19 reazioni) |

Precisione: la calibrazione degli assi e' esatta sui tick; l'incertezza
e' lo spessore della linea, circa 40-60 K sulle temperature e circa 0.01
decadi sul tempo. Per rigenerare i csv: `python3 extract-paper-curves.py`.

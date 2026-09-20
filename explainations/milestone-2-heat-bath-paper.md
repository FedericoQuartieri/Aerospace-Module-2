# Milestone 2 — Heat bath del paper: programmi 0D, solver e confronto con le curve

**Data:** 17 settembre 2026
**Obiettivo:** riprodurre tutti i test 0D del capitolo 3 di Casseau et al., *A
Two-Temperature Open-Source CFD Model for Hypersonic Reacting Flows, Part One*,
Aerospace 2016, 3, 34 (`applications/test/nonEqTTv/aerospace-03-00034-1.pdf`),
su due strade: programmi 0D standalone con Mutation++ (tutte le figure) e solver
`shockThermo` con la libreria caricata a runtime (fig. 3a, 3b, 4, 5, 7), e
confrontare tutto con le curve estratte dal pdf.

---

## 1. Il risultato principale: solver e programma 0D fanno la stessa fisica

Ogni test che esiste su entrambe le strade è confrontato in automatico
(`compare.py`): stesse condizioni iniziali, stessi dati Mutation++, stesso passo di
1 ns. Il solver risolve le equazioni conservative in `e` ed `eve` con le sorgenti
di `mutationSources.H`, il programma 0D integra le stesse sorgenti a mano.

| Figura | caso | Ttr: scarto max solver vs 0D | Tv: scarto max | composizione |
|---|---|---|---|---|
| 3b | N2, raffreddamento | 0.1 K (0.00 %) | 0.3 K (0.00 %) | — |
| 3a | N2, riscaldamento | 2.6 K (0.03 %) | 6.8 K (0.21 %) | — |
| 4 senza E_el | N2 da 30000 K | 34 K (0.14 %) | 86 K (0.69 %) | — |
| 4 con E_el | N2 da 30000 K | 112 K (0.51 %) | 101 K (1.18 %) | — |
| 5 | N2 + N | 4.3 K (0.02 %) | 18 K (0.13 %) | — |
| 7 | N2 + N reagente, Park 0.7 | 5.8 K (0.03 %) | 26 K (0.18 %) | 0.03 % (0.0001 decadi) |

L'accordo resta allo 0.03-0.2 % anche con la chimica accesa (specie, calore di
reazione e accoppiamento chimica-vibrazione): è la verifica che l'integrazione del
modello nel solver è corretta. Gli scarti più grandi (fig. 4) sono nei primi 10 ns
del rilassamento a 30000 K, dove il solver fa due correttori PIMPLE per passo e il
programma 0D un solo Eulero esplicito.

---

## 2. Confronto con il paper, figura per figura

Le curve di hy2Foam sono state estratte dal pdf (`paper-data/`, vedi il README lì):
le figure sono vettoriali, la calibrazione è esatta sui tick degli assi e
l'incertezza è lo spessore della linea (circa 40-60 K, 0.01 decadi in tempo).
"T finale" è il valore alla fine della curva del paper; dove il paper dà un numero
nel testo o in figura è indicato fra parentesi. Le T di equilibrio "attese" le
calcola il programma 0D dalla conservazione dell'energia.

| Figura | condizioni | T finale paper | T finale 0D | scarto max nel transitorio (0D vs paper) |
|---|---|---|---|---|
| 3b | N2, 3000 / 10000 K, 1 atm | 4971.6 K | 4973.0 K | Ttr 2.8 %, Tv 5.1 % |
| 3a | N2, 10000 / 1000 K, 1 atm | 7623.3 K (testo: 7623.3) | 7623.3 K | Ttr 0.3 %, Tv 1.9 % |
| 4 senza E_el | N2, 30000 / 1000 K | 21936 K (fig.: 21.9 kK) | 21932 K | Ttr 0.05 %, Tv 11 % (solo t < 5 ns) |
| 4 con E_el | idem, livelli elettronici | 17659 K (fig.: 17.7 kK) | 17650 K | Ttr 0.3 %, Tv 17 % (solo t < 5 ns) |
| 5 | N2 + N, 5e22 m⁻³ ciascuno, 30000 / 1000 K | 24356 K (fig.: 24.4 kK) | 24353 K | Ttr 2.5 %, Tv 20 % |
| 6 senza V-V | N2 + O2, 5000 / 30000 K, 1 atm | 12129 K (fig.: 12.1 kK) | 12133 K | Ttr 2.7 %, Tv,N2 2.7 %, Tv,O2 6.4 % |
| 6 con V-V | idem | 12129 K | 12133 K | Ttr 2.0 %, Tv,N2 2.5 %, Tv,O2 3.4 % |
| 7 | N2 + N reagente, 30000 / 1000 K | 8292 / 8281 K; n_N2/n0 0.139, n_N/n0 1.222 | 8309 / 8312 K; 0.138, 1.225 | Ttr 2.3 %, Tv 27 %; n: 0.006 decadi |
| 8 | idem da 30000 / 30000 K | 12101 / 12043 K; 0.091, 1.319 | 12121 / 12141 K; 0.089, 1.322 | Ttr 2.0 %, Tv 10 %; n: 0.01 decadi |
| 9 | aria 5 specie, 0.063 atm, 10000 K, una T | 5074 K; vedi §3.6 | 5096 K | T 0.7 %; n: 0.004-0.08 decadi |

Le temperature finali tornano ovunque entro l'incertezza di lettura. Gli scarti nel
transitorio si concentrano dove le curve sono ripide e sono quasi tutti uno
spostamento nel tempo del rilassamento vibrazionale (vedi §3.4). Grafici in
`applications/test/nonEqTTv/output/figX.png` (e `figX-n.png` per le densità).

---

## 3. Scelte di modello e perché

### 3.1 Dati termodinamici del paper

`mutation-data/thermo/species.xml` contiene le tabelle A1 e A2 del paper: θ_v
(N2 3371 K, O2 2256 K, NO 2719 K), entalpie di formazione, livelli elettronici
convertiti da K a cm⁻¹ (÷ 1.438777). I livelli di N2 coincidono con i default di
Mutation++: la vecchia nota "i dati dell'appendice non riproducono la fig. 4" era
un errore di unità nel file precedente.

### 3.2 Energia elettronica: accesa solo per la fig. 4 "with E_el"

Il paper la dichiara solo per la fig. 4. Per le fig. 5, 6, 7, 8 le temperature
finali tornano **solo senza** livelli elettronici: con i livelli di N2 e N la T di
equilibrio della fig. 5 sarebbe 18308 K invece di 24.4 kK, quella della fig. 6
19241 K invece di 12.1 kK. Il testo della fig. 5 cita il modo elettronico di N, ma
i numeri della figura dicono il contrario; si è seguita la figura. Per questo ci
sono due cartelle dati: `mutation-data` e `mutation-data-noElectronic` (solo il
livello fondamentale, link simbolici per il resto), scelte con
`MPP_DATA_DIRECTORY`.

### 3.3 Forza motrice del V-T: e_ve = e_v + e_el come nell'eq. 8

L'`OmegaVT` di Mutation++ usa come forza motrice solo l'energia dell'oscillatore
armonico; il paper usa e_ve = e_v + e_el della molecola che rilassa. Per la fig. 4
con E_el la differenza è enorme: con Mutation++ gli scarti erano 28 % su Ttr e
56 % su Tv, con l'eq. 8 scritta come nel paper (`sourceVT` in `mutationSources.H`,
τ sempre quello di Mutation++) scendono a 0.3 % e 4 % (17 % solo nei primi 4 ns).
Senza livelli elettronici le due forme coincidono, quindi le altre figure non
cambiano. È una scelta di modello, applicata ovunque, non una correzione locale.

### 3.4 Tempi di rilassamento: quelli di Mutation++

τ_VT viene da `MillikanWhiteModel` di Mutation++ (costanti A, B di Park 1993 dal
suo `VT.xml`, correzione di Park con σ = 3e-21 m²). Per N2 puro coincide con le
eq. 9-17 del paper. Per le miscele Mutation++ media τ_MW aritmeticamente sui
partner e applica la correzione di Park una sola volta con la densità numerica
del vibratore, mentre il paper (eq. 9-10, e il sorgente di hy2Foam) fa la media
armonica coppia per coppia con la densità della coppia n_m + n_s. Da qui gli
scarti su Tv nel tratto ripido delle fig. 5 (20 %), 7 (27 %) e 8 (10 %): un
ritardo del rilassamento vibrazionale, con le T finali giuste. Si è tenuto
Mutation++ per non riscrivere il modello; la convenzione del paper è documentata
qui se in futuro servisse.

### 3.5 Chimica: temperatura di Park con esponente 0.7 e Q_C-V non preferenziale

Mutation++ valuta le dissociazioni a `sqrt(T Tv)` (esponente 0.5) e non è
configurabile; il paper usa `T^0.7 Tv^0.3` (eq. 29). Senza toccare la libreria,
`productionRates` le passa lo stato `{T_P, T_P}` prima di chiedere le velocità e
poi ripristina `{T, Tv}`. Vale per meccanismi di sola dissociazione (`N2_Park`,
tabella 2 del paper, A in cm³/mol/s). Effetto dell'esponente (programma 0D):

| Figura | esponente | Ttr | Tv | n_N2 | n_N |
|---|---|---|---|---|---|
| 7 | 0.5 (Mutation++) | 6.6 % | 27 % | 0.025 decadi (6 %) | 0.038 decadi (9 %) |
| 7 | 0.7 (paper) | 2.3 % | 27 % | 0.006 decadi (1.4 %) | 0.006 decadi (1.5 %) |
| 8 | 0.5 (Mutation++) | 3.3 % | 9.5 % | 0.019 decadi (4 %) | 0.010 decadi (2 %) |
| 8 | 0.7 (paper) | 2.0 % | 10 % | 0.010 decadi (2 %) | 0.004 decadi (1 %) |

Con 0.7 la composizione della fig. 7 passa dal 6-9 % all'1.5 %; Tv non cambia
perché il suo scarto viene dal V-T (§3.4). L'accoppiamento chimica-vibrazione è il
modello non preferenziale (eq. 30-31, `sourceCV`), l'unico di Mutation++; il paper
non dichiara quale usa e il tutorial di hy2Foam usa il preferenziale con α = 0.3,
ma con il non preferenziale la composizione torna entro l'1.5 %, quindi non si è
provato altro.

### 3.6 Fig. 9: costanti QK, non Park

Il paper usa 19 reazioni irreversibili (15 dissociazioni, 4 scambi) con le costanti
QK di Scanlon et al. 2015, che cita ma non riporta (la tabella 2 ha solo N2 + N2).
Si sono provati i due set:

| meccanismo | T (max) | N2 | O2 | NO | N | O |
|---|---|---|---|---|---|---|
| `air5_Park` di Mutation++ (5 reazioni reversibili, Park) | 8.9 % | 0.009 dec | 0.36 dec | 0.57 dec | 0.41 dec | 0.55 dec |
| `air5_QK` da hyStrath (19 irreversibili, QK) | 0.7 % | 0.004 dec | 0.06 dec | 0.62 dec (*) | 0.05 dec | 0.08 dec |

(*) solo nel primo punto della curva di NO, che entra nel grafico dal bordo
inferiore (1e-5); dopo, NO è entro 0.05 decadi. Con Park le densità di N, O, NO
sono lontane di un fattore 2.5-3.7 nel transitorio: la fig. 9 è fatta con le
costanti QK. Il set QK è in `hystrath-data/` con README (URL, commit, licenza
GPL-3); Mutation++ non accetta nello stesso meccanismo due reazioni irreversibili
una inversa dell'altra, quindi i due scambi inversi stanno in un secondo file e
`Test-air5` somma le velocità dei due.

### 3.7 Altre scelte

- Scambio V-V (fig. 6): eq. 18 del paper con P = 0.01; la sezione d'urto N2-O2
  non è nel paper, si è preso il valore di hy2Foam (2.667e-19 m²), citato nel
  commento. L'effetto del V-V è quello atteso: Tv,O2 e Tv,N2 si avvicinano e il
  rilassamento è più rapido.
- Caso solver: thermo base `rrho`, quello del template del corso (la libreria lo
  istanzia accanto ai thermo del core in `highEnthalpyMulticomponentThermos.C`;
  `rrhoThermo` è una copia di `janafThermo`, e con `janaf` i numeri sono
  identici). La termodinamica a due temperature non viene dai polinomi di
  OpenFOAM ma dal ponte con Mutation++: il thermo base serve solo a costruire i
  campi. `T` riletta dal file perché il thermo base la limitava a 20000 K; tolto
  `limitTemperature` che tagliava l'energia sopra 20000 K.
- Nessun campo `e`, `eve`, `T` viene più toccato dal thermo base: tutto passa da
  Mutation++ (`setState` con energie, `vars = 0`).

---

## 4. Limiti noti

- Il rilassamento vibrazionale nelle miscele è più lento di quello del paper
  (§3.4): fino al 20-27 % su Tv nel tratto ripido di fig. 5 e 7.
- Nei primi nanosecondi delle fig. 4 gli scarti relativi su Tv sono alti (11-17 %)
  perché Tv è piccola e le curve del paper partono da t = 1 ns con la loro
  incertezza di lettura.
- Il solver risolve `eve` senza trasporto: va bene per l'heat bath, non ancora per
  un caso con flusso.
- Il thermo base rrho emette avvisi sopra 20000 K (innocui: T e psi le fa
  Mutation++).
- Un milione di passi nel solver a cella singola costano circa 28 minuti (contro
  2 secondi del programma 0D): l'overhead è di OpenFOAM, non di Mutation++.
- Le energie elettroniche degli atomi stanno nel serbatoio vibro-elettronico come
  in Mutation++ (in hy2Foam seguono Tve ma pesano sull'energia traslazionale):
  irrilevante per le figure riprodotte, tutte senza livelli elettronici degli atomi.

---

## 5. Come rilanciare

```sh
# programmi 0D, tutte le figure, con confronti e grafici
applications/test/nonEqTTv/Allwmake
applications/test/nonEqTTv/Allrun

# solver, una figura per volta (fig3a fig3b fig4-noEl fig4-el fig5 fig7)
applications/test/nonEqTTv/solverHeatBath/Allrun fig7        # esponente 0.7
applications/test/nonEqTTv/solverHeatBath/Allrun fig7 0.5    # esponente 0.5
python3 applications/test/nonEqTTv/compare.py fig7

# curve del paper (servono solo se si vogliono rigenerare)
python3 applications/test/nonEqTTv/paper-data/extract-paper-curves.py
```

Struttura: `common/heatBath.H` (energie ed equilibrio per i test 0D),
`src/.../mutationSources.H` (sorgenti condivise con il solver), `N2/`, `N2N/`,
`N2O2/`, `air5/` (un programma per miscela), `paper-data/` (curve estratte),
`hystrath-data/` (set QK), `solverHeatBath/` (caso a una cella, condizioni scelte
da `Allrun <figura>`), `output/` (csv e grafici, campionati a non più di 10000
righe).

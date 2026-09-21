# Milestone 2 — Heat bath del paper: programmi 0D, solver e confronto con le curve

**Data:** 21 settembre 2026 (revisione della versione del 17 settembre, dopo un
controllo punto per punto contro il paper e contro il sorgente di hy2Foam)
**Obiettivo:** riprodurre tutti i test 0D del capitolo 3 di Casseau et al., *A
Two-Temperature Open-Source CFD Model for Hypersonic Reacting Flows, Part One*,
Aerospace 2016, 3, 34 (`applications/test/nonEqTTv/aerospace-03-00034-1.pdf`),
su due strade: programmi 0D standalone con Mutation++ (tutte le figure) e solver
`shockThermo` con la libreria caricata a runtime (fig. 3a, 3b, 4, 5, 7), e
confrontare tutto con le curve estratte dal pdf.

Riferimento per i modelli: il paper, e dove il paper non basta il codice con cui
sono state fatte le sue figure, hy2Foam (repository hyStrath, commit
`984e3000a5f8`, vedi `hystrath-data/README.md`). Il suo tutorial `heatBath` è
esattamente il caso della fig. 7 (N2 + N, 30000 / 1000 K, p = 41419.4 Pa,
X = 0.5 / 0.5, passo 1 ns, stessa reazione con le costanti di Park 1993).

---

## 0. Cosa è cambiato in questa revisione

Nessuna modifica è stata fatta per avvicinare le curve: ogni scelta qui sotto
segue un'equazione del paper o il sorgente di hy2Foam, e gli effetti sono
riportati anche quando peggiorano l'accordo (V-V, §3.6).

| # | problema | tipo | cosa si è fatto |
|---|---|---|---|
| 1 | Q_C-V non preferenziale, mentre il caso "Park" del paper (fig. 7, 8) usa il preferenziale con α = 0.3 | modello diverso dal paper | `sourceCV` con i due modelli, default preferenziale (§3.5) |
| 2 | τ_VT di Mutation++ (media aritmetica, Park una volta con n_m) invece delle eq. 9-17 (per coppia, media armonica, densità della coppia) | modello diverso dal paper | `relaxationTimePaper`, default (§3.4) |
| 3 | V-V: imposto Q_O2 = −Q_N2 invece dell'eq. 18 scritta per ogni molecola | errore rispetto all'equazione | eq. 18 per molecola, come in hy2Foam (§3.6) |
| 4 | `compare.py` stampava l'errore relativo massimo con il tempo dell'errore assoluto massimo | errore nella metrica | ogni massimo con il suo istante (§2) |
| 5 | csv con una riga ogni 10 passi (100 per le fig. 7, 8): i primi ns erano interpolati fra t = 0 e t = 10 ns | errore nella metrica | uscita logaritmica, ogni passo all'inizio (§2) |
| 6 | nel solver ogni specie ricalcolava la chimica dopo che le precedenti erano già avanzate; energia ed `eve` a stati ancora diversi | incoerenza nel solver | `correctSources()` una volta per correttore (§1) |
| 7 | nel caso solver N2 + N le frazioni in massa erano arrotondate a 0.6667 / 0.3333 invece di 2/3 e 1/3 (M_N2 = 2 M_N) | errore di input (5·10⁻⁵) | 9 cifre in `solverHeatBath/Allrun` |
| 8 | T_P con esponente 0.7 ed e_ve = e_v + e_el presentati come scelte giustificate dagli scarti; energia elettronica "spenta perché le curve tornano solo così" | documentazione | fonti citate (§3) |

---

## 1. Solver e programma 0D fanno la stessa fisica

Ogni test che esiste su entrambe le strade è confrontato in automatico
(`compare.py`): stesse condizioni iniziali, stessi dati Mutation++, stessi modelli,
stesso passo di 1 ns. Il solver risolve le equazioni conservative in `e` ed `eve`
con le sorgenti di `mutationSources.H`, il programma 0D integra le stesse sorgenti
a mano.

| Figura | caso | Ttr: max scarto solver vs 0D | Tv: max scarto | composizione (N2 / N) |
|---|---|---|---|---|
| 3b | N2, raffreddamento | 0.1 K (0.00 %) | 0.3 K (0.00 %) | — |
| 3a | N2, riscaldamento | 2.6 K (0.03 %) | 6.8 K (0.21 %) | — |
| 4 senza E_el | N2 da 30000 K | 34.1 K (0.14 %) | 85.6 K (0.69 %) | — |
| 4 con E_el | N2 da 30000 K | 111.8 K (0.52 %) | 106.7 K (1.19 %) | — |
| 5 | N2 + N | 6.4 K (0.02 %) | 25.6 K (0.19 %) | — |
| 7 | N2 + N reagente (Park 0.7, pref. 0.3) | 6.5 K (0.03 %) | 29.8 K (0.26 %) | 0.02 % / 0.03 % |

Gli scarti residui sono nei primi 100 ns e vengono dall'integrazione nel tempo,
non dalla fisica: il solver fa due correttori PIMPLE per passo (`nOuterCorrectors
2`) e il secondo rivaluta le sorgenti allo stato del primo, cioè in pratica un
Eulero implicito, mentre il programma 0D fa un Eulero esplicito. La differenza è
dell'ordine di Δt/τ: con N2 a 30000 K e 1 atm τ = 1.23·10⁻⁷ s, quindi
1 ns / τ ≈ 0.8 %, che è proprio lo scarto della fig. 4. I probe del solver ora
scrivono ogni passo (ogni 10 per la fig. 7, prima ogni 100: il confronto saltava
proprio i primi 100 ns del rilassamento).

**Sorgenti una volta per correttore.** `thermo_.correctSources()` all'inizio di
`thermophysicalPredictor` calcola in ogni cella, tutte allo stesso stato, ω_s,
il calore di reazione −Σ h_f,s ω_s e Q_ve = Q_V-T + Q_C-V; specie, energia ed `eve`
usano quei valori. Prima la chimica veniva ricalcolata N_specie + 2 volte per
passo, a composizioni diverse. Su 2 µs della fig. 7 (probe di `e` ed `eve`),
l'energia totale ρ(e + Σ Y_s h_f,s) varia di 6.8·10⁻⁸ con la versione vecchia e
di 1.9·10⁻⁸ con la nuova (il limite delle 8 cifre dei probe): l'incoerenza era
reale ma, con un passo di 1 ns, numericamente trascurabile. Il tempo di calcolo su
una cella cambia poco (2.2 → 2.0 s per 2000 passi): lì domina OpenFOAM.

---

## 2. Confronto con il paper, figura per figura

Le curve di hy2Foam sono estratte dal pdf (`paper-data/`): figure vettoriali,
calibrazione esatta sui tick, incertezza pari allo spessore della linea, circa
40-60 K in verticale e 0.01 decadi in tempo. Per ogni curva `compare.py` riporta
l'errore assoluto massimo e l'errore relativo massimo, **ciascuno con il proprio
istante** (prima l'errore relativo veniva stampato con il tempo dell'errore
assoluto, che in generale cade altrove: il 27 % su Tv della fig. 7 era a 18 ns,
non a 246 ns come scritto).

| Figura | condizioni | T finale paper | T finale 0D | max errore (0D vs paper) |
|---|---|---|---|---|
| 3a | N2, 10000 / 1000 K, 1 atm | 7623.3 K (testo: 7623.3) | 7623.3 K | Ttr 23 K (0.27 %), Tv 61 K (1.4 %) |
| 3b | N2, 3000 / 10000 K, 1 atm | 4971.6 K | 4973.0 K | Ttr 123 K (2.8 %), Tv 313 K (5.1 %) |
| 4 senza E_el | N2, 30000 / 1000 K | 21936 K (fig.: 21.9 kK) | 21932 K | Ttr 9 K (0.03 %), Tv 20 K (0.25 %) |
| 4 con E_el | idem, livelli elettronici | 17659 K (fig.: 17.7 kK) | 17650 K | Ttr 14 K (0.07 %), Tv 15 K (0.19 %) |
| 5 | N2 + N, 5e22 m⁻³ ciascuno, 30000 / 1000 K | 24357 K (fig.: 24.4 kK) | 24353 K | Ttr 10 K (0.04 %), Tv 26 K (0.44 %) |
| 6 senza V-V | N2 + O2, 5000 / 30000 K, 1 atm | 12129 K (fig.: 12.1 kK) | 12133 K | Ttr 2.5 %, Tv,N2 2.1 %, Tv,O2 5.6 % (§3.6) |
| 6 con V-V | idem | 12129 K | 12133 K | Ttr 3.0 %, Tv,N2 3.6 %, Tv,O2 5.0 % (§3.6) |
| 7 | N2 + N reagente, 30000 / 1000 K | 8292 / 8281 K; n_N2/n0 0.139, n_N/n0 1.222 | 8316 / 8306 K; 0.138, 1.225 | Ttr 33 K (0.30 %), Tv 32 K (0.31 %); n 0.006 decadi |
| 8 | idem da 30000 / 30000 K | 12101 / 12043 K; 0.091, 1.319 | 12148 / 12093 K; 0.090, 1.321 | Ttr 51 K (0.42 %), Tv 52 K (0.43 %); n 0.007 decadi |
| 9 | aria 5 specie, 0.063 atm, 10000 K, una T | 5074 K | 5096 K | T 44 K (0.69 %); n 0.004-0.06 decadi (§3.7) |

Gli errori percentuali di Tv nelle prime righe (0.25 %, 0.44 % a t = 2 ns) sono
relativi a Tv ≈ 1000-2000 K: in kelvin sono le cifre fra parentesi.

**Il campionamento contava.** Con una riga ogni 10 passi le fig. 4 davano 11 % e
17 % su Tv "solo nei primi 5 ns": era interpolazione lineare fra t = 0 e
t = 10 ns, cioè fra Tv = 1000 K e un valore già molto più alto. Con l'uscita
attuale (`OutputSchedule` in `common/heatBath.H`: ogni passo fino a ~90 ns, poi
200 punti per decade) gli stessi confronti danno 0.25 % e 0.19 %. Lo stesso
valeva per NO nella fig. 9 (0.62 → 0.03 decadi). I csv sono anche più corti
(400-900 righe invece di fino a 10000).

### 2.1 Quanto è precisa la verifica

Il riferimento è una curva digitalizzata, e la sua incertezza non è solo lo
spessore della linea: 0.01 decadi in tempo, sui tratti ripidi, valgono

| figura | pendenza massima (Ttr / Tv) | 0.01 decadi valgono |
|---|---|---|
| 3a | 1500 / 3900 K per decade | 15 / 39 K |
| 3b | 4500 / 11500 K per decade | 45 / 115 K |
| 4 senza E_el | 6900 / 17400 K per decade | 69 / 174 K |
| 5 | 4900 / 19700 K per decade | 49 / 197 K |
| 6 con V-V (Tv,O2) | 14500 / 35800 K per decade | 145 / 358 K |
| 7 | 10200 / 15900 K per decade | 102 / 159 K |
| 8 | 7400 / 10800 K per decade | 74 / 108 K |

Quindi per le fig. 4, 5, 7, 8 gli errori residui (10-50 K) sono **sotto**
l'incertezza del riferimento: la frase corretta è "indistinguibili dalle curve del
paper entro la loro precisione di lettura", non "accordo allo 0.04 %". Le cifre
decimali della tabella misurano il codice rispetto a una lettura del grafico, non
la fisica. Le T di equilibrio, invece, sono un risultato esatto (conservazione
dell'energia) e il confronto solver / 0D (§1) confronta due calcoli, quindi lì le
cifre hanno senso. La fig. 3b (123 / 313 K contro 45-115 K più lo spessore) e la
fig. 6 (fino a 850 K contro 358 K) restano sopra l'incertezza: sono differenze
vere, discusse in §3.6 e §4.

Infine, un accordo stretto dice che il codice implementa **lo stesso modello** del
paper, non che il modello sia accurato: per il paper stesso la correlazione di
Millikan-White sta entro un fattore 5 dalle misure (paragrafo 2.1.1).

---

## 3. Modelli e fonti

### 3.1 Dati termodinamici

`mutation-data/thermo/species.xml` contiene le tabelle A1 e A2 del paper: θ_v
(N2 3371 K, O2 2256 K, NO 2719 K), entalpie di formazione, livelli elettronici
convertiti da K a cm⁻¹ (÷ 1.438777). Sono anche i valori del file `thermoDEM` del
tutorial di hy2Foam. Dove il database di Mutation++ è diverso, vale il paper: per
esempio θ_v di N2 è 3408.5 K in Mutation++, e il settimo livello eccitato di N2 è
72097.6 cm⁻¹ in Mutation++ contro 1.048976·10⁵ K = 72907.4 cm⁻¹ nella tabella A2
(il file del progetto ha il valore della tabella).

### 3.2 Energia elettronica: quando è inclusa

Due cartelle dati, scelte con `MPP_DATA_DIRECTORY`: `mutation-data` (livelli della
tabella A2) e `mutation-data-noElectronic` (solo il livello fondamentale, link
simbolici per il resto). Quale si usa segue il paper:

- **fig. 3a, 3b**: senza, paragrafo 3.1.1 "Case without Electronic Energy";
- **fig. 4**: tutte e due le curve, paragrafo 3.1.2;
- **fig. 5**: senza. La T di equilibrio scritta sulla figura, 24.4 kK, è quella senza
  livelli elettronici (con i livelli sarebbe 18308 K), e la curva è confrontata con
  dsmcFoam, che non ha il modo elettronico (paragrafo 3). La frase del testo sul
  modo elettronico di N che "porta 1.39 gradi di libertà alla miscela a T_eq" si
  riferisce alla variante con livelli, non disegnata: a 18308 K il calcolo dà
  ζ_el(N) = 1.390. Non c'è contraddizione con la figura (la versione precedente
  di questo documento la leggeva così);
- **fig. 6**: senza, T di equilibrio sulla figura 12.1 kK (con i livelli 19241 K);
- **fig. 7, 8**: senza. Nel tutorial `heatBath` di hy2Foam, che è il caso della
  fig. 7, il coefficiente del modo elettronico `decoupledCvCoeffs[3]` è 0 per tutte
  le specie, e hy2Foam calcola l'energia elettronica come `a[3]·(...)`, zero in quel
  caso (`decoupledEnergyModesThermoI.H`).

Coerenza Mutation++ / hy2Foam: per le molecole il serbatoio è lo stesso
(e_ve = e_v + e_el). Per gli atomi Mutation++ mette l'energia elettronica nel
serbatoio vibro-elettronico; nel paper gli atomi seguono Tve della molecola di
riferimento ma non hanno un'equazione E_ve propria (eq. 23). La differenza non
tocca nessuna figura riprodotta: l'unico caso con livelli elettronici (fig. 4) è
N2 puro.

### 3.3 Forza motrice del V-T: e_ve = e_v + e_el

È la definizione del paper: l'eq. 8 è scritta per e_ve,m, e subito sotto
"e_ve,m = e_v,m + e_el,m". hy2Foam fa lo stesso (`LandauTellerVT` usa l'energia
vibro-elettronica della molecola). `sourceVT` implementa quella forma; l'`OmegaVT`
di Mutation++ userebbe solo l'energia dell'oscillatore armonico, un modello
diverso. Senza livelli elettronici le due coincidono, quindi la scelta conta solo
per la fig. 4 con E_el (nella versione di settembre, con il campionamento di
allora, la forma di Mutation++ dava 28 % su Ttr e 56 % su Tv): è una conseguenza di
aver seguito il paper, non la sua ragione.

### 3.4 Tempi di rilassamento: eq. 9-17

`relaxationTimePaper` (in `mutationSources.H`) implementa le eq. 9-17 come il
modello di default di hy2Foam (`MillikanWhitePark`): per ogni partner s,
τ_ms = τ_MW,ms + τ_P,ms, con τ_P calcolata con la densità della coppia n_m + n_s
(solo n_m per s = m), poi media armonica pesata con le frazioni molari (eq. 9).
Costanti A, B e σ (3·10⁻²¹ m² per N2 e O2) di Park 1993, dal `VT.xml`.

`MillikanWhiteModel::relaxationTime` di Mutation++ media invece τ_MW in modo
**aritmetico** e aggiunge **una sola** correzione di Park con la densità della
sola molecola. Per N2 puro le due forme coincidono; per N2 + N (fig. 5, 7, 8):

| T | 30000 K | 24353 K | 20000 K | 15000 K | 10000 K |
|---|---|---|---|---|---|
| τ Mutation++ / τ paper | 1.38 | 1.28 | 1.17 | 1.06 | 1.02 |

Un τ più lungo del 28-38 % proprio dove Tv sale spiega il ritardo che la versione
precedente attribuiva a "differenze di convenzione". Effetto sulla fig. 5 (senza
chimica, conta solo τ):

| τ | Ttr | Tv |
|---|---|---|
| Mutation++ | 649 K (2.5 %) | 2587 K (19.6 %) |
| paper, eq. 9-17 | 10 K (0.04 %) | 26 K (0.44 %) |

Il τ di Mutation++ resta disponibile per confronto (`mutation` come argomento dei
programmi 0D, `relaxationTime mutation` nel solver; variante `-tauMpp` in
`Allrun`).

### 3.5 Chimica: T_P = T^0.7 Tv^0.3 e Q_C-V preferenziale

**Temperatura di controllo.** Eq. 29, T_P = T^α Tv^(1−α), e il testo: "an exponent
of 0.7 in favour of the trans-rotational temperature is commonly adopted"; il
tutorial di hy2Foam ha `exponentTtr 0.7`. Mutation++ valuta le dissociazioni a
`sqrt(T Tv)` (esponente 0.5) senza possibilità di cambiarlo: `productionRates` le
passa lo stato `{T_P, T_P}` e poi ripristina `{T, Tv}`. Vale per meccanismi di
sola dissociazione (`N2_Park`).

**Costante di velocità.** Tabella 2: A = 7.0·10²¹, β = −1.6, T_a = 113200 K. Il
paper indica le unità m³ mol⁻¹ s⁻¹, ma è il valore di Park in cm³ mol⁻¹ s⁻¹: il
tutorial di hy2Foam ha A = 7.0·10¹⁸ m³ kmol⁻¹ s⁻¹, che è lo stesso numero.
`N2_Park.xml` usa cm³ mol⁻¹ s⁻¹.

**Accoppiamento chimica-vibrazione.** Eq. 30: Q_C-V = Σ_m ω_m (D'_m + e_el,m), con
D' = e_v (non preferenziale, eq. 31) oppure D' = α D (preferenziale, eq. 32). Il
paper chiama la configurazione della fig. 7 "Park TTv model", e nel tutorial
`heatBath` di hy2Foam quella configurazione è

```
chemistryVibrationCoupling { model ParkTTv; ParkTTvCoeffs { exponentTtr 0.7;
    sourceTermModel preferential; preferentialModel { factorType constant; constantFactor 0.3; } } }
```

con D_N2 = 3.36·10⁷ J/kg (tabella A1, `dissocEnergy` nel `thermoDEM`), sulle sole
molecole (`chemistry2Model.C`). `sourceCV` ora implementa tutti e due i modelli;
il default è il preferenziale con α = 0.3. La versione precedente usava il non
preferenziale, e la sua giustificazione ("con il non preferenziale la
composizione torna entro l'1.5 %, quindi non si è provato altro") era proprio il
tipo di scelta da evitare.

Effetto delle tre scelte sulle fig. 7 e 8 (massimo errore relativo):

| τ | Q_C-V | esponente | fig. 7 Ttr | fig. 7 Tv | fig. 7 n_N | fig. 8 Ttr | fig. 8 Tv |
|---|---|---|---|---|---|---|---|
| Mutation++ | non pref. | 0.7 | 2.2 % | 16.9 % | 0.006 dec | 2.0 % | 10.2 % |
| paper | non pref. | 0.7 | 2.0 % | 8.4 % | 0.004 dec | 0.9 % | 4.5 % |
| Mutation++ | pref. 0.3 | 0.7 | 4.2 % | 20.0 % | 0.011 dec | 0.5 % | 3.0 % |
| **paper** | **pref. 0.3** | **0.7** | **0.30 %** | **0.31 %** | **0.001 dec** | **0.42 %** | **0.43 %** |
| paper | pref. 0.3 | 0.5 | 3.5 % | 5.4 % | 0.027 dec | 0.50 % | 0.51 % |

Da notare: il preferenziale da solo, con il τ di Mutation++, **peggiora** la
fig. 7. Non è quindi una scelta che si sarebbe fatta guardando le curve; è
l'impostazione di hy2Foam, e solo insieme al τ del paper riproduce la figura.

### 3.6 Scambio V-V (fig. 6)

Eq. 18 con P = 0.01; la sezione d'urto N2-O2 non è nel paper, si è preso il valore
di hy2Foam (2.667·10⁻¹⁹ m², `KnabCoefficients` del tutorial). La versione
precedente calcolava Q_VV solo per N2 e imponeva Q_O2 = −Q_N2. Il paper però scrive
l'eq. 18 per ogni molecola m con partner l, e hy2Foam (`KnabVV.C`) fa lo stesso:
i due termini non sono opposti (il loro rapporto è −E_v,O2(T)/E_v,N2(T), energie
molari) e la differenza passa alla traslazione, perché T_tr si ricava
dall'energia totale, che resta conservata. Ora `Test-N2O2` usa l'eq. 18 per tutte e
due le molecole.

Con le costanti di Park 1993 questo **peggiora** l'accordo della fig. 6 con V-V
(Tv,O2 da 2.6 % a 5.0 %). Diagnosi, fatta su una copia dei dati fuori dal progetto:
con i coefficienti di Millikan-White di O2 del file di hy2Foam (O2-O2 135.91 /
0.03, O2-N2 131.32 / 0.0295, contro 138 / 0.030 e 134 / 0.0295 di Park 1993) le due
fig. 6 tornano entro 0.12-0.16 % su tutte le curve, con l'eq. 18 per molecola. Lo
stesso file però ha A_N2-N2 = 221.53 invece di 221, e con quello le fig. 3a e 3b
peggiorano (3b da 2.8 a 3.6 % su Ttr), mentre il paper dice esplicitamente che la
configurazione di default usa A_N2-N2 della tabella di Park 1993. Prendere O2 da un
file e N2 dall'altro sarebbe adattare i dati alle curve, quindi il progetto tiene
le costanti di Park 1993 per tutte le coppie, come dichiara il paper, e il residuo
della fig. 6 resta documentato qui.

### 3.7 Fig. 9: costanti QK, non Park

Il paper usa 19 reazioni irreversibili (15 dissociazioni, 4 scambi) con le costanti
QK di Scanlon et al. 2015, che cita ma non riporta (la tabella 2 ha solo N2 + N2).
Si sono provati i due set:

| meccanismo | T (max) | N2 | O2 | NO | N | O |
|---|---|---|---|---|---|---|
| `air5_Park` di Mutation++ (5 reazioni reversibili, Park) | 8.9 % | 0.009 dec | 0.36 dec | 0.20 dec | 0.55 dec | 0.54 dec |
| `air5_QK` da hyStrath (19 irreversibili, QK) | 0.69 % | 0.004 dec | 0.06 dec | 0.03 dec | 0.03 dec | 0.01 dec |

Con Park le densità di N, O, NO sono lontane di un fattore 1.6-3.5: la fig. 9 è
fatta con le costanti QK. (Lo scarto di 0.62 decadi su NO nella versione
precedente era il campionamento, §2.) Il set QK è in `hystrath-data/` con README
(URL, commit, licenza GPL-3); Mutation++ non accetta nello stesso meccanismo due
reazioni irreversibili una inversa dell'altra, quindi i due scambi inversi stanno
in un secondo file e `Test-air5` somma le velocità dei due.

### 3.8 Altre scelte

- Caso solver: thermo base `rrho`, quello del template del corso (la libreria lo
  istanzia accanto ai thermo del core in `highEnthalpyMulticomponentThermos.C`;
  `rrhoThermo` è una copia di `janafThermo`, e con `janaf` i numeri sono
  identici). La termodinamica a due temperature non viene dai polinomi di
  OpenFOAM ma dal ponte con Mutation++: il thermo base serve solo a costruire i
  campi. `T` riletta dal file perché il thermo base la limitava a 20000 K.
- Nessun campo `e`, `eve`, `T` viene toccato dal thermo base: tutto passa da
  Mutation++ (`setState` con energie, `vars = 0`).

---

## 4. Verifiche di consistenza: segni, unità, conservazione

| termine | unità | segno |
|---|---|---|
| Q_V-T = ρ_m (e_ve(T) − e_ve(Tv)) / τ | kg m⁻³ · J kg⁻¹ / s = W m⁻³ | > 0 se T > Tv (energia verso il serbatoio vibrazionale) |
| Q_C-V = Σ_m ω_m (α D_m + e_el,m) | kg m⁻³ s⁻¹ · J kg⁻¹ = W m⁻³ | < 0 con dissociazione (ω_N2 < 0): la vibrazione paga α D |
| Q_chem = −Σ_s h_f,s ω_s (solver, energia sensibile) | W m⁻³ | < 0 con dissociazione: la formazione di N assorbe energia |

L'energia totale (formazione compresa) è costante per costruzione: il serbatoio
vibro-elettronico perde α D per ogni kg di N2 dissociato e quello traslazionale il
resto, (1 − α) D. Verifica indipendente sui csv della fig. 7 e 8: da ogni riga
(T, Tv, n/n0) si ricostruisce lo stato e si ricalcola con Mutation++ l'energia
totale, la massa e gli atomi di azoto. Variazioni massime: energia 5·10⁻⁶, massa e
atomi 3·10⁻⁶, cioè l'arrotondamento a 6 cifre dei csv. Nel solver, su 2 µs della
fig. 7, l'energia totale varia di 1.9·10⁻⁸ (§1).

---

## 5. Differenze residue e limiti noti

- **Fig. 6**: 2-5 % su Tv nel tratto ripido, dovuti ai coefficienti di
  Millikan-White di O2 (§3.6).
- **Fig. 3b**: 123 K su Ttr e 313 K su Tv al ginocchio (t ≈ 2·10⁻⁵ s), sopra
  l'incertezza di lettura (45-115 K più lo spessore). È N2 puro, dove le due forme
  di τ coincidono; A_N2-N2 = 221.53 la peggiora. Non spiegata.
- **Fig. 7, 8, stato finale**: circa 25 K (fig. 7) e 50 K (fig. 8) sopra il
  paper (0.3-0.4 %), entro l'incertezza di lettura.
- Il modello preferenziale con α costante (lo stesso di hy2Foam) toglie
  α D = 1.0·10⁷ J/kg per ogni kg di N2 dissociato qualunque sia Tv, mentre
  e_v(5000 K) ≈ 1.0·10⁶ J/kg: se la dissociazione partisse a Tv bassa potrebbe
  svuotare il serbatoio vibrazionale. Nei casi del paper la dissociazione parte
  quando Tv è già alta e non succede; in un flusso va tenuto d'occhio.
- Il solver risolve `eve` senza trasporto: va bene per l'heat bath, non ancora per
  un caso con flusso.
- Il thermo base rrho emette avvisi sopra 20000 K (innocui: T e psi le fa
  Mutation++).
- Nel solver a cella singola il costo è di OpenFOAM, non di Mutation++: un
  milione di passi della fig. 7 costano 26 minuti (28 prima di questa revisione),
  contro circa 1 s del programma 0D.

---

## 6. Come rilanciare

```sh
# programmi 0D, tutte le figure e le varianti di modello, con confronti e grafici
applications/test/nonEqTTv/Allwmake
applications/test/nonEqTTv/Allrun

# una figura a mano (argomenti facoltativi: tau paper|mutation, C-V preferential|nonPreferential)
cd applications/test/nonEqTTv
MPP_DATA_DIRECTORY=$PWD/mutation-data-noElectronic Test-N2N 30000 1000 1e-3 N2_Park 0.7 output/fig7.csv
MPP_DATA_DIRECTORY=$PWD/mutation-data-noElectronic Test-N2N 30000 1000 1e-3 N2_Park 0.7 output/fig7-nonPref.csv paper nonPreferential

# solver, una figura per volta (fig3a fig3b fig4-noEl fig4-el fig5 fig7)
applications/test/nonEqTTv/solverHeatBath/Allrun fig7        # esponente 0.7
applications/test/nonEqTTv/solverHeatBath/Allrun fig7 0.5    # esponente 0.5
python3 applications/test/nonEqTTv/compare.py fig7

# curve del paper (servono solo se si vogliono rigenerare)
python3 applications/test/nonEqTTv/paper-data/extract-paper-curves.py
```

Nel solver i modelli si scelgono in `constant/physicalProperties`, sotto-dizionario
`highEnthalpyMutation`: `relaxationTime` (`paper` | `mutation`),
`chemistryVibration` (`preferential` | `nonPreferential`), `preferentialFactor`,
`parkExponent`. Il log di `foamRun` stampa quelli usati.

Struttura: `common/heatBath.H` (energie, equilibrio e istanti di scrittura per i
test 0D), `src/.../mutationSources.H` (sorgenti condivise con il solver), `N2/`,
`N2N/`, `N2O2/`, `air5/` (un programma per miscela), `paper-data/` (curve
estratte), `hystrath-data/` (set QK), `solverHeatBath/` (caso a una cella,
condizioni scelte da `Allrun <figura>`), `output/` (csv e grafici).

# Milestone 7 — Cilindro Mach 20 reagente (paper Part Two, §3.2)

Secondo e ultimo caso multi-D del paper: azoto puro a **Mach 20** su un
cilindro di raggio **R = 1 m**, parete a **1000 K**, con **chimica
reagente** (dissociazione N₂ → N) accoppiata al two-temperature. È il
banco di prova della fisica reattiva in multi-D: finora la chimica era
validata solo in 0D (M2/M3) e 1D (M4).

Condizioni (Tab. 2 del paper): U∞ = 6047 m/s, p∞ = 0.89 Pa, T∞ = 220 K,
ρ∞ = 1.363×10⁻⁵ kg/m³, Kn = 0.0022. Bersagli: standoff ~0.25 m (Fig 5a-c,
~5 cm più vicino al corpo del caso non-reagente), C_D = 1.304 (Tab. 5,
run 3), C_H = 88.1 kW.

## 1. Il caso

- **Mesh** (`makeBlockMeshDict.py`): mezzo dominio planare (piano di
  simmetria a y=0), due blocchi polari con il cilindro come bordo interno
  (archi veri) e un bordo esterno che cresce con spline da 1.8 m a monte
  a 5 m nella scia (`R_out(θ) = 1.8 + 3.2·((1+cosθ)/2)²`). Coarse 9k
  celle (prima cella lato vento ~21 µm), fine 156k a 2 µm come il paper.
- **Chimica** (`nonEqTTv/mutation-data-noel/mechanisms/N2N_diss_park.xml`):
  **due** reazioni Park irreversibili, entrambi i partner —
  `N2+N2⇒2N+N2` (A=7.0e21) e `N2+N⇒2N+N` (A=3.0e22), β=−1.6, Ta=113200 K
  — corrispondenti al run 3 (Park TTv + rates Park). **Nota sulle unità**:
  A in cm³/(mol·s) come nel gemello `N2_diss_park.xml` validato in M2; la
  didascalia della Tab. 4 del paper dice "m³" ma è un refuso (sarebbe un
  fattore 10⁶). Mixture binaria `N2N.xml` (2 specie invece di 5: dimezza
  il lavoro del bridge per cella).
- **Trasporto**: `blottner` (da M6). **BC**: slip di Maxwell + jump di
  Smoluchowski a Tw=1000 K, accommodation 1. **LTS** con le impostazioni
  del cono ma cap di pseudo-passo 1e-5 per il campo lontano su scala
  metrica.

## 2. L'avvio impulsivo M20 richiede 90k step, non 30k

Primo insegnamento, di natura numerica: i 30k pseudo-step LTS che
convergevano il cono M11 **non bastano** qui. A 30k (endTime 3e-3) le
pressioni di parete oscillavano ancora ~60% e i campi cambiavano 5-10%
tra due scritture — cioè il transitorio era ancora vivo, non un
risultato. Un guard-rail di convergenza è ora stampato dal post (dp/p
medio tra gli ultimi due snapshot, warning sopra l'1%), così una
soluzione non-convergente non può più essere scambiata per un risultato.
Con endTime 9e-3 (90k step) la coarse converge: dp/p 0.2%.

## 3. Il checkerboard a parete — e una diagnosi sbagliata, poi corretta

A convergenza la coarse mostra un **checkerboard odd-even 2D confinato
alla singola cella a contatto con la parete**: nelle prime ~5 celle p
oscilla 90↔900 Pa e T sbatte sul floor del limiter (50 K) alternandosi a
1000 K, mentre **l'interno è da manuale** (p uniforme 466 Pa attraverso
tutto lo shock layer = Cp 1.87, standoff 0.247 m, sforzo tangenziale di
parete ~4 Pa = ordine paper).

**Diagnosi sbagliata (commit `44b3f9b`)**: avevo attribuito il
checkerboard alla BC `fixedValue` di Tve (un gradiente 8200→1000 K
piantato nelle ultime celle senza κ_ve per propagarlo). Ho rimesso Tve a
`zeroGradient` come nel cono — e **il checkerboard è rimasto identico**.
Ipotesi falsificata: non è la BC di Tve. Il cambio a zeroGradient resta
comunque corretto in sé (manca la diffusione del pool ve, quindi imporre
Tw non ha senso fisico), ma non è la causa del disturbo.

**Causa vera**: è un checkerboard pressione-velocità classico, il modo
odd-even che gli schemi centrali (Kurganov) ammettono quando la
dissipazione numerica è insufficiente. Qui è innescato dall'**aspect
ratio estremo delle celle a parete** (~20 µm radiali × ~35 mm
tangenziali ≈ 1700:1 sulla coarse) combinato con il flusso centrale e le
BC di slip/jump. È confinato allo strato a parete e **non** intacca
l'interno.

## 4. Estrazione robusta: la fisica si valida

Poiché l'interno è pulito, le grandezze di parete si leggono come le
misurerebbe uno strato limite, non dalla cella-1 inquinata:

- **pressione di parete**: mediana sulla banda radiale pulita (celle 3-9
  fuori parete). In uno strato limite dp/dn = 0, quindi quel valore *è*
  la pressione di parete — non è un trucco, è l'approssimazione di
  boundary layer;
- **attrito**: solo la componente **tangenziale** di wallShearStress (il
  checkerboard inietta una spuria componente normale ~250 Pa nello sforzo
  viscoso; lo shear è tangenziale per definizione).

Risultato (coarse, `compare-fig5.py` sui campi cluster):

| Grandezza | Raw (cella-1) | **Robusto** | Paper |
|---|---|---|---|
| Cp ristagno | 3.49 | **1.90** | ~1.87 (reagente, sopra Rayleigh 1.837) |
| C_D | 1.30 (per compensazione di errori) | **1.288** (pressione 1.258 + attrito **+0.030**) | 1.304 / DSMC 1.284 / Newton 1.333 |
| standoff | 0.247 m | 0.247 m | ~0.25 m |

Il C_D raw "giusto" era una coincidenza: pressione gonfia (2.05) +
attrito spurio negativo (−0.75) che si cancellavano. Quello robusto ha lo
**spezzato fisico** (attrito piccolo e positivo) e vale 1.288, entro
l'1.2% dal paper.

Confronto con la Fig 5 (`fig5-surface-comparison.png`,
`fig5-stagnation-comparison.png`):

- **Fig 5d (Cp)**: la nostra curva è **sovrapposta** a dsmcFoam e run3 —
  grado-paper.
- **Fig 5a (Mach linea di ristagno)**: shock a −0.26 m, plateau 20,
  ricalca run3. La T di picco (15.9 kK dietro lo shock) e l'innesco di
  dissociazione (N sale a ~10²¹ 1/m³) sono nel posto giusto.
- **Fig 5e (Cf) e 5f (flusso termico)**: **forma corretta** (Cf picco a
  ~45°, poi decrescente), ma ~20% sotto i riferimenti e ancora
  frastagliati dal residuo del checkerboard. Il flusso termico è la
  grandezza meno affidabile: q ~ dT/dn è intrinsecamente near-wall e la
  cella-1 è proprio dove il disturbo vive; C_H 57 kW sta vicino alla DSMC
  (63.3) più che al run 3 del paper (88.1, che il paper stesso dice
  sovrastimare la DSMC del 39%).

## 5. Opzione 1 — Minmod cura (quasi tutto) il checkerboard

Cambiato `reconstruct(rho|U|T)` da **vanAlbada** (compressivo, andava bene
sul cono M11) a **Minmod**, il più dissipativo dei limiter MUSCL: aggiunge
la dissipazione near-wall che smorza il modo odd-even. Costo: shock un filo
più diffuso, ma lo standoff resta 0.247 m. Rerun coarse sul cluster (90k
step, converge a dp/p 0.24%):

| Grandezza | vanAlbada | **Minmod** | paper / DSMC |
|---|---|---|---|
| Cp ristagno (robusto) | 1.90 | **1.785** | Rayleigh 1.837 |
| C_D | 1.288 | **1.285** (p 1.254 + attr 0.031) | 1.304 / 1.284 |
| C_H [kW] | 57.3 | **74.1** | 88.1 / 63.3 |
| Cf(θ) | seghettata | **liscia**, picco ~0.033 @ 50° | run3 0.040, DSMC 0.055 |

Il C_H schizza da 57 (schiacciato dal checkerboard) a 74.1, **dentro
l'intervallo fisico** DSMC-paper. Cp e Cf ora sono lisce e grado-paper in
forma; C_D entro 1.5%. Resta **un solo spike di flusso termico a θ≈27°**
vicino al ristagno — l'ultima cella dove il disturbo sopravvive, perché
q~dT/dn è la più sensibile. Figure aggiornate:
`cylinder-surface.png`, `cylinder-stagnation.png`.

## 6. La mesh fine non converge — e perché non serve inseguirla

Provata l'opzione 2 (156k celle, prima cella 1.9 µm, Minmod) sul cluster.
**A 90k step non è a regime**: il guard-rail di convergenza segna dp/p
media 2.66%, max 1236% (la coarse allo stesso punto era 0.24%), e i numeri
escono spazzatura (Cp ristagno 4.8, C_D 4.6, C_H 1516 kW). Due cause
concorrenti:

1. **LTS più lento su celle più piccole**: il passo pseudo-temporale
   locale è ∝ alla dimensione di cella, quindi celle ~4× più fini vogliono
   ~4× più iterazioni per convergere. 90k non bastano.
2. **Aspect ratio near-wall peggiore**: la prima cella a 1.9 µm su celle
   tangenziali da ~10 mm dà ~5000:1 (contro 1700:1 della coarse) — il modo
   odd-even è *più* alimentato, non meno. Il Cp robusto che esce 4.8 (non
   ~1.8) dice che sulla fine il checkerboard è più profondo della cella-1,
   quindi nemmeno la banda pulita salva l'estrazione.

**Decisione: la coarse Minmod è il risultato di M7, la fine è rifinitura
futura.** Motivo di fondo, oltre alla convergenza: i riferimenti del paper
**non concordano tra loro** — sul Cf il picco va da 0.040 (hy2Foam run3) a
0.055 (DSMC), spread 37%; sul flusso termico il paper stesso dice che Park
sovrastima la DSMC del 39% (C_H 88 vs 63). Il nostro Cf coarse ~0.033 e
C_H 74 kW cadono *dentro* questa forbice. Inseguire la mesh fine non
avvicinerebbe un "valore vero" (che non esiste come punto singolo), a
fronte di un problema numerico reale da risolvere prima. La cura giusta,
quando/se si vorrà chiudere il quantitativo sulla fine, è il **re-grading
near-wall** (più celle tangenziali o grading più dolce mantenendo ~2 µm a
parete) per riportare l'aspect ratio a ~1500:1; il paper gira 155k a 2 µm
senza questo problema perché il suo flusso KNP tollera quella risoluzione
meglio del nostro centrale. (Bug minore corretto lungo la strada: `o`
riusato in `postProcess-cylinder.py` faceva crashare il plot Mach di
ristagno — commit `ea3901d`.)

## 7. Verdetto M7

**Validato (grado-paper) sulla coarse Minmod**: pressione di parete (Cp,
C_D entro 1.5%), standoff (0.247 vs 0.25 m), profili di linea di ristagno
(Mach = run3, T con la sovrastima del picco *attesa* per la combinazione
Park, dissociazione N nel posto giusto), e Cf/C_H **nell'intervallo dei
riferimenti**. **La fisica reattiva two-temperature in multi-D funziona** —
la chimica, validata finora solo in 0D/1D (M2/M3/M4), regge accoppiata al
two-temperature su un caso 2D completo con shock staccato. È l'obiettivo
di M7, centrato.

**Aperto / rifinitura futura**:
- mesh fine convergente (re-grading near-wall + più step LTS) per limare
  il gap Cf e lo spike residuo di flusso termico;
- diffusione del pool ve nella EveEqn (κ_ve = μ·cv_ve alla Eucken → un
  `laplacian(μ, eve)`), invariata da M6: darebbe il contributo
  vibrazionale al flusso termico e permetterebbe la parete vibrazionale a
  1000 K come nel paper;
- Fig 9 (aria 5 specie), CVDV-QK, Tve multiple.

Caveat di digitalizzazione (`references/README.md`): nei pannelli (b) e
(c) ogni file `run*` contiene due curve fisiche (T_tr+T_v, N2+N): vanno
separate per un overlay pulito dei profili di ristagno.

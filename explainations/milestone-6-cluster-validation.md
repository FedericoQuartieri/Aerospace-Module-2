# Milestone 6 — Campagna cluster e risoluzione del giallo del Cp

**Data:** 9-10 luglio 2026
**Branch:** `milestone-6-cluster-validation`, che si innesta sulla coda di `milestone-5-blunted-cone`. Storia ripulita con un rebase interattivo: gli **8 commit** della trafila di debug del cluster (compreso un vicolo cieco — lo spostamento di `set -e`, poi superato) sono stati **fusi in un unico commit** `5dd942a` *"cluster: PBS job setup for the Mach 11.3 cone"*, che questo documento spiega nel §1. Restano separati, perché sono cose diverse: `05e418d` (il fix vero al clamp `minTemperature`, §3-4) e i commit M5 precedenti.
**Obiettivo:** portare il caso del cono sul cluster PBS di Federico (4 job × 28 core), eseguire il run fine che la M5 aveva definito come *test diagnostico* per il Cp di ristagno basso, e chiudere il giallo.

---

## 1. La trafila del cluster (tutto ciò che NON era ovvio)

Portare il codice su un cluster mai usato prima ha richiesto sei fix successivi, ognuno scoperto empiricamente. In ordine, perché serviranno per ogni cluster futuro:

1. **Sintassi PBS**: il cluster usa PBS Professional/OpenPBS, non Torque classico. La sintassi legacy `-l nodes=1:ppn=28` viene rifiutata con `Job violates queue and/or server resource limits`; serve `-l select=1:ncpus=28:mpiprocs=28`. Inoltre la coda va indicata esplicitamente (`-q cpu`: max 28 ncpus, walltime 48h, `max_user_run=4` — verificato con `qstat -Qf`); senza, il job viene instradato su code con limiti incompatibili (es. `scalability`, walltime max 30 min).
2. **Lo stage-out dell'output PBS non funziona**: il file `-o` non compare mai, né con path relativo né assoluto, nemmeno per un job banale `hostname; sleep`. E il cluster non tiene job history (`qstat -xf` inutilizzabile). Soluzione: ogni job scrive il proprio log con `exec > "$LOG" 2>&1` come prima riga eseguibile — bypass completo di PBS.
3. **OpenFOAM-13 è preinstallato** in `/opt/openfoam13` ma nessun modulefile lo espone (`module avail` vuoto) e il `bashrc` del nostro repo NON lo sorgenta: presuppone un ambiente OpenFOAM già attivo (in locale è nel `~/.bashrc` personale). Senza sourcing esplicito, `foamRun` non è in PATH e i job muoiono all'istante — che combinato col punto 2 (nessun log) dava il sintomo iniziale: job in stato `E` a 00:00, silenzio totale. Diagnosi fatta con una sessione interattiva (`qsub -I -q interactive`), l'unico modo per vedere output live.
4. **Bug bash `pop_var_context`**: nei job batch (ma non in sessione interattiva!) il sourcing diretto del bashrc di OpenFOAM termina la shell con `pop_var_context: head of shell_variables not a function context`. Causa esatta mai individuata (SHELLOPTS ereditato escluso: non contiene errexit; lo spostamento di `set -e` dopo i source NON è bastato). Workaround robusto: sourcing in un sottoprocesso pulito — `bash -c "source ... && export -p" > dump` poi `source dump` — si buttano le funzioni shell interne di OpenFOAM (inutili allo script) e si tengono solo le variabili d'ambiente.
5. **Mutation++ da GitHub, non da rsync**: `thirdParty/Mutationpp` è gitignorata, ma il checkout locale è un clone pulito del mirror pubblico (`github.com/mutationpp/Mutationpp`, commit `117df0b`, zero patch — verificato). Quindi sul cluster si clona direttamente (dal **login node**: i nodi di calcolo potrebbero non avere internet) e si pinna il commit. Niente trasferimento da ~1 GB.
6. **`setup-cluster.sh`**: job una tantum (coda cpu, 3h) che compila Mutation++ + librerie custom del progetto, con verifica dei prerequisiti e log auto-gestito. Completato con successo al primo colpo dopo i fix 1-5.

Bonus di processo: il walltime del job fine è stato portato da 12h (stima teorica) a 24h, ancorato al dato reale del run coarse (4548 s per 7200 celle × 4 core → ~3h proiettate per 120k × 28, con margine per lo scaling mai verificato oltre 4 core).

## 2. Il run fine: il test diagnostico dà il suo verdetto

La M5 si era chiusa con il Cp di ristagno coarse a **1.304** (29% sotto il Rayleigh-Pitot ideale 1.833) e due diagnosi già smentite o in verifica:
- *non-ortogonalità di mesh* → *smentita* (le facce critiche sono lontane dal ristagno, e la mesh fine ha identica non-ortogonalità);
- *sotto-risoluzione dello shock* → il run fine doveva essere il giudice: 5× risoluzione nello strato d'urto, "se il Cp sale verso 1.8 era risoluzione".

**Verdetto del run fine (120k celle, 28 core): Cp = 1.335.** Invariato al 2% con 17× più celle. **Anche la seconda diagnosi era sbagliata**: errore sistematico, indipendente dalla discretizzazione.

## 3. La causa vera: il clamp che riscaldava il free-stream

Con la mesh scagionata, la rianalisi dei campi convergiuti (coarse, locale) ha mostrato il colpevole in una riga di dict:

**Il bridge Mutation++ ha un clamp di sicurezza `minTemperature 200` (default nato nelle milestone 0D, dove tutto era caldo), ma il free-stream del paper è a T∞ = 144.4 K.** Il clamp riscaldava silenziosamente tutto il flusso indisturbato a 200.00 K esatti — verificato: `T upstream = 200.00 K` in ogni cella a monte dello shock. Conseguenze a catena:

- T∞ effettiva 200 K → ρ∞ effettiva 3.69×10⁻⁴ (non 5.113×10⁻⁴) → **Mach effettivo 9.59, non 11.3**
- il postProcess confrontava un flusso a M9.6 con la teoria a M11.3 → "29% di errore"

**La prova che il solver era esatto da sempre**: rifacendo il confronto alle condizioni *effettivamente simulate* (M=9.59), il Rayleigh-Pitot dà Cp=1.831 e il CFD misurava 1.806 — **accordo all'1.4%**. Tre diagnosi (mesh, schema, clamp) per arrivare a scoprire che non c'era nessun errore da spiegare: solo un flusso diverso da quello richiesto.

## 4. Fix e conferma

Fix: `minTemperature 50` nel `physicalProperties` del caso (il valore è configurabile da dict, nessuna modifica al codice). Rerun coarse locale (179 min su 4 core — più lento dei run precedenti: il free-stream a 144 K ha scale temporali più piccole):

| Grandezza | Prima (clamp attivo) | Dopo il fix | Teorico |
|---|---|---|---|
| T∞ effettiva | 200.0 K | **144.4 K** | 144.4 K |
| Mach | 9.59 | **11.29** | 11.3 |
| Cp di ristagno | 1.304 | **1.802** | 1.833 (**errore 1.7%**) |
| Standoff | 0.159 Rn | **0.147 Rn** | 0.1–0.15 Rn (centrato) |

## 5. Il guard-rail permanente

`postProcess-cone.py` ora **misura il free-stream effettivamente simulato** (T, p, M mediati nelle celle a monte dello shock), lo stampa accanto ai valori nominali, **avvisa esplicitamente** se differiscono oltre 2 K, e valuta ogni confronto teorico (Rayleigh-Pitot, normalizzazioni di Cp) **al Mach effettivo**, mai a quello nominale.

**La lezione di metodo, in una riga**: mai confrontare col nominale senza verificare cosa è stato davvero simulato. Tre diagnosi consecutive sbagliate — tutte tecnicamente sofisticate (non-ortogonalità, dissipazione dello schema) — perché l'assunzione più banale ("il free-stream è quello imposto") non era mai stata controllata. Il segnale c'era già nel primissimo dump cella-per-cella (T upstream = 200.5 K, visibile a occhio), ma nessuno lo stava guardando.

## 6. Il run fine definitivo: validazione a grado-paper

Eseguito sul cluster (11 lug 2026, coda `cpu`, 28 core, mesh 120k celle con prima cella 2.3 µm) dal branch M6 con il fix `minTemperature 50`. Il guard-rail conferma per primo che stavolta il flusso simulato è quello giusto:

```
free-stream EFFETTIVO: T=144.4 K (nominale 144.4), p=21.91 Pa, M=11.29 (nominale 11.3)
```

| Grandezza | Coarse (fix) | **Fine (fix)** | Teorico |
|---|---|---|---|
| Cp di ristagno | 1.802 | **1.843** | 1.833 (**errore 0.5%**) |
| Standoff | 0.147 Rn | **0.150 Rn** | 0.1–0.15 Rn |

La progressione completa del Cp racconta l'intera vicenda: **1.30 (clamp attivo, −29%) → 1.802 (coarse col fix, −1.7%) → 1.843 (fine col fix, +0.5%)**. Il leggero eccesso rispetto al Rayleigh-Pitot ideale è fisicamente atteso: la formula assume γ=1.4 esatto, ma a ~3700 K post-shock la vibrazione dell'N₂ inizia a eccitarsi (γ effettivo appena sotto 1.4) → compressione al ristagno leggermente maggiore. La differenza coarse→fine (~2%) è la normale convergenza di griglia sulla cattura dello shock.

**La pressione di parete al punto di ristagno è validata a grado-paper.**

## 7. I riferimenti della Fig 2: estrazione vettoriale, non ricalco manuale

Per il confronto quantitativo sui profili serve la Fig 2 del paper in forma numerica. Invece del pattern Engauge (ricalco manuale sul raster, già usato per la Fig 5 in M2), qui si è sfruttato il fatto che le figure MDPI sono **grafica vettoriale**: le curve esistono nel PDF come comandi di disegno con coordinate esatte. `references/digitize-fig2.py` le estrae con PyMuPDF (`get_drawings`), calibra gli assi dai sei riquadri dei pannelli (range noti dalle etichette) e classifica le serie per proprietà geometriche del tracciato:

- **colore blu** → hy2Foam prima cella 10 µm; **path neri connessi e lunghi** (contiguità >80%) → hy2Foam solida; **tratti neri staccati** → CFD Michigan (dash-dot); **simboli pieni** → triangoli DSMC (centroide); **piccoli path gambo+cap** → barre d'errore CUBRC (punto = centro barra).

L'accuratezza è limitata solo dalla calibrazione degli assi. Le trappole trovate (documentate nel README di `references/`): i **tick degli assi** vanno filtrati per segmento, non per path (gnuplot emette i tick speculari di due bordi opposti in un unico path), e la distinzione solida/tratteggiata va fatta sulla **contiguità dei segmenti**, non sulla lunghezza. QA obbligatorio: `fig2-overlay.png` sovrappone i punti estratti al render della pagina. Risultato: 21 CSV (`references/fig2{a-f}-{serie}.csv`).

## 8. Primo confronto quantitativo (mesh coarse)

`compare-fig2.py` calcola le grandezze di parete del paper (eq. 30–32: Cp, Cf, St normalizzati sul free-stream nominale) leggendo `p` nelle celle owner della parete e i campi generati da `foamPostProcess -solver shockThermo -func wallHeatFlux` / `-func wallShearStress` (il flag `-solver` è necessario: il foamPostProcess liscio non costruisce il thermophysicalTransport). Sul **coarse** locale (prima cella 24.5 µm), scarti rispetto alla hy2Foam digitalizzata su 0.5–4 cm:

| Grandezza | Scarto medio | Scarto max | Lettura |
|---|---|---|---|
| Cp | **2.0%** | 3.7% | sovrapposto a paper/DSMC/esperimenti |
| Cf | 15.9% | 18.1% | sotto le curve paper sul fianco: attesa, è la grandezza più sensibile alla prima cella (24.5 µm vs 2–10 µm del paper) |
| St | 10.6% | 14.5% | forma giusta; nel plateau più vicino agli **esperimenti CUBRC** che alle CFD del paper |

Sulla linea di ristagno (Fig 2a) il profilo T/T∞ segue la hy2Foam (picco ~25.9 vs 25.8) e la Tve ricalca i punti Michigan; lo shock DSMC è più diffuso, com'è fisicamente atteso a Kn~0.002. Figure: `fig2-comparison.png`, `fig2-stagnation-comparison.png` nel caso.

## 9. Il confronto fine: Cp chiuso, e il gap Cf/St smascherato

Confronto rifatto sui campi del run fine (cluster, `job-post-fine.sh`, prima cella 2.3 µm — comparabile ai 2–10 µm del paper):

| Grandezza | Coarse (24.5 µm) | **Fine (2.3 µm)** | Lettura |
|---|---|---|---|
| Cp | 2.0% / 3.7% | **0.7% / 1.7%** | convergenza di griglia confermata, sovrapposto a paper ed esperimenti |
| Cf | 15.9% | **16.3%** | **invariato con 10× di risoluzione a parete** |
| St | 10.6% | **11.6%** | invariato; nel plateau sta *tra* le CFD del paper e gli esperimenti CUBRC (che sono ~15% sotto le CFD) |

L'invarianza di Cf/St rispetto alla prima cella scagiona definitivamente la mesh: è un errore di **modello di trasporto**. Il caso usa `transport sutherland` con coefficienti fittati sull'intervallo caldo (1000–5000 K: entro ±10% dalla Blottner usata da hy2Foam, eccellente a 5000 K), ma un fit a 2 parametri non può coprire 144–6000 K: a **297 K — la temperatura di parete, dove si valutano τ_w = μ∂u/∂y e q_w = κ∂T/∂y — la Sutherland dà μ il 21% più basso della Blottner** (e −47% a T∞=144 K). Il deficit del 16%/12% su Cf/St è tutto lì; il Cp non dipende da μ e infatti converge perfettamente.

Morale gemella di quella del clamp: anche stavolta la grandezza sbagliata era "di contorno" (prima il floor di temperatura, ora il fit di viscosità), non la fisica two-temperature. Fix possibile: trasporto Blottner per specie + Eucken per κ (la miscela Wilke c'è già: `coefficientWilkeMulticomponentMixture`) — è il modello del paper. Figure fine: `fig2-comparison-fine.png`, `fig2-stagnation-comparison-fine.png`.

## 10. Il trasporto Blottner: il gap si chiude

Implementato `blottnerTransport` in `src/thermophysicalModels/specie/transport/blottner/` — il modello di viscosità del paper (hy2Foam/LeMANS):

    mu = 0.1*exp((A*lnT + B)*lnT + C)      [Blottner 1971]
    kappa = mu*Cv*(1.32 + 1.77*R/Cv)       [Eucken modificata, come sutherland]

La classe ricalca `sutherlandTransport` (stessa API: `mu(p,T)`, `kappa(p,T)`, operatori di miscela pesati in massa). **Non serve registrarla in nessuna type table statica**: il thermo del caso è compilato al volo via dynamicCode — il template stock istanzia `${transport}Transport<...>` per nome, quindi basta (a) l'header nell'include path (lnInclude della specie POLIMI, già davanti a quello stock) e (b) il nome `blottner` nella whitelist `etc/codeTemplates/dynamicCode/fluidMulticomponentThermo` (stesso meccanismo con cui era stato aggiunto `rrho`). Gotcha noto: la whitelist è **cacheata** in `~/.OpenFOAM/13` — la copia avviene al source di `etc/bashrc`, che i job cluster fanno già a ogni run.

Coefficienti Blottner (tab. Gupta 1990) aggiunti alle 5 specie attive in `speciesThermo.janaf` (chiavi `A`/`B`/`C` accanto ad `As`/`Ts`, che restano per tornare indietro) e `transport blottner` nel `physicalProperties` del cono. Risultato sul **coarse** (stessa mesh, stesso tutto, solo μ(T) diverso):

| Grandezza | Sutherland | **Blottner** | Verifica incrociata |
|---|---|---|---|
| Cp | 2.0% / 3.7% | **2.0% / 3.7%** | invariato al centesimo, e Cp ristagno 1.802 identico: il trasporto non tocca la pressione — se il thermo fosse sbagliato si vedrebbe qui |
| Cf | 15.9% | **7.3%** | plateau 0.041→0.048 (+17% netto da +27% di μ a parete: lo strato limite si ispessisce e ∂u/∂y cala, risposta fisica corretta) |
| St | 10.6% | **1.8%** | chiuso |

Il residuo Cf ~7% sulla coarse ha due candidati: la prima cella (24.5 µm vs 2.3 µm del fine — ora che il modello è giusto la risoluzione può tornare a contare) e lo spessore della linea tratteggiata digitalizzata. Verdetto al rerun fine sul cluster.

## 11. Stato e prossimi passi (scope M6)

- **Fatto**: cluster operativo end-to-end, causa del Cp trovata e corretta, guard-rail nel postProcess, run fine validato (Cp ristagno 0.5%), riferimenti Fig 2 digitalizzati per via vettoriale, confronto Fig 2 coarse+fine, gap Cf/St diagnosticato (fit Sutherland a parete) **e chiuso con il trasporto Blottner (coarse: Cf 16→7%, St 11→2%)**.
- **Prossimo passo**: rerun fine sul cluster con Blottner (`qsub job-cone-fine.sh` poi `qsub job-post-fine.sh` dopo `git pull`) per il verdetto a grado-paper; poi cilindro Mach 20 reagente (secondo caso del paper Part Two).
- Rimandati noti invariati: diffusione del pool ve nella EveEqn (κ_ve), Fig 9, CVDV-QK, Tve multiple.

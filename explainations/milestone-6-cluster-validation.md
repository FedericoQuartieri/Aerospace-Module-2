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

## 6. Stato e prossimi passi (scope M6)

- **Fatto**: cluster operativo end-to-end (setup + job + log), run fine diagnostico, causa del Cp trovata e corretta, coarse validato all'1.7% dal teorico, guard-rail nel postProcess.
- **Prossimo passo immediato**: run fine definitivo sul cluster alle condizioni corrette — `git pull` + `qsub job-cone-fine.sh`; nel log verificare la riga `free-stream EFFETTIVO: T=144.4 K ... M=11.29` e attendersi Cp ≈ 1.80–1.83.
- **Poi**: digitalizzazione dei riferimenti del paper (Fig 2: linea di ristagno Wang&Boyd/MONACO, Cp/Cf/St con esperimenti CUBRC run 31) per il confronto quantitativo completo — pattern Engauge già usato per la Fig 5.
- **Più avanti nello scope M6 o successivo**: cilindro Mach 20 reagente (secondo caso del paper Part Two); rimandati noti invariati (diffusione del pool ve nella EveEqn, Fig 9, CVDV-QK, Tve multiple).

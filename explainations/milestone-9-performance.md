# Milestone 9 — Analisi di performance (scaling MPI del solver multi-D)

Milestone di **rifinitura opzionale**: la validazione fisica del paper
(M1–M8) è già completa. M9 copre l'asse "High Performance" del corso — che
finora avevamo usato (MPI sul cluster) ma mai *misurato* in modo
sistematico.

## 1. Perché la nostra storia HPC è diversa (e più solida)

Un altro gruppo del corso ha fatto la sua analisi di performance su un
**kernel 0D parallelizzato in OpenMP**: N celle *indipendenti* di heat bath,
ognuna che evolve da sola, distribuite sui thread (embarrassingly
parallel, nessuno scambio dati). Scaling misurato: 2.81× su 8 thread
(efficienza 0.35) su un laptop, con nCells=256.

Noi misuriamo una cosa qualitativamente diversa: lo **strong scaling MPI
del solver 2D reagente vero** (`shockThermo`), con:
- **decomposizione di dominio** reale (scambio di halo tra i rank a ogni
  step, non celle indipendenti);
- **accoppiamento fisico** completo (flusso, chimica, due temperature)
  risolto in modo distribuito;
- il **bridge Mutation++ per cella**, una `Mutation::Mixture` per rank.

È distributed-memory su un problema accoppiato — HPC reale, non repliche di
un kernel 0D.

## 2. Il prerequisito di correttezza (già in M5)

L'analogo del loro "thread-local ownership model" (una Mixture per thread
per la thread-safety) da noi è: **una `Mutation::Mixture` per rank MPI**,
verificata in M5 — `relaxing` seriale vs `mpirun -np 4` dà campi
**identici a ~1e-5**. Il parallelismo non cambia il risultato fisico:
prerequisito indispensabile prima di parlare di velocità.

## 3. Metodo di misura

`perf/run-scaling.sh` gira la stessa mesh su un elenco di rank su una
**finestra di pseudo-tempo fissa** dall'avvio impulsivo, e registra
`ExecutionTime` e il numero di step. Con LTS il passo locale è per-cella
(indipendente dalla decomposizione), quindi il numero di step per
raggiungere `endTime` è ~costante fra i rank → `ExecutionTime` è una
metrica di strong scaling pulita. Normalizziamo per step
(`tstep = ExecutionTime/nSteps`) per togliere ogni residuo jitter di
conteggio. Speedup S(n)=t(1)/t(n), efficienza E=S/n.

## 4. Costo per cella del bridge (dato già raccolto)

Il costo dominante è la chiamata Mutation++ per cella. Misure lungo il
progetto:

| Caso | µs/cella/step | Note |
|---|---|---|
| 1D shockTube reagente (M4) | 34–42 | rilassamento N₂, laptop |
| 2D cono M11 (M7 smoke, cluster) | ~55 | non-reagente (solo thermo) |
| 2D cilindro M20 (M7 smoke, cluster) | ~84 | reagente = ~1.5× il non-reagente |
| 2D cilindro M20 (M9 locale, np=1) | 60.5 | reagente, laptop, avvio impulsivo |

La chimica reagente aggiunge ~50% di costo per cella rispetto al solo
thermo — coerente col fatto che Mutation++ valuta anche le velocità di
reazione a due temperature.

## 5. Strong scaling — preliminare locale

Harness validato in locale (laptop 8 core, cilindro coarse 9000 celle):

| ranks | ms/step | speedup | efficienza |
|---|---|---|---|
| 1 | 544.5 | 1.00 | 1.00 |
| 2 | 392.3 | 1.39 | 0.69 |
| 4 | 295.1 | 1.84 | 0.46 |

Sub-lineare: 9000 celle su 4 rank fa ~2250 celle/rank, e su un laptop il
collo di bottiglia è la banda di memoria condivisa più che la comunicazione.
È il regime "caso piccolo" — utile a validare l'harness, non
rappresentativo di un nodo HPC.

## 6. Strong scaling — sweep rappresentativo sul cluster (da lanciare)

`perf/job-scaling.sh` fa lo sweep su **1 nodo da 28 core**, rank
{1,2,4,7,14,28}, su **due mesh**:
- **coarse** (9k celle): a 28 rank ~320 celle/rank → *communication-bound*,
  lo scaling satura — lo stesso regime del kernel 0D del gruppo 1, ma su un
  solver accoppiato vero;
- **fine** (156k celle): ~5600 celle/rank a 28 → *compute-bound*, scala
  bene.

Sovrapporre le due curve mostra il compromesso compute/comunicazione —
un'analisi di strong scaling completa. Sul cluster:

```
qsub perf/job-scaling.sh          # ~1-2h su 28 core
python3 perf/plot-scaling.py      # tabella + perf/scaling.png (coarse+fine)
```

_(Risultati cluster da inserire qui dopo il run: tabella coarse/fine +
`perf/scaling.png`.)_

## 7. Cosa manca ancora (rifinitura, non obbligatorio)

- weak scaling (celle/rank costante) — richiede mesh scalate coi rank;
- multi-nodo (oltre le 28 core di un nodo) per vedere il costo
  dell'interconnessione;
- profiling fine (frazione del tempo in Mutation++ vs flusso vs solutore
  lineare) con gprof/perf.

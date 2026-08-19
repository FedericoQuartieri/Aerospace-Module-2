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

## 6. Strong scaling sul cluster (1 nodo, 28 core) — RISULTATI

`perf/job-scaling.sh` ha fatto lo sweep su rank {1,2,4,7,14,28} su due mesh.
Figura: `perf/scaling.png`.

**Coarse (9000 celle)** — a 28 rank fa ~320 celle/rank:

| ranks | ms/step | speedup | efficienza |
|---|---|---|---|
| 1 | 383.8 | 1.00 | 1.00 |
| 2 | 200.0 | 1.92 | 0.96 |
| 4 | 115.9 | 3.31 | 0.83 |
| 7 | 85.0 | 4.51 | 0.64 |
| 14 | 54.0 | 7.11 | 0.51 |
| 28 | 48.0 | **8.00** | **0.29** |

**Fine (156000 celle)** — a 28 rank fa ~5600 celle/rank:

| ranks | ms/step | speedup | efficienza |
|---|---|---|---|
| 1 | 11257 | 1.00 | 1.00 |
| 2 | 5698 | 1.98 | 0.99 |
| 4 | 2949 | 3.82 | 0.95 |
| 7 | 1735 | 6.49 | 0.93 |
| 14 | 920 | 12.24 | 0.87 |
| 28 | 566 | **19.90** | **0.71** |

**Il compromesso compute/comunicazione, da manuale.** La coarse segue
l'ideale fino a ~4 rank poi piega e satura (8× a 28, efficienza 0.29):
troppo poche celle per rank, il costo di scambio degli halo domina — lo
stesso regime del kernel 0D di Group 1, ma qui su un solver accoppiato
vero. La fine invece resta **quasi lineare fino a 14 rank** (efficienza
0.87) e arriva a **19.9× su 28 core (efficienza 0.71)**: con ~5600 celle
per rank il costo dominante è il calcolo per cella (bridge Mutation++), non
la comunicazione → compute-bound, scala.

**Confronto diretto con Group 1**: loro OpenMP su kernel 0D di celle
indipendenti → **2.81× su 8 thread** (efficienza 0.35, satura). Noi MPI su
solver 2D reagente accoppiato → **19.9× su 28 rank** (efficienza 0.71),
near-linear fino a 14. È scaling distributed-memory reale su un problema
con scambio di halo, non repliche di un kernel embarrassingly-parallel.

## 7. Verdetto M9

L'asse "High Performance" del corso è coperto, e con una storia più forte:
strong scaling MPI di un solver reagente multi-D accoppiato, correttezza
seriale-vs-parallelo già garantita (M5, 1e-5), efficienza 0.71 su 28 core
sul caso rappresentativo (fine 156k). Il progetto è ora completo su
entrambi gli assi — **validazione fisica del paper (M1–M8) + performance
HPC (M9)**.

## 8. Cosa manca ancora (rifinitura, non obbligatorio)

- weak scaling (celle/rank costante) — richiede mesh scalate coi rank;
- multi-nodo (oltre le 28 core di un nodo) per vedere il costo
  dell'interconnessione;
- profiling fine (frazione del tempo in Mutation++ vs flusso vs solutore
  lineare) con gprof/perf.

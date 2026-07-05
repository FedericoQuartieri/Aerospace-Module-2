# Run sul cluster — blunted cone Mach 11.3

Il caso fine (600×200, ~120k celle, prima cella a parete ~2 µm come il
paper) è pensato per un job da **28 core**. Il bridge Mutation++ è
per-cella e gira sotto MPI con la decomposizione di dominio standard di
OpenFOAM (verificato: run seriale e parallelo coincidono a ~1e-5
relativo sul caso `shockTube1D/relaxing`).

## Setup una tantum sul cluster

Verificato passo-passo con una sessione interattiva (`qsub -I -q
interactive -l ncpus=1 -l walltime=00:05:00`) sul cluster di Federico
il 9 luglio 2026:

- **OpenFOAM-13 è preinstallato** in `/opt/openfoam13` (nessuna build
  necessaria per OpenFOAM stesso). Il `bashrc` del repo *non* lo
  sorgenta da solo: presuppone un ambiente OpenFOAM già attivo (in
  locale Federico ce l'ha nel suo `~/.bashrc` personale). Sul cluster
  va sorgentato esplicitamente **prima** — omesso, `foamRun` non è in
  `PATH` e i job batch muoiono all'istante senza log leggibile (è
  esattamente il sintomo diagnosticato: job in stato `E` a 00:00, nessun
  file di output, riprodotto anche con un job banale `hostname; sleep`).
- `gcc 13.3.0`, `OpenMPI 4.1.6`, `cmake 3.28.3` disponibili di sistema.
- **Nessuna coda di modulefile per OpenFOAM** (`module avail` non lo
  elenca) — il path fisso `/opt/openfoam13` è la fonte di verità.

Sequenza:

```sh
# 1. repo (già fatto: /work/u10806848/Aerospace-Module-2)

# 2. Mutation++ vendorizzata: thirdParty/ è gitignorata. È un mirror
#    pubblico su GitHub (github.com/mutationpp/Mutationpp), nessuna
#    patch locale (verificato con `git status`/`git remote -v` sul
#    checkout del laptop) — si può clonare direttamente, niente
#    trasferimento a mano di ~1 GB. Farlo QUI, sul nodo di login (i
#    nodi di calcolo spesso non hanno accesso internet):
cd ~/Aerospace-Module-2/thirdParty
git clone https://github.com/mutationpp/Mutationpp Mutationpp
cd Mutationpp
git checkout 117df0b14baf41715c9cf04be9fc7f438fa4513f   # commit verificato in dev
cd ../..

# 3. build di Mutation++ + librerie custom del progetto: UN SOLO job,
#    non richiede sessione interattiva (walltime fino a 3h, coda cpu)
cd applications/test/bluntedCone2D/cluster
qsub setup-cluster.sh
qstat -u $USER                              # attendere lo stato C
cat setup-cluster.log                       # deve terminare con
                                             # "=== setup completato ==="
```

`setup-cluster.sh` sorgenta `/opt/openfoam13/etc/bashrc` + il bashrc
del progetto, poi lancia `thirdParty/makeMutationpp`,
`src/thermophysicalModels/Allwmake` e `applications/modules/Allwmake`
in sequenza. Se il passo 2 è già stato fatto (come sopra), lo script lo
rileva e salta il clone; se non lo è e il nodo di calcolo ha comunque
accesso internet prova a clonare da solo — ma è più affidabile farlo a
mano prima, come sopra.

**Nota dynamicCode**: il thermo `highEnthalpyThermo` viene compilato
just-in-time nella directory del caso al primo run (serve il compilatore
sui nodi, o almeno sul nodo di lancio). Se i nodi di calcolo non hanno
g++, fare un run brevissimo sul nodo di login per generare
`dynamicCode/` e poi lanciare il job.

## Lancio

```sh
cd applications/test/bluntedCone2D/cluster
qsub job-cone-fine.sh
qstat -u $USER                 # monitoraggio
```

Configurazione verificata con `qstat -Qf` sul cluster di Federico
(coda `cpu`: `resources_max.ncpus=28`, walltime max 48h, `max_user_run=4`).
Il job usa la sintassi moderna `-l select=1:ncpus=28:mpiprocs=28` (PBS
Professional/OpenPBS) invece della vecchia `-l nodes=N:ppn=M` di Torque
classico, che su questo cluster non si traduce correttamente e dà
`qsub: Job violates queue and/or server resource limits` — probabile
causa: senza `-q cpu` esplicito il job viene instradato altrove (es.
`scalability`, che ha walltime max 30 minuti, molto meno dei 12h
richiesti). Se il tuo cluster usa un altro PBS, ricontrollare `qstat -Qf`
e adattare `select=...`/`-q` di conseguenza.

Se `mpirun` non vede i 28 core assegnati dal job (dipende da come è
compilato OpenMPI sul cluster), usare la machinefile di PBS:
`mpirun -np 28 -machinefile $PBS_NODEFILE foamRun -parallel`
(nel qual caso conviene copiare la riga dentro l'Allrun o lanciare i
passi a mano).

Stima aggiornata con dati reali (non più teorica): il run coarse
(7200 celle, 4 core) ha impiegato 4548 s (ExecutionTime OpenFOAM) per
arrivare a convergenza reale — verificato confrontando i campi tra due
istanti separati, non solo per timeout. Scalando linearmente a 120k
celle / 28 core: ~4548×(120000/7200)×(4/28) ≈ 10800 s ≈ 3h. Walltime
del job impostato a 24h per margine ampio, dato che lo scaling
parallelo a 28 core non è mai stato verificato (solo fino a 4 core) e
resta ben dentro il tetto di 48h della coda. Il paper: 2.8h su 24 core
per il caso equivalente — coerente con la stima.

## Cosa riportare indietro

- `cone-stagnation.png`, `cone-surface.png` e l'output testuale del
  postProcess (standoff, Cp di ristagno)
- `log.foamRun` (per i tempi e la convergenza)
- la directory dell'ultimo time step se si vuole rifare il post in locale

I 4 job disponibili si possono usare per varianti (mesh, maxCo,
accommodation) in parallelo — ogni caso è una copia della directory.

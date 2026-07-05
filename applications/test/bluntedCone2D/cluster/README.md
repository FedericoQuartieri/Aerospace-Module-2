# Run sul cluster — blunted cone Mach 11.3

Il caso fine (600×200, ~120k celle, prima cella a parete ~2 µm come il
paper) è pensato per un job da **28 core**. Il bridge Mutation++ è
per-cella e gira sotto MPI con la decomposizione di dominio standard di
OpenFOAM (verificato: run seriale e parallelo coincidono a ~1e-5
relativo sul caso `shockTube1D/relaxing`).

## Setup una tantum sul cluster

```sh
# 1. repo (niente push: copiare/clonare come preferisci)
git clone <repo> ~/Aerospace-Module-2
cd ~/Aerospace-Module-2

# 2. OpenFOAM-13 (module load openfoam/13 se disponibile, altrimenti build)

# 3. Mutation++ vendorizzata (thirdParty/ è gitignorata: copiare la
#    directory thirdParty/Mutationpp dalla macchina locale, poi:)
cd thirdParty && ./makeMutationpp && cd ..

# 4. build della catena
source etc/bashrc
src/thermophysicalModels/Allwmake
applications/modules/Allwmake
applications/test/nonEqTTv/Allwmake     # opzionale, per i test 0D
```

**Nota importante — thirdParty gitignorata**: `thirdParty/Mutationpp`
(sorgenti + install) non è nel repo. Va copiata a mano sul cluster
(rsync della directory) prima del punto 3, oppure ri-scaricata e
compilata con lo stesso commit/versione.

**Nota dynamicCode**: il thermo `highEnthalpyThermo` viene compilato
just-in-time nella directory del caso al primo run (serve il compilatore
sui nodi, o almeno sul nodo di lancio). Se i nodi di calcolo non hanno
g++, fare un run brevissimo sul nodo di login per generare
`dynamicCode/` e poi lanciare il job.

## Lancio

```sh
cd applications/test/bluntedCone2D/cluster
sbatch job-cone-fine.sh        # adattare partition/account/module load
```

Stima: ~120k celle × 2 correct × ~40 µs / 28 core ≈ 0.35 s/iterazione
→ 10000 iterazioni ≈ 1 h; con margine LTS/convergenza budget 6–12 h
(il paper: 2.8 h su 24 core per il caso equivalente).

## Cosa riportare indietro

- `cone-stagnation.png`, `cone-surface.png` e l'output testuale del
  postProcess (standoff, Cp di ristagno)
- `log.foamRun` (per i tempi e la convergenza)
- la directory dell'ultimo time step se si vuole rifare il post in locale

I 4 job disponibili si possono usare per varianti (mesh, maxCo,
accommodation) in parallelo — ogni caso è una copia della directory.

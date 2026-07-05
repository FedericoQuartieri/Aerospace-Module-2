#!/bin/bash
# ---------------------------------------------------------------------------
# Blunted cone Mach 11.3, paper-grade mesh (600x200, ~120k cells), 28 cores.
#
# Template OpenPBS/PBS Professional (qsub), coda "cpu" del cluster di
# Federico (resources_max.ncpus=28, walltime max 48h, max_user_run=4 —
# verificato con `qstat -Qf`). Usa la sintassi moderna "-l select=..."
# (non "-l nodes=N:ppn=M", che su questo PBS non si traduce correttamente
# e fa scattare "violates queue and/or server resource limits").
#
# Prerequisiti sul cluster (una tantum: qsub setup-cluster.sh, vedi
# cluster/README.md):
#   1. repo copiata sul cluster
#   2. thirdParty/Mutationpp copiata a mano (e' gitignorata, vedi README)
#   3. setup-cluster.sh eseguito con successo (build Mutation++ +
#      librerie custom del progetto)
#
# Verificato sul cluster di Federico: OpenFOAM-13 e' preinstallato in
# /opt/openfoam13 (va sorgentato esplicitamente PRIMA del bashrc del
# progetto - quest'ultimo aggiunge solo le variabili POLIMI_*/Mutation++
# sopra un ambiente OpenFOAM gia' attivo, non lo sostituisce. Senza
# questo passaggio foamRun/blockMesh non sono in PATH e il job muore
# subito senza log utile).
#
# Uso:  qsub job-cone-fine.sh
#
# NOTA: il meccanismo di stage-out dell'output di PBS (-o/-j oe) su
# questo cluster non e' affidabile (verificato anche con job banali).
# Il job scrive quindi il proprio log direttamente su disco con `exec`.
# ---------------------------------------------------------------------------
#PBS -N cone-m11
#PBS -q cpu
#PBS -l select=1:ncpus=28:mpiprocs=28
#PBS -l walltime=24:00:00
#PBS -j oe
#PBS -o cone-fine.log

# Stima walltime: il run coarse (7200 celle, 4 core) ha impiegato
# 4548 s (ExecutionTime OpenFOAM) per arrivare a convergenza reale
# (verificato: campi stabili <0.5% oltre quel punto). Scalando
# linearmente a 120k celle / 28 core: ~4548*(120000/7200)*(4/28) =
# ~10800 s = ~3h. 24h da margine ampio (scaling parallelo a 28 core
# mai testato finora, solo fino a 4) restando ben dentro il tetto di
# 48h della coda cpu. Se comunque non bastasse: il job puo' essere
# rilanciato (startFrom latestTime riprende dall'ultimo checkpoint
# LTS, verificato funzionante dopo lo spegnimento PC del 6 lug 2026).

# --- ambiente (ADATTARE se il path della repo e' diverso) ------------------
REPO="/work/u10806848/Aerospace-Module-2"
CASE_DIR="$REPO/applications/test/bluntedCone2D"

exec > "$CASE_DIR/cluster/cone-fine.log" 2>&1
echo "=== job avviato $(date) su $(hostname) ==="

# Sourcing isolato in un sottoprocesso bash pulito: nei job batch (ma
# non in una sessione interattiva) il sourcing diretto del bashrc di
# OpenFOAM scatena un bug noto di bash ("pop_var_context: head of
# shell_variables not a function context"), verosimilmente per opzioni
# ereditate dall'ambiente PBS. Vedi setup-cluster.sh per i dettagli.
ENV_DUMP=$(mktemp)
bash -c "source /opt/openfoam13/etc/bashrc && source '$REPO/etc/bashrc' && export -p" > "$ENV_DUMP"
if [ ! -s "$ENV_DUMP" ]; then
    echo "ERRORE: il sourcing isolato non ha prodotto variabili d'ambiente."
    exit 1
fi
source "$ENV_DUMP"
rm -f "$ENV_DUMP"

# PBS parte dalla home: spostarsi nel caso
cd "$CASE_DIR" || exit 1

# --- run --------------------------------------------------------------------
# Con OpenMPI compilato con supporto PBS, mpirun rileva i core dal job.
# In caso contrario, decommenta la variante con la machinefile:
#   export OMPI_MCA_orte_default_hostfile=$PBS_NODEFILE
./Allclean
./Allrun -fine 28

# Il postProcess (plot + metriche) gira alla fine dentro Allrun.
# Risultati: cone-stagnation.png, cone-surface.png, log.foamRun
# ---------------------------------------------------------------------------

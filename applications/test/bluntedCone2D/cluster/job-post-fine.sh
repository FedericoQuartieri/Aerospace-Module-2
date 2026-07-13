#!/bin/bash
# ---------------------------------------------------------------------------
# Post-processing del run fine: grandezze di parete + confronto con la Fig 2
# del paper (riferimenti digitalizzati in references/).
#
# Gira in serie su 1 core e dura pochi minuti, ma va comunque sottomesso
# come job: sul nodo di login un processo in foreground muore al primo calo
# di connessione (SIGHUP), portandosi via l'output.
#
# Prerequisiti: il run fine (job-cone-fine.sh) deve essere terminato e
# ricostruito — Allrun chiama reconstructPar da solo, quindi le directory
# temporali stanno nel caso, non nei processor*/.
#
# Uso:  qsub job-post-fine.sh
# Log:  cluster/cone-post.log   (contiene gli scarti Cp/Cf/St stampati da
#       compare-fig2.py; le figure sono fig2-comparison.png e
#       fig2-stagnation-comparison.png nella dir del caso)
# ---------------------------------------------------------------------------
#PBS -N cone-post
#PBS -q cpu
#PBS -l select=1:ncpus=1
#PBS -l walltime=00:30:00
#PBS -j oe
#PBS -o cone-post.log

REPO="/work/u10806848/Aerospace-Module-2"
CASE_DIR="$REPO/applications/test/bluntedCone2D"

# Lo stage-out di PBS (-o) su questo cluster non consegna mai il file:
# il job scrive il proprio log direttamente su disco. Vedi job-cone-fine.sh.
exec > "$CASE_DIR/cluster/cone-post.log" 2>&1
echo "=== post-processing avviato $(date) su $(hostname) ==="

# Sourcing isolato in un sottoprocesso pulito: nei job batch il sourcing
# diretto del bashrc di OpenFOAM scatena il bug bash "pop_var_context".
ENV_DUMP=$(mktemp)
bash -c "source /opt/openfoam13/etc/bashrc && source '$REPO/etc/bashrc' && export -p" > "$ENV_DUMP"
if [ ! -s "$ENV_DUMP" ]; then
    echo "ERRORE: il sourcing isolato non ha prodotto variabili d'ambiente."
    exit 1
fi
source "$ENV_DUMP"
rm -f "$ENV_DUMP"

cd "$CASE_DIR" || exit 1

echo "--- tempi disponibili ---"
foamListTimes

echo "--- wallHeatFlux ---"
foamPostProcess -solver shockThermo -func wallHeatFlux -latestTime \
    || { echo "FALLITO wallHeatFlux"; exit 1; }

echo "--- wallShearStress ---"
foamPostProcess -solver shockThermo -func wallShearStress -latestTime \
    || { echo "FALLITO wallShearStress"; exit 1; }

echo "--- confronto con la Fig 2 del paper ---"
python3 compare-fig2.py || { echo "FALLITO compare-fig2.py"; exit 1; }

echo "=== post-processing completato $(date) ==="
# ---------------------------------------------------------------------------

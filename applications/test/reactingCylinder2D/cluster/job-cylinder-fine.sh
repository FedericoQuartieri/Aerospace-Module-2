#!/bin/bash
# ---------------------------------------------------------------------------
# Mach 20 reacting cylinder, paper-grade mesh (2 x 300 x 260 = 156k celle,
# prima cella 2 um lato vento), 28 core.
#
# Stesso schema del job del cono (job-cone-fine.sh, vedi i commenti li' per
# la trafila PBS di questo cluster: -q cpu esplicita, select=...,
# stage-out -o inaffidabile -> log auto-gestito, sourcing isolato per il
# bug bash pop_var_context nei job batch).
#
# Walltime: il cono fine (120k celle, non reagente) stava largamente in
# 24h; la chimica N2/N con Mutation++ costa ~2-3x per cella e le celle
# sono 1.3x -> 36h di margine (tetto coda: 48h). Il caso riparte dai
# checkpoint LTS (startFrom latestTime) se serve rilanciare.
#
# Uso:  qsub job-cylinder-fine.sh
# ---------------------------------------------------------------------------
#PBS -N cyl-m20
#PBS -q cpu
#PBS -l select=1:ncpus=28:mpiprocs=28
#PBS -l walltime=36:00:00
#PBS -j oe
#PBS -o cylinder-fine.log

REPO="/work/u10806848/Aerospace-Module-2"
CASE_DIR="$REPO/applications/test/reactingCylinder2D"

exec > "$CASE_DIR/cluster/cylinder-fine.log" 2>&1
echo "=== job avviato $(date) su $(hostname) ==="

ENV_DUMP=$(mktemp)
bash -c "source /opt/openfoam13/etc/bashrc && source '$REPO/etc/bashrc' && export -p" > "$ENV_DUMP"
if [ ! -s "$ENV_DUMP" ]; then
    echo "ERRORE: il sourcing isolato non ha prodotto variabili d'ambiente."
    exit 1
fi
source "$ENV_DUMP"
rm -f "$ENV_DUMP"

cd "$CASE_DIR" || exit 1

./Allclean
./Allrun -fine 28

# Allrun chiude gia' con foamPostProcess (wallHeatFlux/wallShearStress) e
# postProcess-cylinder.py: standoff, C_D, C_H e i due PNG stanno in coda
# a questo log e nella dir del caso.
# ---------------------------------------------------------------------------

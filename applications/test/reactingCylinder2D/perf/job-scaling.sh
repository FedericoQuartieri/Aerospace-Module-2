#!/bin/bash
# ---------------------------------------------------------------------------
# M9 strong-scaling sweep of the Mach 20 reacting cylinder on one 28-core
# node. Runs the same coarse case on 1,2,4,7,14,28 MPI ranks over a fixed
# short pseudo-time window and records the per-step wall cost.
#
# Same PBS trafila as the other jobs on this cluster (vedi job-cone-fine.sh):
# -q cpu esplicita, select=..., stage-out -o inaffidabile -> log auto-gestito,
# sourcing isolato per il bug bash pop_var_context nei job batch.
#
# Uso:  qsub job-scaling.sh
# Log:  perf/scaling.log  (contiene la tabella np/ExecTime/us-per-cell-step;
#       il CSV perf/results-scaling.csv poi lo dai a plot-scaling.py)
# ---------------------------------------------------------------------------
#PBS -N cyl-scaling
#PBS -q cpu
#PBS -l select=1:ncpus=28:mpiprocs=28
#PBS -l walltime=03:00:00
#PBS -j oe
#PBS -o scaling.log

REPO="/work/u10806848/Aerospace-Module-2"
CASE_DIR="$REPO/applications/test/reactingCylinder2D"

exec > "$CASE_DIR/perf/scaling.log" 2>&1
echo "=== scaling sweep avviato $(date) su $(hostname) ==="

ENV_DUMP=$(mktemp)
bash -c "source /opt/openfoam13/etc/bashrc && source '$REPO/etc/bashrc' && export -p" > "$ENV_DUMP"
if [ ! -s "$ENV_DUMP" ]; then
    echo "ERRORE: sourcing isolato vuoto."; exit 1
fi
source "$ENV_DUMP"
rm -f "$ENV_DUMP"

cd "$CASE_DIR" || exit 1

# strong scaling: mesh fissa, rank crescenti. 28 core = 1 nodo pieno.
# 7 e 14 sono divisori di 28 per una griglia di punti regolare.
#
# Due mesh, per mostrare il compromesso compute/comunicazione:
#  - coarse (9k celle): a 28 rank fa ~320 celle/rank -> communication-bound,
#    lo scaling satura (come il kernel 0D del gruppo 1, ma su solver vero);
#  - fine (156k celle): ~5600 celle/rank a 28 -> compute-bound, scala bene.
# La finestra di timing e' piu' corta sul fine (celle piccole = piu' step
# per pseudo-tempo, e ogni step costa di piu').
echo "########## STRONG SCALING - COARSE ##########"
./perf/run-scaling.sh "1 2 4 7 14 28" coarse 2e-5
echo "########## STRONG SCALING - FINE ##########"
./perf/run-scaling.sh "1 2 4 7 14 28" fine 3e-6

echo "=== completato $(date) ==="
echo "CSV: perf/results-scaling-coarse.csv e perf/results-scaling-fine.csv"
echo "grafico: python3 perf/plot-scaling.py"
# ---------------------------------------------------------------------------

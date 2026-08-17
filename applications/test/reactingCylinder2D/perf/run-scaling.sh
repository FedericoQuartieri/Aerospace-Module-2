#!/bin/bash
# ---------------------------------------------------------------------------
# Strong-scaling harness for the Mach 20 reacting cylinder (M9).
#
# Unlike a 0D kernel replicated over independent cells (embarrassingly
# parallel), this measures the REAL distributed-memory solver: MPI domain
# decomposition of a coupled 2D reacting two-temperature simulation, with
# halo exchange and a per-cell Mutation++ bridge (one Mixture per rank).
#
# Method: run a FIXED short pseudo-time window from the impulsive start on
# each rank count, record ExecutionTime and the step count. LTS local dt
# is per-cell (decomposition-independent), so the step count to reach
# endTime is ~constant across ranks -> ExecutionTime is a clean strong-
# scaling metric. Per-step wall time = ExecutionTime/nSteps normalises out
# any residual step-count jitter.
#
# Usage:  ./perf/run-scaling.sh "1 2 4 8"  [coarse|fine]  [timeWindow]
# Output: perf/results-scaling.csv  (np, cells, execTime, nSteps, ...)
# ---------------------------------------------------------------------------
set -u
cd "$(dirname "$0")/.." || exit 1          # -> case dir
CASE=$(pwd)
NPROCS="${1:-1 2 4 8}"
MESH="${2:-coarse}"
TWIN="${3:-2e-5}"                          # timing window [pseudo-s]

export MPP_DATA_DIRECTORY="$(cd ../nonEqTTv && pwd)/mutation-data-noel"
mkdir -p perf
RES="perf/results-scaling-${MESH}.csv"
echo "np,cells,execTime_s,nSteps,tstep_ms,us_per_cell_step" > "$RES"

# one mesh, shared by every rank count
python3 makeBlockMeshDict.py "$MESH" > perf/log.mesh 2>&1
blockMesh > perf/log.blockMesh 2>&1 || { echo "blockMesh FAILED"; tail -5 perf/log.blockMesh; exit 1; }
# makeBlockMeshDict.py stampa "... cells = 9000, first ..."
NCELLS=$(grep -oE "cells = [0-9]+" perf/log.mesh | grep -oE "[0-9]+" | head -1)
[ -z "${NCELLS:-}" ] && NCELLS=$(grep -oE "nCells:?[[:space:]]+[0-9]+" perf/log.blockMesh | grep -oE "[0-9]+" | head -1)
[ -z "${NCELLS:-}" ] && NCELLS=1   # fallback: non rompere il calcolo
echo "mesh $MESH: $NCELLS celle, finestra di timing endTime=$TWIN"

for NP in $NPROCS; do
    rm -rf 0 processor*
    foamListTimes -rm > /dev/null 2>&1
    cp -r 0.bak 0
    foamDictionary -entry endTime -set "$TWIN" system/controlDict > /dev/null 2>&1
    LOG="perf/log.np${NP}"
    if [ "$NP" = "1" ]; then
        foamRun > "$LOG" 2>&1
    else
        foamDictionary -entry numberOfSubdomains -set "$NP" \
            system/decomposeParDict > /dev/null 2>&1
        decomposePar -force > "perf/log.decompose.np${NP}" 2>&1
        # --oversubscribe: innocuo sul cluster (slot PBS bastano), serve
        # solo per testare in locale piu' rank dei core fisici disponibili
        mpirun --oversubscribe -np "$NP" foamRun -parallel > "$LOG" 2>&1
    fi
    ET=$(grep "ExecutionTime" "$LOG" | tail -1 | awk '{print $3}')
    NS=$(grep -c "^Time = " "$LOG")
    if [ -z "${ET:-}" ] || [ "$NS" = "0" ]; then
        echo "np=$NP FALLITO:"; tail -6 "$LOG"; continue
    fi
    TSTEP=$(python3 -c "print(f'{$ET/$NS*1000:.4f}')")
    UPC=$(python3 -c "print(f'{$ET/$NS/$NCELLS*1e6:.3f}')")
    echo "$NP,$NCELLS,$ET,$NS,$TSTEP,$UPC" >> "$RES"
    printf "np=%-3s cells=%s ExecTime=%ss steps=%s  %sms/step  %sus/cell/step\n" \
        "$NP" "$NCELLS" "$ET" "$NS" "$TSTEP" "$UPC"
done

# restore the committed run settings
foamDictionary -entry endTime -set 9e-3 system/controlDict > /dev/null 2>&1
foamDictionary -entry numberOfSubdomains -set 28 system/decomposeParDict > /dev/null 2>&1
rm -rf 0 processor*
echo "risultati in $RES"

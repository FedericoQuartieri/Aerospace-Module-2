#!/bin/bash
# ---------------------------------------------------------------------------
# Blunted cone Mach 11.3, paper-grade mesh (600x200, ~120k cells), 28 cores.
#
# SLURM template - ADATTARE alle direttive del tuo cluster (partition,
# account, module load). Se il cluster usa PBS/altro, tradurre l'header.
#
# Prerequisiti sul cluster (una tantum, vedi cluster/README.md):
#   1. repo clonata e OpenFOAM-13 disponibile (module o build)
#   2. thirdParty/makeMutationpp completato
#   3. source etc/bashrc && src/thermophysicalModels/Allwmake
#      && applications/modules/Allwmake
#
# Uso:  sbatch job-cone-fine.sh
# ---------------------------------------------------------------------------
#SBATCH --job-name=cone-m11
#SBATCH --nodes=1
#SBATCH --ntasks=28
#SBATCH --time=12:00:00
#SBATCH --output=cone-fine-%j.log

# --- ambiente (ADATTARE) ---------------------------------------------------
# module load openfoam/13  # oppure:
REPO="$HOME/Aerospace-Module-2"          # ADATTARE al path della repo
source "$REPO/etc/bashrc"

cd "$REPO/applications/test/bluntedCone2D" || exit 1

# --- run --------------------------------------------------------------------
./Allclean
./Allrun -fine 28

# Il postProcess (plot + metriche) gira alla fine dentro Allrun.
# Risultati: cone-stagnation.png, cone-surface.png, log.foamRun
# ---------------------------------------------------------------------------

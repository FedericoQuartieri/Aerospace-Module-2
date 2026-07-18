#!/bin/bash
# ---------------------------------------------------------------------------
# Mach 20 reacting cylinder, mesh COARSE (9k celle), 28 core, da zero.
#
# L'avvio impulsivo a M20 e' molto piu' violento del cono M11: i 30k
# pseudo-step LTS che bastavano li' NON bastano qui (verificato in locale:
# a 30k step le pressioni di parete oscillano ancora, in decadimento
# lento). endTime 9e-3 = 90k step; su 28 core ~1h. Il postProcess in coda
# stampa anche la variazione dei campi tra gli ultimi due snapshot: se
# non e' scesa sotto ~1% il run va prolungato (endTime e' runTimeModifiable,
# il job si puo' risottomettere e riparte da latestTime... ma Allclean
# in testa rifa' tutto: per riprendere, commentare Allclean e lanciare
# a mano i passi di Allrun senza rigenerare mesh/0).
#
# Schema PBS identico agli altri job (vedi job-cone-fine.sh per i
# dettagli della trafila su questo cluster).
#
# Uso:  qsub job-cylinder-coarse.sh
# ---------------------------------------------------------------------------
#PBS -N cyl-m20-coarse
#PBS -q cpu
#PBS -l select=1:ncpus=28:mpiprocs=28
#PBS -l walltime=06:00:00
#PBS -j oe
#PBS -o cylinder-coarse.log

REPO="/work/u10806848/Aerospace-Module-2"
CASE_DIR="$REPO/applications/test/reactingCylinder2D"

exec > "$CASE_DIR/cluster/cylinder-coarse.log" 2>&1
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
./Allrun -parallel 28

# Il postProcess (standoff, C_D, C_H, check di convergenza, PNG) gira
# dentro Allrun; il riassunto sta in coda a questo log.
# ---------------------------------------------------------------------------

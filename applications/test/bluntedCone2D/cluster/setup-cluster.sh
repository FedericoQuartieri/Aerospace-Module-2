#!/bin/bash
# ---------------------------------------------------------------------------
# Setup una tantum sul cluster: build di Mutation++ e delle librerie custom
# del progetto (thermo bridge + shockThermo). Va lanciato UNA VOLTA prima
# di job-cone-fine.sh (o di qualunque altro caso).
#
# thirdParty/Mutationpp e' gitignorata: se manca, questo script la clona
# direttamente da GitHub (e' un mirror pubblico, nessuna patch locale -
# verificato con `git status`/`git remote -v` sul checkout del laptop),
# pinnata al commit usato in sviluppo. Molto piu' comodo di un rsync di
# ~1 GB dal laptop.
#
# Verificato sul cluster di Federico:
#   - OpenFOAM-13 e' preinstallato in /opt/openfoam13 (nessun modulo,
#     nessuna build necessaria per OpenFOAM stesso - va solo sorgentato)
#   - gcc 13.3.0, OpenMPI 4.1.6, cmake 3.28.3 presenti di sistema
#
# Uso:  qsub setup-cluster.sh
#       (poi controllare setup-cluster.log per l'esito)
#
# NOTA: il meccanismo di stage-out dell'output di PBS (-o/-j oe) su
# questo cluster non e' affidabile - il file non compare mai, ne' col
# path relativo ne' assoluto, anche per job banali (verificato). Lo
# script scrive quindi il proprio log direttamente su disco con `exec`,
# bypassando PBS.
# ---------------------------------------------------------------------------
#PBS -N cone-setup
#PBS -q cpu
#PBS -l select=1:ncpus=8:mpiprocs=8
#PBS -l walltime=03:00:00
#PBS -j oe
#PBS -o setup-cluster.log

REPO="/work/u10806848/Aerospace-Module-2"          # ADATTARE se serve
LOG="$REPO/applications/test/bluntedCone2D/cluster/setup-cluster.log"

exec > "$LOG" 2>&1
echo "=== job avviato $(date) su $(hostname) ==="

# OpenFOAM-13 di sistema PRIMA del bashrc del progetto: il bashrc del
# progetto aggiunge solo le variabili POLIMI_*/Mutation++ sopra un
# ambiente OpenFOAM gia' attivo, non lo sostituisce.
#
# NOTA: sourcing isolato in un sottoprocesso bash pulito. OpenFOAM usa
# variabili "local" dentro funzioni nel suo bashrc; in questo job batch
# (ma non in una sessione interattiva sullo stesso cluster) il sourcing
# diretto scatena un bug noto di bash ("pop_var_context: head of
# shell_variables not a function context") che termina la shell -
# verosimilmente per opzioni ereditate dall'ambiente PBS (es. SHELLOPTS)
# che non si presentano in una sessione interattiva pulita. Il
# sottoprocesso isola il problema: si butta via lo stato delle funzioni
# interne di OpenFOAM e si importano solo le variabili d'ambiente
# risultanti, che sono tutto cio' che serve al resto dello script.
ENV_DUMP=$(mktemp)
bash -c "source /opt/openfoam13/etc/bashrc && source '$REPO/etc/bashrc' && export -p" > "$ENV_DUMP"
if [ ! -s "$ENV_DUMP" ]; then
    echo "ERRORE: il sourcing isolato non ha prodotto variabili d'ambiente."
    echo "Il sottoprocesso e' probabilmente fallito - vedi sopra per errori."
    exit 1
fi
source "$ENV_DUMP"
rm -f "$ENV_DUMP"

set -e

cd "$REPO" || exit 1

echo "=== verifica/recupero Mutation++ ==="
MPP_COMMIT="117df0b14baf41715c9cf04be9fc7f438fa4513f"   # master, verificato in dev
if [ ! -d thirdParty/Mutationpp ]; then
    # I nodi di calcolo spesso NON hanno accesso internet (solo il login
    # node). Se questo clone fallisce/appende, e' quello il motivo:
    # fare `git clone` a mano sul login node PRIMA di sottomettere il
    # job (vedi README.md), poi rilanciare.
    echo "thirdParty/Mutationpp non trovata: provo a clonare da GitHub..."
    echo "(se il nodo di calcolo non ha internet, questo fallisce: clonare"
    echo " a mano sul login node prima - vedi README.md)"
    git clone https://github.com/mutationpp/Mutationpp thirdParty/Mutationpp
    (cd thirdParty/Mutationpp && git checkout "$MPP_COMMIT")
fi
echo "WM_PROJECT_DIR=$WM_PROJECT_DIR"
echo "Mutation++ commit: $(cd thirdParty/Mutationpp && git rev-parse HEAD)"

echo "=== 1/3: build Mutation++ ==="
cd thirdParty
./makeMutationpp
cd ..

echo "=== 2/3: build thermophysicalModels custom ==="
src/thermophysicalModels/Allwmake

echo "=== 3/3: build applications/modules (shockThermo) ==="
applications/modules/Allwmake

echo "=== verifica finale ==="
which foamRun
ls "$FOAM_USER_LIBBIN" 2>/dev/null | grep -E "highEnthalpy|shockThermo|shockFluid" \
    || ls platforms/linux64GccDPInt32Opt/lib | grep -E "highEnthalpy|shockThermo|shockFluid"

echo "=== setup completato ==="
# ---------------------------------------------------------------------------

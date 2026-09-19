#!/bin/bash
# Prepara l'ambiente e passa il comando. Ordine obbligato: prima OpenFOAM
# (definisce WM_*, FOAM_USER_LIBBIN, WM_THIRD_PARTY_DIR), poi il bashrc del
# progetto (definisce POLIMI_*), che si appoggia alle variabili del primo.

set -e

source /opt/openfoam13/etc/bashrc
source /project/etc/bashrc

# Mutation++ sta nell'immagine, non in thirdParty/: etc/config.sh/mutationpp ha
# appena puntato MPP_DIRECTORY al bind mount, lo riportiamo su /opt
export MPP_DIRECTORY=/opt/Mutationpp
export MPP_EIGEN=$MPP_DIRECTORY/thirdparty/eigen
export MPP_DATA_DIRECTORY=$MPP_DIRECTORY/data
export PATH=$MPP_DIRECTORY/install/bin:$PATH
export LD_LIBRARY_PATH=$MPP_DIRECTORY/install/lib:$LD_LIBRARY_PATH

exec "$@"

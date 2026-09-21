#!/bin/bash
# Sets up the environment and runs the command. The order is mandatory: first
# OpenFOAM (defines WM_*, FOAM_USER_LIBBIN, WM_THIRD_PARTY_DIR), then the
# project bashrc (defines POLIMI_*), which relies on the variables of the first.

set -e

source /opt/openfoam13/etc/bashrc
source /project/etc/bashrc

# Mutation++ lives in the image, not in thirdParty/: etc/config.sh/mutationpp has
# just pointed MPP_DIRECTORY to the bind mount, so we point it back to /opt
export MPP_DIRECTORY=/opt/Mutationpp
export MPP_EIGEN=$MPP_DIRECTORY/thirdparty/eigen
export MPP_DATA_DIRECTORY=$MPP_DIRECTORY/data
export PATH=$MPP_DIRECTORY/install/bin:$PATH
export LD_LIBRARY_PATH=$MPP_DIRECTORY/install/lib:$LD_LIBRARY_PATH

exec "$@"

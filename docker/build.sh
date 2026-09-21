#!/bin/sh
# Builds the development image (OpenFOAM-13 + Mutation++).
# Only needs to be rerun if the Dockerfile changes: the project code is mounted,
# not copied, so changes to the sources do not require a rebuild.
cd ${0%/*} || exit 1
exec docker build -t aero-m2:dev .

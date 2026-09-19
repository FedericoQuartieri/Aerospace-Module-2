#!/bin/sh
# Costruisce l'immagine di sviluppo (OpenFOAM-13 + Mutation++ + OpenMP).
# Da rifare solo se cambia il Dockerfile: il codice del progetto e' montato,
# non copiato, quindi le modifiche ai sorgenti non richiedono un rebuild.
cd ${0%/*} || exit 1
exec docker build -t aero-m2:dev .

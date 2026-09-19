#!/bin/sh
# Apre una shell nel container con il progetto montato su /project.
#
#   ./docker/run.sh                      shell interattiva
#   ./docker/run.sh ./Allwmake           compila e esce
#   OMP_NUM_THREADS=4 ./docker/run.sh    con 4 thread invece di tutti
#
# Il progetto e' un bind mount: quello che si compila dentro finisce in
# platforms/ sul Mac, e le modifiche fatte sul Mac si vedono subito dentro.
#
# /root sta invece in un volume Docker, non sul bind mount: e' li' che OpenFOAM
# installa le librerie dell'utente ($FOAM_USER_LIBBIN), e senza il volume
# sparirebbero a ogni --rm, costringendo a ricompilare. Il volume e' anche un
# filesystem nativo della VM, molto piu' veloce di quello condiviso con macOS.
cd ${0%/*}/.. || exit 1

if [ -t 0 ]; then tty="-it"; else tty=""; fi

exec docker run --rm $tty \
    -v "$(pwd)":/project \
    -v aero-m2-home:/root \
    -e OMP_NUM_THREADS \
    aero-m2:dev "${@:-bash}"

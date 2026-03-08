#!/bin/sh
set -eu

cd "${0%/*}" || exit 1

probeRoot="postProcessing/probeCenter"

if [ ! -d "$probeRoot" ]; then
    echo "Errore: directory $probeRoot non trovata. Esegui prima foamRun." >&2
    exit 1
fi

probeTimeDir=$(find "$probeRoot" -mindepth 1 -maxdepth 1 -type d | sort -V | tail -n 1)

if [ -z "$probeTimeDir" ]; then
    echo "Errore: nessun output probe in $probeRoot." >&2
    exit 1
fi

for field in p T Tve rho; do
    if [ ! -f "$probeTimeDir/$field" ]; then
        echo "Errore: file $probeTimeDir/$field non trovato." >&2
        exit 1
    fi
done

tmpP=$(mktemp)
tmpT=$(mktemp)
tmpTve=$(mktemp)
tmpRho=$(mktemp)
trap 'rm -f "$tmpP" "$tmpT" "$tmpTve" "$tmpRho"' EXIT

awk '!/^#/ && NF >= 2 {v=$2; gsub(/[()]/, "", v); print $1, v}' "$probeTimeDir/p" > "$tmpP"
awk '!/^#/ && NF >= 2 {v=$2; gsub(/[()]/, "", v); print $1, v}' "$probeTimeDir/T" > "$tmpT"
awk '!/^#/ && NF >= 2 {v=$2; gsub(/[()]/, "", v); print $1, v}' "$probeTimeDir/Tve" > "$tmpTve"
awk '!/^#/ && NF >= 2 {v=$2; gsub(/[()]/, "", v); print $1, v}' "$probeTimeDir/rho" > "$tmpRho"

{
    echo "time,p,T,Tve,rho"
    paste "$tmpP" "$tmpT" "$tmpTve" "$tmpRho" | awk '{print $1 "," $2 "," $4 "," $6 "," $8}'
} > results.csv

echo "Creato results.csv da $probeTimeDir"

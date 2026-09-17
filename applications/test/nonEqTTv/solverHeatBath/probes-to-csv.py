#!/usr/bin/env python3
# raccoglie i probe del solver (una colonna per grandezza, una riga per tempo)
# nel file ../output/<figura>-solver.csv con le stesse colonne del programma 0D
#
# uso: python3 probes-to-csv.py <figura>

import os
import sys
import numpy as np

fig = sys.argv[1]
here = os.path.dirname(os.path.abspath(__file__))
probes = os.path.join(here, "postProcessing", "probes", "0")


def load(field):
    # file probe: "# commenti" poi "tempo valore"
    data = np.loadtxt(os.path.join(probes, field), comments="#", ndmin=2)
    return data[:, 0], data[:, 1]


t, T = load("T")
t, Tve = load("Tve")
columns = [("t", t), ("Ttr", T), ("Tv", Tve)]

# con la chimica servono anche le densita' numeriche normalizzate n/n0
if fig.startswith("fig7"):
    t, rho = load("rho")
    NA = 6.02214076e23
    n = {}
    for name, M in (("N2", 28.0134e-3), ("N", 14.0067e-3)):
        t, Y = load(name)
        n[name] = rho * Y / M * NA
    n0 = n["N2"][0] + n["N"][0]
    columns += [("N2", n["N2"] / n0), ("N", n["N"] / n0)]

out = os.path.join(here, "..", "output", fig + "-solver.csv")
with open(out, "w") as f:
    f.write(",".join(name for name, _ in columns) + "\n")
    for i in range(len(t)):
        f.write(",".join("%.8g" % col[i] for _, col in columns) + "\n")
print("scritto " + os.path.relpath(out, here))

#!/usr/bin/env python3
# collects the solver probes (one column per quantity, one row per time)
# into the file ../output/<figure>-solver.csv with the same columns as the 0D program
#
# usage: python3 probes-to-csv.py <figure>

import os
import sys
import numpy as np

fig = sys.argv[1]
here = os.path.dirname(os.path.abspath(__file__))
probes = os.path.join(here, "postProcessing", "probes", "0")


def load(field):
    # probe file: "# comments" then "time value"
    data = np.loadtxt(os.path.join(probes, field), comments="#", ndmin=2)
    return data[:, 0], data[:, 1]


t, T = load("T")
t, Tve = load("Tve")
columns = [("t", t), ("Ttr", T), ("Tv", Tve)]

# with chemistry the normalised number densities n/n0 are also needed
if fig.startswith("fig7"):
    t, rho = load("rho")
    NA = 6.02214076e23
    n = {}
    for name, M in (("N2", 28.0134e-3), ("N", 14.0067e-3)):
        t, Y = load(name)
        n[name] = rho * Y / M * NA
    n0 = n["N2"][0] + n["N"][0]
    columns += [("N2", n["N2"] / n0), ("N", n["N"] / n0)]

# rows to keep: every step at the beginning, then 200 points per decade (same
# rule as OutputSchedule in the 0D programs), always the last one
keep = []
nextStep = 1.0
for i, ti in enumerate(t):
    step = round(ti / 1.0e-9)
    if step >= nextStep or i == len(t) - 1:
        keep.append(i)
        nextStep = max(nextStep + 1.0, nextStep * 10.0 ** (1.0 / 200.0), step + 1.0)

out = os.path.join(here, "..", "output", fig + "-solver.csv")
with open(out, "w") as f:
    f.write(",".join(name for name, _ in columns) + "\n")
    for i in keep:
        f.write(",".join("%.8g" % col[i] for _, col in columns) + "\n")
print("written " + os.path.relpath(out, here))

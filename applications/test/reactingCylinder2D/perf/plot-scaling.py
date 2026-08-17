"""Strong-scaling analysis for the Mach 20 reacting cylinder (M9).

Reads perf/results-scaling-<mesh>.csv (from run-scaling.sh) for every mesh
that exists (coarse, fine) and produces:
  - a table of speedup S(n)=t(1)/t(n) and parallel efficiency E=S/n per mesh
  - perf/scaling.png (wall time + speedup vs ranks, ideal line, both meshes)

Speedup uses per-step wall time (execTime/nSteps), decomposition-independent
in step count, so it isolates the compute+communication cost. Overlaying the
coarse (few cells/rank -> communication-bound, saturates) and fine (many
cells/rank -> compute-bound, scales) meshes shows the classic strong-scaling
trade-off on a REAL coupled solver - not a 0D kernel of independent cells.

Usage: python3 perf/plot-scaling.py
"""

import os
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.ticker as mticker
import matplotlib.pyplot as plt

HERE = os.path.dirname(os.path.abspath(__file__))
MESHES = [("coarse", "crimson", "o"), ("fine", "navy", "s")]

fig, ax = plt.subplots(1, 2, figsize=(12, 4.6))
allranks = set()
found = False

for mesh, col, mk in MESHES:
    csv = os.path.join(HERE, f"results-scaling-{mesh}.csv")
    if not os.path.exists(csv):
        continue
    found = True
    d = np.genfromtxt(csv, delimiter=",", names=True)
    ranks = np.atleast_1d(d["np"]).astype(int)
    tstep = np.atleast_1d(d["tstep_ms"])
    cells = int(np.atleast_1d(d["cells"])[0])
    o = np.argsort(ranks)
    ranks, tstep = ranks[o], tstep[o]
    allranks.update(ranks.tolist())

    t1 = tstep[ranks == 1][0] if (ranks == 1).any() else tstep[0]
    speedup = t1/tstep
    eff = speedup/ranks

    print(f"\nstrong scaling - {mesh} ({cells} celle)")
    print(f"{'ranks':>6} {'ms/step':>10} {'speedup':>9} {'efficienza':>11}")
    for n, t, s, e in zip(ranks, tstep, speedup, eff):
        print(f"{n:>6} {t:>10.3f} {s:>9.2f} {e:>10.2f}")

    lbl = f"{mesh} ({cells} celle)"
    ax[0].plot(ranks, tstep, mk + "-", c=col, label=lbl)
    ax[1].plot(ranks, speedup, mk + "-", c=col, label=lbl)

if not found:
    raise SystemExit("nessun perf/results-scaling-*.csv trovato: gira prima "
                     "run-scaling.sh o job-scaling.sh")

xr = sorted(allranks)
ax[0].set_xlabel("MPI ranks"); ax[0].set_ylabel("wall time per step [ms]")
ax[0].set_title("costo per step"); ax[0].grid(alpha=0.3, which="both")
ax[0].set_xscale("log", base=2); ax[0].set_yscale("log")
ax[0].set_xticks(xr); ax[0].get_xaxis().set_major_formatter(mticker.ScalarFormatter())
ax[0].legend(fontsize=8)

ax[1].plot(xr, xr, "--", c="gray", label="ideale (lineare)")
ax[1].set_xlabel("MPI ranks"); ax[1].set_ylabel("speedup  t(1)/t(n)")
ax[1].set_title("strong scaling"); ax[1].grid(alpha=0.3)
ax[1].legend(fontsize=8)

fig.suptitle("Mach 20 reacting cylinder - MPI strong scaling")
fig.tight_layout()
out = os.path.join(HERE, "scaling.png")
fig.savefig(out, dpi=150)
print(f"\nsalvato {out}")

"""Plots for the reacting N2-N heat bath (paper Figs. 7 and 8).

Overlays the two integration paths of Test-N2N-reacting:
  - mpp-park05:    Mutation++ end-to-end (dissociation at Tf=sqrt(T*Tv))
  - manual-park07: manual chemistry at the paper exponent Tf=T^0.7*Tv^0.3

Usage: python3 plot-temps-reacting.py <T_tr_init> <T_ve_init>
Reads  output/results-N2N-reacting-<Ttr>-<Tve>-{mpp-park05,manual-park07}.csv
Writes output/temp-curves-N2N-reacting-<Ttr>-<Tve>.png
       output/density-curves-N2N-reacting-<Ttr>-<Tve>.png
"""

import sys
import numpy as np
import matplotlib.pyplot as plt

Ttr, Tve = int(sys.argv[1]), int(sys.argv[2])
base = f"output/results-N2N-reacting-{Ttr}-{Tve}-"

curves = {
    "Mutation++ (q=0.5)": ("-", np.loadtxt(base + "mpp-park05.csv",
                                           delimiter=",", skiprows=1)),
    "manual (q=0.7, paper)": ("--", np.loadtxt(base + "manual-park07.csv",
                                               delimiter=",", skiprows=1)),
}

# --- temperatures (paper Fig. 7a / 8a) ---
plt.figure()
for label, (ls, d) in curves.items():
    m = d[:, 0] > 0
    plt.plot(d[m, 0], d[m, 1] / 1e3, ls, label=f"T_tr, {label}")
    plt.plot(d[m, 0], d[m, 2] / 1e3, ls, label=f"T_ve, {label}")

dlast = curves["Mutation++ (q=0.5)"][1]
plt.axhline(dlast[-1, 1] / 1e3, color="red", linestyle=":", linewidth=0.8)
plt.text(dlast[-1, 0], dlast[-1, 1] / 1e3 * 1.03, f"{dlast[-1, 1]:.0f} K",
         ha="right")

plt.xscale("log")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K x 10^3)")
plt.title(f"reacting N2-N, T_tr(0)={Ttr} K, T_ve(0)={Tve} K")
plt.legend(fontsize=8)
plt.tight_layout()
out = f"output/temp-curves-N2N-reacting-{Ttr}-{Tve}.png"
plt.savefig(out)
print(f"Plot saved to {out}")

# --- normalised number densities (paper Fig. 7b / 8b) ---
plt.figure()
for label, (ls, d) in curves.items():
    m = d[:, 0] > 0
    plt.plot(d[m, 0], d[m, 3], ls, label=f"N2, {label}")
    plt.plot(d[m, 0], d[m, 4], ls, label=f"N, {label}")

plt.xscale("log")
plt.xlabel("Time (s)")
plt.ylabel("Normalised number density, n / n0")
plt.title(f"reacting N2-N, T_tr(0)={Ttr} K, T_ve(0)={Tve} K")
plt.legend(fontsize=8)
plt.tight_layout()
out = f"output/density-curves-N2N-reacting-{Ttr}-{Tve}.png"
plt.savefig(out)
print(f"Plot saved to {out}")

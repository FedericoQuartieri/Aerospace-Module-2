#!/usr/bin/env python3
# compares, for one figure of the paper, the curves of:
#   - paper       paper-data/<fig>-<curve>.csv   (extracted from the pdf, reference)
#   - standalone  output/<fig>.csv               (0D program with Mutation++)
#   - solver      output/<fig>-solver.csv        (OpenFOAM case, if it exists)
# prints the maximum error and the final values, and saves the plots output/<fig>.png
# (temperatures) and output/<fig>-n.png (number densities, where the paper gives them)
#
# usage: python3 compare.py <fig>      (e.g. fig3a)

import os
import sys
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# equilibrium temperatures stated in the text or in the figures of the paper
PAPER_TEQ = {
    "fig3a": 7623.3,
    "fig4-noEl": 21900.0,
    "fig4-el": 17700.0,
    "fig5": 24400.0,
    "fig6-noVV": 12100.0,
    "fig6-VV": 12100.0,
}

fig = sys.argv[1]
here = os.path.dirname(os.path.abspath(__file__))


def load_csv(path):
    # returns {column name: array}; rows starting with # are comments
    with open(path) as f:
        lines = [l for l in f if not l.startswith("#")]
    names = lines[0].strip().split(",")
    data = np.loadtxt(lines[1:], delimiter=",", ndmin=2)
    return {n: data[:, i] for i, n in enumerate(names)}


def read_Teq(path):
    # the 0D program writes at the top "# T_eq = ... K (energy conservation)"
    with open(path) as f:
        first = f.readline()
    if first.startswith("# T_eq"):
        return float(first.split("=")[1].split()[0])
    return None


std = load_csv(os.path.join(here, "output", fig + ".csv"))
Teq = read_Teq(os.path.join(here, "output", fig + ".csv"))

solver_path = os.path.join(here, "output", fig + "-solver.csv")
solver = load_csv(solver_path) if os.path.exists(solver_path) else None

# curves of the paper: files fig-<curve>.csv, or figa-/figb- for the figures
# with two panels (temperatures and number densities); the model variants
# (-park05, -park, -QK, -tauMpp, -nonPref, -vibOnly) are compared with the same curves
base = fig
for suffix in ("-park05", "-park", "-QK", "-tauMpp", "-nonPref", "-vibOnly"):
    base = base.replace(suffix, "")
paper = {}
for prefix in (base, base + "a", base + "b"):
    for path in sorted(glob.glob(os.path.join(here, "paper-data", prefix + "-*.csv"))):
        name = os.path.basename(path)[len(prefix) + 1:-4]
        paper[name] = load_csv(path)

curves = [n for n in std if n != "t"]

print("=== %s ===" % fig)
if Teq is not None:
    print("T_eq from energy conservation: %.1f K" % Teq, end="")
    if base in PAPER_TEQ:
        print("   (paper: %.1f K, deviation %.1f K)" % (PAPER_TEQ[base], Teq - PAPER_TEQ[base]))
    else:
        print("   (the paper does not give a value)")

# ---- errors with respect to the paper and final values, curve by curve
for name in curves:
    # temperatures are called T..., the rest are number densities n/n0
    is_temperature = name.startswith("T")
    unit = "K" if is_temperature else ""
    if name not in paper:
        print("%-6s no curve in the paper" % name)
        continue
    tp, yp = paper[name]["t"], paper[name][name]
    # comparison at the times of the paper (within the simulated interval)
    m = (tp >= std["t"][0]) & (tp <= std["t"][-1]) & (tp > 0)
    ys = np.interp(tp[m], std["t"], std[name])
    err = np.abs(ys - yp[m])
    if is_temperature:
        # maximum absolute error and maximum relative error, each with its own
        # instant: in general they do not occur at the same point
        rel = err / yp[m]
        print("%-6s standalone vs paper: max |err| = %7.1f K at t = %.2e s;"
              " max rel err = %.2f %% at t = %.2e s;"
              "  final (t = %.1e s): paper %.1f, standalone %.1f"
              % (name, err.max(), tp[m][err.argmax()], 100 * rel.max(), tp[m][rel.argmax()],
                 tp[m][-1], yp[m][-1], ys[-1]))
    else:
        # densities on a logarithmic scale: the deviation is measured in decades,
        # |log10(standalone/paper)|, as it is read on the plot
        ok = (ys > 0) & (yp[m] > 0)
        dec = np.abs(np.log10(ys[ok] / yp[m][ok]))
        print("%-6s standalone vs paper: max deviation = %.3f decades (factor %.2f) at t = %.2e s;"
              "  final (t = %.1e s): paper %.4g, standalone %.4g"
              % (name, dec.max(), 10 ** dec.max(), tp[m][ok][dec.argmax()],
                 tp[m][-1], yp[m][-1], ys[-1]))

    if solver is not None and name in solver:
        # solver vs standalone, at the solver times
        ts = solver["t"]
        m2 = (ts >= std["t"][0]) & (ts <= std["t"][-1]) & (ts > 0)
        ys2 = np.interp(ts[m2], std["t"], std[name])
        err2 = np.abs(solver[name][m2] - ys2)
        if is_temperature:
            rel2 = err2 / np.abs(ys2)
            print("%-6s solver vs standalone: max |err| = %7.1f K at t = %.2e s;"
                  " max rel err = %.2f %% at t = %.2e s;  final solver %.1f K"
                  % ("", err2.max(), ts[m2][err2.argmax()], 100 * rel2.max(),
                     ts[m2][rel2.argmax()], solver[name][-1]))
        else:
            ok2 = (ys2 > 0) & (solver[name][m2] > 0)
            dec2 = np.abs(np.log10(solver[name][m2][ok2] / ys2[ok2]))
            print("%-6s solver vs standalone: max deviation = %.4f decades (%.2f %%) at t = %.2e s;"
                  "  final solver %.4g"
                  % ("", dec2.max(), 100 * (10 ** dec2.max() - 1), ts[m2][ok2][dec2.argmax()],
                     solver[name][-1]))
# ---- end of errors

# ---- plots: one for the temperatures, one for the number densities
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
for is_temperature in (True, False):
    names = [n for n in curves if n.startswith("T") == is_temperature]
    # densities are plotted only if the paper reports them (fig. 7, 8, 9)
    if not names or (not is_temperature and not any(n in paper for n in names)):
        continue
    plt.figure(figsize=(7, 4.5))
    for i, name in enumerate(names):
        c = colors[i % len(colors)]
        if name in paper:
            plt.plot(paper[name]["t"], paper[name][name], "-", color=c, lw=2, alpha=0.4,
                     label="paper " + name)
        m = std["t"] > 0
        plt.plot(std["t"][m], std[name][m], "--", color=c, label="standalone " + name)
        if solver is not None and name in solver:
            m = solver["t"] > 0
            plt.plot(solver["t"][m][::5], solver[name][m][::5], "o", color=c, ms=3,
                     label="solver " + name)
    if is_temperature and Teq is not None:
        plt.axhline(Teq, color="k", ls=":", lw=1)
        plt.text(std["t"][1], Teq, " T_eq = %.0f K" % Teq, va="bottom", fontsize=8)
    plt.xscale("log")
    if not is_temperature:
        plt.yscale("log")
    plt.xlabel("t [s]")
    plt.ylabel("T [K]" if is_temperature else "n / n0")
    plt.title(fig if is_temperature else fig + " (number densities)")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    out = os.path.join(here, "output", fig + ("" if is_temperature else "-n") + ".png")
    plt.savefig(out, dpi=130)
    print("plot: " + os.path.relpath(out, here))
# ---- end of plots

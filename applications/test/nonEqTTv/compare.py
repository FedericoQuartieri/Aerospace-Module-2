#!/usr/bin/env python3
# confronta per una figura del paper le curve di:
#   - paper       paper-data/<fig>-<curva>.csv   (estratte dal pdf, riferimento)
#   - standalone  output/<fig>.csv               (programma 0D con Mutation++)
#   - solver      output/<fig>-solver.csv        (caso OpenFOAM, se esiste)
# stampa errore massimo e valori finali, e salva i grafici output/<fig>.png
# (temperature) e output/<fig>-n.png (densita' numeriche, dove il paper le da')
#
# uso: python3 compare.py <fig>      (es. fig3a)

import os
import sys
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# temperature di equilibrio dichiarate nel testo o nelle figure del paper
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
    # ritorna {nome colonna: array}; le righe con # sono commenti
    with open(path) as f:
        lines = [l for l in f if not l.startswith("#")]
    names = lines[0].strip().split(",")
    data = np.loadtxt(lines[1:], delimiter=",", ndmin=2)
    return {n: data[:, i] for i, n in enumerate(names)}


def read_Teq(path):
    # il programma 0D scrive in testa "# T_eq = ... K (conservazione dell'energia)"
    with open(path) as f:
        first = f.readline()
    if first.startswith("# T_eq"):
        return float(first.split("=")[1].split()[0])
    return None


std = load_csv(os.path.join(here, "output", fig + ".csv"))
Teq = read_Teq(os.path.join(here, "output", fig + ".csv"))

solver_path = os.path.join(here, "output", fig + "-solver.csv")
solver = load_csv(solver_path) if os.path.exists(solver_path) else None

# curve del paper: file fig-<curva>.csv, oppure figa-/figb- per le figure
# con due pannelli (temperature e densita' numeriche); le varianti di modello
# (-park05, -park, -QK) si confrontano con le stesse curve
base = fig
for suffix in ("-park05", "-park", "-QK"):
    base = base.replace(suffix, "")
paper = {}
for prefix in (base, base + "a", base + "b"):
    for path in sorted(glob.glob(os.path.join(here, "paper-data", prefix + "-*.csv"))):
        name = os.path.basename(path)[len(prefix) + 1:-4]
        paper[name] = load_csv(path)

curves = [n for n in std if n != "t"]

print("=== %s ===" % fig)
if Teq is not None:
    print("T_eq dalla conservazione dell'energia: %.1f K" % Teq, end="")
    if fig in PAPER_TEQ:
        print("   (paper: %.1f K, scarto %.1f K)" % (PAPER_TEQ[fig], Teq - PAPER_TEQ[fig]))
    else:
        print("   (il paper non da' un valore)")

# ---- errori rispetto al paper e valori finali, curva per curva
for name in curves:
    # le temperature si chiamano T..., il resto sono densita' numeriche n/n0
    is_temperature = name.startswith("T")
    unit = "K" if is_temperature else ""
    if name not in paper:
        print("%-6s nessuna curva del paper" % name)
        continue
    tp, yp = paper[name]["t"], paper[name][name]
    # confronto sui tempi del paper (dentro l'intervallo simulato)
    m = (tp >= std["t"][0]) & (tp <= std["t"][-1]) & (tp > 0)
    ys = np.interp(tp[m], std["t"], std[name])
    err = np.abs(ys - yp[m])
    if is_temperature:
        print("%-6s standalone vs paper: max |err| = %7.1f K (%.2f %%) a t = %.2e s;"
              "  finale (t = %.1e s): paper %.1f, standalone %.1f"
              % (name, err.max(), 100 * (err / yp[m]).max(), tp[m][err.argmax()],
                 tp[m][-1], yp[m][-1], ys[-1]))
    else:
        # densita' in scala logaritmica: lo scarto si misura in decadi,
        # |log10(standalone/paper)|, come si legge sul grafico
        ok = (ys > 0) & (yp[m] > 0)
        dec = np.abs(np.log10(ys[ok] / yp[m][ok]))
        print("%-6s standalone vs paper: max scarto = %.3f decadi (fattore %.2f) a t = %.2e s;"
              "  finale (t = %.1e s): paper %.4g, standalone %.4g"
              % (name, dec.max(), 10 ** dec.max(), tp[m][ok][dec.argmax()],
                 tp[m][-1], yp[m][-1], ys[-1]))

    if solver is not None and name in solver:
        # solver vs standalone, sui tempi del solver
        ts = solver["t"]
        m2 = (ts >= std["t"][0]) & (ts <= std["t"][-1]) & (ts > 0)
        ys2 = np.interp(ts[m2], std["t"], std[name])
        err2 = np.abs(solver[name][m2] - ys2)
        if is_temperature:
            print("%-6s solver vs standalone: max |err| = %7.1f K (%.2f %%);  finale solver %.1f K"
                  % ("", err2.max(), 100 * (err2 / np.abs(ys2)).max(), solver[name][-1]))
        else:
            ok2 = (ys2 > 0) & (solver[name][m2] > 0)
            dec2 = np.abs(np.log10(solver[name][m2][ok2] / ys2[ok2]))
            print("%-6s solver vs standalone: max scarto = %.4f decadi (%.2f %%);  finale solver %.4g"
                  % ("", dec2.max(), 100 * (10 ** dec2.max() - 1), solver[name][-1]))
# ---- fine errori

# ---- grafici: uno per le temperature, uno per le densita' numeriche
colors = plt.rcParams["axes.prop_cycle"].by_key()["color"]
for is_temperature in (True, False):
    names = [n for n in curves if n.startswith("T") == is_temperature]
    # le densita' si disegnano solo se il paper le riporta (fig. 7, 8, 9)
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
    plt.title(fig if is_temperature else fig + " (densita' numeriche)")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    out = os.path.join(here, "output", fig + ("" if is_temperature else "-n") + ".png")
    plt.savefig(out, dpi=130)
    print("grafico: " + os.path.relpath(out, here))
# ---- fine grafici

#!/usr/bin/env python3
# grafici della relazione: stesse curve di applications/test/nonEqTTv/compare.py
# (riferimento estratto dal pdf, programma 0D, solver), ma senza titolo, con
# la legenda in italiano e con nomi che non richiamano la numerazione del paper
#
# uso: python3 genera-grafici.py      (scrive i png in questa cartella)

import os
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

qui = os.path.dirname(os.path.abspath(__file__))
test = os.path.join(qui, "..", "..", "applications", "test", "nonEqTTv")

# nome del grafico -> (risultato del progetto, curve del paper, cosa disegnare)
GRAFICI = {
    "azoto-riscaldamento":   ("fig3a",     "fig3a", "T"),
    "azoto-raffreddamento":  ("fig3b",     "fig3b", "T"),
    "azoto-30000K-senza-el": ("fig4-noEl", "fig4-noEl", "T"),
    "azoto-30000K-con-el":   ("fig4-el",   "fig4-el", "T"),
    "azoto-atomico":         ("fig5",      "fig5", "T"),
    "azoto-ossigeno-VV":     ("fig6-VV",   "fig6-VV", "T"),
    "reagente-temperature":  ("fig7",      "fig7", "T"),
    "reagente-densita":      ("fig7",      "fig7", "n"),
    "aria-temperatura":      ("fig9-QK",   "fig9", "T"),
    "aria-densita":          ("fig9-QK",   "fig9", "n"),
}

NOMI = {
    "Ttr": r"$T_{tr}$", "Tv": r"$T_v$", "T": r"$T$",
    "Tv_N2": r"$T_{v,\mathrm{N_2}}$", "Tv_O2": r"$T_{v,\mathrm{O_2}}$",
    "N2": r"$\mathrm{N_2}$", "O2": r"$\mathrm{O_2}$", "NO": r"$\mathrm{NO}$",
    "N": r"$\mathrm{N}$", "O": r"$\mathrm{O}$",
}


def carica(path):
    # {colonna: valori}; le righe con # sono commenti
    with open(path) as f:
        righe = [r for r in f if not r.startswith("#")]
    nomi = righe[0].strip().split(",")
    dati = np.loadtxt(righe[1:], delimiter=",", ndmin=2)
    return {n: dati[:, i] for i, n in enumerate(nomi)}


def T_equilibrio(path):
    with open(path) as f:
        prima = f.readline()
    return float(prima.split("=")[1].split()[0]) if prima.startswith("# T_eq") else None


colori = plt.rcParams["axes.prop_cycle"].by_key()["color"]

for nome, (ris, base, tipo) in GRAFICI.items():
    std = carica(os.path.join(test, "output", ris + ".csv"))
    p_solver = os.path.join(test, "output", ris + "-solver.csv")
    solver = carica(p_solver) if os.path.exists(p_solver) else None
    paper = {}
    for pre in (base, base + "a", base + "b"):
        for p in sorted(glob.glob(os.path.join(test, "paper-data", pre + "-*.csv"))):
            paper[os.path.basename(p)[len(pre) + 1:-4]] = carica(p)

    curve = [c for c in std if c != "t" and c.startswith("T") == (tipo == "T")]
    plt.figure(figsize=(7, 4.5))
    for i, c in enumerate(curve):
        col = colori[i % len(colori)]
        if c in paper:
            plt.plot(paper[c]["t"], paper[c][c], "-", color=col, lw=2, alpha=0.4,
                     label="riferimento " + NOMI[c])
        m = std["t"] > 0
        plt.plot(std["t"][m], std[c][m], "--", color=col, label="programma 0D " + NOMI[c])
        if solver is not None and c in solver:
            m = solver["t"] > 0
            plt.plot(solver["t"][m][::5], solver[c][m][::5], "o", color=col, ms=3,
                     label="solver " + NOMI[c])
    Teq = T_equilibrio(os.path.join(test, "output", ris + ".csv"))
    if tipo == "T" and Teq is not None:
        plt.axhline(Teq, color="k", ls=":", lw=1)
        plt.text(std["t"][1], Teq, " equilibrio: %.0f K" % Teq, va="bottom", fontsize=8)
    plt.xscale("log")
    if tipo == "n":
        plt.yscale("log")
    plt.xlabel("t [s]")
    plt.ylabel("T [K]" if tipo == "T" else r"$n/n_0$")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(qui, nome + ".png"), dpi=130)
    plt.close()
    print("scritto " + nome + ".png")

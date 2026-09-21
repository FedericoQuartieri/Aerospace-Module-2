#!/usr/bin/env python3
# grafici della relazione: stesse curve di applications/test/nonEqTTv/compare.py
# (riferimento estratto dal pdf, programma 0D, solver), ma senza titolo, con
# la legenda in italiano e con nomi che non richiamano la numerazione del paper
#
# uso: python3 genera-grafici.py        (scrive i png in questa cartella)
#      python3 genera-grafici.py --en   (legenda in inglese, png in ../images-en,
#                                        per la versione inglese main.tex)

import os
import sys
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

qui = os.path.dirname(os.path.abspath(__file__))
test = os.path.join(qui, "..", "..", "applications", "test", "nonEqTTv")

# lingua della legenda e cartella di uscita
inglese = "--en" in sys.argv[1:]
uscita = os.path.join(qui, "..", "images-en") if inglese else qui
os.makedirs(uscita, exist_ok=True)
TESTI = {
    "riferimento": "reference" if inglese else "riferimento",
    "programma":   "0D program" if inglese else "programma 0D",
    "solver":      "solver",
    "equilibrio":  "equilibrium" if inglese else "equilibrio",
}

# etichette delle varianti con le due formule del tempo di rilassamento V-T
TAU_PAPER = r"$\tau$ reference study" if inglese else r"$\tau$ del riferimento"
TAU_MPP = r"$\tau$ Mutation++" if inglese else r"$\tau$ di Mutation++"

# nome del grafico -> (risultato, curve del paper, cosa disegnare, variante)
# la variante, se c'e', e' (altro risultato, etichetta del primo, etichetta della
# variante): le due curve finiscono nello stesso grafico e il solver non si disegna
GRAFICI = {
    "nitrogen-heating":   ("fig3a",     "fig3a", "T"),
    "nitrogen-cooling":  ("fig3b",     "fig3b", "T"),
    "nitrogen-hot-no-electronic": ("fig4-noEl", "fig4-noEl", "T"),
    "nitrogen-hot-electronic":   ("fig4-el",   "fig4-el", "T"),
    "nitrogen-atomic":         ("fig5",      "fig5", "T"),
    "nitrogen-oxygen-no-vv":  ("fig6-noVV", "fig6-noVV", "T"),
    "nitrogen-oxygen-vv":     ("fig6-VV",   "fig6-VV", "T"),
    "reacting-noneq-temperatures": ("fig7", "fig7", "T"),
    "reacting-noneq-densities":    ("fig7", "fig7", "n"),
    "reacting-eq-temperatures":    ("fig8", "fig8", "T"),
    "reacting-eq-densities":       ("fig8", "fig8", "n"),
    "air-temperature":      ("fig9-QK",   "fig9", "T"),
    "air-densities":          ("fig9-QK",   "fig9", "n"),
    "park-exponent-densities": ("fig7", "fig7", "n",
                                ("fig7-park05", "$a = 0.7$", "$a = 0.5$")),
    "air-densities-park-rates": ("fig9-park", "fig9", "n"),
    "relaxation-time-nitrogen-atomic": ("fig5", "fig5", "T",
                                        ("fig5-tauMpp", TAU_PAPER, TAU_MPP)),
    "relaxation-time-reacting": ("fig7", "fig7", "T",
                                 ("fig7-tauMpp", TAU_PAPER, TAU_MPP)),
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

for nome, voce in GRAFICI.items():
    ris, base, tipo = voce[:3]
    variante = voce[3] if len(voce) > 3 else None
    if variante and not inglese:
        # virgola decimale nelle etichette della versione italiana
        variante = (variante[0], variante[1].replace(".", "{,}"),
                    variante[2].replace(".", "{,}"))
    std = carica(os.path.join(test, "output", ris + ".csv"))
    alt = carica(os.path.join(test, "output", variante[0] + ".csv")) if variante else None
    p_solver = os.path.join(test, "output", ris + "-solver.csv")
    solver = carica(p_solver) if os.path.exists(p_solver) and not variante else None
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
                     label=TESTI["riferimento"] + " " + NOMI[c])
        suff = " (" + variante[1] + ")" if variante else ""
        m = std["t"] > 0
        plt.plot(std["t"][m], std[c][m], "--", color=col,
                 label=TESTI["programma"] + " " + NOMI[c] + suff)
        if alt is not None and c in alt:
            m = alt["t"] > 0
            plt.plot(alt["t"][m], alt[c][m], ":", color=col,
                     label=TESTI["programma"] + " " + NOMI[c] + " (" + variante[2] + ")")
        if solver is not None and c in solver:
            m = solver["t"] > 0
            plt.plot(solver["t"][m][::5], solver[c][m][::5], "o", color=col, ms=3,
                     label=TESTI["solver"] + " " + NOMI[c])
    Teq = T_equilibrio(os.path.join(test, "output", ris + ".csv"))
    if tipo == "T" and Teq is not None:
        plt.axhline(Teq, color="k", ls=":", lw=1)
        plt.text(std["t"][1], Teq, " %s: %.0f K" % (TESTI["equilibrio"], Teq),
                 va="bottom", fontsize=8)
    plt.xscale("log")
    if tipo == "n":
        plt.yscale("log")
    plt.xlabel("t [s]")
    plt.ylabel("T [K]" if tipo == "T" else r"$n/n_0$")
    plt.legend(fontsize=8)
    plt.grid(alpha=0.3)
    plt.tight_layout()
    plt.savefig(os.path.join(uscita, nome + ".png"), dpi=130)
    plt.close()
    print("scritto " + nome + ".png")

"""Quantitative comparison against Figure 2 of the Part Two paper.

Prerequisites (run in the case dir, on a converged solution):
    foamPostProcess -solver shockThermo -func wallHeatFlux   -time <t>
    foamPostProcess -solver shockThermo -func wallShearStress -time <t>
    references/fig2*.csv from references/digitize-fig2.py

Wall quantities (paper eqs. 30-32, all with NOMINAL free-stream values -
the guard-rail in postProcess-cone.py already checks that the simulated
free-stream matches them):
    Cp = (p - p_inf)/(0.5 rho_inf U_inf^2)
    Cf = |tau_w|   /(0.5 rho_inf U_inf^2)
    St = |q_w|     /(0.5 rho_inf U_inf^3)

Outputs: fig2-comparison.png (3 panels: Cp, Cf, St vs axial distance)
         fig2-stagnation-comparison.png (panel a: T/Tinf profiles)
         printed deviation metrics vs the digitised hy2Foam curves.
"""

import os
import re
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RHO_INF = 5.113e-4
U_INF = 2764.5
P_INF = 21.9139
T_INF = 144.4
RN = 6.35e-3
R_N2 = 8.31446261815324/28.0134e-3
Q_INF = 0.5*RHO_INF*U_INF**2

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, "references")

times = sorted(
    (d for d in os.listdir(HERE)
     if re.fullmatch(r"[0-9.e+-]+", d) and float(d) > 0),
    key=float
)
tdir = times[-1]
print(f"using time {tdir}")

# ------------------------------------------------------- mesh: wall faces --
def strip_comments(txt):
    txt = re.sub(r"/\*.*?\*/", "", txt, flags=re.S)
    return re.sub(r"//.*", "", txt)

pm = os.path.join(HERE, "constant", "polyMesh")

btxt = strip_comments(open(os.path.join(pm, "boundary")).read())
m = re.search(r"wall\s*\{[^}]*?nFaces\s+(\d+);[^}]*?startFace\s+(\d+);",
              btxt, re.S)
nF, startF = int(m.group(1)), int(m.group(2))

ptxt = strip_comments(open(os.path.join(pm, "points")).read())
m = re.search(r"\n(\d+)\s*\(\s*(.*?)\)\s*$", ptxt, re.S)
pts = re.findall(r"\(([^)]+)\)", m.group(2))
points = np.array([[float(q) for q in p.split()] for p in pts])

ftxt = strip_comments(open(os.path.join(pm, "faces")).read())
m = re.search(r"\n(\d+)\s*\(\s*(.*?)\)\s*$", ftxt, re.S)
faceItems = re.findall(r"\d+\(([^)]+)\)", m.group(2))
faces = [list(map(int, it.split())) for it in faceItems]

otxt = strip_comments(open(os.path.join(pm, "owner")).read())
m = re.search(r"\n(\d+)\s*\(\s*(.*?)\)\s*$", otxt, re.S)
owner = np.array([int(x) for x in m.group(2).split()])

wallCentres = np.array([
    points[faces[startF + i]].mean(axis=0) for i in range(nF)
])
wallOwner = owner[startF:startF + nF]

# axial distance from the stagnation point [cm]
xAx = (wallCentres[:, 0] + RN)*100.0
order = np.argsort(xAx)

# ----------------------------------------------------------- field reads --
def read_internal(name):
    txt = open(os.path.join(HERE, tdir, name)).read()
    m = re.search(r"internalField\s+nonuniform\s+List<(scalar|vector)>\s*"
                  r"\d+\s*\((.*?)\)\s*;", txt, re.S)
    if m.group(1) == "vector":
        rows = re.findall(r"\(([^)]+)\)", m.group(2))
        return np.array([[float(q) for q in r.split()] for r in rows])
    return np.array([float(q) for q in m.group(2).split()])

def read_wall_boundary(name):
    txt = open(os.path.join(HERE, tdir, name)).read()
    m = re.search(r"boundaryField.*?wall\s*\{(.*?)\n    \}", txt, re.S)
    blk = m.group(1)
    mv = re.search(r"value\s+nonuniform\s+List<(scalar|vector)>\s*"
                   r"\d+\s*\((.*?)\)\s*;", blk, re.S)
    if mv is None:
        mu = re.search(r"value\s+uniform\s+(\(?[^;]+?\)?)\s*;", blk)
        val = mu.group(1)
        if "(" in val:
            v = [float(q) for q in val.strip("()").split()]
            return np.tile(v, (nF, 1))
        return np.full(nF, float(val))
    if mv.group(1) == "vector":
        rows = re.findall(r"\(([^)]+)\)", mv.group(2))
        return np.array([[float(q) for q in r.split()] for r in rows])
    return np.array([float(q) for q in mv.group(2).split()])

p_cells = read_internal("p")
pWall = p_cells[wallOwner]              # dp/dn = 0 at the wall
qWall = read_wall_boundary("wallHeatFlux")
tauWall = read_wall_boundary("wallShearStress")

Cp = (pWall - P_INF)/Q_INF
Cf = np.linalg.norm(np.atleast_2d(tauWall), axis=-1)/Q_INF
St = np.abs(qWall)/(0.5*RHO_INF*U_INF**3)

# ------------------------------------------------------------- references --
def ref(name):
    path = os.path.join(REF, name)
    if not os.path.exists(path):
        return None
    return np.loadtxt(path, delimiter=",", skiprows=1)

def plot_refs(ax, panel):
    styles = {
        "hy2foam": dict(c="k", lw=1.2, ls="-", label="hy2Foam (paper)"),
        "hy2foam10um": dict(c="gray", lw=1.0, ls="-",
                            label="hy2Foam 10um (paper)"),
        "michigan": dict(c="g", marker=".", ls="none", ms=3,
                         label="CFD Michigan (paper)"),
        "dsmc": dict(c="b", marker="^", ls="none", ms=4, mfc="none",
                     label="DSMC MONACO (paper)"),
        "experiments": dict(c="m", marker="s", ls="none", ms=5, mfc="none",
                            label="Experiments CUBRC (paper)"),
    }
    for name, st in styles.items():
        d = ref(f"fig2{panel}-{name}.csv")
        if d is None:
            continue
        if name in ("hy2foam", "hy2foam10um"):
            o = np.argsort(d[:, 0])
            ax.plot(d[o, 0], d[o, 1], **st)
        else:
            ax.plot(d[:, 0], d[:, 1], **st)

# ------------------------------------------------------- surface figure ---
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
ours = dict(c="r", lw=1.8, ls="--", label="shockThermo (this work)")

for ax, (panel, data, ylab, ymax) in zip(axes, [
    ("d", Cp, "pressure coefficient", 1.0),
    ("e", Cf, "friction coefficient", 0.15),
    ("f", St, "Stanton number", 0.2),
]):
    plot_refs(ax, panel)
    ax.plot(xAx[order], data[order], **ours)
    ax.set_xlim(0, 4)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("axial distance from stagnation point [cm]")
    ax.set_ylabel(ylab)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)

fig.suptitle(f"Mach 11.3 blunted cone vs paper Fig. 2 (t = {tdir})")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig2-comparison.png"), dpi=150)
print("saved fig2-comparison.png")

# --------------------------------------------------- deviation metrics ----
print("\nscarti vs hy2Foam digitalizzata (0.5-4 cm, zone regolari):")
for panel, data, label in (("d", Cp, "Cp"), ("e", Cf, "Cf"), ("f", St, "St")):
    d = ref(f"fig2{panel}-hy2foam.csv")
    if d is None:
        continue
    o = np.argsort(d[:, 0])
    xr, yr = d[o, 0], d[o, 1]
    mask = (xAx[order] > 0.5) & (xAx[order] < min(4.0, xr.max()))
    yi = np.interp(xAx[order][mask], xr, yr)
    yo = data[order][mask]
    err = np.abs(yo - yi)/np.maximum(np.abs(yi), 1e-12)
    print(f"  {label}: media {err.mean()*100:.1f}%, max {err.max()*100:.1f}%")

# ------------------------------------------------ stagnation line (a) -----
def read_vec_or_scal(name):
    return read_internal(name)

import subprocess
if not os.path.exists(os.path.join(HERE, tdir, "C")):
    subprocess.run(["foamPostProcess", "-func", "writeCellCentres",
                    "-time", tdir], cwd=HERE, capture_output=True)

C = read_internal("C")
T = read_internal("T")
Tve = read_internal("Tve")
x, y = C[:, 0], C[:, 1]

axisBand = np.unique(y[np.abs(x) < 2.9*RN])
yTol = np.sort(np.abs(axisBand))[0]*3
line = (np.abs(y) <= max(yTol, 1e-4)) & (x < -RN*0.99)
xs = (x[line] + RN)*1000.0  # stagnation line position [mm], wall at 0
o = np.argsort(xs)

fig2, ax = plt.subplots(figsize=(8, 6))
plot_refs(ax, "a")
ax.plot(xs[o], T[line][o]/T_INF, "r--", lw=1.8, label="shockThermo T")
ax.plot(xs[o], Tve[line][o]/T_INF, "r:", lw=1.8, label="shockThermo Tve")
ax.set_xlim(-3, 0)
ax.set_ylim(0, 35)
ax.set_xlabel("stagnation line position [mm]")
ax.set_ylabel("normalised temperature T/T_inf")
ax.grid(alpha=0.3)
ax.legend(fontsize=8)
fig2.suptitle("Stagnation line vs paper Fig. 2a")
fig2.tight_layout()
fig2.savefig(os.path.join(HERE, "fig2-stagnation-comparison.png"), dpi=150)
print("saved fig2-stagnation-comparison.png")

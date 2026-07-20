"""Quantitative comparison against Figure 5 of the Part Two paper
(Mach 20 reacting cylinder). Overlays our run-3-equivalent solution
(Park TTv + Park rates) on the digitised references in references/.

ROBUST wall extraction (see the M7 note in the milestone doc): the
solver develops a bounded odd-even checkerboard in the SINGLE wall-
adjacent cell (extreme near-wall aspect ratio ~1700:1 + central flux +
slip/jump walls). The interior is clean. Surface quantities are therefore
read the way a boundary-layer measurement would:
  - surface pressure: median over a clean radial band a few cells off the
    wall (dp/dn = 0 through a boundary layer, so this IS the wall value);
  - skin friction: tangential component of wallShearStress only (the
    checkerboard leaks a spurious wall-NORMAL component into the viscous
    stress; true shear is tangential by definition).
Heat flux stays the shakiest quantity (q ~ dT/dn is genuinely near-wall).

Data dir defaults to this case; override with argv[1] to point at a
reconstructed run copied elsewhere (e.g. the cluster tarball).

Outputs: fig5-surface-comparison.png, fig5-stagnation-comparison.png
         + printed C_D / C_H and deviation notes.
"""

import os
import re
import sys
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

RHO_INF = 1.363e-5
U_INF = 6047.0
P_INF = 0.89
T_INF = 220.0
R_CYL = 1.0
R_N2 = 8.31446261815324/28.0134e-3
GAMMA = 1.4
Q_INF = 0.5*RHO_INF*U_INF**2
NA = 6.02214076e23
M_N2 = 28.0134e-3
M_N = 14.0067e-3

HERE = os.path.dirname(os.path.abspath(__file__))
REF = os.path.join(HERE, "references")
DATA = sys.argv[1] if len(sys.argv) > 1 else HERE

times = sorted(
    (d for d in os.listdir(DATA)
     if re.fullmatch(r"[0-9.e+-]+", d) and float(d) > 0),
    key=float
)
tdir = times[-1]
print(f"dati da {DATA}, tempo {tdir}")

# ------------------------------------------------------------ mesh parsing --
def strip_comments(txt):
    txt = re.sub(r"/\*.*?\*/", "", txt, flags=re.S)
    return re.sub(r"//.*", "", txt)

pm = os.path.join(DATA, "constant", "polyMesh")

btxt = strip_comments(open(os.path.join(pm, "boundary")).read())
m = re.search(r"wall\s*\{[^}]*?nFaces\s+(\d+);[^}]*?startFace\s+(\d+);",
              btxt, re.S)
nF, startF = int(m.group(1)), int(m.group(2))

ptxt = strip_comments(open(os.path.join(pm, "points")).read())
m = re.search(r"\n(\d+)\s*\(\s*(.*?)\)\s*$", ptxt, re.S)
points = np.array([[float(q) for q in p.split()]
                   for p in re.findall(r"\(([^)]+)\)", m.group(2))])

ftxt = strip_comments(open(os.path.join(pm, "faces")).read())
m = re.search(r"\n(\d+)\s*\(\s*(.*?)\)\s*$", ftxt, re.S)
faces = [list(map(int, it.split()))
         for it in re.findall(r"\d+\(([^)]+)\)", m.group(2))]

otxt = strip_comments(open(os.path.join(pm, "owner")).read())
m = re.search(r"\n(\d+)\s*\(\s*(.*?)\)\s*$", otxt, re.S)
owner = np.array([int(x) for x in m.group(2).split()])

wallFaces = [faces[startF + i] for i in range(nF)]
wallCentres = np.array([points[f].mean(axis=0) for f in wallFaces])
wallOwner = owner[startF:startF + nF]

def areaVec(f):
    P = points[f]
    n = np.zeros(3)
    for k in range(len(P)):
        n += np.cross(P[k], P[(k + 1) % len(P)])
    return 0.5*n

wallA = np.array([areaVec(f) for f in wallFaces])
depth = 0.1
# theta from the stagnation direction (0 = front, paper convention)
thetaW = np.degrees(np.arctan2(wallCentres[:, 1], -wallCentres[:, 0]))

# ------------------------------------------------------------- field reads --
def internal(name):
    txt = open(os.path.join(DATA, tdir, name)).read()
    m = re.search(r"internalField\s+nonuniform\s+List<(scalar|vector)>\s*"
                  r"\d+\s*\((.*?)\)\s*;", txt, re.S)
    if m.group(1) == "vector":
        return np.array([[float(q) for q in r.split()]
                         for r in re.findall(r"\(([^)]+)\)", m.group(2))])
    return np.array([float(q) for q in m.group(2).split()])

def wall_bnd(name):
    txt = open(os.path.join(DATA, tdir, name)).read()
    m = re.search(r"boundaryField.*?wall\s*\{(.*?)\n    \}", txt, re.S)
    mv = re.search(r"value\s+nonuniform\s+List<(scalar|vector)>\s*"
                   r"\d+\s*\((.*?)\)\s*;", m.group(1), re.S)
    if mv.group(1) == "vector":
        return np.array([[float(q) for q in r.split()]
                         for r in re.findall(r"\(([^)]+)\)", mv.group(2))])
    return np.array([float(q) for q in mv.group(2).split()])

C = internal("C")
p = internal("p")
T = internal("T")
Tve = internal("Tve")
U = internal("U")
YN = internal("N")
rho = internal("rho")
x, y = C[:, 0], C[:, 1]
rC = np.sqrt(x**2 + y**2)
thC = np.degrees(np.arctan2(y, -x))

# -------------------------------------------- robust surface extraction ----
CpS = np.zeros(nF)
qS = np.zeros(nF)
kappa_air = 0.0  # not used; heat flux taken from wallHeatFlux field
for i in range(nF):
    sel = np.abs(thC - thetaW[i]) < 0.9
    o = np.argsort(rC[sel])
    band = p[sel][o][2:9]           # clean cells 3..9 off the wall
    CpS[i] = (np.median(band) - P_INF)/Q_INF

tau = wall_bnd("wallShearStress")
nhat = wallCentres.copy()
nhat[:, 2] = 0
nhat /= np.linalg.norm(nhat[:, :2], axis=1, keepdims=True)
tanShear = tau - np.sum(tau*nhat, axis=1)[:, None]*nhat
CfS = np.linalg.norm(tanShear, axis=1)/Q_INF

qw = wall_bnd("wallHeatFlux")
qS = np.abs(qw)/1e4  # W/cm^2

order = np.argsort(thetaW)

# --------------------------------------------------------------- C_D / C_H --
CdP = 2.0*np.sum(CpS*Q_INF*wallA[:, 0])/(Q_INF*2*R_CYL*depth)
Amag = np.linalg.norm(wallA, axis=1)
CdF = -2.0*np.sum(tanShear[:, 0]*Amag)/(Q_INF*2*R_CYL*depth)
Ch = np.sum(np.abs(qw)*Amag)/depth/1000.0
print(f"C_D = {CdP + CdF:.3f} (pressione {CdP:.3f} + attrito {CdF:.3f}; "
      f"paper 1.304, DSMC 1.284, Newton 1.333)")
print(f"C_H (mezzo cilindro, per m) = {Ch:.1f} kW "
      f"(paper run3 88.1, DSMC 63.3; la meno affidabile: q~dT/dn e' near-wall)")

# ------------------------------------------------------------- references ---
def ref(name):
    path = os.path.join(REF, name)
    return np.loadtxt(path, delimiter=",", skiprows=1) if os.path.exists(path) else None

def plot_refs(ax, panel):
    for nm, st in (
        ("run3", dict(c="b", lw=1.3, ls="-", label="hy2Foam run3 (paper)")),
        ("dsmc", dict(c="k", marker="+", ls="none", ms=5,
                      label="dsmcFoam (paper)")),
    ):
        d = ref(f"fig5{panel}-{nm}.csv")
        if d is None:
            continue
        if nm == "run3":
            o = np.argsort(d[:, 0])
            ax.plot(d[o, 0], d[o, 1], **st)
        else:
            ax.plot(d[:, 0], d[:, 1], **st)

# ------------------------------------------------------- surface figure ----
fig, axes = plt.subplots(1, 3, figsize=(16, 5))
ours = dict(c="r", lw=1.8, ls="--", label="shockThermo (this work)")
for ax, (panel, data, ylab, ymax) in zip(axes, [
    ("d", CpS, "pressure coefficient", 2.0),
    ("e", CfS, "skin-friction coefficient", 0.06),
    ("f", qS, "surface heat flux [W/cm2]", 15.0),
]):
    plot_refs(ax, panel)
    ax.plot(thetaW[order], data[order], **ours)
    ax.set_xlim(0, 180)
    ax.set_ylim(0, ymax)
    ax.set_xlabel("theta [deg]")
    ax.set_ylabel(ylab)
    ax.grid(alpha=0.3)
    ax.legend(fontsize=7)
fig.suptitle(f"Mach 20 reacting cylinder vs paper Fig. 5d-f (t = {tdir})")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "fig5-surface-comparison.png"), dpi=150)
print("saved fig5-surface-comparison.png")

# ---------------------------------------------------- stagnation figure ----
dTh = thC[(x < -R_CYL) & (thC < 10) & (thC > 0)].min()*1.001
line = (np.abs(thC) <= 1.5*dTh) & (x < -R_CYL*0.999)
xs = x[line] + R_CYL
o = np.argsort(xs)
Mline = np.linalg.norm(U[line][o], axis=1)/np.sqrt(GAMMA*R_N2*T[line][o])
nN2 = np.maximum(rho*(1 - YN)/M_N2*NA, 1e10)
nN = np.maximum(rho*YN/M_N*NA, 1e10)

fig2, axes2 = plt.subplots(1, 3, figsize=(16, 5))
axes2[0].plot(xs[o], Mline, "r--", lw=1.8, label="shockThermo")
d = ref("fig5a-run3.csv")
if d is not None:
    axes2[0].plot(d[:, 0] + R_CYL, d[:, 1], "b-", lw=1.2,
                  label="hy2Foam run3")
axes2[0].set_ylabel("Mach"); axes2[0].set_ylim(0, 21)
axes2[1].plot(xs[o], T[line][o]/1000, "r--", lw=1.8, label="T")
axes2[1].plot(xs[o], Tve[line][o]/1000, "r:", lw=1.8, label="Tve")
db = ref("fig5b-run3.csv")
if db is not None:
    axes2[1].plot(db[:, 0] + R_CYL, db[:, 1], "b-", lw=1.0,
                  label="hy2Foam run3")
axes2[1].set_ylabel("temperature [kK]"); axes2[1].set_ylim(0, 16)
axes2[2].semilogy(xs[o], nN2[line][o], "r--", lw=1.8, label="N2")
axes2[2].semilogy(xs[o], nN[line][o], "g--", lw=1.8, label="N")
dc = ref("fig5c-run3.csv")
if dc is not None:
    axes2[2].semilogy(dc[:, 0] + R_CYL, dc[:, 1], "b-", lw=1.0,
                      label="hy2Foam N2")
axes2[2].set_ylabel("number density [1/m3]"); axes2[2].set_ylim(1e18, 1e22)
for ax in axes2:
    ax.set_xlim(-0.5, 0)
    ax.set_xlabel("stagnation line x+R [m]")
    ax.grid(alpha=0.3)
    ax.legend(fontsize=8)
fig2.suptitle(f"Mach 20 reacting cylinder vs paper Fig. 5a-c (t = {tdir})")
fig2.tight_layout()
fig2.savefig(os.path.join(HERE, "fig5-stagnation-comparison.png"), dpi=150)
print("saved fig5-stagnation-comparison.png")

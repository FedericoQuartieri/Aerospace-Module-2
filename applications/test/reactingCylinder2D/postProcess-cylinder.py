"""Post-processing for the Mach 20 reacting cylinder (Part Two, Sec. 3.2).

Paper targets:
  - shock standoff ~ 0.25 m on the stagnation line (Sec. 3.2, Fig. 5a-c;
    ~5 cm closer to the body than the non-reacting case)
  - drag coefficient C_D = 1.304 (Table 5, run 3; Newtonian ~ 4/3)
  - integrated heat flux C_H = 88.1 kW (Table 5, run 3 - the paper itself
    notes the Park combination overpredicts dsmcFoam, 63.3 kW, by 39%)

Guard-rail (M6 lesson): measure the EFFECTIVE simulated free stream and
compare theory at the effective Mach, warning on nominal mismatch.

Needs wallHeatFlux/wallShearStress fields for C_D/C_H (Allrun generates
them via foamPostProcess -solver shockThermo); the script degrades
gracefully to pressure-only C_D if they are missing.

Outputs: cylinder-stagnation.png, cylinder-surface.png + printed summary.
"""

import os
import re
import subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# free stream, Table 2 of the paper
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

times = sorted(
    (d for d in os.listdir(HERE)
     if re.fullmatch(r"[0-9.e+-]+", d) and float(d) > 0),
    key=float
)
tdir = times[-1]
nt = len(times)
print(f"t = {tdir} s ({nt} snapshots)")

# ------------------------------------------------- convergence guard-rail --
# LTS verso lo stazionario: la soluzione e' credibile solo se i campi non
# cambiano piu' tra gli ultimi due snapshot (lezione del run locale: a
# 30k step il transitorio M20 e' ancora vivo e Cp/Cf oscillano).
def _read_scalar(td, name):
    txt = open(os.path.join(HERE, td, name)).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\d+\s*"
                  r"\((.*?)\)\s*;", txt, re.S)
    return np.array([float(q) for q in m.group(1).split()])

if nt >= 2:
    pPrev = _read_scalar(times[-2], "p")
    pLast = _read_scalar(times[-1], "p")
    dp = np.abs(pLast - pPrev)/np.maximum(np.abs(pPrev), 1e-3)
    print(f"convergenza: dp/p tra {times[-2]} e {times[-1]}: "
          f"media {dp.mean()*100:.2f}%, max {dp.max()*100:.1f}%")
    if dp.mean() > 0.01:
        print("  *** WARNING: NON a regime (media > 1%): prolungare il run "
              "prima di fidarsi di Cp/Cf/C_D/C_H ***")

# ------------------------------------------------------------ mesh parsing --
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

wallFaces = [faces[startF + i] for i in range(nF)]
wallCentres = np.array([points[f].mean(axis=0) for f in wallFaces])
wallOwner = owner[startF:startF + nF]

def faceAreaVector(f):
    # Newell's formula: area vector pointing out of the owner cell,
    # i.e. INTO the cylinder for wall faces
    p = points[f]
    n = np.zeros(3)
    for i in range(len(p)):
        a, b = p[i], p[(i + 1) % len(p)]
        n += np.cross(a, b)
    return 0.5*n

wallAreas = np.array([faceAreaVector(f) for f in wallFaces])
depth = 0.1  # planar mesh thickness (makeBlockMeshDict.py)

# theta from the stagnation direction (paper convention: 0 = front)
thetaGeo = np.degrees(np.arctan2(wallCentres[:, 1], wallCentres[:, 0]))
thetaPaper = 180.0 - thetaGeo
orderTh = np.argsort(thetaPaper)

# ------------------------------------------------------------- field reads --
def read_internal(name):
    txt = open(os.path.join(HERE, tdir, name)).read()
    m = re.search(r"internalField\s+nonuniform\s+List<(scalar|vector)>\s*"
                  r"\d+\s*\((.*?)\)\s*;", txt, re.S)
    if m is None:
        mu = re.search(r"internalField\s+uniform\s+(\(?[^;]+?\)?)\s*;", txt)
        val = mu.group(1)
        if "(" in val:
            return np.tile([float(q) for q in val.strip("()").split()],
                           (1, 1))
        return np.array([float(val)])
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

if not os.path.exists(os.path.join(HERE, tdir, "C")):
    subprocess.run(["foamPostProcess", "-func", "writeCellCentres",
                    "-time", tdir], cwd=HERE, capture_output=True)

C = read_internal("C")
T = read_internal("T")
Tve = read_internal("Tve")
p = read_internal("p")
rho = read_internal("rho") if os.path.exists(
    os.path.join(HERE, tdir, "rho")) else p/(R_N2*T)
U = read_internal("U")
YN2 = read_internal("N2")
YN = read_internal("N")
x, y = C[:, 0], C[:, 1]

# ------------------------------------------- effective free stream (M6!) ---
fs = x < -1.6
T_eff = T[fs].mean()
p_eff = p[fs].mean()
U_eff = np.linalg.norm(U[fs], axis=1).mean()
M_eff = U_eff/np.sqrt(GAMMA*R_N2*T_eff)
print(f"free-stream EFFETTIVO: T={T_eff:.1f} K (nominale {T_INF}), "
      f"p={p_eff:.3f} Pa, M={M_eff:.2f} (nominale 20.0)")
if abs(T_eff - T_INF) > 2:
    print("  *** WARNING: free-stream simulato != nominale (clamp? BC?) ***")

# ------------------------------------------------------- stagnation line ---
# angular selection: the first theta-row of cells next to the -x axis.
# A fixed |y| band fails on a polar mesh (it runs out of cells at small
# radius and the "standoff" becomes the band edge, not the shock).
thC = np.degrees(np.arctan2(np.abs(y), -x))   # 0 on the upstream axis
dTh = thC[(x < -R_CYL) & (thC < 10)].min()*1.001
line = (thC <= 1.5*dTh) & (x < -R_CYL*0.999)
xs = x[line] + R_CYL  # 0 at the wall, negative upstream
o = np.argsort(xs)
xso, Tso, Tveo = xs[o], T[line][o], Tve[line][o]

# standoff: half-rise of T between free stream and peak
Tmax = Tso.max()
half = T_INF + 0.5*(Tmax - T_INF)
ix = np.where(Tso > half)[0]
standoff = -xso[ix[0]] if len(ix) else float("nan")
print(f"shock standoff (linea di ristagno) = {standoff:.3f} m "
      f"(paper: ~0.25 m, entrambi i codici)")

# ------------------------------------------------------------- C_D / C_H ---
# ROBUST (vedi il commit "Tve wall" e la doc M7): il solver sviluppa un
# checkerboard odd-even confinato alla SOLA cella a parete (aspect ratio
# ~1700:1 + flusso centrale + BC slip/jump); l'interno e' pulito. Le
# grandezze di parete si leggono come le misurerebbe uno strato limite:
# pressione dalla banda radiale pulita 3-9 celle fuori parete (dp/dn=0
# nello strato limite -> e' il valore di parete), attrito dalla sola
# componente tangenziale (il checkerboard inietta una spuria componente
# NORMALE nello sforzo viscoso).
rCC = np.sqrt(x**2 + y**2)
thCC = np.degrees(np.arctan2(y, -x))
thetaWgeo = np.degrees(np.arctan2(wallCentres[:, 1], -wallCentres[:, 0]))
CpS = np.zeros(nF)
for i in range(nF):
    sel = np.abs(thCC - thetaWgeo[i]) < 0.9
    ro = np.argsort(rCC[sel])          # NON riusare 'o': e' della linea di ristagno
    CpS[i] = (np.median(p[sel][ro][2:9]) - P_INF)/Q_INF

# stagnation Cp (robusto: faccia a theta minimo) + Rayleigh a M effettivo
CpStag = CpS[np.argmin(thetaWgeo)]
g = GAMMA
Mr = M_eff
pratio = ((g + 1)**2*Mr**2/(4*g*Mr**2 - 2*(g - 1)))**(g/(g - 1)) \
    *(1 - g + 2*g*Mr**2)/(g + 1)
CpRay = (pratio - 1)*P_INF/Q_INF
print(f"stagnation Cp = {CpStag:.3f}  (Rayleigh pitot ideale frozen a "
      f"M={Mr:.2f}: {CpRay:.3f}; reagente atteso leggermente diverso)")

Dp = np.sum(CpS*Q_INF*wallAreas[:, 0])
CdP = 2.0*Dp/(Q_INF*2*R_CYL*depth)             # x2: mezzo dominio

CdF = 0.0
Ch = float("nan")
tanShear = None
try:
    tau = read_wall_boundary("wallShearStress")
    nh = wallCentres.copy(); nh[:, 2] = 0
    nh /= np.linalg.norm(nh[:, :2], axis=1, keepdims=True)
    tanShear = tau - np.sum(tau*nh, axis=1)[:, None]*nh
    Df = -np.sum(tanShear[:, 0]*np.linalg.norm(wallAreas, axis=1))
    CdF = 2.0*Df/(Q_INF*2*R_CYL*depth)
except FileNotFoundError:
    print("  (wallShearStress assente: C_D solo pressione)")
try:
    qw = read_wall_boundary("wallHeatFlux")
    Qint = np.sum(np.abs(qw)*np.linalg.norm(wallAreas, axis=1))/depth
    Ch = Qint/1000.0  # kW per metro di profondita', mezzo cilindro
except FileNotFoundError:
    print("  (wallHeatFlux assente: C_H non calcolato)")

print(f"C_D = {CdP + CdF:.3f} (pressione {CdP:.3f} + attrito {CdF:.3f}; "
      f"paper run 3: 1.304, DSMC 1.284, Newtoniano 1.333)")
if np.isfinite(Ch):
    print(f"C_H (mezzo cilindro, per m di profondita') = {Ch:.1f} kW "
          f"(paper run 3: 88.1, DSMC 63.3)")

# ------------------------------------------------------------------ plots ---
nN2 = rho*YN2/M_N2*NA
nN = np.maximum(rho*YN/M_N*NA, 1e10)

fig, axes = plt.subplots(1, 3, figsize=(15, 4.6))
axes[0].plot(xso - 0*R_CYL, Tso/1000, "r-", label="T")
axes[0].plot(xso, Tveo/1000, "b--", label="Tve")
axes[0].set_xlabel("stagnation line position x+R [m]")
axes[0].set_ylabel("T [kK]")
axes[0].set_xlim(-0.8, 0)
axes[0].legend(); axes[0].grid(alpha=0.3)
axes[0].set_title(f"standoff {standoff:.3f} m (paper ~0.25)")

Mline = np.linalg.norm(U[line][o], axis=1)/np.sqrt(GAMMA*R_N2*Tso)
axes[1].plot(xso, Mline, "k-")
axes[1].set_xlabel("stagnation line position x+R [m]")
axes[1].set_ylabel("Mach")
axes[1].set_xlim(-0.8, 0); axes[1].grid(alpha=0.3)

axes[2].semilogy(xso, nN2[line][o], "k-", label="N2")
axes[2].semilogy(xso, nN[line][o], "g--", label="N")
axes[2].set_xlabel("stagnation line position x+R [m]")
axes[2].set_ylabel("number density [1/m3]")
axes[2].set_xlim(-0.8, 0); axes[2].set_ylim(1e18, 1e22)
axes[2].legend(); axes[2].grid(alpha=0.3)

fig.suptitle(f"Mach 20 reacting cylinder - stagnation line (t = {tdir})")
fig.tight_layout()
fig.savefig(os.path.join(HERE, "cylinder-stagnation.png"), dpi=140)
print("Plot saved to cylinder-stagnation.png")

fig2, axes2 = plt.subplots(1, 3, figsize=(15, 4.6))
# CpS robusto (banda pulita); thetaWgeo = theta geometrico (0 = ristagno)
oS = np.argsort(thetaWgeo)
axes2[0].plot(thetaWgeo[oS], CpS[oS], "r-")
axes2[0].set_ylabel("pressure coefficient")
if tanShear is not None:
    CfS = np.linalg.norm(tanShear, axis=1)/Q_INF
    axes2[1].plot(thetaWgeo[oS], CfS[oS], "r-")
axes2[1].set_ylabel("friction coefficient")
try:
    axes2[2].plot(thetaWgeo[oS], np.abs(qw[oS])/1e4, "r-")
except NameError:
    pass
axes2[2].set_ylabel("surface heat flux [W/cm2]")
for ax in axes2:
    ax.set_xlabel("theta [deg]")
    ax.set_xlim(0, 180)
    ax.grid(alpha=0.3)
fig2.suptitle(f"Mach 20 reacting cylinder - surface (t = {tdir}); "
              "paper Fig. 5d-f")
fig2.tight_layout()
fig2.savefig(os.path.join(HERE, "cylinder-surface.png"), dpi=140)
print("Plot saved to cylinder-surface.png")

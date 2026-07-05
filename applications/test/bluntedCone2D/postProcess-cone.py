"""Post-processing for the Mach 11.3 blunted cone (paper Part Two, Fig. 2).

Reads the latest time dir plus cell centres (written by foamPostProcess
-func writeCellCentres, invoked here if missing) and produces:
  - cone-stagnation.png : T, Tve, rho/rho_inf along the stagnation line
  - cone-surface.png    : surface pressure coefficient (first-cell row)
  - printed metrics     : shock standoff, stagnation-point Cp vs the
                          modified-Newtonian estimate

Reference curves from the paper (Wang & Boyd CFD/MONACO, CUBRC run 31)
can be digitised and overlaid later; this script establishes the pipeline.
"""

import os
import re
import subprocess
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

R_N2 = 8.31446261815324/28.0134e-3
RHO_INF = 5.113e-4
U_INF = 2764.5
P_INF = 21.9139
T_INF = 144.4
RN = 6.35e-3
DELTA = np.radians(25.0)


def read_scalar(tdir, name):
    txt = open(os.path.join(tdir, name)).read()
    m = re.search(r"internalField\s+nonuniform\s+List<scalar>\s*\d+\s*"
                  r"\((.*?)\)\s*;", txt, re.S)
    if m:
        return np.array([float(q) for q in m.group(1).split()])
    mu = re.search(r"internalField\s+uniform\s+([^;]+);", txt)
    return float(mu.group(1))


def read_vector(tdir, name):
    txt = open(os.path.join(tdir, name)).read()
    m = re.search(r"internalField\s+nonuniform\s+List<vector>\s*\d+\s*"
                  r"\((.*?)\)\s*;", txt, re.S)
    rows = re.findall(r"\(([^)]+)\)", m.group(1))
    return np.array([[float(q) for q in r.split()] for r in rows])


times = sorted(
    (d for d in os.listdir(".")
     if re.fullmatch(r"[0-9.e+-]+", d) and float(d) > 0),
    key=float
)
tdir = times[-1]

if not os.path.exists(os.path.join(tdir, "C")):
    subprocess.run(
        ["foamPostProcess", "-func", "writeCellCentres",
         "-time", tdir],
        check=True, capture_output=True
    )

C = read_vector(tdir, "C")
T = read_scalar(tdir, "T")
Tve = read_scalar(tdir, "Tve")
p = read_scalar(tdir, "p")

x, y = C[:, 0], C[:, 1]
rho = p/(R_N2*T)

# --- stagnation line: first cell row along the axis -----------------------
# block A cells adjacent to the axis have the smallest |y| for x < -Rn
axisBand = np.unique(y[np.abs(x) < 2.9*RN])
yTol = np.sort(np.abs(axisBand))[0]*3
line = (np.abs(y) <= max(yTol, 1e-4)) & (x < -RN*0.99)
xs = x[line]
order = np.argsort(xs)
xs = xs[order]
Ts = T[line][order]
Tvs = Tve[line][order]
rhos = rho[line][order]

# shock standoff from the TEMPERATURE half-rise: the density gradient is
# useless here (the cold-wall boundary layer density spike dominates it)
Tpk = Ts.max()
i_sh = np.argmax(Ts > 0.5*(Tpk + T_INF))
standoff = -RN - xs[i_sh]
print(f"t = {tdir} s ({len(times)} snapshots)")
print(f"stagnation-line cells: {line.sum()}")
print(f"shock standoff (axis) = {standoff*1e3:.2f} mm "
      f"({standoff/RN:.3f} Rn; strong-shock correlations ~0.1-0.15 Rn)")

# stagnation Cp vs modified-Newtonian
i_w = np.argmax(xs)  # cell closest to the wall
cp_stag = (p[line][order][i_w] - P_INF)/(0.5*RHO_INF*U_INF**2)
M_inf = U_INF/np.sqrt(1.4*R_N2*T_INF)
cp_mn = 2.0/(1.4*M_inf**2)*(
    ((1.4 + 1)**2*M_inf**2/(4*1.4*M_inf**2 - 2*(1.4 - 1)))**(1.4/(1.4 - 1))
    *((1 - 1.4 + 2*1.4*M_inf**2)/(1.4 + 1)) - 1
)
print(f"stagnation Cp = {cp_stag:.3f}  (Rayleigh pitot ideal: {cp_mn:.3f})")

plt.figure(figsize=(7, 5))
plt.plot(xs*1e3, Ts/T_INF, "r-", label="T / T_inf")
plt.plot(xs*1e3, Tvs/T_INF, "b--", label="Tve / T_inf")
plt.plot(xs*1e3, rhos/RHO_INF, "k-", label="rho / rho_inf")
plt.axvline(-RN*1e3, color="gray", ls=":", lw=0.8)
plt.xlabel("x [mm] (wall at -6.35)")
plt.ylabel("normalised")
plt.title("Stagnation line, Mach 11.3 blunted cone")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("cone-stagnation.png", dpi=150)
print("Plot saved to cone-stagnation.png")

# --- surface pressure coefficient (first-cell row along the body) ---------
# nose: cells with r ~ Rn + half first cell; cone: distance from cone line
r = np.sqrt(x**2 + y**2)
onNose = (x < -RN*np.cos(np.pi/2 - DELTA)) & (r < RN*1.02) & (r > RN*0.9)
dCone = (y - (x + RN*np.cos(np.pi/2 - DELTA))*np.tan(DELTA)
         - RN*np.sin(np.pi/2 - DELTA))*np.cos(DELTA)
onCone = (x >= -RN*np.cos(np.pi/2 - DELTA)) & (dCone > 0) & (dCone < 1e-3)
surf = onNose | onCone
xw = x[surf]
cpw = (p[surf] - P_INF)/(0.5*RHO_INF*U_INF**2)
order = np.argsort(xw)

plt.figure(figsize=(7, 5))
plt.plot((xw[order] + RN)*1e3, cpw[order], "r.-", ms=3,
         label="shockThermo (first-cell row)")
plt.xlabel("axial distance from stagnation point [mm]")
plt.ylabel("pressure coefficient")
plt.title("Surface pressure, Mach 11.3 blunted cone")
plt.legend()
plt.grid(alpha=0.3)
plt.tight_layout()
plt.savefig("cone-surface.png", dpi=150)
print("Plot saved to cone-surface.png")

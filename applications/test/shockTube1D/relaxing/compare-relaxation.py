"""Compare the CFD post-shock relaxation zone against the steady ODE
reference (Test-postShockRelax).

Steps:
 1. locate the shock front in the last two snapshots (max |d rho/dx|)
    -> shock speed us;
 2. the caller (Allrun) runs Test-postShockRelax with the measured us,
    producing reference.csv (x, T, Tve, rho, u_sf, p, Y_N2, Y_N in the
    shock-attached frame, x = 0 at the front);
 3. this script (second invocation with --compare) extracts the CFD
    profile behind the front, maps it to the shock frame and overlays
    it with the reference; metrics over the relaxation zone.

Usage: python3 compare-relaxation.py --speed        (prints us only)
       python3 compare-relaxation.py --compare      (full comparison)
"""

import os
import re
import sys
import numpy as np

XMIN, XMAX, NCELL = -0.3, 0.7, 2000
R_N2 = 8.31446261815324/28.0134e-3


def read_field(tdir, name):
    txt = open(os.path.join(tdir, name)).read()
    m = re.search(r"internalField\s+nonuniform\s+List<(scalar|vector)>\s*"
                  r"\d+\s*\((.*?)\)\s*;", txt, re.S)
    if not m:
        mu = re.search(r"internalField\s+uniform\s+(\(?[^;]+?\)?)\s*;", txt)
        val = mu.group(1)
        if "(" in val:
            v = [float(q) for q in val.strip("()").split()]
            return np.tile(v, (NCELL, 1))
        return np.full(NCELL, float(val))
    body = m.group(2)
    if m.group(1) == "vector":
        rows = re.findall(r"\(([^)]+)\)", body)
        return np.array([[float(q) for q in r.split()] for r in rows])
    return np.array([float(q) for q in body.split()])


def cell_centres():
    dx = (XMAX - XMIN)/NCELL
    return np.linspace(XMIN, XMAX, NCELL, endpoint=False) + dx/2


def shock_position(tdir):
    T = read_field(tdir, "T")
    p = read_field(tdir, "p")
    rho = p/(R_N2*T)  # good enough to locate the front (pure N2 dominant)
    x = cell_centres()
    # search in the right half (ahead of the contact)
    mask = x > 0.02
    i = np.argmax(-np.diff(rho[mask]))
    return x[mask][i]


times = sorted(
    (d for d in os.listdir(".")
     if re.fullmatch(r"[0-9.e+-]+", d) and float(d) > 0),
    key=float
)

if "--speed" in sys.argv:
    t1, t2 = times[-2], times[-1]
    x1, x2 = shock_position(t1), shock_position(t2)
    us = (x2 - x1)/(float(t2) - float(t1))
    # integer m/s: keeps the shell and the C++ lround() file naming in sync
    print(round(us))
    sys.exit(0)

# ------------------------------------------------------------- compare ---
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

tdir = times[-1]
t = float(tdir)

T = read_field(tdir, "T")
Tve = read_field(tdir, "Tve")
p = read_field(tdir, "p")
U = read_field(tdir, "U")[:, 0]
YN2 = read_field(tdir, "N2")
YN = read_field(tdir, "N")
x = cell_centres()

t1, t2 = times[-2], times[-1]
us = (shock_position(t2) - shock_position(t1))/(float(t2) - float(t1))
x_sh = shock_position(tdir)

ref = np.loadtxt("reference.csv", delimiter=",", skiprows=1)
x_ref = ref[:, 0]
T_ref, Tve_ref = ref[:, 1], ref[:, 2]
YN2_ref, YN_ref = ref[:, 6], ref[:, 7]

# CFD profile in the shock frame: xi = distance behind the front
behind = (x < x_sh) & (x > x_sh - 0.06)
xi = x_sh - x[behind]
order = np.argsort(xi)
xi = xi[order]
T_c = T[behind][order]
Tve_c = Tve[behind][order]
YN2_c = YN2[behind][order]
YN_c = YN[behind][order]

# Metrics over the relaxation zone proper: from behind the smeared front
# (4 cells) to where the ODE reference reaches thermal equilibrium
# ((T - Tve)/T < 2%). Beyond that lies the contact surface, which belongs
# to the neighbouring wave, not to the relaxation physics.
dx = (XMAX - XMIN)/NCELL
i_eq = np.argmax((T_ref - Tve_ref)/T_ref < 0.02)
x_eq = x_ref[i_eq] if i_eq > 0 else x_ref[-1]
zone = (xi > 4*dx) & (xi < min(x_eq, xi.max()))
Ti = np.interp(xi[zone], x_ref, T_ref)
Tvei = np.interp(xi[zone], x_ref, Tve_ref)

errT = np.abs(T_c[zone] - Ti)/Ti
errTve = np.abs(Tve_c[zone] - Tvei)/np.maximum(Tvei, 1.0)

print(f"t = {t:.2e} s, measured us = {us:.1f} m/s, front at x = {x_sh:.4f} m")
print(f"relaxation zone ({(4*dx)*1e3:.1f} .. "
      f"{min(x_eq, xi.max())*1e3:.1f} mm behind front, "
      f"{zone.sum()} cells):")
print(f"  T   : max rel err = {errT.max()*100:.2f} %, "
      f"mean = {errT.mean()*100:.2f} %")
print(f"  Tve : max rel err = {errTve.max()*100:.2f} %, "
      f"mean = {errTve.mean()*100:.2f} %")

# --------------------------------------------------------------- plots ---
fig, (ax1, ax2) = plt.subplots(1, 2, figsize=(12, 5))

ax1.plot(x_ref*1e3, T_ref, "k-", lw=1.2, label="ODE reference T")
ax1.plot(x_ref*1e3, Tve_ref, "k--", lw=1.2, label="ODE reference Tve")
ax1.plot(xi*1e3, T_c, "ro", ms=3, label="shockThermo T")
ax1.plot(xi*1e3, Tve_c, "bs", ms=3, label="shockThermo Tve")
ax1.set_xlim(0, min(0.06, xi.max())*1e3)
ax1.set_xlabel("distance behind shock front [mm]")
ax1.set_ylabel("Temperature [K]")
ax1.grid(alpha=0.3)
ax1.legend()

ax2.plot(x, T, "r-", lw=1, label="T")
ax2.plot(x, Tve, "b-", lw=1, label="Tve")
ax2.axvline(x_sh, color="k", ls=":", lw=0.8)
ax2.set_xlabel("x [m]")
ax2.set_ylabel("Temperature [K]")
ax2.set_title(f"full domain at t = {t*1e6:.0f} us")
ax2.grid(alpha=0.3)
ax2.legend()

fig.suptitle(
    f"Relaxing shock tube, N2: post-shock zone vs steady ODE "
    f"(us = {us:.0f} m/s)"
)
fig.tight_layout()
fig.savefig("relaxation-comparison.png", dpi=150)
print("Plot saved to relaxation-comparison.png")

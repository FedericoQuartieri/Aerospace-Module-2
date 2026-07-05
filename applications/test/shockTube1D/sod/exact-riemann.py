"""Exact Riemann solution for the Sod tube in cold N2 (gamma = 1.4) and
quantitative comparison with the shockThermo solution.

Standard exact solver (Toro, Ch. 4) for an ideal gas with constant gamma.
The CFD profile is read from the last written time directory (fields T, p,
U along x). Comparison metrics: L1 errors of rho, u, p over the domain and
positions of the shock/contact.

Usage: python3 exact-riemann.py [caseDir=.]
Writes sod-comparison.png and prints the metrics.
"""

import os
import re
import sys
import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

GAMMA = 1.4
R_N2 = 8.31446261815324/28.0134e-3

# Initial states (must match system/setFieldsDict)
P_L, T_L = 1.0e5, 348.432
P_R, T_R = 1.0e4, 278.746
RHO_L = P_L/(R_N2*T_L)
RHO_R = P_R/(R_N2*T_R)
X0 = 0.0
XMIN, XMAX, NCELL = -5.0, 5.0, 1000

case = sys.argv[1] if len(sys.argv) > 1 else "."


# ---------------------------------------------------------------- exact ---
def exact_riemann(rhoL, uL, pL, rhoR, uR, pR, g):
    """Returns p*, u* of the star region (Toro's two-rarefaction start +
    Newton on the pressure function)."""
    aL = np.sqrt(g*pL/rhoL)
    aR = np.sqrt(g*pR/rhoR)

    def f_side(p, rho, ps, a):
        if p > ps:  # shock
            A = 2.0/((g + 1)*rho)
            B = (g - 1)/(g + 1)*ps
            return (p - ps)*np.sqrt(A/(p + B))
        # rarefaction
        return 2*a/(g - 1)*((p/ps)**((g - 1)/(2*g)) - 1)

    def df_side(p, rho, ps, a):
        if p > ps:
            A = 2.0/((g + 1)*rho)
            B = (g - 1)/(g + 1)*ps
            return np.sqrt(A/(B + p))*(1 - (p - ps)/(2*(B + p)))
        return 1.0/(rho*a)*(p/ps)**(-(g + 1)/(2*g))

    # initial guess (two-rarefaction)
    p = ((aL + aR - 0.5*(g - 1)*(uR - uL)) /
         (aL/pL**((g - 1)/(2*g)) + aR/pR**((g - 1)/(2*g))))**(2*g/(g - 1))

    for _ in range(50):
        f = f_side(p, rhoL, pL, aL) + f_side(p, rhoR, pR, aR) + (uR - uL)
        df = df_side(p, rhoL, pL, aL) + df_side(p, rhoR, pR, aR)
        dp = -f/df
        p = max(p + dp, 1e-8)
        if abs(dp) < 1e-12*p:
            break

    u = 0.5*(uL + uR) + 0.5*(f_side(p, rhoR, pR, aR)
                             - f_side(p, rhoL, pL, aL))
    return p, u


def sample(x, t, rhoL, uL, pL, rhoR, uR, pR, g):
    """Samples the exact solution at positions x for time t."""
    aL = np.sqrt(g*pL/rhoL)
    aR = np.sqrt(g*pR/rhoR)
    ps, us = exact_riemann(rhoL, uL, pL, rhoR, uR, pR, g)

    rho = np.zeros_like(x)
    u = np.zeros_like(x)
    p = np.zeros_like(x)

    s = (x - X0)/t

    # left rarefaction fan bounds
    shL = uL - aL
    aLs = aL*(ps/pL)**((g - 1)/(2*g))
    stL = us - aLs

    # right shock speed (Sod: right wave is a shock)
    rhoRs = rhoR*((ps/pR + (g - 1)/(g + 1)) /
                  ((g - 1)/(g + 1)*ps/pR + 1))
    sR = uR + aR*np.sqrt((g + 1)/(2*g)*ps/pR + (g - 1)/(2*g))

    rhoLs = rhoL*(ps/pL)**(1/g)

    for i, si in enumerate(s):
        if si < shL:                      # undisturbed left
            rho[i], u[i], p[i] = rhoL, uL, pL
        elif si < stL:                    # rarefaction fan
            u[i] = 2/(g + 1)*(aL + (g - 1)/2*uL + si)
            a = aL - (g - 1)/2*(u[i] - uL)
            rho[i] = rhoL*(a/aL)**(2/(g - 1))
            p[i] = pL*(a/aL)**(2*g/(g - 1))
        elif si < us:                     # left star region
            rho[i], u[i], p[i] = rhoLs, us, ps
        elif si < sR:                     # right star region
            rho[i], u[i], p[i] = rhoRs, us, ps
        else:                             # undisturbed right
            rho[i], u[i], p[i] = rhoR, uR, pR

    return rho, u, p, us, sR


# ------------------------------------------------------------------ CFD ---
def read_field(tdir, name, ncell):
    txt = open(os.path.join(tdir, name)).read()
    m = re.search(r"internalField\s+nonuniform\s+List<(scalar|vector)>\s*"
                  r"\d+\s*\((.*?)\)\s*;", txt, re.S)
    if not m:
        mu = re.search(r"internalField\s+uniform\s+"
                       r"(\(?[^;]+?\)?)\s*;", txt)
        val = mu.group(1)
        if "(" in val:
            v = [float(q) for q in val.strip("()").split()]
            return np.tile(v, (ncell, 1))
        return np.full(ncell, float(val))
    body = m.group(2)
    if m.group(1) == "vector":
        rows = re.findall(r"\(([^)]+)\)", body)
        return np.array([[float(q) for q in r.split()] for r in rows])
    return np.array([float(q) for q in body.split()])


times = sorted(
    (d for d in os.listdir(case)
     if re.fullmatch(r"[0-9.e+-]+", d) and float(d) > 0),
    key=float
)
tdir = times[-1]
t = float(tdir)

T = read_field(os.path.join(case, tdir), "T", NCELL)
p = read_field(os.path.join(case, tdir), "p", NCELL)
U = read_field(os.path.join(case, tdir), "U", NCELL)
u = U[:, 0]
rho = p/(R_N2*T)

x = np.linspace(XMIN, XMAX, NCELL, endpoint=False) + (XMAX - XMIN)/NCELL/2

rho_e, u_e, p_e, us, sR = sample(x, t, RHO_L, 0, P_L, RHO_R, 0, P_R, GAMMA)

# ------------------------------------------------------------- metrics ---
def l1(a, b, norm):
    return np.mean(np.abs(a - b))/norm

print(f"t = {t} s, exact shock position x = {sR*t:.3f} m, "
      f"contact x = {us*t:.3f} m")
print(f"L1(rho)/rho_L = {l1(rho, rho_e, RHO_L):.4f}")
print(f"L1(u)/u*      = {l1(u, u_e, abs(us)):.4f}")
print(f"L1(p)/p_L     = {l1(p, p_e, P_L):.4f}")

# CFD shock position: steepest rho gradient right of the contact
mask = x > us*t + 0.2
i_sh = np.argmax(-np.diff(rho[mask]))
x_sh_cfd = x[mask][i_sh]
print(f"shock position: CFD = {x_sh_cfd:.3f} m, exact = {sR*t:.3f} m "
      f"(err = {abs(x_sh_cfd - sR*t):.3f} m)")

# ---------------------------------------------------------------- plots ---
fig, axes = plt.subplots(2, 2, figsize=(11, 7), sharex=True)
for ax, (name, num, ex) in zip(
    axes.flat,
    [("rho [kg/m^3]", rho, rho_e), ("u [m/s]", u, u_e),
     ("p [Pa]", p, p_e), ("T [K]", T, p_e/(R_N2*rho_e))]
):
    ax.plot(x, ex, "k-", lw=1, label="exact Riemann")
    ax.plot(x[::5], num[::5], "ro", ms=2.5, label="shockThermo")
    ax.set_ylabel(name)
    ax.grid(alpha=0.3)
axes[0, 0].legend()
for ax in axes[1]:
    ax.set_xlabel("x [m]")
fig.suptitle(f"Sod shock tube, pure N2 (frozen), t = {t*1e3:.1f} ms")
fig.tight_layout()
fig.savefig("sod-comparison.png", dpi=150)
print("Plot saved to sod-comparison.png")

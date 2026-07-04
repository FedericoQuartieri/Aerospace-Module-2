"""Compare the reacting solver heat bath against the standalone 0D test.

Reference: output/results-N2N-reacting-30000-1000-mpp-park05.csv from
Test-N2N-reacting (Mutation++ end-to-end path: same kinetics, VT and CV
sources the solver bridge uses). Columns: t, T_tr, T_ve, nN2/n0, nN/n0.

Solver: probe histories of T, Tve, rho, N2 (=Y_N2), N (=Y_N).
Number densities: n_i = rho*Y_i/Mw_i*NA, normalised by n0 = 1.0e23 m^-3.
"""

import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

NA = 6.02214076e23
MW_N2 = 28.0134e-3   # kg/mol
MW_N = 14.0067e-3
N0 = 1.0e23          # total initial number density [m^-3]

ref = np.loadtxt("reference.csv", delimiter=",", skiprows=1)
t_ref = ref[:, 0]
Ttr_ref, Tve_ref = ref[:, 1], ref[:, 2]
nN2_ref, nN_ref = ref[:, 3], ref[:, 4]


def load_probe(field):
    path = glob.glob(f"postProcessing/probes/0/{field}")[0]
    data = np.loadtxt(path, comments="#")
    return data[:, 0], data[:, 1]


t_T, T = load_probe("T")
t_Tve, Tve = load_probe("Tve")
t_rho, rho = load_probe("rho")
t_YN2, YN2 = load_probe("N2")
t_YN, YN = load_probe("N")

nN2 = rho*YN2/MW_N2*NA/N0
nN = rho*YN/MW_N*NA/N0

# Interpolate the solver histories onto the reference times for the metrics
mask = (t_ref >= max(t_T[0], 1e-12)) & (t_ref <= t_T[-1])
tm = t_ref[mask]

T_i = np.interp(tm, t_T, T)
Tve_i = np.interp(tm, t_Tve, Tve)
nN2_i = np.interp(tm, t_rho, nN2)
nN_i = np.interp(tm, t_rho, nN)

err_T = np.abs(T_i - Ttr_ref[mask])
err_Tve = np.abs(Tve_i - Tve_ref[mask])
err_nN2 = np.abs(nN2_i - nN2_ref[mask])
err_nN = np.abs(nN_i - nN_ref[mask])

print(f"T_tr   : max |err| = {err_T.max():8.2f} K   "
      f"({100*(err_T/Ttr_ref[mask]).max():.2f} %)")
print(f"T_ve   : max |err| = {err_Tve.max():8.2f} K   "
      f"({100*(err_Tve/np.maximum(Tve_ref[mask], 1.0)).max():.2f} %)")
print(f"nN2/n0 : max |err| = {err_nN2.max():.4f}")
print(f"nN /n0 : max |err| = {err_nN.max():.4f}")
print(f"final @ {t_T[-1]:.2e} s: solver T={T[-1]:.1f} Tve={Tve[-1]:.1f} "
      f"nN2/n0={nN2[-1]:.4f} nN/n0={nN[-1]:.4f}")
i_end = np.searchsorted(t_ref, t_T[-1])
i_end = min(i_end, len(t_ref) - 1)
print(f"reference        : T={Ttr_ref[i_end]:.1f} Tve={Tve_ref[i_end]:.1f} "
      f"nN2/n0={nN2_ref[i_end]:.4f} nN/n0={nN_ref[i_end]:.4f}")

# --- temperature plot ---
plt.figure()
m = t_ref > 0
plt.plot(t_ref[m], Ttr_ref[m], "k-", label="reference T_tr (0D mpp)")
plt.plot(t_ref[m], Tve_ref[m], "k--", label="reference T_ve (0D mpp)")
m = t_T > 0
plt.plot(t_T[m][::30], T[m][::30], "ro", ms=3, label="shockThermo T")
plt.plot(t_Tve[m][::30], Tve[m][::30], "bs", ms=3, label="shockThermo Tve")
plt.xscale("log")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.title("Reacting N2-N heat bath: solver vs standalone (Park, q=0.5)")
plt.legend()
plt.tight_layout()
plt.savefig("heatbath-reacting-comparison.png", dpi=150)
print("Plot saved to heatbath-reacting-comparison.png")

# --- composition plot ---
plt.figure()
m = t_ref > 0
plt.plot(t_ref[m], nN2_ref[m], "k-", label="reference N2 (0D mpp)")
plt.plot(t_ref[m], nN_ref[m], "k--", label="reference N (0D mpp)")
m = t_rho > 0
plt.plot(t_rho[m][::30], nN2[m][::30], "ro", ms=3, label="shockThermo N2")
plt.plot(t_rho[m][::30], nN[m][::30], "bs", ms=3, label="shockThermo N")
plt.xscale("log")
plt.xlabel("Time (s)")
plt.ylabel("Normalised number density, n / n0")
plt.title("Reacting N2-N heat bath: composition, solver vs standalone")
plt.legend()
plt.tight_layout()
plt.savefig("heatbath-reacting-composition.png", dpi=150)
print("Plot saved to heatbath-reacting-composition.png")

import glob
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Reference: standalone Mutation++ heat bath (t, T_tr, T_ve)
ref = np.loadtxt("reference.csv", delimiter=",", skiprows=1)
t_ref, Ttr_ref, Tve_ref = ref[:, 0], ref[:, 1], ref[:, 2]

# Solver probes: "time value" per line, one probe
def load_probe(field):
    path = glob.glob(f"postProcessing/probes/0/{field}")[0]
    data = np.loadtxt(path, comments="#")
    return data[:, 0], data[:, 1]

t_T, T = load_probe("T")
t_Tve, Tve = load_probe("Tve")

# Interpolate the solver histories onto the reference times for the metrics
mask = (t_ref >= max(t_T[0], 1e-12)) & (t_ref <= t_T[-1])
T_i = np.interp(t_ref[mask], t_T, T)
Tve_i = np.interp(t_ref[mask], t_Tve, Tve)

err_T = np.abs(T_i - Ttr_ref[mask])
err_Tve = np.abs(Tve_i - Tve_ref[mask])

print(f"T_tr : max |err| = {err_T.max():8.2f} K   "
      f"({100 * (err_T / Ttr_ref[mask]).max():.2f} %)")
print(f"T_ve : max |err| = {err_Tve.max():8.2f} K   "
      f"({100 * (err_Tve / np.maximum(Tve_ref[mask], 1.0)).max():.2f} %)")
print(f"final: solver T = {T[-1]:.1f} / Tve = {Tve[-1]:.1f}   "
      f"reference T = {Ttr_ref[-1]:.1f} / Tve = {Tve_ref[-1]:.1f}")

plt.figure()

m = t_ref > 0
plt.plot(t_ref[m], Ttr_ref[m], "k-", label="reference T_tr (Mutation++ 0D)")
plt.plot(t_ref[m], Tve_ref[m], "k--", label="reference T_ve (Mutation++ 0D)")
m = t_T > 0
plt.plot(t_T[m][::20], T[m][::20], "ro", ms=3, label="shockThermo T")
plt.plot(t_Tve[m][::20], Tve[m][::20], "bs", ms=3, label="shockThermo Tve")

plt.xscale("log")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.title("N2 heat bath: solver vs standalone Mutation++ (VT relaxation)")
plt.legend()
plt.tight_layout()
plt.savefig("heatbath-comparison.png", dpi=150)
print("Plot saved to heatbath-comparison.png")

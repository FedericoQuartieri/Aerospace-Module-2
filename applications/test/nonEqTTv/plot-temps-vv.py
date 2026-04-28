import sys
import numpy as np
import matplotlib.pyplot as plt

# Load data
data = np.loadtxt(sys.argv[1], delimiter=",", skiprows=1)

# Columns: t, T_tr, T_ve_N2, T_ve_O2
t = data[:, 0]
T_tr = data[:, 1]
T_vN2 = data[:, 2]
T_vO2 = data[:, 3]

# Remove t = 0 for log scale
mask = t > 0
t = t[mask]
T_tr = T_tr[mask]
T_vN2 = T_vN2[mask]
T_vO2 = T_vO2[mask]

# Reference equilibrium line
plt.axhline(T_tr[-1], color="black", linestyle="--", alpha=0.5)
plt.text(t[-1] * 1.05, T_tr[-1] * 1.02, f"{T_tr[-1]:.1f}", color="black", va="bottom", ha="right")

plt.plot(t, T_tr, label="T_tr", linewidth=2)
plt.plot(t, T_vN2, label="T_ve (N2)", linestyle="-.")
plt.plot(t, T_vO2, label="T_ve (O2)", linestyle=":")

plt.xscale("log")
plt.xlabel("Time (s)")
plt.ylabel("Temperature (K)")
plt.legend()

plt.tight_layout()
plt.savefig("temp-curves-vv.png")

print(f"Plot saved to temp-curves-vv.png")

import sys
import numpy as np
import matplotlib.pyplot as plt

# Load data (skip header if present)
data = np.loadtxt("results.csv", delimiter=",", skiprows=1)

# Columns
t = data[:, 0]
T_tr = data[:, 1]
T_ve = data[:, 2]

# Remove t = 0 for log scale
mask = t > 0
t = t[mask]
T_tr = T_tr[mask]
T_ve = T_ve[mask]

# Plot
plt.figure()

plt.axhline(T_tr[-1], color="red", linestyle="--")

plt.plot(t, T_tr, label="T_tr")
plt.plot(t, T_ve, label="T_ve")

plt.xscale("log")  # log scale on x-axis
plt.xlabel("Time (s)")
plt.ylabel("Temperature")
plt.legend()

plt.tight_layout()
plt.savefig("temp-curves.png")

print("Plot saved to temp-curves.png")

import matplotlib.pyplot as plt
import pandas as pd
import numpy as np

# Constants
NA = 6.02214076e23
Mw_N = 14.0067e-3  # kg/mol
Mw_N2 = 28.0134e-3 # kg/mol
n0 = 5.0e22        # Initial number density for normalization

# Read Data
try:
    data = pd.read_csv('results_refactor.csv')
except ImportError:
    print("Pandas not found, using numpy")
    data_np = np.genfromtxt('results_refactor.csv', delimiter=',', names=True)
    data = pd.DataFrame(data_np)

# Rename 'time[s]' to 'Time' for consistency with original script's plotting logic
data.rename(columns={'time[s]': 'Time'}, inplace=True)

# Calculate Number Densities (m^-3)
# n = (rho / Mw) * NA
data['n_N'] = (data['rho_N[kg/m3]'] / Mw_N) * NA
data['n_N2'] = (data['rho_N2[kg/m3]'] / Mw_N2) * NA

# Normalize
data['n_N_norm'] = data['n_N'] / n0
data['n_N2_norm'] = data['n_N2'] / n0

# --- Plot 1: Temperatures ---
plt.figure(figsize=(10, 6))
plt.semilogx(data['Time'], data['T_tr[K]'], label='T_tr (hy2Foam)', color='black', linewidth=2)
plt.semilogx(data['Time'], data['T_vib[K]'], label='T_vib (hy2Foam)', color='black', linestyle='--', linewidth=2)

plt.title('Reacting N2-N Heat Bath: Temperatures (Park Model)')
plt.xlabel('Time (s)')
plt.ylabel('Temperature (K)')
plt.grid(True, which="both", ls="-", alpha=0.4)
plt.legend()
plt.xlim(0, 1e-5)
plt.ylim(4000, 10000)
plt.savefig('temperature_plot_refactor.png')
print("Saved temperature_plot_refactor.png")

# --- Plot 2: Normalized Number Densities ---
plt.figure(figsize=(10, 6))
plt.semilogx(data['Time'], data['n_N_norm'], label='N (hy2Foam)', color='black', marker='o', markevery=0.1, linestyle='None')
plt.semilogx(data['Time'], data['n_N2_norm'], label='N2 (hy2Foam)', color='black', marker='+', markevery=0.1, linestyle='None')

# Plot lines for clarity as well
plt.semilogx(data['Time'], data['n_N_norm'], color='black', alpha=0.5)
plt.semilogx(data['Time'], data['n_N2_norm'], color='black', alpha=0.5)

plt.title('Reacting N2-N Heat Bath: Densities (Park Model)')
plt.xlabel('Time (s)')
plt.ylabel('Normalized Number Density (n/n0)')
plt.grid(True, which="both", ls="-", alpha=0.4)
plt.legend()
plt.xlim(1e-6, 2e-4)
plt.ylim(0, 1.4)
plt.savefig('density_plot_refactor.png')
print("Saved density_plot_refactor.png")

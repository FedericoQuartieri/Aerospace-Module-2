import matplotlib
matplotlib.use('Agg')  # Use non-interactive backend
import matplotlib.pyplot as plt
import numpy as np

# Constants
NA = 6.02214076e23
Mw_N = 14.0067e-3  # kg/mol
Mw_N2 = 28.0134e-3 # kg/mol
n0 = 5.0e22        # Initial number density for normalization

# Read Data using numpy
print("Loading hy2Foam data...")
data_np = np.genfromtxt('results.csv', delimiter=',', names=True)

# Create dict-like access  
data = {}
data['Time'] = data_np['time']
data['T_tr[K]'] = data_np['T_tr']
data['T_vib[K]'] = data_np['T_vib']
data['rho_N[kg/m3]'] = data_np['rho_N']
data['rho_N2[kg/m3]'] = data_np['rho_N2']

# Calculate Number Densities (m^-3)
# n = (rho / Mw) * NA
data['n_N'] = (data['rho_N[kg/m3]'] / Mw_N) * NA
data['n_N2'] = (data['rho_N2[kg/m3]'] / Mw_N2) * NA

# Normalize
data['n_N_norm'] = data['n_N'] / n0
data['n_N2_norm'] = data['n_N2'] / n0

# Read reference data from Engauge (handle European decimal format with quotes)
print("Loading reference data...")
import csv
ref_data_raw = []
with open('Engauge_results_module2.csv', 'r') as f:
    reader = csv.reader(f)
    next(reader)  # Skip header
    for row in reader:
        # First column has European format with comma as decimal separator
        time_str = row[0].replace(',', '.')
        time = float(time_str)
        t_tr = float(row[1])
        t_vib = float(row[2])
        ref_data_raw.append([time, t_tr, t_vib])

ref_data = np.array(ref_data_raw)
ref_time = ref_data[:, 0]
ref_T_tr = ref_data[:, 1]
ref_T_vib = ref_data[:, 2]
print(f"  Loaded {len(ref_time)} reference data points")
print(f"  Time range: [{ref_time.min():.3e}, {ref_time.max():.3e}] s")
print(f"  T_tr range: [{ref_T_tr.min():.1f}, {ref_T_tr.max():.1f}] K")

# --- Plot 1: Temperatures ---
plt.figure(figsize=(12, 7))

# Add reference data FIRST (so it appears behind)
plt.semilogx(ref_time, ref_T_tr, label='T_tr (Reference)', color='red', linewidth=3, linestyle='-', marker='o', markersize=4, markevery=2)
plt.semilogx(ref_time, ref_T_vib, label='T_vib (Reference)', color='blue', linewidth=3, linestyle='--', marker='s', markersize=4, markevery=2)

# Then hy2Foam data
plt.semilogx(data['Time'], data['T_tr[K]'], label='T_tr (hy2Foam)', color='black', linewidth=2.5, alpha=0.8)
plt.semilogx(data['Time'], data['T_vib[K]'], label='T_vib (hy2Foam)', color='darkgreen', linestyle='--', linewidth=2.5, alpha=0.8)

plt.title('Reacting N2-N Heat Bath: Temperatures (Park Model)', fontsize=14, fontweight='bold')
plt.xlabel('Time (s)', fontsize=12)
plt.ylabel('Temperature (K)', fontsize=12)
plt.grid(True, which="both", ls="-", alpha=0.3)
plt.legend(loc='best', fontsize=11)
plt.tight_layout()
plt.savefig('temperature_plot_refactor.png', dpi=200, bbox_inches='tight')
print("Saved temperature_plot_refactor.png")
print(f"  hy2Foam data range: T_tr=[{data['T_tr[K]'].min():.1f}, {data['T_tr[K]'].max():.1f}] K, time=[{data['Time'].min():.2e}, {data['Time'].max():.2e}] s")
print(f"  Reference data range: T_tr=[{ref_T_tr.min():.1f}, {ref_T_tr.max():.1f}] K, time=[{ref_time.min():.2e}, {ref_time.max():.2e}] s")

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

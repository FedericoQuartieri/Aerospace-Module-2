import matplotlib.pyplot as plt
import numpy as np

# Read CSV
data = np.genfromtxt("results.csv", delimiter=",", names=True)

# Log-scale needs strictly positive time
mask = data["time"] > 0
time = data["time"][mask]

def save_semilogx_single(y: np.ndarray, label: str, ylabel: str, marker: str, out_png: str):
    plt.figure(figsize=(7.2, 4.8))
    plt.semilogx(time, y, marker=marker, linestyle="None", label=label)
    plt.xlabel("Time, t [s]")
    plt.ylabel(ylabel)
    plt.legend(frameon=True, loc="best")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

def save_temperature_plot(out_png: str):
    plt.figure(figsize=(7.2, 4.8))
    T = data["T"][mask] / 1000.0
    plt.semilogx(time, T, marker="+", linestyle="-", linewidth=1.0, markersize=4, label="T_tr")

    if "Tve" in data.dtype.names:
        Tve = data["Tve"][mask] / 1000.0
        plt.semilogx(
            time,
            Tve,
            marker="o",
            markerfacecolor="none",
            linestyle="-",
            linewidth=1.0,
            markersize=4,
            label="T_v"
        )

    plt.xlabel("Time, t [s]")
    plt.ylabel(r"Temperature [$10^3$ K]")
    plt.xlim(1e-9, 1e-4)
    plt.legend(frameon=True, loc="best")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

# Same “paper-ish” feel: log-x + markers only
save_semilogx_single(data["p"][mask], "p", "Pressure, p [Pa]", marker="v", out_png="p_vs_time.png")
save_temperature_plot(out_png="T_Tve_vs_time.png")
save_semilogx_single(data["rho"][mask], "rho", r"Density, $\rho$ [kg/m$^3$]", marker="o", out_png="rho_vs_time.png")

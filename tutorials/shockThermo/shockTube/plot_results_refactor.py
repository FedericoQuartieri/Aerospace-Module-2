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
    plt.semilogx(time, T, markersize=4, label="T_tr")

    if "Tve" in data.dtype.names:
        Tve = data["Tve"][mask] / 1000.0
        plt.semilogx(
            time,
            Tve,
            markersize=4,
            label="T_ve"
        )

    plt.xlabel("Time, t [s]")
    plt.ylabel(r"Temperature [$10^3$ K]")
    plt.xlim(1e-9, 2.5e-5)
    plt.legend(frameon=True, loc="best")
    plt.tight_layout()
    plt.savefig(out_png, dpi=300)
    plt.close()

save_temperature_plot(out_png="Ttr_Tve.png")

from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from collections.abc import Sequence

from analyte import Analyte
from electrolyte import Electrolyte
from solver import solve


PZC_SIO2 = 2.0  # (pKa + pKb) / 2


def ph_sweep(analyte: Analyte,electrolyte: Electrolyte,ph_min: float = 1.0,ph_max: float = 12.0,n: int = 200 ) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    #Solve for psi_0 across a range of bulk pH values.
    ph_values = np.linspace(ph_min, ph_max, n)
    psi_values = np.array([solve(analyte, electrolyte, ph) for ph in ph_values])
    return ph_values, psi_values

def slope(ph_values, psi_values):
    return np.gradient(psi_values, ph_values)

def plot_ph_sweep(ph_values: NDArray[np.float64],psi_values: NDArray[np.float64], out_path: Path, pzc: float | None = PZC_SIO2, ) -> None:
    #Plot psi_0 against pH and save it.
    fig, ax = plt.subplots(figsize=(6, 4))
    ax.plot(ph_values, psi_values * 1e3, linewidth=2)

    if pzc is not None:
        ax.axhline(0.0, linestyle="--", linewidth=1, color="grey")
        ax.axvline(pzc, linestyle="--", linewidth=1, color="grey")
        ax.annotate(
            f"pzc = {pzc:g}",
            xy=(pzc, 0.0),
            xytext=(6, 6),
            textcoords="offset points",
            fontsize=9,
            color="grey",
        )

    ax.set_xlabel("bulk pH")
    ax.set_ylabel("surface potential $\\psi_0$ (mV)")
    ax.set_title("SiO$_2$ ISFET response")
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)



#Test
def plot_salt_sweep(
    analyte: Analyte,
    electrolytes: Sequence[Electrolyte],
    labels: Sequence[str],
    out_path: Path,
    pzc: float | None = PZC_SIO2,
) -> None:
    """Plot psi_0 vs pH at several ionic strengths on one axes.
 
    The point of the figure is that the curves flatten as salt rises:
    more salt means a shorter Debye length, a bigger C_DL, and a larger
    share of the divider taken by the double layer, so less of the
    surface charge shows up as potential.
    """
    fig, ax = plt.subplots(figsize=(6, 4))
 
    for electrolyte, label in zip(electrolytes, labels):
        ph_values, psi_values = ph_sweep(analyte, electrolyte)
        peak = np.max(np.abs(slope(ph_values, psi_values))) * 1e3
        ax.plot(
            ph_values,
            psi_values * 1e3,
            linewidth=2,
            label=f"{label}  ({peak:.0f} mV/pH)",
        )
 
    if pzc is not None:
        ax.axhline(0.0, linestyle="--", linewidth=1, color="grey")
        ax.axvline(pzc, linestyle="--", linewidth=1, color="grey")
 
    ax.set_xlabel("bulk pH")
    ax.set_ylabel("surface potential $\\psi_0$ (mV)")
    ax.set_title("Screening: sensitivity falls as salt rises")
    ax.legend(title="ionic strength", fontsize=9)
    ax.grid(alpha=0.3)
 
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)
 


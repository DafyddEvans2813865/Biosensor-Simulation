from pathlib import Path

import matplotlib.pyplot as plt
import numpy as np
from numpy.typing import NDArray
from collections.abc import Sequence

from analyte import Analyte
from electrolyte import Electrolyte
from solver import solve
from surface import AmphotericSurface, ChargedSurface, beta_int
from charge import surface_pH

PZC_SIO2 = 2.0  # (pKa + pKb) / 2


def ph_sweep(surface: ChargedSurface, electrolyte: Electrolyte, ph_min: float = 1.0, ph_max: float = 12.0, n: int = 200) -> tuple[NDArray[np.float64], NDArray[np.float64]]:
    ph_values = np.linspace(ph_min, ph_max, n)
    psi_values = np.array([solve(surface, electrolyte, ph) for ph in ph_values])
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


def plot_analyte_sweep_comparison(
    analytes: Sequence[Analyte],
    labels: Sequence[str],
    electrolyte: Electrolyte,
    out_path: Path,
) -> None:
    if len(analytes) != len(labels):
        raise ValueError("analytes and labels must have the same length")

    fig, ax = plt.subplots(figsize=(6, 4))
    styles = ("-", "--", "-", "--")
    for index, (analyte, label) in enumerate(zip(analytes, labels)):
        ph_values, psi_values = ph_sweep(analyte, electrolyte)
        ax.plot(
            ph_values,
            psi_values * 1e3,
            linewidth=2.5 if index % 2 == 0 else 1.8,
            linestyle=styles[index % len(styles)],
            label=label,
            zorder=index + 2,
        )

    ax.set_xlabel("bulk pH")
    ax.set_ylabel("surface potential $\\psi_0$ (mV)")
    ax.set_title("Baseline and TPSA-blocked peptide pH response")
    ax.legend()
    ax.grid(alpha=0.3)
    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)


def plot_tpsa_comparison(
    baseline: Analyte,
    ytsf_tpsa: Analyte,
    yaaf_tpsa: Analyte,
    ytsf_tpsa_value: float,
    yaaf_tpsa_value: float,
    electrolyte: Electrolyte,
    out_path: Path,
) -> None:
    ph_values, baseline_psi = ph_sweep(baseline, electrolyte)
    _, ytsf_psi = ph_sweep(ytsf_tpsa, electrolyte)
    _, yaaf_psi = ph_sweep(yaaf_tpsa, electrolyte)

    fig, (ax_response, ax_difference) = plt.subplots(1, 2, figsize=(11, 4), sharex=True)

    ax_response.plot(
        ph_values,
        baseline_psi * 1e3,
        color="0.35",
        linestyle="--",
        linewidth=2,
        label="Baseline (YTSF = YAAF)",
    )
    ax_response.plot(
        ph_values,
        ytsf_psi * 1e3,
        color="tab:blue",
        linewidth=2,
        label=f"YTSF, TPSA = {ytsf_tpsa_value:.0f} Å²",
    )
    ax_response.plot(
        ph_values,
        yaaf_psi * 1e3,
        color="tab:orange",
        linewidth=2,
        label=f"YAAF, TPSA = {yaaf_tpsa_value:.0f} Å²",
    )

    ax_response.set_xlabel("bulk pH")
    ax_response.set_ylabel("surface potential $\\psi_0$ (mV)")
    ax_response.set_title("Surface potential")
    ax_response.legend(fontsize=9)
    ax_response.grid(alpha=0.3)

    ax_difference.plot(
        ph_values,
        (ytsf_psi - yaaf_psi) * 1e3,
        color="black",
        linewidth=2,
        label="Δψ₀ (YTSF − YAAF)",
    )
    ax_difference.axhline(0.0, color="0.35", linestyle="--", linewidth=1)
    ax_difference.set_xlabel("bulk pH")
    ax_difference.set_ylabel("surface potential difference (mV)")
    ax_difference.set_title("TPSA response difference")
    ax_difference.legend(fontsize=9)
    ax_difference.grid(alpha=0.3)

    fig.tight_layout()
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=150)
    plt.close(fig)



#Test
def plot_salt_sweep(analyte: Analyte,electrolytes: Sequence[Electrolyte],labels: Sequence[str],out_path: Path,pzc: float | None = PZC_SIO2,) -> None:
    #Plot psi_0 vs pH at several ionic strengths on one axes.
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


def plot_van_hal(surfaces: Sequence[AmphotericSurface], labels: Sequence[str], electrolyte: Electrolyte, out_path: Path) -> None:
    if len(surfaces) != len(labels):
        raise ValueError("surfaces and labels must have the same length")

    vt = electrolyte.thermal_voltage
    fig, (ax_beta, ax_cdif, ax_alpha) = plt.subplots(1, 3, figsize=(14, 4.5), sharex=True, constrained_layout=True)
    colors = plt.get_cmap("tab10").colors

    for index, (surface, label) in enumerate(zip(surfaces, labels)):
        ph_values, psi_values = ph_sweep(surface, electrolyte, surface.pzc - 4.0, surface.pzc + 4.0)
        delta_ph = ph_values - surface.pzc
        ph_surface = surface_pH(ph_values, psi_values, vt)

        betas = np.array([beta_int(surface, ph) for ph in ph_surface])

        # C_dif is evaluated at the diffuse-layer potential, after the Stern drop
        sigma_0 = np.array([surface.charge(ph) for ph in ph_surface])
        psi_d = psi_values - sigma_0 / electrolyte.c_stern if electrolyte.c_stern is not None else psi_values
        c_difs = np.array([electrolyte.c_dif(p) for p in psi_d])

        alphas = -slope(ph_values, psi_values) / (np.log(10.0) * vt)
        color = colors[index % len(colors)]

        ax_beta.plot(delta_ph, betas, color=color, linewidth=2, label=label)
        ax_cdif.plot(delta_ph, c_difs, color=color, linewidth=2, label=label)
        ax_alpha.plot(delta_ph, alphas, color=color, linewidth=2, label=label)

    ax_beta.set_yscale("log")
    ax_beta.set_ylim(1e15, 1e19)
    ax_beta.set_ylabel(r"$\beta_{int}$ (groups/m$^2$)")
    ax_beta.set_title("Intrinsic buffer capacity (Fig. 1)")

    ax_cdif.set_ylim(0.10, 0.20)
    ax_cdif.set_ylabel(r"$C_{dif}$ (F/m$^2$)")
    ax_cdif.set_title("Differential capacitance (Fig. 2)")

    ax_alpha.set_ylim(0.0, 1.05)
    ax_alpha.set_ylabel(r"$\alpha$")
    ax_alpha.set_title("Sensitivity parameter (Fig. 3)")

    for axis in (ax_beta, ax_cdif, ax_alpha):
        axis.set_xlabel(r"$\Delta$pH")
        axis.set_xlim(-4.0, 4.0)
        axis.axvline(0.0, color="0.45", linestyle="-", linewidth=1)
        axis.grid(alpha=0.25)
        axis.legend(frameon=False, fontsize=9)

    fig.suptitle("Reproduction of van Hal et al. (1995), 0.1 M, $C_{Stern}$ = 0.2 F/m$^2$", fontsize=13)
    out_path.parent.mkdir(parents=True, exist_ok=True)
    fig.savefig(out_path, dpi=180)
    plt.close(fig)
from pathlib import Path

from electrolyte import Electrolyte
from surface import AmphotericSurface
from sweeps import plot_van_hal


OXIDES = [
    AmphotericSurface(pka=6.0, pkb=-2.0, site_density=5e18),
    AmphotericSurface(pka=10.0, pkb=6.0, site_density=8e18),
    AmphotericSurface(pka=4.0, pkb=2.0, site_density=10e18),
]
LABELS = ["SiO$_2$", "Al$_2$O$_3$", "Ta$_2$O$_5$"]


if __name__ == "__main__":
    plot_van_hal(
        OXIDES,
        LABELS,
        Electrolyte(ionic_strength=100.0, c_stern=0.2),
        Path("figs/van_hal_repro.png"),
    )

from math import isclose
from pathlib import Path

from electrolyte import Electrolyte
from surface import AmphotericSurface
from sweeps import plot_van_hal


def test_amphoteric_surface_charge_is_zero_at_pzc() -> None:
    surface = AmphotericSurface(pka=6.0, pkb=-2.0, site_density=5e18)

    assert isclose(surface.charge(2.0), 0.0, abs_tol=1e-12)


def test_amphoteric_surface_charge_is_negative_above_pzc() -> None:
    surface = AmphotericSurface(pka=6.0, pkb=-2.0, site_density=5e18)

    assert surface.charge(8.0) < 0.0


def test_plot_van_hal_writes_a_figure(tmp_path: Path) -> None:
    surfaces = [AmphotericSurface(pka=6.0, pkb=-2.0, site_density=5e18)]
    output = tmp_path / "nested" / "van_hal.png"

    plot_van_hal(surfaces, ["SiO$_2$"], Electrolyte(ionic_strength=100.0, c_stern=0.2), output)

    assert output.is_file()
    assert output.stat().st_size > 0

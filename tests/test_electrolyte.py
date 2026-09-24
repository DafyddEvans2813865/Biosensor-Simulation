from math import isclose
from pathlib import Path

import numpy as np
from scipy.constants import epsilon_0 as EPS_0


from electrolyte import Electrolyte
from charge import * 
from solver import solve
from analyte import *
from sweeps import *

SIO2 = Analyte(
    name="SiO2",
    sites=[Site(-2.0, SiteKind.ACIDIC), Site(6.0, SiteKind.BASIC)],
    site_density=1e18,  # m^-2, from 1e14 cm^-2
)

SALTS_MM = [1, 10, 100, 1000]


ELECTROLYTE_1_MM = Electrolyte(
	ionic_strength=1.0,  # 1 mM converted to 1 mol/m^3
	permittivity=80.0,
)

ELECTROLYTE_10_MM = Electrolyte(
	ionic_strength=10.0,  # 10 mM converted to 10 mol/m^3
	permittivity=80.0,
)

ELECTROLYTE_100_MM = Electrolyte(
	ionic_strength=100.0,  # 100 mM converted to 100 mol/m^3
	permittivity=80.0,
)


def test_thermal_voltage() -> None:
	assert isclose(ELECTROLYTE_10_MM.thermal_voltage, 0.0257, rel_tol=0.01)


def test_debye_length() -> None:
	assert isclose(ELECTROLYTE_10_MM.debye_length, 3.1e-9, rel_tol=0.05)

def test_debye_length_scales() -> None:
    assert isclose(ELECTROLYTE_1_MM.debye_length, ELECTROLYTE_100_MM.debye_length * 10.0, rel_tol=1e-9)


def test_surface_ph_shifts_one_unit_per_59_mv() -> None:
    v_t = ELECTROLYTE_10_MM.thermal_voltage
    assert isclose(surface_pH(7.0, 0.0592, v_t), 8.0, rel_tol=1e-3)

def test_psi_zero_at_pzc() -> None:
    #Should be 0 
    assert isclose(solve(SIO2, ELECTROLYTE_10_MM, 2.0), 0.0, abs_tol=1e-9)

def test_psi_is_positive_below_pzc() -> None:
    # Below pH 2 the surface is net protonated, positive.
    assert solve(SIO2, ELECTROLYTE_10_MM, 1.0) > 0.0


def test_psi_is_negative_above_pzc() -> None:
    # Above pH 2 the goes negative.
    assert solve(SIO2, ELECTROLYTE_10_MM, 7.0) < 0.0

def test_psi_decreases_with_ph() -> None:
    # More protons at low pH means a more positive surface, monotonically.
    values = [solve(SIO2, ELECTROLYTE_10_MM, ph) for ph in range(1, 13)]
    assert all(a > b for a, b in zip(values, values[1:]))
    
def test_slope_never_exceeds_nernst():
    ph, psi = ph_sweep(SIO2, ELECTROLYTE_10_MM)
    assert np.max(np.abs(slope(ph, psi))) < 0.059

def test_slope_is_sub_nernstian():
    ph, psi = ph_sweep(SIO2, ELECTROLYTE_10_MM)
    peak = np.max(np.abs(slope(ph, psi)))
    assert 0.030 < peak < 0.045

def test_q0_matches_linear_capacitance():
    e = ELECTROLYTE_10_MM
    c_dl = e.permittivity * EPS_0 / e.debye_length
    assert isclose(e.q0, 2 * c_dl * e.thermal_voltage, rel_tol=1e-9)

def test_sensitivity_falls_with_salt():
    peaks = []
    for e in (ELECTROLYTE_1_MM, ELECTROLYTE_10_MM, ELECTROLYTE_100_MM):
        ph, psi = ph_sweep(SIO2, e)
        peaks.append(np.max(np.abs(slope(ph, psi))))
    assert peaks[0] > peaks[1] > peaks[2]


def test_sensitivity_gap_is_much_larger_without_stern() -> None:
    no_stern = [Electrolyte(ionic_strength=i, c_stern=None) for i in SALTS_MM]
    with_stern = [Electrolyte(ionic_strength=i, c_stern=0.8) for i in SALTS_MM]

    def sensitivity(electrolytes):
        peaks = []
        for electrolyte in electrolytes:
            ph, psi = ph_sweep(SIO2, electrolyte)
            peaks.append(np.max(np.abs(slope(ph, psi))))
        return peaks

    no_stern_sens = sensitivity(no_stern)
    with_stern_sens = sensitivity(with_stern)

    no_stern_gap = no_stern_sens[0] - no_stern_sens[-1]
    with_stern_gap = with_stern_sens[0] - with_stern_sens[-1]

    assert no_stern_gap > 3.0 * with_stern_gap


for c_stern, tag in [(None, "no_stern"), (0.8, "stern")]:
    electrolytes = [Electrolyte(ionic_strength=i, c_stern=c_stern) for i in SALTS_MM]
    labels = [f"{i} mM" for i in SALTS_MM]
    plot_salt_sweep(SIO2, electrolytes, labels, Path(f"figs/salt_sweep_{tag}.png"))
from math import isclose

from electrolyte import Electrolyte




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

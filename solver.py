from math import cosh

from scipy.optimize import brentq

from charge import *
from surface import ChargedSurface

def diffuse_potential(electrolyte, psi_0: float, sigma_0: float) -> float:
    # Potential at the start of the diffuse layer
    if electrolyte.c_stern is None:
        return psi_0
    return psi_0 - sigma_0 / electrolyte.c_stern

def solve(surface, electrolyte, ph_bulk, boltzmann_ph: bool = True, stern: bool = True):
    vt = electrolyte.thermal_voltage

    def residual(psi_0):
        ph = surface_pH(ph_bulk, psi_0, vt) if boltzmann_ph else ph_bulk
        if hasattr(surface, "charge"):
            sigma_0 = surface.charge(ph)
        else:
            sigma_0 = surface_charge(surface, ph)
        psi_d = diffuse_potential(electrolyte, psi_0, sigma_0) if stern else psi_0
        return sigma_0 + electrolyte.diffuse_charge(psi_d)

    return brentq(residual, -2.0, 2.0, xtol=1e-12)
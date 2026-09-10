
from scipy.optimize import brentq
from charge import * 

def solve(analyte, electrolyte, ph_bulk):
    def residual(psi_0):
        ph_s = surface_pH(ph_bulk, psi_0, electrolyte.thermal_voltage)
        return surface_charge(analyte, ph_s) + electrolyte.diffuse_charge(psi_0)

    return brentq(residual, -0.5, 0.5)

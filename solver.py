
from scipy.optimize import brentq
from charge import * 

def solve(analyte, electrolyte, ph_bulk):
    def residual(psi_0):
        ph_s = surface_pH(ph_bulk, psi_0, electrolyte.thermal_voltage)
        sigma_0 = surface_charge(analyte,ph_s)

        if electrolyte.c_stern is not None:
            psi_d = psi_0 - sigma_0 / electrolyte.c_stern
        else:
            psi_d = psi_0

        return sigma_0 + electrolyte.diffuse_charge(psi_d)

    return brentq(residual, -0.5, 0.5)

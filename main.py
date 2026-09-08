from analyte import *
from charge import *

SIO2 = Analyte(
    name="SiO2",
    sites=[Site(-2.0, SiteKind.ACIDIC), Site(6.0, SiteKind.BASIC)],
    site_density=1e18,  # m^-2, from 1e14 cm^-2
)

charge = surface_charge(SIO2, 2.0)

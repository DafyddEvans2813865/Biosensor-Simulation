from analyte import Analyte, SiteKind
from constants import Q


def surface_charge(analyte: Analyte, ph_surface: float) -> float:
    # Net surface charge density in C/m^2 at a given surface pH.
    # A "site" is a titratable OH group on the oxide, holds or drops a proton.

    # Fraction of acidic sites that have deprotonated
    acidic_charge = sum(-1.0 / (1.0 + 10.0 ** (site.pk - ph_surface)) for site in analyte.sites if site.kind is SiteKind.ACIDIC)

    # Fraction of basic sites still holding a proton
    basic_charge = sum( 1.0 / (1.0 + 10.0 ** (ph_surface - site.pk)) for site in analyte.sites if site.kind is SiteKind.BASIC)

    return Q * analyte.site_density * (acidic_charge + basic_charge)

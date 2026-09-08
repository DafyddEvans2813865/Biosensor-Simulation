from dataclasses import dataclass
from enum import Enum


class SiteKind(Enum):
	ACIDIC = "acidic"
	BASIC = "basic"


@dataclass(frozen=True)
class Site:
	pk: float
	kind: SiteKind


@dataclass(frozen=True)
class Analyte:
	name: str
	sites: list[Site]
	site_density: float

SIO2 = Analyte(
    name="SiO2",
    sites=[Site(-2.0, SiteKind.ACIDIC), Site(6.0, SiteKind.BASIC)],
    site_density=1e18,  # m^-2, from 1e14 cm^-2
)
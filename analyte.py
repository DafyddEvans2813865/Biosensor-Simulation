from dataclasses import dataclass
from enum import Enum

from rdkit import Chem
from rdkit.Chem import Descriptors


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


def _sequence_tpsa(sequence: str) -> float:
	molecule = Chem.MolFromSequence(sequence)
	if molecule is None:
		raise ValueError(f"Could not parse peptide sequence: {sequence}")
	return Descriptors.TPSA(molecule)


def _tpsa_site_density(tpsa: float) -> float:
	# One molecule per TPSA footprint; 1 A^2 = 1e-20 m^2.
	return 1e20 / tpsa


SIO2 = Analyte(
	name="SiO2",
	sites=[Site(-2.0, SiteKind.ACIDIC), Site(6.0, SiteKind.BASIC)],
	site_density=1e18,
)

YTSF = Analyte("YTSF", [
	Site(2.2, SiteKind.ACIDIC),
	Site(9.1, SiteKind.BASIC),
	Site(10.1, SiteKind.ACIDIC),
], site_density=1e18)

YAAF = Analyte("YAAF", [
	Site(2.2, SiteKind.ACIDIC),
	Site(9.1, SiteKind.BASIC),
	Site(10.1, SiteKind.ACIDIC),
], site_density=1e18)

YTSF_TPSA = _sequence_tpsa("YTSF")
YAAF_TPSA = _sequence_tpsa("YAAF")

ytsf_TPSA = Analyte("YTSF_TPSA", YTSF.sites, _tpsa_site_density(YTSF_TPSA))
yaaf_TPSA = Analyte("YAAF_TPSA", YAAF.sites, _tpsa_site_density(YAAF_TPSA))
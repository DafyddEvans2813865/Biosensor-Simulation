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
	
	# N0 = 1e18 m^-2 (A0 = 100 Å^2);
	N0, A0 = 1e18, 100.0
	return N0 * min(1.0, A0 / tpsa)


SIO2 = Analyte(
	name="SiO2",
	sites=[Site(-2.0, SiteKind.ACIDIC), Site(6.0, SiteKind.BASIC)],
	site_density=1e18,
)

YTSF = Analyte("YTSF", [
	Site(9.21, SiteKind.BASIC),    # N-terminal amine (CRC)
	Site(10.46, SiteKind.ACIDIC),  # Tyr phenol (CRC)
], site_density=1e18)

YAAF = Analyte("YAAF", [
	Site(9.21, SiteKind.BASIC),    # N-terminal amine (CRC)
	Site(10.46, SiteKind.ACIDIC),  # Tyr phenol (CRC)
], site_density=1e18)

YTSF_TPSA = _sequence_tpsa("YTSF")
YAAF_TPSA = _sequence_tpsa("YAAF")

yaaf_TPSA = Analyte("YAAF_TPSA", [
	Site(7.15, SiteKind.BASIC),   # N-terminal amine (MolGpKa)
	Site(9.89, SiteKind.ACIDIC),  # Tyr phenol (MolGpKa)
], _tpsa_site_density(YAAF_TPSA))

ytsf_TPSA = Analyte("YTSF_TPSA", [
	Site(7.06, SiteKind.BASIC),   # N-terminal amine (MolGpKa)
	Site(9.82, SiteKind.ACIDIC),  # Tyr phenol (MolGpKa)
], _tpsa_site_density(YTSF_TPSA))
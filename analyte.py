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
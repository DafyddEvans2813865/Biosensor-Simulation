from analyte import Analyte, Site, SiteKind
from tpsa import sequence_tpsa, smiles_tpsa, tpsa_site_density


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


YAAF_TPSA = sequence_tpsa("YAAF")
yaaf_TPSA = Analyte("YAAF_TPSA", [
	Site(7.15, SiteKind.BASIC),   # N-terminal amine (MolGpKa)
	Site(9.89, SiteKind.ACIDIC),  # Tyr phenol (MolGpKa)
], tpsa_site_density(YAAF_TPSA))

YTSF_TPSA = sequence_tpsa("YTSF")
ytsf_TPSA = Analyte("YTSF_TPSA", [
	Site(7.06, SiteKind.BASIC),   # N-terminal amine (MolGpKa)
	Site(9.82, SiteKind.ACIDIC),  # Tyr phenol (MolGpKa)
], tpsa_site_density(YTSF_TPSA))


YTSF_NBOC_TPSA    = smiles_tpsa("CC(C)(C)OC(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](CO)C(=O)N[C@@H](Cc1ccccc1)C(N)=O")
ytsf_nboc = Analyte("YTSF_NBoc", [
	Site(9.95, SiteKind.ACIDIC),  # Tyr phenol (MolGpKa); amine blocked by Boc
], tpsa_site_density(YTSF_NBOC_TPSA))

YTSF_PS_TPSA      = smiles_tpsa("N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](COP(=O)(O)O)C(=O)N[C@@H](Cc1ccccc1)C(N)=O")
ytsf_ps = Analyte("YTSF_pS", [
	Site(7.06, SiteKind.BASIC),   # N-terminal amine (MolGpKa)
	Site(1.87, SiteKind.ACIDIC),  # phosphate, 1st step (MolGpKa)
	Site(7.2, SiteKind.ACIDIC),   # phosphate, 2nd step (CRC, as in the paper)
	Site(9.65, SiteKind.ACIDIC),  # Tyr phenol (MolGpKa)
], tpsa_site_density(YTSF_PS_TPSA))

YTSF_OGALNAC_TPSA = smiles_tpsa("N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H]([C@@H](C)OC1OC(CO)C(O)C(O)C1NC(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@@H](Cc1ccccc1)C(N)=O")
ytsf_ogalnac = Analyte("YTSF_OGalNAc", [
	Site(7.07, SiteKind.BASIC),   # N-terminal amine (MolGpKa)
	Site(9.81, SiteKind.ACIDIC),  # Tyr phenol (MolGpKa)
], tpsa_site_density(YTSF_OGALNAC_TPSA))

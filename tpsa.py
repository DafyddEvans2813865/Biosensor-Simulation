from rdkit import Chem
from rdkit.Chem import Descriptors


def sequence_tpsa(sequence: str) -> float:
	molecule = Chem.MolFromSequence(sequence)
	if molecule is None:
		raise ValueError(f"Could not parse peptide sequence: {sequence}")
	return Descriptors.TPSA(molecule)

def smiles_tpsa(smiles: str) -> float:
	molecule = Chem.MolFromSmiles(smiles)
	if molecule is None:
		raise ValueError(f"Could not parse SMILES: {smiles}")
	return Descriptors.TPSA(molecule)


def tpsa_site_density(tpsa: float) -> float:
	# N0 = 1e18 m^-2 (A0 = 100 Å^2);
	N0, A0 = 1e18, 100.0
	return N0 * min(1.0, A0 / tpsa)

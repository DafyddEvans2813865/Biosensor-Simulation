import os.path as osp
import torch

from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize

from utils.ionization_group import get_ionization_aid
from utils.descriptor import mol2vec
from utils.net import GCNNet


root = osp.dirname(osp.abspath(__file__)) 

def load_model(model_file):
    model = GCNNet()
    model.load_state_dict(torch.load(model_file, map_location="cpu"))
    model.eval()
    return model

def model_pred(mol, aid, model):
    with torch.no_grad():
        return float(model(mol2vec(mol, aid))[0][0])  #molecule -> graph -> network -> single pKa number

def predict(smiles):
    mol = Chem.MolFromSmiles(smiles)                                    # parse SMILES into a molecule
    mol = rdMolStandardize.Uncharger().uncharge(mol)                    # neutralise any charges
    mol = AllChem.AddHs(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))     # add acid sites (H's) - RDKit doesn't store H's as atoms by default,

    results = []
    for kind in ["acid", "base"]:                                       
        model = load_model(osp.join(root, f"../models/weight_{kind}.pth"))  
        for aid in get_ionization_aid(mol, acid_or_base=kind):         
            atom = mol.GetAtomWithIdx(aid)                              
            parent = atom.GetNeighbors()[0] if atom.GetSymbol() == "H" else atom  # if it's an H, use the atom it's bonded to
            on_ring = any(n.GetIsAromatic() for n in parent.GetNeighbors())       # next to an aromatic ring?
            results.append((kind, parent.GetSymbol(), on_ring, round(model_pred(mol, aid, model), 2)))  # save (type, element, aromatic?, pKa)
    return results

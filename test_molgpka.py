from rdkit import Chem
from rdkit.Chem import AllChem
from rdkit.Chem.MolStandardize import rdMolStandardize
import os.path as osp
import torch
from utils.ionization_group import get_ionization_aid
from utils.descriptor import mol2vec
from utils.net import GCNNet

root = osp.abspath("")

def load_model(model_file):
    model = GCNNet()
    model.load_state_dict(torch.load(model_file, map_location="cpu"))
    model.eval()
    return model

def model_pred(mol, aid, model):
    with torch.no_grad():
        return float(model(mol2vec(mol, aid))[0][0])

def predict(smiles):
    mol = Chem.MolFromSmiles(smiles)
    mol = rdMolStandardize.Uncharger().uncharge(mol)
    mol = AllChem.AddHs(Chem.MolFromSmiles(Chem.MolToSmiles(mol)))
    results = []
    for kind in ["acid", "base"]:
        model = load_model(osp.join(root, f"../models/weight_{kind}.pth"))
        for aid in get_ionization_aid(mol, acid_or_base=kind):
            results.append((kind, round(model_pred(mol, aid, model), 2)))
    return results

# Authors' own MolGpKa predictions, Pan et al. 2021, Figure 4
paper = {
    "butan-1-amine":          ("CCCCN",           10.72),
    "2-fluorobutan-1-amine":  ("NCC(F)CC",         8.69),
    "3-fluorobutan-1-amine":  ("NCCC(F)C",         9.45),
    "4-fluorobutan-1-amine":  ("NCCCCF",          10.38),
    "pyridine":               ("c1ccncc1",         5.08),
    "3-fluoropyridine":       ("Fc1cccnc1",        2.48),
    "3-acetylpyridine":       ("CC(=O)c1cccnc1",   3.20),
    "3-methoxypyridine":      ("COc1cccnc1",       4.63),
    "3-aminopyridine":        ("Nc1cccnc1",        5.60),
    "3-methylpyridine":       ("Cc1cccnc1",        5.32),
}
for name, (smi, ref) in paper.items():
    print(f"{name:24s} paper {ref:5.2f}  Output {predict(smi)}")
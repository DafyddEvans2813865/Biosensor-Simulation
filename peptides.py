from molgpka import predict
from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

# Authors' own MolGpKa predictions, Pan et al. 2021, Figure 4
PAPER = {
    "butan-1-amine":          ("CCCCN",           10.72),
    "2-fluorobutan-1-amine":  ("NCC(F)CC",         8.69),   # F close to N -> lower pKa
    "3-fluorobutan-1-amine":  ("NCCC(F)C",         9.45),
    "4-fluorobutan-1-amine":  ("NCCCCF",          10.38),   # F far from N -> small effect
    "pyridine":               ("c1ccncc1",         5.08),
    "3-fluoropyridine":       ("Fc1cccnc1",        2.48),   # electron-withdrawing -> lower
    "3-acetylpyridine":       ("CC(=O)c1cccnc1",   3.20),
    "3-methoxypyridine":      ("COc1cccnc1",       4.63),
    "3-aminopyridine":        ("Nc1cccnc1",        5.60),   # electron-donating -> higher
    "3-methylpyridine":       ("Cc1cccnc1",        5.32),
}

# Tyr-Thr-Ser-Phe and variants, all with an amidated C-terminus (C(N)=O)
PEPTIDES = {
    "YTSF":         "N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](CO)C(=O)N[C@@H](Cc1ccccc1)C(N)=O",
    "YAAF":         "N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H](C)C(=O)N[C@@H](C)C(=O)N[C@@H](Cc1ccccc1)C(N)=O",                     # Thr, Ser -> Ala
    "YTSF_NBoc":    "CC(C)(C)OC(=O)N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](CO)C(=O)N[C@@H](Cc1ccccc1)C(N)=O",  # Boc-capped N-terminus
    "YTSF_pS":      "N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H]([C@@H](C)O)C(=O)N[C@@H](COP(=O)(O)O)C(=O)N[C@@H](Cc1ccccc1)C(N)=O",       # phosphorylated Ser
    "YTSF_OGalNAc": "N[C@@H](Cc1ccc(O)cc1)C(=O)N[C@@H]([C@@H](C)OC1OC(CO)C(O)C(O)C1NC(C)=O)C(=O)N[C@@H](CO)C(=O)N[C@@H](Cc1ccccc1)C(N)=O",  # GalNAc sugar on Thr
}

for name, (smi, ref) in PAPER.items():                      
    print(f"{name:24s} Paper {ref:5.2f}  Output {predict(smi)}")


for name, smi in PEPTIDES.items():
    print(name, CalcMolFormula(Chem.MolFromSmiles(smi)))    
    print(predict(smi))                                     
from search.local_catalog import LOCAL_SCIENTIFIC_CATALOG as C
from rdkit import Chem
from chemistry_engine.formula_analyzer import parse_chemical_formula
from collections import Counter
print(len(C))
keys=Counter(e["key"] for e in C); print([k for k,v in keys.items() if v>1])
for e in C:
    s=e.get("canonical_smiles")
    if s and Chem.MolFromSmiles(s) is None: print("BAD SMILES",e["key"],s)
    if not parse_chemical_formula(e["symbol_or_formula"]) and e["entity_type"]!="biomolecule_protein": print("BAD FORMULA",e["key"],e["symbol_or_formula"])
    if e.get("molar_mass") is None: print("NOMASS", e["key"])

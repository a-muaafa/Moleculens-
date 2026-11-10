import io
from rdkit import Chem
from rdkit.Chem import Draw, rdMolDescriptors
from rdkit.Chem.Scaffolds import MurckoScaffold


def validate_smiles(smiles: str):
    if not smiles or not isinstance(smiles, str) or not smiles.strip():
        return None, None, "الرجاء إدخال صيغة SMILES غير فارغة."
    cleaned = smiles.strip()
    mol = Chem.MolFromSmiles(cleaned)
    if mol is None:
        return None, None, "صيغة SMILES غير صحيحة كيميائيا ولا يمكن تحويلها إلى بنية جزيئية."
    canonical = Chem.MolToSmiles(mol, canonical=True)
    return mol, canonical, None


def get_molecular_formula(mol) -> str:
    if mol is None:
        return ""
    return rdMolDescriptors.CalcMolFormula(mol)


def get_murcko_scaffold(mol):
    if mol is None:
        return None, ""
    scaf_smi = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
    if not scaf_smi:
        return mol, ""
    scaf_mol = Chem.MolFromSmiles(scaf_smi)
    return (scaf_mol if scaf_mol is not None else mol), scaf_smi


def render_mol_png(mol, size=(360, 260)) -> bytes:
    if mol is None:
        raise ValueError("لا يمكن رسم بنية جزيئية فارغة.")
    img = Draw.MolToImage(mol, size=size)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


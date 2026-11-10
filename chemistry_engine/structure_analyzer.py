from io import BytesIO
from typing import Any
from rdkit import Chem
from rdkit.Chem import AllChem, Crippen, Descriptors, Draw, Lipinski, rdMolDescriptors
from rdkit.Chem.Draw import rdMolDraw2D

FUNCTIONAL_GROUP_PATTERNS: list[tuple[str, str, str]] = [
    ("carboxylic_acid", "مجموعة كربوكسيل حمضية (-COOH)", "[CX3](=O)[OX2H1]"),
    ("ester", "مجموعة إستر (-COO-)", "[#6][CX3](=O)[OX2H0][#6]"),
    ("aldehyde", "مجموعة ألدهيد (-CHO)", "[CX3H1](=O)[#6,#1]"),
    ("ketone", "مجموعة كيتون / كربونيل (>C=O)", "[#6][CX3](=O)[#6]"),
    ("amide", "مجموعة أميد (-CONH-)", "[CX3](=[OX1])[NX3]"),
    ("alcohol_oh", "مجموعة هيدروكسيل كحولية (-OH)", "[#6][OX2H]"),
    ("phenol", "مجموعة فينول عطرية (Ar-OH)", "[OX2H][cX3]:[c]"),
    ("ether", "مجموعة إيثر (-O-)", "[OD2]([#6])[#6]"),
    ("primary_amine", "مجموعة أمين أولية (-NH2)", "[NX3;H2;!$(NC=O)]"),
    ("secondary_amine", "مجموعة أمين ثانوية (-NH-)", "[NX3;H1;!$(NC=O)]"),
    ("nitrile", "مجموعة نتريل / سيانو (-C≡N)", "[NX1]#[CX2]"),
    ("nitro", "مجموعة نيترو (-NO2)", "[$([NX3](=O)=O),$([NX3+](=O)[O-])]"),
    ("halide", "ذرة هالوجين عضوية (-F / -Cl / -Br / -I)", "[#6][F,Cl,Br,I]"),
    ("aromatic_ring", "حلقة عطرية (Aromatic Ring)", "a1aaaa1"),
    ("benzene_ring", "حلقة بنزين سداسية عطرية (Benzene Ring)", "c1ccccc1"),
    ("alkene", "رابطة ألكين ثنائية كربون-كربون (C=C)", "[CX3]=[CX3]"),
    ("alkyne", "رابطة ألكاين ثلاثية كربون-كربون (C≡C)", "[CX2]#[CX2]"),
    ("peroxide", "رابطة بيروكسيد (-O-O-)", "[OX2][OX2]"),
    ("phosphate", "مجموعة فوسفات (-PO4)", "[PX4](=O)([O,OH])([O,OH])[O,OH]"),
]

CPK_ELEMENT_COLORS: dict[str, str] = {
    "H": "#E2E8F0",
    "C": "#334155",
    "N": "#2563EB",
    "O": "#DC2626",
    "F": "#16A34A",
    "Cl": "#15803D",
    "Br": "#991B1B",
    "I": "#6B21A8",
    "S": "#CA8A04",
    "P": "#EA580C",
    "Na": "#7C3AED",
    "K": "#8B5CF6",
    "Ca": "#0D9488",
    "Mg": "#059669",
    "Fe": "#B45309",
    "Cu": "#D97706",
    "Zn": "#64748B",
}


def validate_and_parse_smiles(smiles: str):
    if not smiles or not isinstance(smiles, str):
        return None, None, "لا توجد صيغة SMILES."
    mol = Chem.MolFromSmiles(smiles.strip())
    if mol is None:
        return None, None, f"صيغة SMILES غير صالحة: {smiles}"
    canon = Chem.MolToSmiles(mol, canonical=True)
    return mol, canon, ""


def extract_rdkit_formula(smiles: str) -> str | None:
    mol, _, _ = validate_and_parse_smiles(smiles)
    if mol is None:
        return None
    return rdMolDescriptors.CalcMolFormula(mol)


def render_molecule_2d_png(smiles: str, size: tuple[int, int] = (420, 280)) -> bytes | None:
    mol, _, _ = validate_and_parse_smiles(smiles)
    if mol is None:
        return None
    AllChem.Compute2DCoords(mol)
    img = Draw.MolToImage(mol, size=size)
    buf = BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


def render_molecule_2d_svg(
    smiles: str,
    width: int = 480,
    height: int = 320,
    show_atom_indices: bool = False,
) -> str | None:
    mol, _, _ = validate_and_parse_smiles(smiles)
    if mol is None:
        return None
    mol_draw = Chem.Mol(mol)
    if mol_draw.GetNumHeavyAtoms() <= 6:
        mol_draw = Chem.AddHs(mol_draw)
    AllChem.Compute2DCoords(mol_draw)
    drawer = rdMolDraw2D.MolDraw2DSVG(int(width), int(height))
    opts = drawer.drawOptions()
    opts.addAtomIndices = bool(show_atom_indices)
    opts.addStereoAnnotation = True
    opts.bondLineWidth = 2.2
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol_draw)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def analyze_molecular_structure(smiles: str) -> dict[str, Any] | None:
    mol, canon, _ = validate_and_parse_smiles(smiles)
    if mol is None:
        return None

    mol_with_h = Chem.AddHs(mol)
    single_bonds = 0
    double_bonds = 0
    triple_bonds = 0
    aromatic_bonds = 0

    for bond in mol_with_h.GetBonds():
        btype = bond.GetBondType()
        if bond.GetIsAromatic() or btype == Chem.BondType.AROMATIC:
            aromatic_bonds += 1
        elif btype == Chem.BondType.SINGLE:
            single_bonds += 1
        elif btype == Chem.BondType.DOUBLE:
            double_bonds += 1
        elif btype == Chem.BondType.TRIPLE:
            triple_bonds += 1

    ring_info = mol.GetRingInfo()
    atom_rings = ring_info.AtomRings()
    ring_sizes = [len(r) for r in atom_rings]
    aromatic_rings = int(rdMolDescriptors.CalcNumAromaticRings(mol))
    aliphatic_rings = int(rdMolDescriptors.CalcNumAliphaticRings(mol))

    detected_groups: list[dict[str, Any]] = []
    for fg_key, fg_label_ar, smarts in FUNCTIONAL_GROUP_PATTERNS:
        pat = Chem.MolFromSmarts(smarts)
        if pat is not None:
            matches = mol.GetSubstructMatches(pat)
            if matches:
                if fg_key == "alcohol_oh":
                    cooh_pat = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
                    cooh_atoms = {idx for m in (mol.GetSubstructMatches(cooh_pat) if cooh_pat else []) for idx in m}
                    filtered = [m for m in matches if not any(idx in cooh_atoms for idx in m)]
                    if not filtered:
                        continue
                    matches = filtered
                detected_groups.append(
                    {
                        "key": fg_key,
                        "label_ar": fg_label_ar,
                        "count": len(matches),
                    }
                )

    chiral_centers = Chem.FindMolChiralCenters(mol, includeUnassigned=True)
    stereo_bonds = []
    for b in mol.GetBonds():
        if b.GetBondType() == Chem.BondType.DOUBLE and b.GetStereo() not in (Chem.BondStereo.STEREONONE, Chem.BondStereo.STEREOANY):
            stereo_bonds.append(str(b.GetStereo()))

    has_stereo = len(chiral_centers) > 0 or len(stereo_bonds) > 0
    if chiral_centers:
        centers_str = ", ".join(f"ذرة رقم {idx} ({cfg})" for idx, cfg in chiral_centers)
        stereo_summary_ar = f"يحتوي الجزيء على {len(chiral_centers)} مركز كيرالي فراغي: {centers_str}."
    elif stereo_bonds:
        stereo_summary_ar = f"يحتوي الجزيء على تماكب هندسي في الروابط الثنائية ({', '.join(stereo_bonds)})."
    else:
        stereo_summary_ar = "جزيء غير كيرالي (لا يحتوي على مراكز كيرالية غير متناظرة)."

    return {
        "canonical_smiles": canon,
        "rdkit_formula": rdMolDescriptors.CalcMolFormula(mol),
        "heavy_atoms_count": int(mol.GetNumHeavyAtoms()),
        "total_atoms_with_h": int(mol_with_h.GetNumAtoms()),
        "bonds": {
            "single": single_bonds,
            "double": double_bonds,
            "triple": triple_bonds,
            "aromatic": aromatic_bonds,
            "total": int(mol_with_h.GetNumBonds()),
            "summary_ar": (
                f"روابط أحادية: {single_bonds} | روابط ثنائية: {double_bonds} | "
                f"روابط ثلاثية: {triple_bonds} | روابط عطرية: {aromatic_bonds}"
            ),
        },
        "rings": {
            "total_rings": len(atom_rings),
            "aromatic_rings": aromatic_rings,
            "aliphatic_rings": aliphatic_rings,
            "ring_sizes": ring_sizes,
            "summary_ar": (
                f"إجمالي الحلقات: {len(atom_rings)} (عطرية: {aromatic_rings}، أليفاتية: {aliphatic_rings})"
                if atom_rings
                else "جزيء مفتوح السلسلة أو لا حلقي (0 حلقات)"
            ),
        },
        "functional_groups": detected_groups,
        "functional_groups_summary_ar": (
            " ، ".join(f"{g['label_ar']} (×{g['count']})" for g in detected_groups)
            if detected_groups
            else "مركب غير عضوي أو هيدروكربوني بسيط بدون مجموعات وظيفية معقدة"
        ),
        "stereochemistry": {
            "has_stereochemistry": has_stereo,
            "chiral_centers_count": len(chiral_centers),
            "chiral_centers": chiral_centers,
            "stereo_double_bonds_count": len(stereo_bonds),
            "summary_ar": stereo_summary_ar,
        },
    }


def generate_3d_conformer_data(smiles: str) -> dict[str, Any] | None:
    mol, canon, _ = validate_and_parse_smiles(smiles)
    if mol is None:
        return None

    mol3d = Chem.AddHs(mol)
    params = AllChem.ETKDGv3()
    params.randomSeed = 42
    res = AllChem.EmbedMolecule(mol3d, params)
    if res != 0:
        res = AllChem.EmbedMolecule(mol3d, useRandomCoords=True, randomSeed=42)
    if res != 0 or mol3d.GetNumConformers() == 0:
        return None

    try:
        if AllChem.MMFFHasAllMoleculeParams(mol3d):
            AllChem.MMFFOptimizeMolecule(mol3d, maxIters=200)
        else:
            AllChem.UFFOptimizeMolecule(mol3d, maxIters=200)
    except Exception:
        pass

    conf = mol3d.GetConformer()
    atoms_3d = []
    for idx in range(mol3d.GetNumAtoms()):
        atom = mol3d.GetAtomWithIdx(idx)
        pos = conf.GetAtomPosition(idx)
        sym = atom.GetSymbol()
        atoms_3d.append(
            {
                "index": idx,
                "symbol": sym,
                "x": round(float(pos.x), 4),
                "y": round(float(pos.y), 4),
                "z": round(float(pos.z), 4),
                "color": CPK_ELEMENT_COLORS.get(sym, "#0F766E"),
                "radius": 8 if sym == "H" else 14,
            }
        )

    bonds_3d = []
    for bond in mol3d.GetBonds():
        b_idx = bond.GetBeginAtomIdx()
        e_idx = bond.GetEndAtomIdx()
        b_order = int(round(bond.GetBondTypeAsDouble()))
        bonds_3d.append(
            {
                "begin": b_idx,
                "end": e_idx,
                "order": max(1, b_order),
                "is_aromatic": bool(bond.GetIsAromatic()),
            }
        )

    mol_block = Chem.MolToMolBlock(mol3d)
    return {
        "canonical_smiles": canon,
        "mol_block": mol_block,
        "atoms": atoms_3d,
        "bonds": bonds_3d,
        "atom_count": len(atoms_3d),
        "bond_count": len(bonds_3d),
    }


def structure_matches_formula(smiles: str, formula: str) -> bool:
    import re

    from chemistry_engine.formula_analyzer import normalize_chemical_formula

    rd_formula = extract_rdkit_formula(smiles)
    if not rd_formula or not formula:
        return False
    rd_clean = re.sub(r"[+\-]\d*$", "", rd_formula)
    a = normalize_chemical_formula(rd_clean)
    b = normalize_chemical_formula(formula)
    return a is not None and a == b


def compute_rdkit_molecular_metrics(smiles: str) -> dict | None:
    mol, canon, _ = validate_and_parse_smiles(smiles)
    if mol is None:
        return None
    mw = float(Descriptors.MolWt(mol))
    logp = float(Crippen.MolLogP(mol))
    hbd = int(Lipinski.NumHDonors(mol))
    hba = int(Lipinski.NumHAcceptors(mol))
    rot = int(Lipinski.NumRotatableBonds(mol))
    tpsa = float(rdMolDescriptors.CalcTPSA(mol))
    heavy_atoms = int(mol.GetNumHeavyAtoms())
    rings = int(rdMolDescriptors.CalcNumRings(mol))
    formula = rdMolDescriptors.CalcMolFormula(mol)

    lipinski_violations = sum([mw > 500.0, logp > 5.0, hbd > 5, hba > 10])
    return {
        "canonical_smiles": canon,
        "rdkit_formula": formula,
        "molecular_weight": round(mw, 3),
        "logp": round(logp, 2),
        "hbd": hbd,
        "hba": hba,
        "rotatable_bonds": rot,
        "tpsa": round(tpsa, 2),
        "heavy_atoms": heavy_atoms,
        "ring_count": rings,
        "lipinski_violations": lipinski_violations,
        "lipinski_compliant": lipinski_violations <= 1,
    }

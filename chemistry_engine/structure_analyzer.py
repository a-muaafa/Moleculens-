from io import BytesIO
import re
from typing import Any

try:
    from rdkit import Chem
    from rdkit.Chem import AllChem, Crippen, Descriptors, Draw, Lipinski, rdMolDescriptors
    from rdkit.Chem.Draw import rdMolDraw2D
    HAS_RDKIT = True
except ImportError:
    Chem = None  # type: ignore
    AllChem = None  # type: ignore
    Crippen = None  # type: ignore
    Descriptors = None  # type: ignore
    Draw = None  # type: ignore
    Lipinski = None  # type: ignore
    rdMolDescriptors = None  # type: ignore
    rdMolDraw2D = None  # type: ignore
    HAS_RDKIT = False

FUNCTIONAL_GROUP_PATTERNS: list[tuple[str, str, str]] = [
    ("carboxylic_acid", "مجموعة كربوكسيل حمضية (-COOH | Carboxylic Acid)", "[CX3](=O)[OX2H1]"),
    ("ester", "مجموعة إستر (-COO- | Ester)", "[#6][CX3](=O)[OX2H0][#6]"),
    ("aldehyde", "مجموعة ألدهيد (-CHO | Aldehyde)", "[CX3H1](=O)[#6,#1]"),
    ("ketone", "مجموعة كيتون / كربونيل (>C=O | Ketone)", "[#6][CX3](=O)[#6]"),
    ("amide", "مجموعة أميد (-CONH- | Amide)", "[CX3](=[OX1])[NX3]"),
    ("alcohol_oh", "مجموعة هيدروكسيل كحولية (-OH | Alcohol)", "[#6][OX2H]"),
    ("phenol", "مجموعة فينول عطرية (Ar-OH | Phenol)", "[OX2H][cX3]:[c]"),
    ("ether", "مجموعة إيثر (-O- | Ether)", "[OD2]([#6])[#6]"),
    ("primary_amine", "مجموعة أمين أولية (-NH2 | Primary Amine)", "[NX3;H2;!$(NC=O)]"),
    ("secondary_amine", "مجموعة أمين ثانوية (-NH- | Secondary Amine)", "[NX3;H1;!$(NC=O)]"),
    ("nitrile", "مجموعة نتريل / سيانو (-C≡N | Nitrile)", "[NX1]#[CX2]"),
    ("nitro", "مجموعة نيترو (-NO2 | Nitro)", "[$([NX3](=O)=O),$([NX3+](=O)[O-])]"),
    ("halide", "ذرة هالوجين عضوية (-F / -Cl / -Br / -I | Halide)", "[#6][F,Cl,Br,I]"),
    ("aromatic_ring", "حلقة عطرية (Aromatic Ring)", "a1aaaa1"),
    ("benzene_ring", "حلقة بنزين سداسية عطرية (Benzene Ring)", "c1ccccc1"),
    ("alkene", "رابطة ألكين ثنائية كربون-كربون (C=C | Alkene)", "[CX3]=[CX3]"),
    ("alkyne", "رابطة ألكاين ثلاثية كربون-كربون (C≡C | Alkyne)", "[CX2]#[CX2]"),
    ("peroxide", "رابطة بيروكسيد (-O-O- | Peroxide)", "[OX2][OX2]"),
    ("phosphate", "مجموعة فوسفات (-PO4 | Phosphate)", "[PX4](=O)([O,OH])([O,OH])[O,OH]"),
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
    s_clean = smiles.strip()
    if not HAS_RDKIT or Chem is None:
        return None, None, "محرك التحليل البنيوي RDKit غير متوفر للتحقق من البنية الجزيئية."
    mol = Chem.MolFromSmiles(s_clean)
    if mol is None:
        return None, None, f"صيغة SMILES غير صالحة: {smiles}"
    canon = Chem.MolToSmiles(mol, canonical=True)
    return mol, canon, ""


def compute_inchikey(smiles: str) -> str | None:
    if not smiles or not isinstance(smiles, str):
        return None
    if HAS_RDKIT and Chem is not None:
        mol = Chem.MolFromSmiles(smiles.strip())
        if mol is None:
            return None
        try:
            return Chem.MolToInchiKey(mol)
        except Exception:
            return None
    return None


def extract_rdkit_formula(smiles: str) -> str | None:
    if not HAS_RDKIT or rdMolDescriptors is None:
        return None
    mol, _, _ = validate_and_parse_smiles(smiles)
    if mol is None:
        return None
    return rdMolDescriptors.CalcMolFormula(mol)


def structure_matches_formula(smiles: str, formula: str) -> bool:
    from chemistry_engine.formula_analyzer import normalize_chemical_formula

    if not smiles or not formula:
        return False
    if not HAS_RDKIT or rdMolDescriptors is None:
        return False
    rd_formula = extract_rdkit_formula(smiles)
    if not rd_formula:
        return False
    rd_clean = re.sub(r"[+\-]\d*$", "", rd_formula)
    a = normalize_chemical_formula(rd_clean)
    b = normalize_chemical_formula(formula)
    return a is not None and a == b


def verify_substance_structure_identity(
    entity: dict[str, Any],
    expected_query: str = "",
) -> dict[str, Any]:
    smiles = (entity.get("canonical_smiles") or "").strip()
    formula = (entity.get("symbol_or_formula") or "").strip()
    name_ar = entity.get("name_ar", "")
    name_en = entity.get("name_en", "")
    entity_type = entity.get("entity_type", "")

    if entity_type == "chemical_element":
        from chemistry_engine.periodic_table import get_element_basic_info

        el_info = get_element_basic_info(formula)
        if el_info is not None:
            return {
                "is_verified": True,
                "verification_status": "verified",
                "formula_matches": True,
                "canonical_smiles": smiles,
                "inchikey": compute_inchikey(smiles),
                "status_ar": f"تم التحقق من هوية العنصر الكيميائي ({name_ar} — {formula}) عبر الجدول الدوري القياسي (IUPAC).",
            }
        return {
            "is_verified": False,
            "verification_status": "unconfirmed",
            "formula_matches": False,
            "canonical_smiles": "",
            "inchikey": None,
            "status_ar": f"النتيجة غير مؤكدة: تعذر التحقق من رمز العنصر ({formula}) في الجدول الدوري القياسي.",
        }

    if not smiles or not formula:
        return {
            "is_verified": False,
            "verification_status": "unconfirmed",
            "formula_matches": False,
            "canonical_smiles": "",
            "inchikey": None,
            "status_ar": "النتيجة غير مؤكدة (Unconfirmed Identity): لا تتوفر بنية SMILES أو صيغة مرجعية كافية للتحقق البنيوي.",
        }

    if not HAS_RDKIT or Chem is None or rdMolDescriptors is None:
        return {
            "is_verified": False,
            "verification_status": "unconfirmed",
            "formula_matches": False,
            "canonical_smiles": smiles,
            "inchikey": None,
            "status_ar": "النتيجة غير مؤكدة (Unconfirmed Identity): يتعذر التحقق الفعلي من تطابق الصيغة والبنية الجزيئية (SMILES) لعدم توفر محرك RDKit.",
        }

    mol, canon, err = validate_and_parse_smiles(smiles)
    if mol is None:
        return {
            "is_verified": False,
            "verification_status": "mismatch",
            "formula_matches": False,
            "canonical_smiles": "",
            "inchikey": None,
            "status_ar": f"فشل التحقق البنيوي: {err}",
        }

    formula_ok = structure_matches_formula(smiles, formula)
    inchikey_val = compute_inchikey(smiles)
    ids = entity.get("identifiers", {}) or {}
    has_ref_id = bool(
        ids.get("pubchem_cid")
        or entity.get("pubchem_cid")
        or ids.get("chembl_id")
        or entity.get("chembl_id")
        or ids.get("chebi_id")
        or entity.get("chebi_id")
        or ids.get("uniprot_id")
        or entity.get("uniprot_id")
        or entity.get("inchikey")
        or entity.get("key")
    )

    is_verified = bool(formula_ok and has_ref_id)
    if is_verified:
        ver_status = "verified"
        status_ar = (
            f"تحقق مزدوج ناجح (Identity & Structure Verified): البنية الجزيئية (SMILES) تطابق الصيغة الجزيئية ({formula}) "
            f"وتعود فعلياً للمادة الموثقة ({name_ar} | {name_en})."
        )
    elif formula_ok:
        ver_status = "unconfirmed"
        status_ar = (
            f"النتيجة غير مؤكدة بالكامل: تطابق الصيغة الجزيئية ({formula}) متحقق، لكن يلزم معرف مرجعي للتأكد من عدم الخلط مع المصاوغات (Isomers)."
        )
    else:
        ver_status = "mismatch"
        status_ar = (
            f"تنبيه عدم تطابق الهوية (Identity Mismatch Warning): نجاح قراءة SMILES لا يعني عودتها للمادة المطلوبة؛ "
            f"الصيغة المستخرجة من البنية لا تطابق الصيغة المسجلة ({formula})."
        )

    return {
        "is_verified": is_verified,
        "verification_status": ver_status,
        "formula_matches": formula_ok,
        "canonical_smiles": canon or smiles,
        "inchikey": inchikey_val,
        "status_ar": status_ar,
    }


def render_molecule_2d_png(smiles: str, size: tuple[int, int] = (480, 320)) -> bytes | None:
    if not HAS_RDKIT or AllChem is None or Draw is None:
        return None
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
    width: int = 640,
    height: int = 420,
    show_atom_indices: bool = False,
) -> str | None:
    if not HAS_RDKIT or Chem is None or AllChem is None or rdMolDraw2D is None:
        return None
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
    opts.bondLineWidth = 2.4
    opts.padding = 0.08
    rdMolDraw2D.PrepareAndDrawMolecule(drawer, mol_draw)
    drawer.FinishDrawing()
    return drawer.GetDrawingText()


def analyze_molecular_structure(smiles: str) -> dict[str, Any] | None:
    if not HAS_RDKIT or Chem is None or rdMolDescriptors is None:
        return None
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
        stereo_summary_ar = f"يحتوي الجزيء على {len(chiral_centers)} مركز كيرالي فراغي (Chiral Center): {centers_str}."
    elif stereo_bonds:
        stereo_summary_ar = f"يحتوي الجزيء على تماكب هندسي في الروابط الثنائية ({', '.join(stereo_bonds)})."
    else:
        stereo_summary_ar = "جزيء غير كيرالي (Achiral — لا يحتوي على مراكز كيرالية غير متناظرة)."

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
                f"روابط أحادية (Single): {single_bonds} | روابط ثنائية (Double): {double_bonds} | "
                f"روابط ثلاثية (Triple): {triple_bonds} | روابط عطرية (Aromatic): {aromatic_bonds}"
            ),
        },
        "rings": {
            "total_rings": len(atom_rings),
            "aromatic_rings": aromatic_rings,
            "aliphatic_rings": aliphatic_rings,
            "ring_sizes": ring_sizes,
            "summary_ar": (
                f"إجمالي الحلقات (Total Rings): {len(atom_rings)} (عطرية Aromatic: {aromatic_rings}، أليفاتية Aliphatic: {aliphatic_rings})"
                if atom_rings
                else "جزيء مفتوح السلسلة أو لا حلقي (Acyclic — 0 حلقات)"
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
    if not HAS_RDKIT or Chem is None or AllChem is None:
        return None
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


def compute_rdkit_molecular_metrics(smiles: str) -> dict[str, Any] | None:
    if not HAS_RDKIT or Descriptors is None or Crippen is None or Lipinski is None or rdMolDescriptors is None:
        return None
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

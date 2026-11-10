import re
from chemistry_engine.periodic_table import (
    ATOMIC_MASSES,
    ELEMENT_NAMES_AR,
    ELEMENT_NAMES_EN,
    VALID_ELEMENTS,
)

SUBSCRIPT_MAP = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")
SUPERSCRIPT_MAP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻", "0123456789+-")
SINGLE_LETTER_ORGANIC = {"C", "H", "N", "O", "P", "S", "F", "K", "I", "B"}
LOWER_TO_SYMBOL = {sym.lower(): sym for sym in VALID_ELEMENTS}

SPECIAL_FORMULA_CASE_MAP = {
    "cl2": "Cl2", "br2": "Br2", "hcl": "HCl", "hbr": "HBr", "nacl": "NaCl",
    "kcl": "KCl", "naoh": "NaOH", "koh": "KOH", "lioh": "LiOH",
    "ca(oh)2": "Ca(OH)2", "mg(oh)2": "Mg(OH)2", "ba(oh)2": "Ba(OH)2",
    "al(oh)3": "Al(OH)3", "fe(oh)3": "Fe(OH)3", "cu(oh)2": "Cu(OH)2",
    "cao": "CaO", "mgo": "MgO", "cuo": "CuO", "zno": "ZnO", "al2o3": "Al2O3",
    "fe2o3": "Fe2O3", "fe3o4": "Fe3O4", "sio2": "SiO2",
    "caco3": "CaCO3", "na2co3": "Na2CO3", "nahco3": "NaHCO3",
    "cuso4": "CuSO4", "feso4": "FeSO4", "znso4": "ZnSO4", "na2so4": "Na2SO4", "baso4": "BaSO4",
    "agno3": "AgNO3", "agcl": "AgCl", "kno3": "KNO3", "nano3": "NaNO3",
    "cacl2": "CaCl2", "mgcl2": "MgCl2", "fecl2": "FeCl2", "fecl3": "FeCl3", "zncl2": "ZnCl2", "bacl2": "BaCl2", "alcl3": "AlCl3",
    "kmno4": "KMnO4", "k2cr2o7": "K2Cr2O7", "kclo3": "KClO3",
    "nh4cl": "NH4Cl", "nh4oh": "NH4OH", "ch3coona": "CH3COONa", "na3po4": "Na3PO4",
    "cu(no3)2": "Cu(NO3)2", "co(nh2)2": "CO(NH2)2",
}


def replace_subscripts(text: str) -> str:
    if not text:
        return ""
    return text.translate(SUBSCRIPT_MAP)


def clean_formula_token(formula_str: str) -> str:
    if not formula_str:
        return ""
    s = replace_subscripts(formula_str.strip()).replace(" ", "")
    s = re.sub(r"[\^]?[0-9]*[+\-]$", "", s)
    return s


def exact_formula_string_match(formula_a: str, formula_b: str) -> bool:
    if not formula_a or not formula_b:
        return False
    ca = replace_subscripts(formula_a.strip()).replace(" ", "").lower()
    cb = replace_subscripts(formula_b.strip()).replace(" ", "").lower()
    return ca == cb


def _parse_group(formula: str) -> dict[str, int] | None:
    stack: list[dict[str, int]] = [{}]
    i = 0
    n = len(formula)
    while i < n:
        ch = formula[i]
        if ch == "(":
            stack.append({})
            i += 1
        elif ch == ")":
            if len(stack) < 2:
                return None
            top = stack.pop()
            i += 1
            m_start = i
            while i < n and formula[i].isdigit():
                i += 1
            mult = int(formula[m_start:i]) if i > m_start else 1
            if mult <= 0 or mult > 500:
                return None
            for el, cnt in top.items():
                stack[-1][el] = stack[-1].get(el, 0) + cnt * mult
        elif ch.isupper():
            el_start = i
            i += 1
            if i < n and formula[i].islower():
                i += 1
            elem = formula[el_start:i]
            if elem not in VALID_ELEMENTS:
                return None
            m_start = i
            while i < n and formula[i].isdigit():
                i += 1
            cnt = int(formula[m_start:i]) if i > m_start else 1
            if cnt <= 0 or cnt > 1000:
                return None
            stack[-1][elem] = stack[-1].get(elem, 0) + cnt
        else:
            return None

    if len(stack) != 1 or not stack[0]:
        return None
    return stack[0]


def parse_chemical_formula(formula_str: str) -> dict[str, int] | None:
    if not formula_str or not isinstance(formula_str, str):
        return None
    cleaned = clean_formula_token(formula_str)
    if not cleaned or not re.fullmatch(r"[A-Za-z0-9()]+", cleaned):
        return None

    low = cleaned.lower()
    if low in SPECIAL_FORMULA_CASE_MAP:
        return _parse_group(SPECIAL_FORMULA_CASE_MAP[low])

    if low in LOWER_TO_SYMBOL and len(cleaned) <= 2:
        sym = LOWER_TO_SYMBOL[low]
        if cleaned == sym:
            return {sym: 1}

    candidate = cleaned
    if candidate.islower():
        if all(c.isdigit() or c in "()" or c.upper() in SINGLE_LETTER_ORGANIC for c in candidate):
            candidate = candidate.upper()

    return _parse_group(candidate)


def format_hill_formula(counts: dict[str, int]) -> str:
    if not counts:
        return ""
    parts = []
    if "C" in counts:
        c_cnt = counts["C"]
        parts.append("C" if c_cnt == 1 else f"C{c_cnt}")
        if "H" in counts:
            h_cnt = counts["H"]
            parts.append("H" if h_cnt == 1 else f"H{h_cnt}")
        for elem in sorted(k for k in counts if k not in {"C", "H"}):
            cnt = counts[elem]
            parts.append(elem if cnt == 1 else f"{elem}{cnt}")
    else:
        for elem in sorted(counts.keys()):
            cnt = counts[elem]
            parts.append(elem if cnt == 1 else f"{elem}{cnt}")
    return "".join(parts)


def normalize_chemical_formula(formula_str: str) -> str | None:
    counts = parse_chemical_formula(formula_str)
    if not counts:
        return None
    return format_hill_formula(counts)


def formulas_are_equivalent(formula_a: str, formula_b: str) -> bool:
    if not formula_a or not formula_b:
        return False
    if exact_formula_string_match(formula_a, formula_b):
        return True
    counts_a = parse_chemical_formula(formula_a)
    counts_b = parse_chemical_formula(formula_b)
    if counts_a is not None and counts_b is not None:
        return counts_a == counts_b
    return False


def calculate_molar_mass(formula_str: str) -> float | None:
    counts = parse_chemical_formula(formula_str)
    if not counts:
        return None
    total = sum(ATOMIC_MASSES[el] * cnt for el, cnt in counts.items())
    return round(total, 4)


def _format_arabic_atom_count(count: int, el_name_ar: str) -> str:
    if count == 1:
        return f"ذرة {el_name_ar} واحدة"
    if count == 2:
        return f"ذرتان {el_name_ar}"
    if 3 <= count <= 10:
        return f"{count} ذرات {el_name_ar}"
    return f"{count} ذرة {el_name_ar}"


def analyze_formula_composition(formula_str: str) -> dict:
    counts = parse_chemical_formula(formula_str)
    if not counts:
        return {
            "is_valid": False,
            "molar_mass": None,
            "hill_formula": "",
            "total_atoms": 0,
            "elements_count": 0,
            "components": [],
            "summary_ar": "صيغة غير قابلة للتحليل الذري المباشر.",
        }

    total_mass = sum(ATOMIC_MASSES[el] * cnt for el, cnt in counts.items())
    total_atoms = sum(counts.values())
    components = []
    descriptions_ar = []

    for el, cnt in counts.items():
        atom_mass = ATOMIC_MASSES[el]
        sub_mass = atom_mass * cnt
        pct = (sub_mass / total_mass * 100.0) if total_mass > 0 else 0.0
        name_ar = ELEMENT_NAMES_AR.get(el, el)
        name_en = ELEMENT_NAMES_EN.get(el, el)
        desc_ar = _format_arabic_atom_count(cnt, name_ar)
        descriptions_ar.append(desc_ar)
        components.append(
            {
                "symbol": el,
                "name_ar": name_ar,
                "name_en": name_en,
                "count": cnt,
                "description_ar": desc_ar,
                "atomic_mass": round(atom_mass, 4),
                "subtotal_mass": round(sub_mass, 4),
                "mass_percentage": round(pct, 2),
            }
        )

    summary_ar = " + ".join(descriptions_ar)
    return {
        "is_valid": True,
        "molar_mass": round(total_mass, 3),
        "hill_formula": format_hill_formula(counts),
        "total_atoms": total_atoms,
        "elements_count": len(counts),
        "atom_counts": counts,
        "components": components,
        "summary_ar": summary_ar,
    }

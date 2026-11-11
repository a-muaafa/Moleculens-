import re
from chemistry_engine.periodic_table import (
    ATOMIC_MASSES,
    ELEMENT_NAMES_AR,
    ELEMENT_NAMES_EN,
    VALID_ELEMENTS,
)

SUBSCRIPT_MAP = str.maketrans("₀₁₂₃₄₅₆₇₈₉", "0123456789")
TO_SUBSCRIPT_MAP = str.maketrans("0123456789", "₀₁₂₃₄₅₆₇₈₉")
SUPERSCRIPT_MAP = str.maketrans("⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻", "0123456789+-")
TO_SUPERSCRIPT_MAP = str.maketrans("0123456789+-", "⁰¹²³⁴⁵⁶⁷⁸⁹⁺⁻")
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

# Standard states at STP/SATP (25 °C, 1 atm) for chemical species
GASEOUS_FORMULAS_STP = {
    "H2", "He", "N2", "O2", "F2", "Ne", "Cl2", "Ar", "Kr", "Xe", "Rn", "Og",
    "H", "N", "O", "F", "Cl",
    "CO", "CO2", "NO", "NO2", "N2O", "SO2", "H2S", "NH3", "HCl", "HBr", "HF", "HI",
    "CH4", "C2H6", "C3H8", "C4H10", "C2H4", "C2H2", "CH3OCH3", "HCHO", "O3", "CH3Cl",
}

LIQUID_FORMULAS_STP = {
    "H2O", "H2O2", "Br2", "Br", "Hg",
    "CH3COOH", "HCOOH", "HCOOCH3", "H2SO4", "HNO3", "C3H6O3",
    "C6H6", "C7H8", "C6H5NH2", "C6H5NO2", "CH3OH", "C2H5OH",
    "(CH3)2CHOH", "CH3CH2CH2OH", "C3H8O3", "CH3CHO", "C6H5CHO",
    "CH3COCH3", "CH3COCH2CH3", "CH3COOCH2CH3", "CH3COOCH3", "C2H4Br2",
}


def replace_subscripts(text: str) -> str:
    if not text:
        return ""
    return text.translate(SUBSCRIPT_MAP)


def format_chemical_formula_display(formula_str: str) -> str:
    """Formats chemical formulas using proper Unicode subscripts and ionic superscripts (e.g., H2O -> H₂O, SO4^2- -> SO₄²⁻)."""
    if not formula_str or not isinstance(formula_str, str):
        return ""
    raw = formula_str.strip()
    if "Protein" in raw or "aa)" in raw:
        return raw
    plain = replace_subscripts(raw)
    # Handle explicit ionic charge patterns at the end like ^2-, ^3-, 2+, +, -
    charge_suffix = ""
    m_caret = re.search(r"\^(\d*[+\-])$", plain)
    if m_caret:
        charge_suffix = m_caret.group(1).translate(TO_SUPERSCRIPT_MAP)
        plain = plain[: m_caret.start()]
    else:
        m_ion = re.search(r"([+\-])$", plain)
        if m_ion:
            # Check if like Ca2+ vs H3O+ vs NH4+ vs SO4^2-
            if plain in {"Ca2+", "Mg2+", "Fe2+", "Fe3+", "Cu2+", "Zn2+", "Ba2+", "Al3+"}:
                charge_suffix = plain[-2:].translate(TO_SUPERSCRIPT_MAP)
                plain = plain[:-2]
            else:
                charge_suffix = plain[-1].translate(TO_SUPERSCRIPT_MAP)
                plain = plain[:-1]

    # Convert digits following letters or closing parenthesis to subscripts
    def _sub_repl(match: re.Match) -> str:
        prefix = match.group(1)
        digits = match.group(2).translate(TO_SUBSCRIPT_MAP)
        return f"{prefix}{digits}"

    formatted = re.sub(r"([A-Za-z\)])(\d+)", _sub_repl, plain)
    return f"{formatted}{charge_suffix}"


def format_chemical_equation_display(equation_text: str) -> str:
    """Formats a balanced chemical equation with Unicode subscripts and right arrow (→) for LTR display."""
    if not equation_text:
        return ""
    normalized_eq = equation_text.replace("->", "→")
    sides = normalized_eq.split("→")
    formatted_sides = []
    for side in sides:
        terms = [t.strip() for t in side.split("+") if t.strip()]
        out_terms = []
        for term in terms:
            m = re.match(r"^(\d+)(.+)$", term)
            if m:
                coeff, form = m.group(1), m.group(2).strip()
                out_terms.append(f"{coeff}{format_chemical_formula_display(form)}")
            else:
                out_terms.append(format_chemical_formula_display(term))
        formatted_sides.append(" + ".join(out_terms))
    return " → ".join(formatted_sides)


def get_standard_phase_code(formula_str: str, physical_state_text: str = "") -> str:
    """Returns the standard state phase code ('g', 'l', 's', 'aq') for a substance at STP/SATP."""
    cleaned = replace_subscripts((formula_str or "").strip()).replace(" ", "")
    if cleaned in GASEOUS_FORMULAS_STP:
        return "g"
    if cleaned in LIQUID_FORMULAS_STP:
        return "l"
    txt = (physical_state_text or "").strip()
    if txt.startswith("غاز") or "غاز عديم" in txt or "غاز سام" in txt or "غاز أصفر" in txt or "غاز بني" in txt or "غاز خانق" in txt:
        return "g"
    if txt.startswith("سائل"):
        return "l"
    if "محلول" in txt or "أيون" in txt:
        return "aq"
    return "s"


def is_gaseous_at_stp(formula_str: str, substance: dict | None = None, rxn_state: str = "") -> bool:
    """Checks whether a substance is a gas eligible for STP molar volume conversion (V / 22.414 L/mol)."""
    if rxn_state:
        if rxn_state == "g":
            return True
        if rxn_state in {"s", "l", "aq"}:
            return False
    cleaned = replace_subscripts((formula_str or "").strip()).replace(" ", "")
    if cleaned in GASEOUS_FORMULAS_STP:
        return True
    if cleaned in LIQUID_FORMULAS_STP:
        return False
    if substance and isinstance(substance, dict):
        phys = substance.get("physical_properties", {}) or {}
        p_state = phys.get("physical_state", "")
        phase = get_standard_phase_code(cleaned, p_state)
        return phase == "g"
    return False


def clean_formula_token(formula_str: str) -> str:
    if not formula_str:
        return ""
    s = replace_subscripts(formula_str.strip()).replace(" ", "")
    s = re.sub(r"\^[0-9]*[+\-]?", "", s)
    s = re.sub(r"\([0-9]*[+\-]\)", "", s)
    s = re.sub(r"(^[A-Z][a-z]?)\d*[+\-]\d*$", r"\1", s)
    s = re.sub(r"[+\-]\d*$", "", s)
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
            "formatted_formula": format_chemical_formula_display(formula_str),
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
        "formatted_formula": format_chemical_formula_display(formula_str),
        "total_atoms": total_atoms,
        "elements_count": len(counts),
        "atom_counts": counts,
        "components": components,
        "summary_ar": summary_ar,
    }

"""
Reaction Prediction Engine — محرك توقع التفاعلات (المختبر الافتراضي).

المبدأ العلمي الحاكم في هذا المحرك:
    1. لا يُختلق أي ناتج كيميائي. كل ناتج مقترح إمّا موثّق في سجل التفاعلات المعتمد،
       أو مُشتق من قاعدة كيميائية عامة معلنة صراحةً (مع ذكر القاعدة ومصدرها وحدودها).
    2. عند غياب دليل كافٍ تُعرض رسالة "لا تتوفر بيانات كافية" بدل نتيجة مضللة.
    3. كل ناتج مستنتج يحمل دليله (Evidence) ومستوى ثقته، ويُفصل بوضوح عن التفاعلات الموثقة مخبرياً.

مستويات الدليل (Evidence Tiers):
    - "documented"         : تفاعل موثّق في سجل التفاعلات المعتمد (مرجع علمي محدّد).
    - "rule_inferred"      : استنتاج من قاعدة كيميائية عامة (General Chemistry Rule) مع ذكر المصدر.
    - "insufficient_evidence": لا يوجد أساس كيميائي موثوق كافٍ → لا تُعرض نواتج.
"""

from __future__ import annotations

import re
from typing import Any

from chemistry_engine.formula_analyzer import (
    format_chemical_equation_display,
    get_standard_phase_code,
    parse_chemical_formula,
)
from chemistry_engine.periodic_table import ELEMENT_NAMES_AR, ELEMENT_NAMES_EN
from reaction_engine.equation_balancer import (
    balance_reaction_algebraically,
    format_balanced_equation,
    verify_atom_conservation,
)
from reaction_engine.stoichiometry import compute_reaction_stoichiometry
from reaction_engine.yield_calculator import compute_theoretical_yields

# --------------------------------------------------------------------------------------
# 1) Knowledge base — قواعد كيميائية عامة
# --------------------------------------------------------------------------------------

NONMETAL_ELEMENTS: set[str] = {
    "H", "C", "N", "O", "F", "P", "S", "Cl", "Se", "Br", "I", "Si", "B", "As", "Te",
}

ALKALI_METALS: set[str] = {"Li", "Na", "K", "Rb", "Cs"}
ALKALINE_EARTH_METALS: set[str] = {"Be", "Mg", "Ca", "Sr", "Ba"}

# سلسلة النشاط الفلزي (الأقوى نشاطاً أولاً) — أساس قواعد الإحلال المفرد
ACTIVITY_SERIES: list[str] = [
    "K", "Ba", "Ca", "Na", "Mg", "Al", "Mn", "Zn", "Cr", "Fe", "Cd", "Co", "Ni",
    "Sn", "Pb", "H", "Cu", "Hg", "Ag", "Au",
]

# معادن لا تتفاعل مع الأحماض المخففة العادية وتُعتبر خاملة تجاه الماء
NOBLE_METALS: set[str] = {"Cu", "Hg", "Ag", "Au", "Pt", "Pd"}

# أكاسيد المعادن الشائعة والأوكسدة المعتمدة لها (لصياغة نواتج التفاعل مع الأكسجين)
METAL_OXIDES: dict[str, str] = {
    "Li": "Li2O", "Na": "Na2O", "K": "K2O", "Mg": "MgO", "Ca": "CaO", "Al": "Al2O3",
    "Zn": "ZnO", "Fe": "Fe2O3", "Cu": "CuO", "Pb": "PbO", "Ag": "Ag2O", "Sn": "SnO2",
    "Mn": "MnO2", "Co": "Co3O4", "Ni": "NiO", "Cr": "Cr2O3", "Ba": "BaO", "Sr": "SrO",
}

# حالة الأكسدة الشائعة للفلز في أملاحه (تُستخدم لتوليد صيغة الملح الصحيح)
METAL_ION_CHARGES: dict[str, int] = {
    "Li": 1, "Na": 1, "K": 1, "Rb": 1, "Cs": 1, "Ag": 1,
    "Mg": 2, "Ca": 2, "Sr": 2, "Ba": 2, "Zn": 2, "Fe": 2, "Cu": 2, "Ni": 2,
    "Co": 2, "Mn": 2, "Sn": 2, "Pb": 2, "Cd": 2, "Ti": 4,
    "Al": 3, "Cr": 3,
}

# أيونات/شجرات معقدة معرّفة بمكوّناتها الذرية وشحنتها
POLYATOMIC_GROUPS: list[dict[str, Any]] = [
    {"key": "NH4", "counts": {"N": 1, "H": 4}, "charge": 1, "is_cation": True, "ar": "الأمونيوم", "en": "Ammonium"},
    {"key": "OH", "counts": {"O": 1, "H": 1}, "charge": -1, "is_cation": False, "ar": "الهيدروكسيد", "en": "Hydroxide"},
    {"key": "NO3", "counts": {"N": 1, "O": 3}, "charge": -1, "is_cation": False, "ar": "النترات", "en": "Nitrate"},
    {"key": "HCO3", "counts": {"H": 1, "C": 1, "O": 3}, "charge": -1, "is_cation": False, "ar": "البيكربونات", "en": "Bicarbonate"},
    {"key": "CO3", "counts": {"C": 1, "O": 3}, "charge": -2, "is_cation": False, "ar": "الكربونات", "en": "Carbonate"},
    {"key": "SO4", "counts": {"S": 1, "O": 4}, "charge": -2, "is_cation": False, "ar": "الكبريتات", "en": "Sulfate"},
    {"key": "PO4", "counts": {"P": 1, "O": 4}, "charge": -3, "is_cation": False, "ar": "الفوسفات", "en": "Phosphate"},
    {"key": "CrO4", "counts": {"Cr": 1, "O": 4}, "charge": -2, "is_cation": False, "ar": "الكرومات", "en": "Chromate"},
    {"key": "Cr2O7", "counts": {"Cr": 2, "O": 7}, "charge": -2, "is_cation": False, "ar": "الدي كرومات", "en": "Dichromate"},
    {"key": "CH3COO", "counts": {"C": 2, "H": 3, "O": 2}, "charge": -1, "is_cation": False, "ar": "الأسيتات", "en": "Acetate"},
    {"key": "Cl", "counts": {"Cl": 1}, "charge": -1, "is_cation": False, "ar": "الكلوريد", "en": "Chloride"},
    {"key": "Br", "counts": {"Br": 1}, "charge": -1, "is_cation": False, "ar": "البروميد", "en": "Bromide"},
    {"key": "I", "counts": {"I": 1}, "charge": -1, "is_cation": False, "ar": "اليوديد", "en": "Iodide"},
    {"key": "F", "counts": {"F": 1}, "charge": -1, "is_cation": False, "ar": "الفلوريد", "en": "Fluoride"},
    {"key": "S", "counts": {"S": 1}, "charge": -2, "is_cation": False, "ar": "الكبريتيد", "en": "Sulfide"},
    {"key": "O", "counts": {"O": 1}, "charge": -2, "is_cation": False, "ar": "الأكسيد", "en": "Oxide"},
]

ANION_LABELS_AR: dict[str, str] = {
    "Cl": "كلوريد", "Br": "بروميد", "I": "يوديد", "F": "فلوريد", "SO4": "كبريتات",
    "NO3": "نترات", "CO3": "كربونات", "HCO3": "بيكربونات", "PO4": "فوسفات",
    "OH": "هيدروكسيد", "CH3COO": "أسيتات", "CrO4": "كرومات", "Cr2O7": "دي كرومات",
    "S": "كبريتيد", "O": "أكسيد", "HCOO": "ميثانوات",
}
ANION_LABELS_EN: dict[str, str] = {
    "Cl": "chloride", "Br": "bromide", "I": "iodide", "F": "fluoride", "SO4": "sulfate",
    "NO3": "nitrate", "CO3": "carbonate", "HCO3": "bicarbonate", "PO4": "phosphate",
    "OH": "hydroxide", "CH3COO": "acetate", "CrO4": "chromate", "Cr2O7": "dichromate",
    "S": "sulfide", "O": "oxide", "HCOO": "methanoate",
}

# شجرات أنيونات الأحماض القوية → الحمض قوي
STRONG_ACID_ANIONS: set[str] = {"Cl", "Br", "I", "NO3"}
STRONG_ACID_FORMULAS: set[str] = {"HClO4", "HBrO4"}

# صيغ حمض ضعيف معروفة (غير قابلة للاشتقاق من شجرات معقدة)
WEAK_ACID_FORMULAS: dict[str, tuple[str, int]] = {
    "HF": ("F", 1), "HCN": ("CN", 1), "H2S": ("S", 1), "HSCN": ("SCN", 1),
    "H2CO3": ("CO3", 2), "H3PO4": ("PO4", 3), "H2SO3": ("SO3", 2), "H2C2O4": ("C2O4", 2),
    "HNO2": ("NO2", 1), "H3BO3": ("BO3", 3), "H3AsO4": ("AsO4", 3), "HClO": ("ClO", 1),
}

# قواعد الذوبانية (Solubility rules) — قاعدة عامة تربوياً
INSOLUBLE_CATIONS_WITH_CHLORIDE: set[str] = {"Ag", "Pb", "Hg"}
INSOLUBLE_CATIONS_WITH_SULFATE: set[str] = {"Ba", "Pb", "Sr", "Ca", "Ag"}
ALWAYS_SOLUBLE_CATIONS: set[str] = {"Na", "K", "NH4", "Li"}
ALWAYS_SOLUBLE_ANIONS: set[str] = {"NO3", "ClO4", "CH3COO"}

# ألوان الرواسب الشائعة (للملاحظات التجريبية المتوقعة)
PRECIPITATE_COLORS_AR: dict[str, str] = {
    "AgCl": "أبيض", "AgBr": "أبيض مائل للصفرة", "AgI": "أصفر", "Ag2CO3": "أصفر",
    "Ag3PO4": "أصفر", "PbI2": "أصفر ذهبي", "PbCrO4": "أصفر", "PbSO4": "أبيض",
    "PbCl2": "أبيض", "BaSO4": "أبيض", "SrSO4": "أبيض", "CaSO4": "أبيض",
    "CaCO3": "أبيض", "MgCO3": "أبيض", "Cu(OH)2": "أزرق سماوي", "Fe(OH)3": "بني محمر",
    "Fe(OH)2": "أخضر شاحب", "Mg(OH)2": "أبيض", "Zn(OH)2": "أبيض", "Al(OH)3": "أبيض هلامي",
    "CuS": "أسود", "FeS": "أسود", "PbS": "أسود", "Ag2S": "أسود", "MnO2": "أسود",
    "BaCO3": "أبيض", "Ca3(PO4)2": "أبيض", "CuCO3": "أخضر مزرق",
}

COMMON_COMPOUND_NAMES_AR: dict[str, str] = {
    "H2O": "الماء", "H2O2": "بيروكسيد الهيدروجين", "CO2": "ثاني أكسيد الكربون",
    "CO": "أول أكسيد الكربون", "H2": "الهيدروجين", "O2": "الأكسجين", "N2": "النيتروجين",
    "NH3": "الأمونيا", "SO2": "ثاني أكسيد الكبريت", "H2S": "كبريتيد الهيدروجين",
    "NO": "أكسيد النيتروجين", "NO2": "ثاني أكسيد النيتروجين", "CH4": "الميثان",
    "NaOH": "هيدروكسيد الصوديوم", "KOH": "هيدروكسيد البوتاسيوم",
    "Ca(OH)2": "هيدروكسيد الكالسيوم", "Ba(OH)2": "هيدروكسيد الباريوم",
    "LiOH": "هيدروكسيد الليثيوم", "HCl": "حمض الهيدروكلوريك",
    "HBr": "حمض الهيدروبروميك", "HI": "حمض الهيدرويويدي", "HF": "حمض الهيدروفلوريك",
    "H2SO4": "حمض الكبريتيك", "HNO3": "حمض النيتريك", "HClO4": "حمض البيركلوريك",
    "H3PO4": "حمض الفوسفوريك", "H2CO3": "حمض الكربونيك", "CH3COOH": "حمض الخليك",
    "HCOOH": "حمض النمل", "H2S": "كبريتيد الهيدروجين", "HCN": "حمض السيانيد",
    "C2H5OH": "الإيثانول", "CH3OH": "الميثانول", "C6H6": "البنزين",
}
COMMON_COMPOUND_NAMES_EN: dict[str, str] = {
    "H2O": "Water", "H2O2": "Hydrogen Peroxide", "CO2": "Carbon Dioxide",
    "CO": "Carbon Monoxide", "H2": "Hydrogen", "O2": "Oxygen", "N2": "Nitrogen",
    "NH3": "Ammonia", "SO2": "Sulfur Dioxide", "H2S": "Hydrogen Sulfide",
    "NO": "Nitric Oxide", "NO2": "Nitrogen Dioxide", "CH4": "Methane",
    "NaOH": "Sodium Hydroxide", "KOH": "Potassium Hydroxide",
    "Ca(OH)2": "Calcium Hydroxide", "Ba(OH)2": "Barium Hydroxide",
    "LiOH": "Lithium Hydroxide", "HCl": "Hydrochloric Acid",
    "HBr": "Hydrobromic Acid", "HI": "Hydroiodic Acid", "HF": "Hydrofluoric Acid",
    "H2SO4": "Sulfuric Acid", "HNO3": "Nitric Acid", "HClO4": "Perchloric Acid",
    "H3PO4": "Phosphoric Acid", "H2CO3": "Carbonic Acid", "CH3COOH": "Acetic Acid",
    "HCOOH": "Formic Acid", "HCN": "Hydrocyanic Acid",
    "C2H5OH": "Ethanol", "CH3OH": "Methanol", "C6H6": "Benzene",
}

# أكاسيد غير فلزية (حموضة) شائعة
ACIDIC_OXIDE_FORMULAS: set[str] = {"SO2", "SO3", "CO2", "CO", "P2O5", "P4O10", "N2O5", "SiO2"}

# مصادر القواعد العامة المستخدمة في الاستنتاج
SOURCE_VOGEL = (
    "Vogel's Textbook of Quantitative Chemical Analysis, 6th Ed., Longman (1999) — "
    "General Inorganic Reaction Rules (Acid–Base, Precipitation, Displacement)"
)
SOURCE_CRC = "CRC Handbook of Chemistry and Physics, 104th Ed. (2023) — Solubility Rules and Reaction Classes"
SOURCE_BROWN = (
    "Brown, T.L., LeMay, H.E. et al., Chemistry: The Central Science, 14th Ed., "
    "Pearson (2018) — Ch. 4 (Aqueous Reactions) & Ch. 5 (Thermodynamics)"
)
SOURCE_CAREY = (
    "Carey, G.A. & Sundberg, R.J., Advanced Organic Chemistry, Part A: Structure and Mechanisms, "
    "Springer (2007) — Nucleophilic Acyl Substitution (Esterification)"
)

# --------------------------------------------------------------------------------------
# 2) Species classification
# --------------------------------------------------------------------------------------

KIND_METAL = "metal"
KIND_NONMETAL_ELEMENT = "nonmetal_element"
KIND_ELEMENTAL_MOLECULE = "elemental_molecule"
KIND_STRONG_ACID = "strong_acid"
KIND_WEAK_ACID = "weak_acid"
KIND_STRONG_BASE = "strong_base"
KIND_WEAK_BASE = "weak_base"
KIND_SALT = "salt"
KIND_OXIDE = "oxide"
KIND_HYDROCARBON = "hydrocarbon"
KIND_ORGANIC_ALCOHOL = "organic_alcohol"
KIND_ORGANIC_OTHER = "organic_compound"
KIND_WATER = "water"
KIND_AMMONIA = "ammonia"
KIND_PEROXIDE = "peroxide"
KIND_UNKNOWN = "unknown"

KIND_LABELS_AR: dict[str, str] = {
    KIND_METAL: "عنصر فلزي",
    KIND_NONMETAL_ELEMENT: "عنصر لافلزي",
    KIND_ELEMENTAL_MOLECULE: "جزيء عنصري",
    KIND_STRONG_ACID: "حمض قوي",
    KIND_WEAK_ACID: "حمض ضعيف",
    KIND_STRONG_BASE: "قاعدة قوية (هيدروكسيد ذائب)",
    KIND_WEAK_BASE: "قاعدة ضعيفة",
    KIND_SALT: "ملح",
    KIND_OXIDE: "أوكسيد",
    KIND_HYDROCARBON: "هيدروكربون (مركب عضوي)",
    KIND_ORGANIC_ALCOHOL: "كحول أحادي الوظيفة أو إيثر",
    KIND_ORGANIC_OTHER: "مركب عضوي",
    KIND_WATER: "ماء",
    KIND_AMMONIA: "أمونيا (قاعدة ضعيفة)",
    KIND_PEROXIDE: "بيروكسيد",
    KIND_UNKNOWN: "مادة غير مصنّفة",
}


def _strip_hydrate(formula: str) -> tuple[str, str]:
    """Removes the hydrate part (e.g. 'CuSO4·5H2O' -> 'CuSO4') and returns a human-readable note."""
    if not formula:
        return "", ""
    cleaned = formula.replace("·", ".").replace("*", ".")
    if "." not in cleaned:
        return cleaned, ""
    base, _, tail = cleaned.partition(".")
    tail_clean = tail.replace("H2O", "").strip()
    note = f" (مادة بلورية مائية: {tail_clean or tail} H₂O)" if tail_clean else ""
    return base, note


def _subtract(base: dict[str, int], part: dict[str, int], times: int) -> dict[str, int] | None:
    out = dict(base)
    for el, cnt in part.items():
        if out.get(el, 0) < cnt * times:
            return None
        out[el] = out.get(el, 0) - cnt * times
        if out[el] == 0:
            del out[el]
    return out


def _fit_group(counts: dict[str, int], group: dict[str, Any]) -> tuple[int, dict[str, int]] | None:
    """Finds the maximum multiplicity of a polyatomic group that fits inside the formula."""
    if not group["counts"]:
        return None
    limit = min(counts.get(el, 0) // cnt for el, cnt in group["counts"].items())
    if limit < 1:
        return None
    remainder = _subtract(counts, group["counts"], limit)
    if remainder is None:
        return None
    return limit, remainder


def _enumerate_ion_candidates(counts: dict[str, int]) -> list[dict[str, Any]]:
    """
    Enumerates every chemically sensible (cation, anion) interpretation of a formula.

    A candidate is accepted only when the charges balance exactly, which avoids the
    false positives produced by naive greedy pattern matching (e.g. NaHCO3 mis-read
    as hydroxide + carbonate, or H2SO4 mis-read as sulfide + oxide).
    """
    candidates: list[dict[str, Any]] = []

    # --- Case 1: the cation is the ammonium ion NH4+ -------------------------------------
    nh4 = next(g for g in POLYATOMIC_GROUPS if g["key"] == "NH4")
    nh4_fit = _fit_group(counts, nh4)
    if nh4_fit:
        nh4_times, remainder = nh4_fit
        for group in POLYATOMIC_GROUPS:
            if group["is_cation"]:
                continue
            anion_fit = _fit_group(remainder, group)
            if not anion_fit:
                continue
            anion_times, rest = anion_fit
            if rest:
                continue
            candidates.append(
                {
                    "cation": "NH4",
                    "cation_charge": 1,
                    "cation_count": nh4_times,
                    "anion": group["key"],
                    "anion_charge": group["charge"],
                    "anion_count": anion_times,
                    "specificity": len(group["counts"]) * 10 + nh4_times,
                }
            )

    # --- Case 2: the cation is a metal (salt) or a proton (acid) -------------------------
    for group in POLYATOMIC_GROUPS:
        if group["is_cation"]:
            continue
        anion_fit = _fit_group(counts, group)
        if not anion_fit:
            continue
        anion_times, remainder = anion_fit

        leftover_nonmetals = {
            el: cnt for el, cnt in remainder.items() if el != "H" and el in NONMETAL_ELEMENTS
        }
        if leftover_nonmetals:
            continue

        hydrogen_count = remainder.get("H", 0)
        metals = {el: cnt for el, cnt in remainder.items() if el != "H"}

        if not metals and hydrogen_count > 0:
            # acid: H+ (hydrogen_count) paired with the anion group
            candidates.append(
                {
                    "cation": "H",
                    "cation_charge": 1,
                    "cation_count": hydrogen_count,
                    "anion": group["key"],
                    "anion_charge": group["charge"],
                    "anion_count": anion_times,
                    "specificity": len(group["counts"]) * 10 + hydrogen_count,
                }
            )
        elif len(metals) == 1 and hydrogen_count == 0:
            symbol, count = next(iter(metals.items()))
            total_anion_charge = abs(group["charge"] * anion_times)
            if count <= 0 or total_anion_charge % count != 0:
                continue
            charge = total_anion_charge // count
            if not 1 <= charge <= 8:
                continue
            candidates.append(
                {
                    "cation": symbol,
                    "cation_charge": charge,
                    "cation_count": count,
                    "anion": group["key"],
                    "anion_charge": group["charge"],
                    "anion_count": anion_times,
                    "specificity": len(group["counts"]) * 10 + anion_times,
                }
            )

    candidates.sort(key=lambda c: -c["specificity"])
    return candidates


def decompose_ions(counts: dict[str, int]) -> dict[str, Any] | None:
    """
    Decomposes an inorganic formula into cation / anion components.

    Returns None when the formula cannot be interpreted as a simple ionic compound.
    """
    candidates = _enumerate_ion_candidates(counts)
    if not candidates:
        return None
    best = candidates[0]
    return {
        "cation": best["cation"],
        "cation_charge": best["cation_charge"],
        "cation_count": best["cation_count"],
        "anion": best["anion"],
        "anion_charge": best["anion_charge"],
        "anion_count": best["anion_count"],
        "candidates_count": len(candidates),
    }


def _salt_formula(
    cation: str,
    cation_charge: int,
    cation_count: int,
    anion: str,
    anion_charge: int,
    anion_count: int = 1,
) -> str:
    """
    Builds the formula of a salt from its ions, balancing the charges.

    Examples: Fe3+ + SO4^2- -> Fe2(SO4)3 ; Mg2+ + OH- -> Mg(OH)2 ; NH4+ + SO4^2- -> (NH4)2SO4.
    """
    if not cation or cation == "H":
        return ""
    z = abs(int(cation_charge))
    q = abs(int(anion_charge))
    if z == 0 or q == 0:
        return ""
    cc = max(int(cation_count), 1)
    ac = max(int(anion_count), 1)
    guard = 0
    while cc * z != ac * q and guard < 24:
        if cc * z < ac * q:
            cc += 1
        else:
            ac += 1
        guard += 1
    if cc > 1 and ac > 1:
        divisor = _gcd(cc, ac)
        if divisor > 1:
            cc //= divisor
            ac //= divisor

    if cation == "NH4":
        cation_part = "NH4"
        if cc > 1:
            cation_part = f"({cation_part}){cc}"
    else:
        cation_part = cation + (str(cc) if cc > 1 else "")

    anion_part = anion
    if ac > 1:
        # only polyatomic groups are wrapped in parentheses (Mg(OH)2, not Mg(OH2))
        if anion in {"Cl", "Br", "I", "F", "S", "O"}:
            anion_part = f"{anion_part}{ac}"
        else:
            anion_part = f"({anion_part}){ac}"
    return (cation_part + anion_part).strip()


def _gcd(a: int, b: int) -> int:
    while b:
        a, b = b, a % b
    return abs(a) or 1


def _with_al(name: str) -> str:
    """Adds the Arabic definite article to element names used inside salt names."""
    if not name:
        return name
    if name.startswith("ال"):
        return name
    return f"ال{name}"


def _compound_name_ar(formula: str, kind: str, info: dict[str, Any]) -> str:
    if formula in COMMON_COMPOUND_NAMES_AR:
        return COMMON_COMPOUND_NAMES_AR[formula]
    cation = info.get("cation")
    anion = info.get("anion")
    if anion and cation and cation != "H":
        anion_ar = ANION_LABELS_AR.get(anion, anion)
        if cation == "NH4":
            return f"{anion_ar} الأمونيوم"
        return f"{anion_ar} {_with_al(ELEMENT_NAMES_AR.get(cation, cation))}"
    if kind == KIND_OXIDE and cation:
        return f"أكسيد {_with_al(ELEMENT_NAMES_AR.get(cation, cation))}"
    return formula


def _compound_name_en(formula: str, kind: str, info: dict[str, Any]) -> str:
    if formula in COMMON_COMPOUND_NAMES_EN:
        return COMMON_COMPOUND_NAMES_EN[formula]
    cation = info.get("cation")
    anion = info.get("anion")
    if anion and cation and cation != "H":
        anion_en = ANION_LABELS_EN.get(anion, anion)
        if cation == "NH4":
            return f"Ammonium {anion_en}"
        return f"{anion_en.capitalize()} of {ELEMENT_NAMES_EN.get(cation, cation)}"
    if kind == KIND_OXIDE and cation:
        return f"{ELEMENT_NAMES_EN.get(cation, cation)} Oxide"
    return formula


def classify_species(formula: str, name_ar: str = "", name_en: str = "") -> dict[str, Any]:
    """Classifies a substance into a chemically meaningful family and extracts its ions."""
    raw_formula = (formula or "").strip()
    base_formula, hydrate_note = _strip_hydrate(raw_formula)
    counts = parse_chemical_formula(base_formula) if base_formula else None

    info: dict[str, Any] = {
        "formula": base_formula or raw_formula,
        "input_formula": raw_formula,
        "name_ar": name_ar or base_formula or raw_formula,
        "name_en": name_en or base_formula or raw_formula,
        "kind": KIND_UNKNOWN,
        "kind_label_ar": KIND_LABELS_AR[KIND_UNKNOWN],
        "counts": counts or {},
        "cation": "",
        "cation_charge": 0,
        "cation_count": 0,
        "anion": "",
        "anion_charge": 0,
        "anion_count": 0,
        "proticity": 0,
        "hydrate_note": hydrate_note,
        "is_valid_formula": counts is not None,
    }

    if not counts:
        info["kind_label_ar"] = "صيغة غير قابلة للتحليل الذري"
        return info

    elements = set(counts)

    # --- water / ammonia / peroxide ---
    if counts == {"H": 2, "O": 1}:
        info["kind"] = KIND_WATER
        info["name_ar"] = name_ar or "الماء"
        info["name_en"] = name_en or "Water"
        info["kind_label_ar"] = KIND_LABELS_AR[KIND_WATER]
        return info
    if counts == {"N": 1, "H": 3}:
        info["kind"] = KIND_AMMONIA
        info["name_ar"] = name_ar or "الأمونيا"
        info["name_en"] = name_en or "Ammonia"
        info["kind_label_ar"] = KIND_LABELS_AR[KIND_AMMONIA]
        return info
    if counts == {"H": 2, "O": 2}:
        info["kind"] = KIND_PEROXIDE
        info["name_ar"] = name_ar or "بيروكسيد الهيدروجين"
        info["name_en"] = name_en or "Hydrogen Peroxide"
        info["kind_label_ar"] = KIND_LABELS_AR[KIND_PEROXIDE]
        return info

    # --- single element ---
    if len(elements) == 1:
        el = next(iter(elements))
        cnt = counts[el]
        if cnt == 1:
            if el in NONMETAL_ELEMENTS:
                info["kind"] = KIND_NONMETAL_ELEMENT
            else:
                info["kind"] = KIND_METAL
                info["cation"] = el
                info["cation_charge"] = 0
                info["cation_count"] = 1
        else:
            info["kind"] = KIND_ELEMENTAL_MOLECULE
        info["name_ar"] = name_ar or ELEMENT_NAMES_AR.get(el, el)
        info["name_en"] = name_en or ELEMENT_NAMES_EN.get(el, el)
        info["kind_label_ar"] = KIND_LABELS_AR[info["kind"]]
        return info

    # --- hydrocarbons (C and H only) ---
    if elements <= {"C", "H"} and "C" in elements:
        info["kind"] = KIND_HYDROCARBON
        info["name_ar"] = name_ar or COMMON_COMPOUND_NAMES_AR.get(base_formula, base_formula)
        info["name_en"] = name_en or COMMON_COMPOUND_NAMES_EN.get(base_formula, base_formula)
        info["kind_label_ar"] = KIND_LABELS_AR[KIND_HYDROCARBON]
        return info

    # --- organic families (formula level) ---
    if elements <= {"C", "H", "O"}:
        o_count = counts.get("O", 0)
        c_count = counts.get("C", 0)
        h_count = counts.get("H", 0)
        if o_count == 2 and c_count >= 1 and h_count >= 2:
            info["kind"] = KIND_WEAK_ACID
            info["anion"] = "CH3COO" if c_count >= 2 else "HCOO"
            info["anion_charge"] = -1
            info["anion_count"] = 1
            info["cation"] = "H"
            info["cation_charge"] = 1
            info["cation_count"] = 1
            info["proticity"] = 1
        elif o_count == 1 and c_count >= 1 and h_count >= 2:
            info["kind"] = KIND_ORGANIC_ALCOHOL
        else:
            info["kind"] = KIND_ORGANIC_OTHER
        info["name_ar"] = name_ar or COMMON_COMPOUND_NAMES_AR.get(base_formula, base_formula)
        info["name_en"] = name_en or COMMON_COMPOUND_NAMES_EN.get(base_formula, base_formula)
        info["kind_label_ar"] = KIND_LABELS_AR[info["kind"]]
        return info

    # --- explicit weak-acid formulas that do not decompose into standard groups ---
    if base_formula in WEAK_ACID_FORMULAS:
        anion_key, prot = WEAK_ACID_FORMULAS[base_formula]
        info["kind"] = KIND_WEAK_ACID
        info["anion"] = anion_key
        info["anion_charge"] = -1
        info["anion_count"] = 1
        info["cation"] = "H"
        info["cation_charge"] = 1
        info["cation_count"] = prot
        info["proticity"] = prot
        info["name_ar"] = name_ar or COMMON_COMPOUND_NAMES_AR.get(base_formula, base_formula)
        info["name_en"] = name_en or COMMON_COMPOUND_NAMES_EN.get(base_formula, base_formula)
        info["kind_label_ar"] = KIND_LABELS_AR[KIND_WEAK_ACID]
        return info

    # --- acidic (non-metal) oxides ---
    if base_formula in ACIDIC_OXIDE_FORMULAS:
        info["kind"] = KIND_OXIDE
        info["is_acidic_oxide"] = True
        info["name_ar"] = name_ar or COMMON_COMPOUND_NAMES_AR.get(base_formula, base_formula)
        info["name_en"] = name_en or COMMON_COMPOUND_NAMES_EN.get(base_formula, base_formula)
        info["kind_label_ar"] = "أوكسيد حمضي (Acid-forming Oxide)"
        return info

    # --- ionic compounds ---
    decomp = decompose_ions(counts)
    if decomp is None:
        info["kind"] = KIND_UNKNOWN
        info["kind_label_ar"] = KIND_LABELS_AR[KIND_UNKNOWN]
        return info

    info["cation"] = decomp["cation"]
    info["cation_charge"] = decomp["cation_charge"]
    info["cation_count"] = decomp["cation_count"]
    info["anion"] = decomp["anion"]
    info["anion_charge"] = decomp["anion_charge"]
    info["anion_count"] = decomp["anion_count"]

    cation = decomp["cation"]
    anion = decomp["anion"]

    if cation == "H":
        # acid
        prot = decomp["cation_count"]
        info["proticity"] = prot
        strong = anion in STRONG_ACID_ANIONS or (anion == "SO4" and prot <= 2) or base_formula in STRONG_ACID_FORMULAS
        info["kind"] = KIND_STRONG_ACID if strong else KIND_WEAK_ACID
        info["name_ar"] = name_ar or _compound_name_ar(base_formula, info["kind"], info)
        info["name_en"] = name_en or _compound_name_en(base_formula, info["kind"], info)
    elif anion == "OH":
        # hydroxide → base or (insoluble) hydroxide salt
        if cation in ALKALI_METALS or cation in {"Ca", "Sr", "Ba"}:
            info["kind"] = KIND_STRONG_BASE
        elif cation == "NH4":
            info["kind"] = KIND_WEAK_BASE
        else:
            info["kind"] = KIND_SALT
        info["name_ar"] = name_ar or _compound_name_ar(base_formula, info["kind"], info)
        info["name_en"] = name_en or _compound_name_en(base_formula, info["kind"], info)
    elif anion == "O":
        info["kind"] = KIND_OXIDE
        info["is_acidic_oxide"] = cation in {"", "C", "S", "P", "N", "Si"}
        info["name_ar"] = name_ar or _compound_name_ar(base_formula, info["kind"], info)
        info["name_en"] = name_en or _compound_name_en(base_formula, info["kind"], info)
    else:
        info["kind"] = KIND_SALT
        info["name_ar"] = name_ar or _compound_name_ar(base_formula, info["kind"], info)
        info["name_en"] = name_en or _compound_name_en(base_formula, info["kind"], info)

    # weak-acid special cases that do not decompose into standard groups
    info["kind_label_ar"] = KIND_LABELS_AR.get(info["kind"], KIND_LABELS_AR[KIND_UNKNOWN])
    return info


# --------------------------------------------------------------------------------------
# 3) Solubility & activity helpers
# --------------------------------------------------------------------------------------


def is_salt_insoluble(cation: str, anion: str) -> bool:
    """Applies the classic solubility rules to decide whether an ionic product precipitates."""
    if not cation or not anion:
        return False
    if cation in ALWAYS_SOLUBLE_CATIONS or cation == "H":
        return False
    if anion in ALWAYS_SOLUBLE_ANIONS:
        return False
    if anion == "Cl" and cation in INSOLUBLE_CATIONS_WITH_CHLORIDE:
        return True
    if anion == "SO4" and cation in INSOLUBLE_CATIONS_WITH_SULFATE:
        return True
    if anion in {"CO3", "PO4", "S", "CrO4"}:
        return True
    if anion == "OH":
        if cation in ALKALI_METALS or cation == "NH4":
            return False
        if cation in {"Ca", "Sr", "Ba"}:
            return False  # slightly soluble hydroxides
        return True
    if anion in {"Br", "I"} and cation in {"Ag", "Pb", "Hg"}:
        return True
    return False


def metal_activity_index(symbol: str) -> int:
    """Position of a metal in the activity series (-1 when not listed)."""
    try:
        return ACTIVITY_SERIES.index(symbol)
    except ValueError:
        return -1


def is_more_active_metal(a: str, b: str) -> bool:
    ia, ib = metal_activity_index(a), metal_activity_index(b)
    if ia < 0 or ib < 0:
        return False
    return ia < ib


# --------------------------------------------------------------------------------------
# 4) Rule-based inference engine
# --------------------------------------------------------------------------------------


def _rule(
    rule_id: str,
    name_ar: str,
    basis_ar: str,
    source: str,
    source_id: str,
    url: str,
    category: str,
    reaction_type_ar: str,
    confidence: str,
) -> dict[str, Any]:
    return {
        "rule_id": rule_id,
        "rule_name_ar": name_ar,
        "basis_ar": basis_ar,
        "source": source,
        "source_id": source_id,
        "reference_url": url,
        "category": category,
        "reaction_type_ar": reaction_type_ar,
        "confidence": confidence,
    }


RULES: dict[str, dict[str, Any]] = {
    "acid_base_neutralization": _rule(
        "acid_base_neutralization",
        "قاعدة: حمض + قاعدة (هيدروكسيد) = ملح + ماء",
        "يتحد أيون الهيدروجين (H⁺) من الحمض مع أيون الهيدروكسيد (OH⁻) من القاعدة بنسبة 1:1 لتكوين جزيء الماء، "
        "بينما يبقى الأيون المعدني مع الأنيون الحمضي في الملح الناتج (المعادلة الصافية: H⁺ + OH⁻ → H₂O).",
        SOURCE_VOGEL,
        "VOGEL-6ed-Ch10-AcidBase-Neutralisation",
        "https://books.google.com/books?id=Jh4gAQAAIAAJ",
        "neutralization",
        "تفاعل تعادل حمض وقاعدة (Acid–Base Neutralization)",
        "high",
    ),
    "acid_carbonate": _rule(
        "acid_carbonate",
        "قاعدة: حمض + كربونات/بيكربونات = ملح + ماء + غاز CO₂",
        "أيونات الكربونات والبيكربونات ضعيفة القاعدية وتتفاعل مع الأحماض فتحرّر ثاني أكسيد الكربون غازاً، "
        "مع ملح المعدن والماء (المعادلة الصافية: HCO₃⁻ + H⁺ → H₂O + CO₂↑).",
        SOURCE_VOGEL,
        "VOGEL-6ed-Ch10-Carbonate-Titration",
        "https://books.google.com/books?id=Jh4gAQAAIAAJ",
        "neutralization",
        "تفاعل حمض مع كربونات ينتج غازاً (Acid–Carbonate Gas Evolution)",
        "high",
    ),
    "acid_oxide": _rule(
        "acid_oxide",
        "قاعدة: حمض + أكسيد قلوي = ملح + ماء",
        "الأوكسيدات القلوية للمعادن (CaO, MgO, CuO …) تتفاعل مع الأحماض لتكوّن ملح المعدن وماء (المعادلة الصافية: O²⁻ + 2H⁺ → H₂O).",
        SOURCE_VOGEL,
        "VOGEL-6ed-Ch10-Metal-Oxides",
        "https://books.google.com/books?id=Jh4gAQAAIAAJ",
        "neutralization",
        "تفاعل حمض مع أكسيد فلزي (Acid–Basic Oxide)",
        "high",
    ),
    "acid_ammonia": _rule(
        "acid_ammonia",
        "قاعدة: حمض + أمونيا = ملح الأمونيوم",
        "الأمونيا قاعدة ضعيفة تتفاعل مع الأحماض مكوِّنةً أيون الأمونيوم (NH₄⁺) في الملح الناتج.",
        SOURCE_VOGEL,
        "VOGEL-6ed-Ch10-Ammonia-Titration",
        "https://books.google.com/books?id=Jh4gAQAAIAAJ",
        "neutralization",
        "تفاعل حمض مع أساس ضعيف (Acid–Weak Base)",
        "high",
    ),
    "acid_active_metal": _rule(
        "acid_active_metal",
        "قاعدة: حمض مخفف + معدن نشط = ملح + غاز الهيدروجين",
        "المعادن التي تسبق الهيدروجين في سلسلة النشاط تُزيح الهيدروجين من الأحماض المخففة (HCl و H₂SO₄ المخفف) "
        "مكوِّنةً ملحها و H₂، بينما لا ينطبق ذلك على الأحماض المؤكسدة مثل HNO₃.",
        SOURCE_BROWN,
        "BROWN-14ed-Ch4-Reactions-of-Metals-with-Acids",
        "https://www.pearson.com/en-us/subject-catalog/P/chemistry-the-central-science/P200000003123",
        "single_displacement",
        "تفاعل إحلال مفرد: حمض ومعدن نشط (Single Displacement — Hydrogen Evolution)",
        "medium",
    ),
    "alkali_water": _rule(
        "alkali_water",
        "قاعدة: معدن قلوي + ماء = هيدروكسيد + غاز الهيدروجين",
        "المعادن القلوية (Li, Na, K) تتفاعل مع الماء تاركةً ملح الهيدروكسيد الذائب ومحررةً غاز الهيدروجين، "
        "وقد تظهر رغوة على سطح الفلز أثناء التفاعل.",
        SOURCE_BROWN,
        "BROWN-14ed-Ch4-Reactions-of-Metals-with-Water",
        "https://www.pearson.com/en-us/subject-catalog/P/chemistry-the-central-science/P200000003123",
        "single_displacement",
        "تفاعل معدن قلوي مع الماء (Alkali Metal–Water Reaction)",
        "high",
    ),
    "combustion": _rule(
        "combustion",
        "قاعدة: احتراق الهيدروكربون = ثاني أكسيد الكربون + ماء",
        "اكتمال الاحتراق في وجود أكسجين زائد (Combustion in excess O₂) ينتج CO₂ و H₂O فقط، "
        "مع اشتعال وتطور حرارة، ومعادلة الاحتراق الكاملة: CₓHᵧ + (x + y/4) O₂ → x CO₂ + (y/2) H₂O.",
        SOURCE_CRC,
        "CRC-104ed-Combustion-Analysis",
        "https://www.routledge.com/CRC-Handbook-of-Chemistry-and-Physics-104th-Edition/p/book/9780367679981",
        "combustion",
        "تفاعل احتراق كامل (Complete Combustion)",
        "high",
    ),
    "metal_oxygen": _rule(
        "metal_oxygen",
        "قاعدة: معدن + أكسجين = أكسيد الفلز",
        "يتأكسد الفلز في وجود الأكسجين مكوّناً أكسيده المطابق لحالة أكسدة العنصر المعتادة.",
        SOURCE_CRC,
        "CRC-104ed-Metal-Oxides-Formation",
        "https://www.routledge.com/CRC-Handbook-of-Chemistry-and-Physics-104th-Edition/p/book/9780367679981",
        "redox_synthesis",
        "تفاعل اتحاد وأكسدة (Oxidation of a Metal by Oxygen)",
        "medium",
    ),
    "metal_halogen": _rule(
        "metal_halogen",
        "قاعدة: معدن + هالوجين = هاليد الفلز",
        "تتفاعل الفلزات مع الهالوجينات (Cl₂, Br₂, I₂, F₂) مكوِّنةً هاليداتها مع تحرير طاقة كبيرة.",
        SOURCE_CRC,
        "CRC-104ed-Metal-Halides",
        "https://www.routledge.com/CRC-Handbook-of-Chemistry-and-Physics-104th-Edition/p/book/9780367679981",
        "redox_synthesis",
        "تفاعل اتحاد مع هالوجين (Metal–Halogen Synthesis)",
        "high",
    ),
    "salt_salt_precipitation": _rule(
        "salt_salt_precipitation",
        "قاعدة: تبادل مزدوج بين ملحين مع راسب",
        "عند مزج محلولين لملحين إذا كانت إحدى النواتج صعبة الذوبان وفق قواعد الذوبانية فإنها تترسّب "
        "وتحصل التفاعل، وإن كانت النواتج جميعها ذائبة فلا يتغير شيء كيميائياً (لا يوجد صافي ردة فعل).",
        SOURCE_CRC,
        "CRC-104ed-Solubility-Rules",
        "https://www.routledge.com/CRC-Handbook-of-Chemistry-and-Physics-104th-Edition/p/book/9780367679981",
        "precipitation",
        "تفاعل ترسيب وتبادل مزدوج (Precipitation & Double Displacement)",
        "high",
    ),
    "salt_metal_displacement": _rule(
        "salt_metal_displacement",
        "قاعدة: معدن + ملح = ملح جديد + معدن منفصل",
        "يزيح المعدن النشط الفلز من ملح فلز أقل نشاطاً منه في سلسلة النشاط، فيترسّب الفلز المنفصل ويظهر في قاع المحلول.",
        SOURCE_BROWN,
        "BROWN-14ed-Ch4-Single-Displacement-Reactions",
        "https://www.pearson.com/en-us/subject-catalog/P/chemistry-the-central-science/P200000003123",
        "single_displacement",
        "تفاعل إحلال مفرد (Single Displacement)",
        "medium",
    ),
    "salt_hydroxide_precipitation": _rule(
        "salt_hydroxide_precipitation",
        "قاعدة: ملح + هيدروكسيد = راسب هيدروكسيد + ملح",
        "إضافة هيدروكسيد قاعدي إلى محلول ملح فلز يسبّب هيدروكسيد الفلز غير الذائب وفق قواعد الذوبانية.",
        SOURCE_CRC,
        "CRC-104ed-Precipitation-of-Metal-Hydroxides",
        "https://www.routledge.com/CRC-Handbook-of-Chemistry-and-Physics-104th-Edition/p/book/9780367679981",
        "precipitation",
        "تفاعل ترسيب هيدروكسيد فلزي (Metal Hydroxide Precipitation)",
        "high",
    ),
    "esterification": _rule(
        "esterification",
        "قاعدة: حمض كربوكسيلي + كحول = إستر + ماء",
        "تكثيف الحمض الكربوكسيلي مع الكحول ينتج إستراً وماءً تحت حفاز حمضي (H₂SO₄) وارتفاع حرارة، "
        "وهي تفاعل عكوس يعتمد على إزالة الماء.",
        SOURCE_CAREY,
        "CAREY-AdvOrgChem-PartA-Ch11-Esterification",
        "https://link.springer.com/book/10.1007/978-0-387-28247-5",
        "organic_reaction",
        "تفاعل استرة (Esterification)",
        "medium",
    ),
}


def _mk_species(formula: str, name_ar: str = "", name_en: str = "", state: str = "") -> dict[str, Any]:
    info = classify_species(formula, name_ar, name_en)
    info["state"] = state or get_standard_phase_code(info["formula"], "")
    info["name_ar"] = info["name_ar"] or info["formula"]
    info["name_en"] = info["name_en"] or info["formula"]
    return info


def _rule_result(
    rule_key: str,
    reactants: list[dict[str, Any]],
    product_formulas: list[str],
    product_kinds: list[str] | None = None,
    extra: dict[str, Any] | None = None,
) -> dict[str, Any]:
    rule = RULES[rule_key]
    products: list[dict[str, Any]] = []
    for idx, f in enumerate(product_formulas):
        kind = product_kinds[idx] if product_kinds and idx < len(product_kinds) else KIND_UNKNOWN
        sp = _mk_species(f)
        if sp["kind"] == KIND_UNKNOWN and kind != KIND_UNKNOWN:
            sp["kind"] = kind
        products.append(
            {
                "formula": sp["formula"],
                "name_ar": sp["name_ar"],
                "name_en": sp["name_en"],
                "coefficient": 1,
                "state": sp["state"],
                "kind": sp["kind"],
                "kind_label_ar": sp["kind_label_ar"],
            }
        )
    payload = {
        "rule": rule,
        "reactants": reactants,
        "products": products,
    }
    if extra:
        payload.update(extra)
    return payload


def _find_species(species: list[dict[str, Any]], kind: str) -> list[dict[str, Any]]:
    return [s for s in species if s["kind"] == kind]


def _infer_by_rules(species: list[dict[str, Any]]) -> dict[str, Any] | None:
    """Applies the documented general-chemistry rules in priority order."""
    if len(species) < 2:
        return None

    by_formula = {s["formula"]: s for s in species}

    # R1 — combustion of a hydrocarbon with oxygen
    hydrocarbons = _find_species(species, KIND_HYDROCARBON)
    oxygen = by_formula.get("O2")
    if hydrocarbons and oxygen:
        hyd = hydrocarbons[0]
        reactants = [hyd, oxygen]
        products = ["CO2", "H2O"]
        return _rule_result(
            "combustion",
            reactants,
            products,
            extra={
                "net_ionic": "لا ينطبق (معادلة جزيئية عامة لتفاعل احتراق) — Reaction is a molecular combustion equation.",
                "observations": [
                    {"type": "heat", "text_ar": "اشتعال وتطور حرارة شديدة (تفاعل طارد للحرارة بشدة).", "basis_ar": "قاعدة الاحتراق الكامل."},
                    {"type": "light", "text_ar": "انبعاث ضوء (لهب) عند الاحتراق.", "basis_ar": "الاحتراق يصدر ضوءاً مرئياً."},
                ],
                "ph_analysis": {"available": False, "summary_ar": "لا ينطبق حساب pH على تفاعل احتراق (لا يوجد محلول مائي).", "detail_ar": ""},
                "conditions_required_ar": "شرارة أو حرارة بدء + أكسجين زائد (الأكسجين في الهواء كافٍ).",
            },
        )

    # R2 — alkali metal with water
    for sp in species:
        if sp["kind"] == KIND_METAL and sp["cation"] in ALKALI_METALS and by_formula.get("H2O"):
            hydroxide = _salt_formula(sp["cation"], sp["cation_charge"], sp["cation_count"], "OH", 1)
            return _rule_result(
                "alkali_water",
                [sp, by_formula["H2O"]],
                [hydroxide, "H2"],
                extra={
                    "net_ionic": "2 Na + 2 H₂O → 2 Na⁺ + 2 OH⁻ + H₂↑ (صافي أيوني مبسّط).",
                    "observations": [
                        {"type": "gas", "text_ar": "تصاعد غاز الهيدروجين H₂ (فقاعات).", "basis_ar": "قاعدة تفاعل المعادن القلوية مع الماء."},
                        {"type": "heat", "text_ar": "ارتفاع محسوس في درجة الحرارة وطفاء المعدن المصهر.", "basis_ar": "تفاعل طارد للحرارة."},
                        {"type": "ph", "text_ar": "المحلول الناتج قاعدي بدرجة عالية.", "basis_ar": "تكون أيون OH⁻ في المحلول."},
                    ],
                    "ph_analysis": {"available": True, "summary_ar": "محلول قاعدي (pH > 7).", "detail_ar": "ينتج هيدروكسيد قابل للذوبان في الماء فيرتفع pH."},
                    "conditions_required_ar": "تفاعل تلقائي مع الماء في درجة الحرارة العادية.",
                },
            )

    # R3 — strong/weak acid with metal (hydrogen evolution) — only non-oxidizing acids
    acids = _find_species(species, KIND_STRONG_ACID) + _find_species(species, KIND_WEAK_ACID)
    for acid in acids:
        if acid["formula"] not in {"HCl", "H2SO4"}:
            continue
        for metal in _find_species(species, KIND_METAL):
            sym = metal["cation"]
            if sym in NOBLE_METALS or sym in ALKALI_METALS:
                continue
            if metal_activity_index(sym) < 0 or not is_more_active_metal(sym, "H"):
                continue
            anion = "Cl" if acid["formula"] == "HCl" else "SO4"
            metal_charge = METAL_ION_CHARGES.get(sym, 2)
            salt = _salt_formula(sym, metal_charge, 1, anion, 1)
            extra: dict[str, Any] = {
                "net_ionic": f"{sym} + 2 H⁺ → {sym}²⁺ + H₂↑ (مع ملح مقابل في المحلول).",
                "observations": [
                    {"type": "gas", "text_ar": "تصاعد فقاعات غاز الهيدروجين H₂.", "basis_ar": "المعدن يزيح الهيدروجين من الحمض."},
                    {"type": "heat", "text_ar": "ارتفاع في درجة حرارة المحلول.", "basis_ar": "تفاعل طارد للحرارة."},
                ],
                "ph_analysis": {"available": True, "summary_ar": "الحمض الزائد يبقى حمضياً (pH < 7)، والمحلول النهائي يكون حمضياً.", "detail_ar": "يتحدد pH النهائي بكمية الحمض المتبقية."},
                "conditions_required_ar": "حمض مخفف (HCl أو H₂SO₄ المخفف) — لا ينطبق على HNO₃ أو H₂SO₄ المركز.",
            }
            return _rule_result("acid_active_metal", [acid, metal], [salt, "H2"], extra=extra)

    # R4 — acid + carbonate / bicarbonate
    for acid in acids:
        for sp in species:
            if sp["kind"] != KIND_SALT or sp["anion"] not in {"CO3", "HCO3"}:
                continue
            if sp["cation"] in {"H", "NH4"} and sp["anion"] == "CO3":
                continue
            salt = _salt_formula(sp["cation"], sp["cation_count"], sp["anion"], 1)
            return _rule_result(
                "acid_carbonate",
                [acid, sp],
                [salt, "H2O", "CO2"],
                extra={
                    "net_ionic": "HCO₃⁻ + H⁺ → H₂O + CO₂↑ (أو CO₃²⁻ + 2 H⁺ → H₂O + CO₂↑).",
                    "observations": [
                        {"type": "gas", "text_ar": "تصاعد غاز ثاني أكسيد الكربون CO₂ (فوران وفقاعات).", "basis_ar": "تحرر CO₂ من الكربونات بالأحماض."},
                        {"type": "heat", "text_ar": "انخفاض حرارة المزيج عادةً (تفاعل بطيء الحرارة في كثير من الحالات).", "basis_ar": "قاعدة تفاعل الأحماض مع الكربونات."},
                    ],
                    "ph_analysis": {"available": True, "summary_ar": "الناتج حمضي أو محايد حسب فائض الحمض.", "detail_ar": "يتحدد pH النهائي بكمية الحمض الزائدة بعد تحرر CO₂."},
                    "conditions_required_ar": "محلول مائي؛ يجب أن يكون الحمض كافياً لتحرر CO₂.",
                },
            )

    # R5 — acid + metal hydroxide (neutralization)
    for acid in acids:
        for sp in species:
            if sp["kind"] not in {KIND_STRONG_BASE, KIND_SALT}:
                continue
            if sp["anion"] != "OH" or sp["cation"] in {"H", "NH4"}:
                continue
            salt = _salt_formula(sp["cation"], sp["cation_charge"], sp["cation_count"], acid["anion"], acid["anion_charge"])
            extra = {
                "net_ionic": "H⁺ + OH⁻ → H₂O (المعادلة الصافية للتعادل).",
                "observations": [
                    {"type": "heat", "text_ar": "ارتفاع في درجة الحرارة (تفاعل تعادل طارد للحرارة، ΔH ≈ −57 kJ/mol للأحماض والقواعد القوية).", "basis_ar": "حرارة التعادل القياسية."},
                    {"type": "indicator", "text_ar": "تحوّل لون مؤشر مناسب نحو نقطة التكافؤ (تُستخدم هذه الخاصية في المعايرة).", "basis_ar": "تكافؤ H⁺ مع OH⁻."},
                ],
                "ph_analysis": {"available": True, "summary_ar": "محايد عند التكافؤ (pH ≈ 7 لحمض قوي وقاعدة قوية).", "detail_ar": "يتحدد pH النهائي بفائض الحمض أو القاعدة."},
                "conditions_required_ar": "محلول مائي (التفاعل أيوني ولا يحدث في وسط غير مائي).",
            }
            if sp["cation"] not in ALKALI_METALS and sp["cation"] not in {"Ca", "Sr", "Ba"}:
                extra["observations"].append(
                    {"type": "precipitate", "text_ar": f"تكوّن راسب من هيدروكسيد {sp['name_ar']} (غير ذائب في الماء).", "basis_ar": "قواعد الذوبانية: هيدروكسيدات معظم الفلزات غير ذائبة."}
                )
            return _rule_result("acid_base_neutralization", [acid, sp], [salt, "H2O"], extra=extra)

    # R6 — acid + ammonia
    for acid in acids:
        for sp in species:
            if sp["kind"] != KIND_AMMONIA:
                continue
            salt = _salt_formula("NH4", 1, 1, acid["anion"], acid["anion_charge"])
            return _rule_result(
                "acid_ammonia",
                [acid, sp],
                [salt],
                extra={
                    "net_ionic": "NH₃ + H⁺ → NH₄⁺ (مع الأيون المقابل في الملح).",
                    "observations": [
                        {"type": "heat", "text_ar": "ارتفاع في درجة الحرارة.", "basis_ar": "تعادل حمض مع أساس ضعيف طارد للحرارة."},
                        {"type": "smoke", "text_ar": "تكوّن سحب بيضاء من كلوريد الأمونيوم عند ملامسة الرطوبة أو خروج الغاز من وعاء مغلق.", "basis_ar": "تكوّن ملح الأمونيوم المتطاير."},
                    ],
                    "ph_analysis": {"available": True, "summary_ar": "محلول حمضي ضعيف (محلول ملح الأمونيوم حمضي) إلا إذا زاد الأمونيا.", "detail_ar": "أيون الأمونيوم يذيب الماء تسامياً (ت أثر أساسي ضعيف)."},
                    "conditions_required_ar": "محلول مائي أو غاز جاف.",
                },
            )

    # R7 — acid + basic metal oxide
    for acid in acids:
        for sp in species:
            if sp["kind"] != KIND_OXIDE or sp["cation"] in {"", "H"}:
                continue
            if sp["cation"] in {"Si", "P", "S", "C", "N"}:
                continue  # acidic oxides follow a different (acid–acid) pathway
            salt = _salt_formula(sp["cation"], sp["cation_charge"], sp["cation_count"], acid["anion"], acid["anion_charge"])
            return _rule_result(
                "acid_oxide",
                [acid, sp],
                [salt, "H2O"],
                extra={
                    "net_ionic": "O²⁻ + 2 H⁺ → H₂O (المعادلة الصافية).",
                    "observations": [
                        {"type": "dissolution", "text_ar": "ذوبان الأوكسيد الصلب تدريجياً في الحمض.", "basis_ar": "الأوكسيدات القلوية تذوب في الأحماض."},
                        {"type": "heat", "text_ar": "ارتفاع في درجة الحرارة.", "basis_ar": "تفاعل أكسيد وحمض طارد للحرارة."},
                    ],
                    "ph_analysis": {"available": True, "summary_ar": "محايد عند التكافؤ أو حمضي بفائض الحمض.", "detail_ar": "يتحدد pH النهائي بالفائض."},
                    "conditions_required_ar": "محلول مائي.",
                },
            )

    # R8 — salt + salt (precipitation or no net reaction)
    salts = [s for s in species if s["kind"] == KIND_SALT and s["anion"] not in {"OH", "O"}]
    if len(salts) >= 2:
        for i in range(len(salts)):
            for j in range(i + 1, len(salts)):
                a, b = salts[i], salts[j]
                if a["cation"] == b["cation"] and a["anion"] == b["anion"]:
                    continue
                prod1 = _salt_formula(a["cation"], a["cation_charge"], a["cation_count"], b["anion"], b["anion_charge"])
                prod2 = _salt_formula(b["cation"], b["cation_charge"], b["cation_count"], a["anion"], a["anion_charge"])
                p1 = _mk_species(prod1)
                p2 = _mk_species(prod2)
                precip = None
                for cand in (p1, p2):
                    if is_salt_insoluble(cand["cation"], cand["anion"]):
                        precip = cand
                        break
                if precip is None:
                    continue
                other = p2 if precip is p1 else p1
                color = PRECIPITATE_COLORS_AR.get(precip["formula"], "")
                obs = [
                    {
                        "type": "precipitate",
                        "text_ar": f"تكوّن راسب من {precip['name_ar']} ({precip['formula']})" + (f" لونه {color}." if color else "."),
                        "basis_ar": "قواعد الذوبانية: الملح الناتج صعب الذوبان في الماء.",
                    }
                ]
                if other["anion"] in {"NO3", "Cl", "Br", "I", "SO4"} and other["cation"] in ALWAYS_SOLUBLE_CATIONS:
                    obs.append({"type": "dissolution", "text_ar": f"يبقى {other['name_ar']} ({other['formula']}) ذائباً في المحلول.", "basis_ar": "ذوبانية الملح المقابل."})
                return _rule_result(
                    "salt_salt_precipitation",
                    [a, b],
                    [precip["formula"], other["formula"]],
                    extra={
                        "net_ionic": f"{precip['cation']}{_charge_text(precip['cation_charge'])} + {other['cation']}{_charge_text(other['cation_charge'])} → {precip['formula']}(↓)",
                        "observations": obs,
                        "ph_analysis": {"available": False, "summary_ar": "لا يتأثر الـ pH جوهرياً في تبادل الأيونات (لا تُستهلك أيونات H⁺ أو OH⁻).", "detail_ar": ""},
                        "conditions_required_ar": "محلول مائي للطرفين.",
                    },
                )

    # R9 — salt + metal (displacement)
    for metal in _find_species(species, KIND_METAL):
        if metal["cation"] in NOBLE_METALS:
            continue
        for salt in salts:
            if salt["cation"] in {"H", "NH4"} or salt["cation"] == metal["cation"]:
                continue
            if not is_more_active_metal(metal["cation"], salt["cation"]):
                continue
            new_salt = _salt_formula(
                metal["cation"],
                METAL_ION_CHARGES.get(metal["cation"], 2),
                1,
                salt["anion"],
                salt["anion_charge"],
            )
            return _rule_result(
                "salt_metal_displacement",
                [metal, salt],
                [new_salt, salt["cation"]],
                extra={
                    "net_ionic": f"{metal['cation']} + {salt['cation']}²⁺ → {metal['cation']}²⁺ + {salt['cation']} (صافي أيوني).",
                    "observations": [
                        {"type": "precipitate", "text_ar": f"ترسّب فلز {salt['cation']} في قاع المحلول (راسب معدني).", "basis_ar": "قاعدة سلسلة النشاط والإحلال المفرد."},
                        {"type": "color", "text_ar": "قد يتغير لون المحلول إذا كان أيون الفلز المنقول ذائناً ملوناً (مثل Cu²⁺ الأزرق).", "basis_ar": "اختلاف ألوان الأيونات الفلزية الذائبة."},
                    ],
                    "ph_analysis": {"available": False, "summary_ar": "لا يتغير pH جوهرياً (تفاعل أكسدة-اختزال لا يستهلك H⁺ أو OH⁻).", "detail_ar": ""},
                    "conditions_required_ar": "محلول مائي للملح.",
                },
            )

    # R10 — salt + hydroxide (precipitation of metal hydroxide)
    for salt in salts:
        for base in species:
            if base["kind"] != KIND_STRONG_BASE or base["anion"] != "OH":
                continue
            if salt["cation"] in {"H", "NH4"}:
                continue
            hydroxide = _salt_formula(salt["cation"], salt["cation_charge"], salt["cation_count"], "OH", 1)
            if not is_salt_insoluble(salt["cation"], "OH"):
                continue
            other_salt = _salt_formula(base["cation"], base["cation_charge"], base["cation_count"], salt["anion"], salt["anion_charge"])
            color = PRECIPITATE_COLORS_AR.get(hydroxide, "")
            return _rule_result(
                "salt_hydroxide_precipitation",
                [salt, base],
                [hydroxide, other_salt],
                extra={
                    "net_ionic": f"{salt['cation']}{_charge_text(salt['cation_charge'])} + 2 OH⁻ → {hydroxide}(↓)",
                    "observations": [
                        {
                            "type": "precipitate",
                            "text_ar": f"تكوّن راسب من {hydroxide}" + (f" لونه {color}." if color else "."),
                            "basis_ar": "قواعد الذوبانية لهيدروكسيدات الفلزات.",
                        },
                        {"type": "ph", "text_ar": "ارتفاع pH مع تكوّن الراسب (استهلاك أيونات OH⁻).", "basis_ar": "التفاعل يستهلك OH⁻ من المحلول."},
                    ],
                    "ph_analysis": {"available": True, "summary_ar": "يتحدد pH النهائي بكمية القاعدة الزائدة.", "detail_ar": "عند التكافؤ يكون المحلول قاعدياً قليلاً بسبب التسرب الحلزوني للفلز."},
                    "conditions_required_ar": "محلول مائي.",
                },
            )

    # R11 — metal + halogen
    for metal in _find_species(species, KIND_METAL):
        for sp in species:
            if sp["kind"] != KIND_ELEMENTAL_MOLECULE or sp["formula"] not in {"Cl2", "Br2", "I2", "F2"}:
                continue
            halogen = sp["formula"][0]
            halide = _salt_formula(metal["cation"], METAL_ION_CHARGES.get(metal["cation"], 2), 1, halogen, 1)
            return _rule_result(
                "metal_halogen",
                [metal, sp],
                [halide],
                extra={
                    "net_ionic": "لا ينطبق (معادلة جزيئية) — Metal–Halogen Synthesis.",
                    "observations": [
                        {"type": "heat", "text_ar": "تطور حرارة كبيرة جداً مع ضوء (احتراق الفلز).", "basis_ar": "تفاعل فلز مع هالوجين طارد للحرارة بشدة."},
                        {"type": "smoke", "text_ar": "تكوّن دخان غبار من هاليد الفلز.", "basis_ar": "معظم الهاليدات الفلزية غبارية ومتطايرة."},
                    ],
                    "ph_analysis": {"available": False, "summary_ar": "لا ينطبق (لا يوجد محلول مائي).", "detail_ar": ""},
                    "conditions_required_ar": "حرارة أو شرارة بدء.",
                },
            )

    # R12 — metal + oxygen
    for metal in _find_species(species, KIND_METAL):
        if "O2" not in by_formula:
            continue
        oxide_formula = METAL_OXIDES.get(metal["cation"])
        if not oxide_formula:
            continue
        return _rule_result(
            "metal_oxygen",
            [metal, by_formula["O2"]],
            [oxide_formula],
            extra={
                "net_ionic": "لا ينطبق (معادلة جزيئية) — Oxidation of a Metal.",
                "observations": [
                    {"type": "heat", "text_ar": "ارتفاع في درجة الحرارة (تكوّن الأكسيد طارد للحرارة).", "basis_ar": "أكسدة الفلز في وجود الأكسجين."},
                ],
                "ph_analysis": {"available": False, "summary_ar": "لا ينطبق (لا يوجد محلول مائي).", "detail_ar": ""},
                "conditions_required_ar": "حرارة أو تسخين (وكثير من الفلزات تتأكسد ببطء في درجة حرارة الغرفة).",
            },
        )

    # R13 — carboxylic acid + alcohol (esterification)
    organic_acids = [s for s in species if s["kind"] == KIND_WEAK_ACID and s["anion"] in {"CH3COO", "HCOO"}]
    alcohols = _find_species(species, KIND_ORGANIC_ALCOHOL)
    if organic_acids and alcohols:
        acid, alcohol = organic_acids[0], alcohols[0]
        ester_counts = dict(acid["counts"])
        for el, cnt in alcohol["counts"].items():
            ester_counts[el] = ester_counts.get(el, 0) + cnt
        ester_counts["H"] = ester_counts.get("H", 0) - 2
        ester_counts["O"] = ester_counts.get("O", 0) - 1
        ester_counts = {el: cnt for el, cnt in ester_counts.items() if cnt > 0}
        ester_formula = _hill_formula(ester_counts)
        return _rule_result(
            "esterification",
            [acid, alcohol],
            [ester_formula, "H2O"],
            extra={
                "net_ionic": "لا ينطبق (تفاعل عضوي تكثيفي غير أيوني) — Condensation (Fischer Esterification).",
                "observations": [
                    {"type": "smell", "text_ar": "رائحة مميزة للإستر الناتج ورائحة خفيفة للكحول أثناء التسخين.", "basis_ar": "الإسترات ذات روائح مميزة."},
                    {"type": "heat", "text_ar": "تسخين لطيف مع حفاز حمضي (H₂SO₄ مركّز) وإزالة الماء.", "basis_ar": "تفاعل تكثيف مُحفَّز حمضياً."},
                ],
                "ph_analysis": {"available": False, "summary_ar": "لا ينطبق (وسط عضوي غير مائي).", "detail_ar": ""},
                "conditions_required_ar": "حمض مركّز + تسخين + إزالة الماء باستمرار.",
                "limitations_ar": [
                    "تم اشتقاق صيغة الإستر على مستوى التركيب الذري (C/H/O) وليس على مستوى البنية الفراغية، "
                    "لذا لا يمكن الجزم بتحديد موضع مجموعة الاستر في الجزيء، وقد تتشارك عدة أيزومرات نفس الصيغة.",
                    "صيغة الكحول C/H/O ذات الأكسجين الواحد قد تكون كحولاً أو إيثراً (مثل CH₃OCH₃)؛ والتفاعل يختلف بين الحالتين.",
                ],
            },
        )

    return None


def _hill_formula(counts: dict[str, int]) -> str:
    from chemistry_engine.formula_analyzer import format_hill_formula

    return format_hill_formula(counts)


def _charge_text(charge: int) -> str:
    """Renders an ionic charge using Unicode superscripts (e.g. 2+, 2−)."""
    if charge == 0:
        return ""
    magnitude = abs(charge)
    sup = {1: "¹", 2: "²", 3: "³", 4: "⁴", 5: "⁵"}.get(magnitude, str(magnitude))
    return f"{sup}+" if charge > 0 else f"{sup}⁻"


# --------------------------------------------------------------------------------------
# 5) Public API
# --------------------------------------------------------------------------------------


def _build_observations_payload(extra: dict[str, Any]) -> list[dict[str, Any]]:
    return list(extra.get("observations", []))


def _insufficient_payload(reason_ar: str, missing_ar: list[str], species: list[dict[str, Any]]) -> dict[str, Any]:
    return {
        "supported": False,
        "mode": "insufficient_evidence",
        "mode_label_ar": "لا تتوفر بيانات كافية (Insufficient Evidence)",
        "verdict_ar": reason_ar,
        "confidence": {"level": "none", "label_ar": "لا يمكن تقديره", "basis_ar": "لا يوجد أساس استدلال موثوق لهذه المجموعة."},
        "reactants": [
            {
                "formula": s["formula"],
                "name_ar": s["name_ar"],
                "name_en": s["name_en"],
                "kind": s["kind"],
                "kind_label_ar": s["kind_label_ar"],
                "state": s["state"],
            }
            for s in species
        ],
        "products": [],
        "observations": [],
        "ph_analysis": {"available": False, "summary_ar": "لا يتوفر حساب pH لعدم وجود مسار تفاعل موثوق.", "detail_ar": ""},
        "evidence": [],
        "missing_information_ar": missing_ar,
        "explanation_ar": (
            "تعذّر تحديد مسار تفاعل موثوق لهذه المواد. "
            "تلتزم منصة MolecuLens بمبدأ عدم اختلاق النواتج: يُعرض غياب الأدلة بدل معادلة مُقدَّرة قد تكون مضللة."
        ),
    }


def predict_reaction(
    reactant_inputs: list[dict[str, Any]],
    domain: str | None = None,
    custom_conditions: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Predicts whether the supplied substances will react and what the products are.

    Evidence tiers:
      1. "documented"          — a reaction recorded in the verified registry.
      2. "rule_inferred"       — derived from a declared general-chemistry rule (with source + limits).
      3. "insufficient_evidence" — no reliable basis; products are NOT invented.
    """
    species: list[dict[str, Any]] = []
    for item in reactant_inputs:
        sub = item.get("substance", {}) or {}
        formula = sub.get("symbol_or_formula") or item.get("formula", "")
        if not formula:
            continue
        sp = _mk_species(formula, sub.get("name_ar", ""), sub.get("name_en", ""))
        sp["input_amount"] = item.get("amount")
        sp["input_unit"] = item.get("unit", "mol")
        species.append(sp)

    if len(species) < 2:
        return _insufficient_payload(
            "يلزم إدخال مادتين على الأقل لاستكشاف التفاعل المحتمل بينهما.",
            ["أضف مادة ثانية (عنصر، مركب، حمض، قاعدة، أو ملح) للمقارنة."],
            species,
        )

    invalid = [s for s in species if not s["is_valid_formula"]]
    if invalid:
        return _insufficient_payload(
            "إحدى المواد المدخلة بصيغة كيميائية غير قابلة للقراءة، لذلك لم يُجرَ التحليل.",
            [f"صحّح الصيغة المدخلة للمادة: {s['input_formula']}" for s in invalid],
            species,
        )

    # ---- Tier 1: documented reactions from the verified registry -----------------------
    from reaction_engine.reaction_registry import find_verified_reactions

    formulas = [s["formula"] for s in species]
    documented = find_verified_reactions(formulas, reaction_type=None, domain=None)
    if documented:
        rxn = documented[0]
        payload = _documented_payload(rxn, species)
        payload["alternatives"] = [
            {"reaction_id": r["reaction_id"], "reaction_name_ar": r["reaction_name_ar"], "category": r.get("category", "")}
            for r in documented[1:]
        ]
        return payload

    # ---- Tier 2: rule-based inference ------------------------------------------------
    inference = _infer_by_rules(species)
    if inference is None:
        return _insufficient_payload(
            "لا توجد قاعدة كيميائية موثقة تسمح بالتنبؤ بتفاعل هذه المواد معاً في المحرك الحالي.",
            [
                "أضف بيانات مرجعية موثقة (مرجع علمي أو تجربة مخبرية) لهذه المجموعة إلى سجل التفاعلات.",
                "أو حدّد نوع التفاعل المتوقع يدوياً واستخدم Reaction Lab الموثّق.",
                "أو استخدم محرك التنبؤ القواعدي المدعوم بقواعد التفاعل العضوي بعد توسيع قاعدة المعرفة.",
            ],
            species,
        )

    return _finalize_inferred(inference, species, reactant_inputs)


def _documented_payload(rxn: dict[str, Any], species: list[dict[str, Any]]) -> dict[str, Any]:
    reactants = [
        {
            "formula": r["formula"],
            "name_ar": r["name_ar"],
            "name_en": r.get("name_en", r["formula"]),
            "coefficient": r.get("coefficient", 1),
            "state": r.get("state", ""),
            "kind": "documented",
            "kind_label_ar": "مادة متفاعلة موثّقة",
        }
        for r in rxn["reactants"]
    ]
    products = [
        {
            "formula": p["formula"],
            "name_ar": p["name_ar"],
            "name_en": p.get("name_en", p["formula"]),
            "coefficient": p.get("coefficient", 1),
            "state": p.get("state", ""),
            "kind": "documented",
            "kind_label_ar": "ناتج موثّق",
        }
        for p in rxn["products"]
    ]
    eq = format_balanced_equation(rxn["reactants"], rxn["products"])
    conservation = verify_atom_conservation(rxn["reactants"], rxn["products"])
    observations, net_ionic, ph_analysis = _enrich_documented_analysis(rxn, products)
    return {
        "supported": True,
        "mode": "documented",
        "mode_label_ar": "تفاعل موثّق في سجل التفاعلات المعتمد (Documented Reaction)",
        "verdict_ar": (
            f"نعم — يوجد تفاعل موثّق مخبرياً لهذه المواد: {rxn['reaction_name_ar']}، "
            f"والنواتج المسجلة هي: " + "، ".join(f"{p['name_ar']} ({p['formula']})" for p in products) + "."
        ),
        "confidence": {
            "level": "high",
            "label_ar": "ثقة عالية — تفاعل موثّق بمرجع علمي محدّد",
            "basis_ar": rxn.get("source", ""),
        },
        "reactants": reactants,
        "products": products,
        "equation_text": eq["text"],
        "equation_unicode": format_chemical_equation_display(eq["text"]),
        "equation_latex": eq["latex"],
        "category": rxn.get("category", ""),
        "category_label_ar": rxn.get("reaction_type_ar", ""),
        "reaction_type_ar": rxn.get("reaction_type_ar", ""),
        "domain": rxn.get("domain", ""),
        "domain_label_ar": rxn.get("domain_label_ar", ""),
        "net_ionic": net_ionic,
        "observations": observations,
        "ph_analysis": ph_analysis,
        "evidence": [
            {
                "rule_id": rxn.get("reaction_id", ""),
                "rule_name_ar": rxn.get("reaction_name_ar", ""),
                "basis_ar": rxn.get("explanation_ar", ""),
                "source": rxn.get("source", ""),
                "source_id": rxn.get("source_id", ""),
                "reference_url": rxn.get("reference_url", ""),
                "confidence": "high",
            }
        ],
        "documentation_type": rxn.get("documentation_type", ""),
        "documentation_label_ar": rxn.get("documentation_label_ar", ""),
        "conditions": dict(rxn.get("conditions", {})),
        "atom_conservation": conservation,
        "explanation_ar": rxn.get("explanation_ar", ""),
        "assumptions_and_limits_ar": rxn.get("assumptions_and_limits_ar", ""),
        "alternatives": [],
    }


def _enrich_documented_analysis(
    rxn: dict[str, Any],
    products: list[dict[str, Any]],
) -> tuple[list[dict[str, Any]], str, dict[str, Any]]:
    """
    Derives the net ionic equation, the expected laboratory observations and the pH behaviour
    of a documented reaction from its category and its product list.
    """
    category = rxn.get("category", "")
    observations: list[dict[str, Any]] = []
    net_ionic = ""
    ph_analysis: dict[str, Any] = {
        "available": False,
        "summary_ar": "لا يُحسب pH إلا في التفاعلات الحمضية-القاعدية.",
        "detail_ar": "",
    }

    product_formulas = [p["formula"] for p in products]
    reactant_formulas = [r["formula"] for r in rxn.get("reactants", [])]

    if category == "neutralization":
        has_carbonate = any("CO3" in f or "HCO3" in f for f in reactant_formulas)
        has_hydroxide = any("OH" in f for f in reactant_formulas)
        if has_hydroxide:
            net_ionic = "H⁺ + OH⁻ → H₂O (المعادلة الصافية للتعادل)."
            observations.append(
                {
                    "type": "heat",
                    "text_ar": "ارتفاع في درجة الحرارة (تفاعل تعادل طارد للحرارة، ΔH ≈ −57 kJ/mol للأحماض والقواعد القوية).",
                    "basis_ar": "حرارة التعادل القياسية في المحاليل المائية المخففة.",
                }
            )
            observations.append(
                {
                    "type": "indicator",
                    "text_ar": "تحوّل لون المؤشر المناسب نحو نقطة التكافؤ (وهي أساس طريقة المعايرة).",
                    "basis_ar": "تكافؤ H⁺ مع OH⁻ بنسبة 1:1.",
                }
            )
            ph_analysis = {
                "available": True,
                "summary_ar": "محايد عند التكافؤ (pH ≈ 7 لحمض قوي وقاعدة قوية)، وحمضي أو قاعدي حسب فائض أحد الطرفين.",
                "detail_ar": "يتحدد pH النهائي بكمية الحمض أو القاعدة الزائدة.",
            }
        elif has_carbonate:
            net_ionic = "HCO₃⁻ + H⁺ → H₂O + CO₂↑ (أو CO₃²⁻ + 2 H⁺ → H₂O + CO₂↑)."
            ph_analysis = {
                "available": True,
                "summary_ar": "الحمض الزائد يبقى حمضياً (pH < 7) بعد تحرر غاز CO₂.",
                "detail_ar": "يتحدد pH النهائي بكمية الحمض المتبقية.",
            }
        else:
            net_ionic = "H⁺ + B (أساس) → BH⁺ (تكافؤ بروتون مع أساس)."
            ph_analysis = {
                "available": True,
                "summary_ar": "يتحدد pH النهائي بفائض الحمض أو الأساس.",
                "detail_ar": "في حالة التكافؤ مع أساس ضعيف يكون المحلول حمضياً قليلاً (تأثر أيوني).",
            }

    if category == "precipitation":
        for p in products:
            info = classify_species(p["formula"])
            if info["cation"] and info["anion"] and is_salt_insoluble(info["cation"], info["anion"]):
                color = PRECIPITATE_COLORS_AR.get(p["formula"], "")
                observations.append(
                    {
                        "type": "precipitate",
                        "text_ar": f"تكوّن راسب من {p['name_ar']} ({p['formula']})" + (f" لونه {color}." if color else "."),
                        "basis_ar": "قواعد الذوبانية المعتمدة في التحليل الكيفي.",
                    }
                )
                net_ionic = (
                    f"{info['cation']}{_charge_text(info['cation_charge'])} + "
                    f"{info['anion']}{_charge_text(info['anion_charge'])} → {p['formula']}(↓)"
                )
                break

    if category == "combustion":
        observations.append(
            {"type": "heat", "text_ar": "اشتعال وتطور حرارة شديدة (احتراق طارد للحرارة بشدة).", "basis_ar": "تفاعل احتراق تام."}
        )
        observations.append(
            {"type": "light", "text_ar": "انبعاث ضوء (لهب) أثناء الاحتراق.", "basis_ar": "الاحتراق يصدر ضوءاً مرئياً."}
        )

    for p in products:
        if p["formula"] in {"CO2", "H2", "O2", "NH3"}:
            observations.append(
                {"type": "gas", "text_ar": f"تصاعد {p['name_ar']} ({p['formula']}) على هيئة فقاعات.", "basis_ar": "ناتج غازي مسجّل في التفاعل الموثّق."}
            )

    if category == "single_displacement" and "H2" in product_formulas:
        net_ionic = "M + 2 H⁺ → M²⁺ + H₂↑ (تفاعل فلز مع حمض مخفف)."
    elif category == "single_displacement":
        metal_products = [p for p in products if classify_species(p["formula"])["kind"] == KIND_METAL]
        salt_products = [p for p in products if classify_species(p["formula"])["kind"] == KIND_SALT]
        if metal_products and salt_products:
            salt_info = classify_species(salt_products[0]["formula"])
            net_ionic = (
                f"M + {salt_info['cation']}{_charge_text(salt_info['cation_charge'])} → "
                f"M²⁺ + {metal_products[0]['formula']} (إحلال مفرد)"
            )

    if category == "redox_synthesis":
        observations.append(
            {"type": "heat", "text_ar": "ارتفاع في درجة الحرارة أثناء التفاعل (تكوّن الأكسيد طارد للحرارة).", "basis_ar": "أكسدة الفلز في وجود الأكسجين."}
        )
    elif category in {"combustion", "redox_synthesis", "organic_reaction"} and not net_ionic:
        net_ionic = "لا ينطبق (معادلة جزيئية عامة وليست أيونية) — Molecular (non-ionic) equation."

    return observations, net_ionic, ph_analysis


def _finalize_inferred(
    inference: dict[str, Any],
    species: list[dict[str, Any]],
    reactant_inputs: list[dict[str, Any]],
) -> dict[str, Any]:
    rule = inference["rule"]
    reactants = inference["reactants"]
    products = inference["products"]

    reactant_formulas = [r["formula"] for r in reactants]
    product_formulas = [p["formula"] for p in products]

    balance = balance_reaction_algebraically(reactant_formulas, product_formulas)
    if balance is not None:
        r_coeffs, p_coeffs = balance
        for idx, c in enumerate(r_coeffs):
            reactants[idx]["coefficient"] = c
        for idx, c in enumerate(p_coeffs):
            products[idx]["coefficient"] = c
    else:
        for r in reactants:
            r.setdefault("coefficient", 1)
        for p in products:
            p.setdefault("coefficient", 1)

    r_specs = [
        {
            "formula": r["formula"],
            "name_ar": r["name_ar"],
            "name_en": r["name_en"],
            "coefficient": r["coefficient"],
            "state": r.get("state", ""),
        }
        for r in reactants
    ]
    p_specs = [
        {
            "formula": p["formula"],
            "name_ar": p["name_ar"],
            "name_en": p["name_en"],
            "coefficient": p["coefficient"],
            "state": p.get("state", ""),
        }
        for p in products
    ]

    eq = format_balanced_equation(r_specs, p_specs)
    conservation = verify_atom_conservation(r_specs, p_specs)

    stoichiometry: dict[str, Any] | None = None
    yields: list[dict[str, Any]] = []
    reactants_table: list[dict[str, Any]] = []
    limiting_reagent = None
    try:
        user_inputs = [
            {
                "substance": {
                    "symbol_or_formula": item.get("substance", {}).get("symbol_or_formula", item.get("formula", "")),
                    "name_ar": item.get("substance", {}).get("name_ar", ""),
                    "name_en": item.get("substance", {}).get("name_en", ""),
                },
                "formula": item.get("formula", ""),
                "amount": item.get("amount", 1.0),
                "unit": item.get("unit", "mol"),
            }
            for item in reactant_inputs
        ]
        stoich = compute_reaction_stoichiometry(r_specs, user_inputs)
        yields = compute_theoretical_yields(p_specs, stoich["reaction_extent_mol"])
        reactants_table = stoich["reactants_breakdown"]
        limiting_reagent = stoich["limiting_reagent"]
        stoichiometry = {
            "reaction_extent_mol": stoich["reaction_extent_mol"],
            "is_exact_stoichiometric_ratio": stoich["is_exact_stoichiometric_ratio"],
            "limiting_reagent": limiting_reagent,
            "excess_reagents": stoich["excess_reagents"],
        }
    except Exception as exc:  # pragma: no cover - defensive: never break the prediction flow
        stoichiometry = {"error_ar": f"تعذّر إتمام الحساب الستوكيومتري: {exc}"}

    extra = inference
    observations = _build_observations_payload(extra)
    limitations = list(extra.get("limitations_ar", []))
    limitations.append(
        "هذا الاستنتاج مبني على قاعدة كيميائية عامة وليس على تجربة مخبرية موثّقة لهذه التركيبة تحديداً؛ "
        "يُنصح بالرجوع إلى مرجع تفاعلي معتمد قبل الاعتماد عليه في تصميم عملية صناعية."
    )

    ph = extra.get("ph_analysis", {"available": False, "summary_ar": "", "detail_ar": ""})
    net_ionic = extra.get("net_ionic", "")

    explanation_ar = _build_explanation_ar(
        species=species,
        reactants=reactants,
        products=products,
        rule=rule,
        equation_unicode=format_chemical_equation_display(eq["text"]),
        observations=observations,
        ph=ph,
        stoichiometry=stoichiometry,
    )

    return {
        "supported": True,
        "mode": "rule_inferred",
        "mode_label_ar": "استنتاج قواعدي — غير موثّق تجريبياً (Rule-Based Inference)",
        "verdict_ar": (
            f"نعم — وفق قاعدة كيميائية عامة ({rule['rule_name_ar']}) يتوقع تفاعل هذه المواد، "
            "والنواتج المتوقعة هي: " + "، ".join(f"{p['name_ar']} ({p['formula']})" for p in products) + "."
        ),
        "confidence": {
            "level": rule["confidence"],
            "label_ar": {
                "high": "ثقة عالية — قاعدة كيميائية أساسية",
                "medium": "ثقة متوسطة — قاعدة صحيحة ضمن شروط محددة",
                "low": "ثقة منخفضة — قاعدة تقريبية تحتاج تحققاً",
            }.get(rule["confidence"], rule["confidence"]),
            "basis_ar": rule["basis_ar"],
        },
        "reactants": reactants,
        "products": products,
        "equation_text": eq["text"],
        "equation_unicode": format_chemical_equation_display(eq["text"]),
        "equation_latex": eq["latex"],
        "category": rule["category"],
        "category_label_ar": rule["reaction_type_ar"],
        "reaction_type_ar": rule["reaction_type_ar"],
        "net_ionic": net_ionic,
        "observations": observations,
        "ph_analysis": ph,
        "evidence": [
            {
                "rule_id": rule["rule_id"],
                "rule_name_ar": rule["rule_name_ar"],
                "basis_ar": rule["basis_ar"],
                "source": rule["source"],
                "source_id": rule["source_id"],
                "reference_url": rule["reference_url"],
                "confidence": rule["confidence"],
            }
        ],
        "conditions_required_ar": extra.get("conditions_required_ar", ""),
        "limitations_ar": limitations,
        "atom_conservation": conservation,
        "stoichiometry": stoichiometry,
        "reactants_table": reactants_table,
        "theoretical_yields": yields,
        "limiting_reagent": limiting_reagent,
        "explanation_ar": explanation_ar,
        "alternatives": [],
    }


def _build_explanation_ar(
    species: list[dict[str, Any]],
    reactants: list[dict[str, Any]],
    products: list[dict[str, Any]],
    rule: dict[str, Any],
    equation_unicode: str,
    observations: list[dict[str, Any]],
    ph: dict[str, Any],
    stoichiometry: dict[str, Any] | None,
) -> str:
    lines: list[str] = []

    lines.append("**الخطوة 1 — تصنيف المواد المدخلة:**")
    for sp in species:
        lines.append(f"- {sp['name_ar']} ({sp['formula']}): {sp['kind_label_ar']}.")
    lines.append("")

    lines.append(f"**الخطوة 2 — القاعدة الكيميائية المطبّقة ({rule['rule_id']}):** {rule['basis_ar']}")
    lines.append("")

    lines.append(f"**الخطوة 3 — المعادلة الموزونة المتوقعة:** {equation_unicode}")
    lines.append("")

    lines.append("**الخطوة 4 — النواتج المتوقعة:** " + "، ".join(f"{p['name_ar']} ({p['formula']})" for p in products) + ".")
    lines.append("")

    if stoichiometry and "error_ar" not in stoichiometry:
        lim = stoichiometry.get("limiting_reagent")
        if lim:
            lines.append(
                f"**الخطوة 5 — المتفاعل المحدد (Limiting Reagent):** {lim['name_ar']} ({lim['formula']}) "
                f"بكمية {lim['initial_mol']:g} مول، والكمية المتبقية بعد التفاعل تُحسب وفق النسبة المولية للتفاعل."
            )
        else:
            lines.append("**الخطوة 5 — النسبة المولية:** الكميات المدخلة متكافئة مع النسبة المولية للتفاعل (لا يوجد متفاعل محدد).")
    else:
        lines.append("**الخطوة 5 — الحساب الستوكيومتري:** لم تُدخل كميات رقمية، لذلك لم يُحسب المتفاعل المحدد أو النواتج النظرية.")
    lines.append("")

    if observations:
        lines.append("**الخطوة 6 — الملاحظات المتوقعة (مبنية على قواعد كيميائية):**")
        for obs in observations:
            lines.append(f"- {obs['text_ar']}")
        lines.append("")

    if ph.get("available"):
        lines.append(f"**الخطوة 7 — التغيّر في درجة الحموضة (pH):** {ph.get('summary_ar', '')} {ph.get('detail_ar', '')}")
    else:
        lines.append("**الخطوة 7 — التغيّر في درجة الحموضة (pH):** لا ينطبق على هذا النوع من التفاعلات أو تتوفر بيانات غير كافية.")

    return "\n".join(lines)

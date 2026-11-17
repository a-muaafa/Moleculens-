"""
Extended acid / base / indicator catalogue for the MolecuLens Titration Lab.

Design rules (scientific honesty):
  * Every species carries its constants (pKa / pKb), classification, solubility,
    properties and a bibliographic source for the constants.
  * Molar masses are NOT typed by hand: they are computed from the molecular formula with
    IUPAC standard atomic weights, so they cannot drift from the formula.
  * Ka / Kb are derived from the pKa / pKb list, so the numbers are internally consistent.
  * Constants are taken from standard compilations (CRC Handbook, Harris QCA, Williams pKa
    table). The citation is at table level; page numbers are deliberately NOT given because
    they differ between editions. `verification_status` states that the numbers were
    transcribed from the tables and still need a spot-check against the physical edition.
  * Species whose constants are not reliably tabulated (e.g. H2S pKa2, whose literature range
    spans many pK units) are intentionally left out instead of guessed.
"""

from __future__ import annotations

import re
from typing import Any

# --------------------------------------------------------------------------------------
# Atomic weights (IUPAC / CIAAW standard atomic weights, abridged to 2-3 decimals)
# --------------------------------------------------------------------------------------
ATOMIC_WEIGHTS: dict[str, float] = {
    "H": 1.008, "Li": 6.94, "B": 10.81, "C": 12.011, "N": 14.007, "O": 15.999,
    "F": 18.998, "Na": 22.990, "Mg": 24.305, "Al": 26.982, "P": 30.974, "S": 32.06,
    "Cl": 35.45, "K": 39.098, "Ca": 40.078, "Br": 79.904, "Rb": 85.468, "Sr": 87.62,
    "I": 126.904, "Cs": 132.905, "Ba": 137.327,
}


def molar_mass_from_formula(formula: str) -> float:
    """Molar mass (g/mol) from a plain molecular formula such as 'Ca(OH)2' or 'C6H8O7'."""
    tokens = re.findall(r"([A-Z][a-z]?|\(|\))(\d*)", formula)
    stack: list[dict[str, int]] = [{}]
    for tok, num in tokens:
        n = int(num) if num else 1
        if tok == "(":
            stack.append({})
        elif tok == ")":
            top = stack.pop()
            for el, c in top.items():
                stack[-1][el] = stack[-1].get(el, 0) + c * n
        else:
            if tok not in ATOMIC_WEIGHTS:
                raise ValueError(f"Unknown element '{tok}' in formula '{formula}'")
            stack[-1][tok] = stack[-1].get(tok, 0) + n
    return round(sum(ATOMIC_WEIGHTS[el] * c for el, c in stack[0].items()), 3)


# --------------------------------------------------------------------------------------
# Sources
# --------------------------------------------------------------------------------------
SRC_CRC_INORG = (
    "CRC Handbook of Chemistry and Physics, 97th ed., table 'Dissociation Constants of Inorganic Acids and Bases'"
)
SRC_CRC_ORG = (
    "CRC Handbook of Chemistry and Physics, 97th ed., table 'Dissociation Constants of Organic Acids and Bases'"
)
SRC_HARRIS = "D.C. Harris, Quantitative Chemical Analysis, 9th ed. (2015), Appendix 'Acid Dissociation Constants'"
SRC_WILLIAMS = "R. Williams, 'pKa Data Compiled by R. Williams' (ACS Organic Division compilation)"
SRC_DERIVED_KW = "Derived: pKb = 14.00 - pKa(conjugate acid) at 25 C (Kw = 1.0e-14); pKa from the source named in 'source'"
SRC_INDICATORS_HARRIS = "D.C. Harris, Quantitative Chemical Analysis, 9th ed. (2015), table 'Common Acid-Base Indicators'"
SRC_INDICATORS_CRC = "CRC Handbook of Chemistry and Physics, 97th ed., table 'Indicators' (pH transition ranges)"
SRC_INDICATORS_VOGEL = "Vogel's Textbook of Quantitative Chemical Analysis, 6th ed., chapter on acid-base indicators"
SRC_ATOMIC = "IUPAC/CIAAW standard atomic weights (molar mass computed from the formula)"

VERIFICATION_NOTE_AR = (
    "ط§ظ„ظ‚ظٹظ… ظ…ظ†ظ‚ظˆظ„ط© ظ…ظ† ط§ظ„ط¬ط¯ط§ظˆظ„ ط§ظ„ظ…ط±ط¬ط¹ظٹط© ط§ظ„ظ…ط°ظƒظˆط±ط© ظˆظ…ط­ط³ظˆط¨ط© ط§ظ„ط§طھط³ط§ظ‚ ط¢ظ„ظٹط§ظ‹ط› ظ„ظ… طھظڈط±ط§ط¬ظژط¹ طµظپط­ط© ط¨طµظپط­ط© "
    "ط¹ظ„ظ‰ ط§ظ„ظ†ط³ط®ط© ط§ظ„ظˆط±ظ‚ظٹط©طŒ ظˆظٹظڈظ†طµط­ ط¨ظ…ط·ط§ط¨ظ‚طھظ‡ط§ ط¹ظ†ط¯ ط§ظ„ط§ط³طھط®ط¯ط§ظ… ط§ظ„ط¨ط­ط«ظٹ."
)
VERIFICATION_STATUS = "transcribed_from_standard_tables_not_page_verified"

CATEGORY_AR = {
    "inorganic_strong": "ط­ظ…ط¶/ظ‚ط§ط¹ط¯ط© ط؛ظٹط± ط¹ط¶ظˆظٹط© ظ‚ظˆظٹط©",
    "inorganic_weak": "ط­ظ…ط¶/ظ‚ط§ط¹ط¯ط© ط؛ظٹط± ط¹ط¶ظˆظٹط© ط¶ط¹ظٹظپط©",
    "organic_carboxylic": "ط­ظ…ط¶ ظƒط±ط¨ظˆظƒط³ظٹظ„ظٹ ط¹ط¶ظˆظٹ",
    "organic_other": "ط­ظ…ط¶ ط¹ط¶ظˆظٹ ط¢ط®ط± (ظپظٹظ†ظˆظ„طŒ ط­ظ…ط¶ ظƒظٹطھظˆظ†ظٹ...)",
    "hydroxide_strong": "ظ‡ظٹط¯ط±ظˆظƒط³ظٹط¯ ظپظ„ط²ظٹ ظ‚ظˆظٹ",
    "hydroxide_sparingly_soluble": "ظ‡ظٹط¯ط±ظˆظƒط³ظٹط¯ ط´ط­ظٹط­ ط§ظ„ط°ظˆط¨ط§ظ†",
    "amine_aliphatic": "ط£ظ…ظٹظ† ط£ظ„ظٹظپط§طھظٹ",
    "amine_aromatic_heterocyclic": "ط£ظ…ظٹظ† ط¹ط·ط±ظٹ / ط­ظ„ظ‚ظٹ ط؛ظٹط± ظ…طھط¬ط§ظ†ط³",
    "ammonia_hydrazine_like": "ط£ظ…ظˆظ†ظٹط§ ظˆظ…ط´طھظ‚ط§طھظ‡ط§ ط؛ظٹط± ط§ظ„ط¹ط¶ظˆظٹط©",
    "salt_basic": "ظ…ظ„ط­ ظ‚ط§ط¹ط¯ظٹ (ط£ظ†ظٹظˆظ† ط¶ط¹ظٹظپ ط§ظ„ط­ظ…ط¶ظٹط© ط§ظ„ظ…ط±ط§ظپظ‚ط©)",
    "buffer_base": "ظ‚ط§ط¹ط¯ط© ظ…ظ†ط¸ظ…ط© (Buffer)",
}


def _fmt_pk(vals: list[float]) -> str:
    return ", ".join(f"{v:g}" for v in vals)


def _acid(
    key: str,
    molecular_formula: str,
    name_ar: str,
    name_en: str,
    pkas: list[float],
    category: str,
    solubility_ar: str,
    properties_ar: str,
    source: str,
) -> dict[str, Any]:
    strong = pkas[0] <= -1.0
    n = len(pkas)
    kind_ar = "ظ‚ظˆظٹ" if strong else "ط¶ط¹ظٹظپ"
    proto_ar = {1: "ط£ط­ط§ط¯ظٹ ط§ظ„ط¨ط±ظˆطھظˆظ†", 2: "ط«ظ†ط§ط¦ظٹ ط§ظ„ط¨ط±ظˆطھظˆظ†", 3: "ط«ظ„ط§ط«ظٹ ط§ظ„ط¨ط±ظˆطھظˆظ†"}[n]
    return {
        "formula": key,
        "molecular_formula": molecular_formula,
        "name_ar": name_ar,
        "name_en": name_en,
        "type": "strong" if strong else "weak",
        "category": category,
        "category_ar": CATEGORY_AR[category],
        "strength_ar": f"ط­ظ…ط¶ {kind_ar} {proto_ar} (pKa = {_fmt_pk(pkas)})",
        "pKa": pkas[0],
        "pKa_list": list(pkas),
        "Ka": 10.0 ** (-pkas[0]),
        "Ka_list": [10.0 ** (-p) for p in pkas],
        "protons": n,
        "molar_mass": molar_mass_from_formula(molecular_formula),
        "solubility_ar": solubility_ar,
        "properties_ar": properties_ar,
        "source": source,
        "source_molar_mass": SRC_ATOMIC,
        "verification_status": VERIFICATION_STATUS,
        "verification_note_ar": VERIFICATION_NOTE_AR,
    }


def _base(
    key: str,
    molecular_formula: str,
    name_ar: str,
    name_en: str,
    pkbs: list[float],
    category: str,
    solubility_ar: str,
    properties_ar: str,
    source: str,
    hydroxides: int | None = None,
    strong: bool = False,
) -> dict[str, Any]:
    n = hydroxides if hydroxides is not None else len(pkbs)
    kind_ar = "ظ‚ظˆظٹط©" if strong else "ط¶ط¹ظٹظپط©"
    proto_ar = {1: "ط£ط­ط§ط¯ظٹط© ط§ظ„طھظƒط§ظپط¤", 2: "ط«ظ†ط§ط¦ظٹط© ط§ظ„طھظƒط§ظپط¤", 3: "ط«ظ„ط§ط«ظٹط© ط§ظ„طھظƒط§ظپط¤"}[n]
    return {
        "formula": key,
        "molecular_formula": molecular_formula,
        "name_ar": name_ar,
        "name_en": name_en,
        "type": "strong" if strong else "weak",
        "category": category,
        "category_ar": CATEGORY_AR[category],
        "strength_ar": f"ظ‚ط§ط¹ط¯ط© {kind_ar} {proto_ar} (pKb = {_fmt_pk(pkbs)})",
        "pKb": pkbs[0],
        "pKb_list": list(pkbs),
        "Kb": 10.0 ** (-pkbs[0]),
        "Kb_list": [10.0 ** (-p) for p in pkbs],
        "hydroxides": n,
        "molar_mass": molar_mass_from_formula(molecular_formula),
        "solubility_ar": solubility_ar,
        "properties_ar": properties_ar,
        "source": source,
        "source_molar_mass": SRC_ATOMIC,
        "verification_status": VERIFICATION_STATUS,
        "verification_note_ar": VERIFICATION_NOTE_AR,
    }


def _pkb_from_conj_pka(pka_conjugate: list[float]) -> list[float]:
    """pKb steps of a polybasic species from the pKa values of its conjugate acid chain (descending)."""
    return [round(14.0 - p, 2) for p in pka_conjugate]


# --------------------------------------------------------------------------------------
# Extra acids (pKa from standard tables)
# --------------------------------------------------------------------------------------
EXTRA_ACIDS: dict[str, dict[str, Any]] = {}


def _add_acid(*args: Any, **kwargs: Any) -> None:
    spec = _acid(*args, **kwargs)
    EXTRA_ACIDS[spec["formula"]] = spec


# Strong mineral acids
_add_acid("HBr", "HBr", "ط­ظ…ط¶ ط§ظ„ظ‡ظٹط¯ط±ظˆط¨ط±ظˆظ…ظٹظƒ", "Hydrobromic Acid", [-9.0], "inorganic_strong",
          "ظٹط°ظˆط¨ ظƒظ„ظٹط§ظ‹ ظپظٹ ط§ظ„ظ…ط§ط، (ط­طھظ‰ ~49% w/w)", "ط­ظ…ط¶ ظ…ط¹ط¯ظ†ظٹ ظ‚ظˆظٹ ط¬ط¯ط§ظ‹ط› ط£ظ‚ظˆظ‰ ظ…ظ† HClط› ط£ط¨ط®ط±ط© ظ…ظ‡ظٹظ‘ط¬ط©", SRC_HARRIS)
_add_acid("HI", "HI", "ط­ظ…ط¶ ط§ظ„ظ‡ظٹط¯ط±ظˆظٹظˆط¯ظٹظƒ", "Hydroiodic Acid", [-10.0], "inorganic_strong",
          "ظٹط°ظˆط¨ ظƒظ„ظٹط§ظ‹ ظپظٹ ط§ظ„ظ…ط§ط،", "ط­ظ…ط¶ ظ…ط¹ط¯ظ†ظٹ ظ‚ظˆظٹ ط¬ط¯ط§ظ‹ط› ظٹطھط£ظƒط³ط¯ طھط¯ط±ظٹط¬ظٹط§ظ‹ ط¨ط§ظ„ظ‡ظˆط§ط، ط¥ظ„ظ‰ ط§ظ„ظٹظˆط¯ (ظ„ظˆظ† ط¨ظ†ظٹ)", SRC_HARRIS)
_add_acid("HClO4", "HClO4", "ط­ظ…ط¶ ط§ظ„ط¨ظٹط±ظƒظ„ظˆط±ظٹظƒ", "Perchloric Acid", [-10.0], "inorganic_strong",
          "ظٹظ…طھط²ط¬ ظƒظ„ظٹط§ظ‹ ظ…ط¹ ط§ظ„ظ…ط§ط،", "ظ…ظ† ط£ظ‚ظˆظ‰ ط§ظ„ط£ط­ظ…ط§ط¶ ط§ظ„ظ…ط¹ط¯ظ†ظٹط©ط› ظ…ط¤ظƒط³ط¯ ط®ط·ط± ط¹ظ†ط¯ ط§ظ„طھط±ظƒظٹط² ط§ظ„ط¹ط§ظ„ظٹ ط£ظˆ ط§ظ„ط³ط§ط®ظ†", SRC_HARRIS)
_add_acid("HClO3", "HClO3", "ط­ظ…ط¶ ط§ظ„ظƒظ„ظˆط±ظٹظƒ", "Chloric Acid", [-1.0], "inorganic_strong",
          "ظٹظˆط¬ط¯ ظپظٹ ط§ظ„ظ…ط­ط§ظ„ظٹظ„ ط§ظ„ظ…ط§ط¦ظٹط© ط­طھظ‰ ~40% ظپظ‚ط·", "ط­ظ…ط¶ ظ‚ظˆظٹ ظˆظ…ط¤ظƒط³ط¯ط› ط؛ظٹط± ظ…ط³طھظ‚ط± ط®ط§ط±ط¬ ط§ظ„ظ…ط­ظ„ظˆظ„", SRC_CRC_INORG)
_add_acid("H2SO3", "H2SO3", "ط­ظ…ط¶ ط§ظ„ظƒط¨ط±ظٹطھظˆط²", "Sulfurous Acid", [1.85, 7.20], "inorganic_weak",
          "ظ„ط§ ظٹظˆط¬ط¯ ظ†ظ‚ظٹط§ظ‹ط› ظٹطھظƒظˆظ‘ظ† ط¨ط¥ط°ط§ط¨ط© SO2 ظپظٹ ط§ظ„ظ…ط§ط،", "ط­ظ…ط¶ ط¶ط¹ظٹظپ ط«ظ†ط§ط¦ظٹ ط§ظ„ط¨ط±ظˆطھظˆظ†ط› ط¹ط§ظ…ظ„ ظ…ط®طھط²ظ„", SRC_CRC_INORG)
# Weak inorganic acids
_add_acid("HNO2", "HNO2", "ط­ظ…ط¶ ط§ظ„ظ†ظٹطھط±ظˆط²", "Nitrous Acid", [3.25], "inorganic_weak",
          "ظٹظˆط¬ط¯ ظپظٹ ط§ظ„ظ…ط­ط§ظ„ظٹظ„ ط§ظ„ظ…ط®ظپظپط© ط§ظ„ط¨ط§ط±ط¯ط© ظپظ‚ط·", "ط؛ظٹط± ظ…ط³طھظ‚ط± ظٹطھظپظƒظƒ ط¥ظ„ظ‰ NO ظˆHNO3ط› ظ…ظ„ط­ظ‡ ط§ظ„ظ†طھط±ظٹطھ", SRC_HARRIS)
_add_acid("HClO2", "HClO2", "ط­ظ…ط¶ ط§ظ„ظƒظ„ظˆط±ظˆط²", "Chlorous Acid", [1.94], "inorganic_weak",
          "ظٹظˆط¬ط¯ ظپظٹ ط§ظ„ظ…ط­ظ„ظˆظ„ ظپظ‚ط·", "ط­ظ…ط¶ ط¶ط¹ظٹظپ ط؛ظٹط± ظ…ط³طھظ‚ط± ظˆظ…ط¤ظƒط³ط¯", SRC_CRC_INORG)
_add_acid("HOCl", "HOCl", "ط­ظ…ط¶ ط§ظ„ظ‡ظٹط¨ظˆظƒظ„ظˆط±ظˆط²", "Hypochlorous Acid", [7.53], "inorganic_weak",
          "ظٹظˆط¬ط¯ ظپظٹ ط§ظ„ظ…ط­ظ„ظˆظ„ ط§ظ„ظ…ط§ط¦ظٹ ظپظ‚ط·", "ط­ظ…ط¶ ط¶ط¹ظٹظپ ظ…ط¤ظƒط³ط¯ ظˆظ…ط·ظ‡ط±ط› ط£ط³ط§ط³ ظ…ط¨ظٹط¶ط§طھ ط§ظ„ظƒظ„ظˆط±", SRC_HARRIS)
_add_acid("HOBr", "HOBr", "ط­ظ…ط¶ ط§ظ„ظ‡ظٹط¨ظˆط¨ط±ظˆظ…ظˆط²", "Hypobromous Acid", [8.63], "inorganic_weak",
          "ظٹظˆط¬ط¯ ظپظٹ ط§ظ„ظ…ط­ظ„ظˆظ„ ط§ظ„ظ…ط§ط¦ظٹ ظپظ‚ط·", "ط­ظ…ط¶ ط¶ط¹ظٹظپ ظ…ط¤ظƒط³ط¯ ط؛ظٹط± ظ…ط³طھظ‚ط±", SRC_CRC_INORG)
_add_acid("HCN", "HCN", "ط­ظ…ط¶ ط§ظ„ظ‡ظٹط¯ط±ظˆط³ظٹط§ظ†ظٹظƒ", "Hydrocyanic Acid", [9.21], "inorganic_weak",
          "ظٹظ…طھط²ط¬ ظƒظ„ظٹط§ظ‹ ظ…ط¹ ط§ظ„ظ…ط§ط،", "ط­ظ…ط¶ ط¶ط¹ظٹظپ ط¬ط¯ط§ظ‹ط› ط´ط¯ظٹط¯ ط§ظ„ط³ظ…ظٹط©ط› ظٹطھط·ظ„ط¨ طھظ‡ظˆظٹط© ظˆط¥ط¬ط±ط§ط،ط§طھ ط³ظ„ط§ظ…ط© طµط§ط±ظ…ط©", SRC_HARRIS)
_add_acid("H2NSO3H", "H3NO3S", "ط­ظ…ط¶ ط§ظ„ط³ظ„ظپط§ظ…ظٹظƒ", "Sulfamic Acid", [0.99], "inorganic_weak",
          "ط°ظˆظ‘ط§ط¨ (~213 g/L ط¹ظ†ط¯ 20 آ°C)", "ط­ظ…ط¶ طµظ„ط¨ ظ…ط³طھظ‚ط± ظٹط³طھط®ط¯ظ… ظ…ط¹ظٹط§ط±ظٹط§ظ‹ ط£ظˆظ„ظٹط§ظ‹ ظپظٹ ط§ظ„ظ…ط¹ط§ظٹط±ط© ظˆط¥ط²ط§ظ„ط© ط§ظ„طھط±ط³ط¨ط§طھ", SRC_CRC_INORG)
# Organic carboxylic and other
_add_acid("C3H6O2", "C3H6O2", "ط­ظ…ط¶ ط§ظ„ط¨ط±ظˆط¨ط§ظ†ظˆظٹظƒ (ط§ظ„ط¨ط±ظˆط¨ظٹظˆظ†ظٹظƒ)", "Propanoic Acid", [4.87], "organic_carboxylic",
          "ظٹظ…طھط²ط¬ ظƒظ„ظٹط§ظ‹ ظ…ط¹ ط§ظ„ظ…ط§ط،", "ط­ظ…ط¶ ط¯ظ‡ظ†ظٹ ظ‚طµظٹط± ط§ظ„ط³ظ„ط³ظ„ط©ط› ط±ط§ط¦ط­ط© ظ†ظپط§ط°ط©", SRC_HARRIS)
_add_acid("C4H8O2", "C4H8O2", "ط­ظ…ط¶ ط§ظ„ط¨ظٹظˆطھط§ظ†ظˆظٹظƒ (ط§ظ„ط¨ظٹظˆطھظٹط±ظٹظƒ)", "Butanoic Acid", [4.82], "organic_carboxylic",
          "ظٹظ…طھط²ط¬ ظƒظ„ظٹط§ظ‹ ظ…ط¹ ط§ظ„ظ…ط§ط،", "ط­ظ…ط¶ ط¯ظ‡ظ†ظٹط› ط±ط§ط¦ط­ط© ظƒط±ظٹظ‡ط© ظ‚ظˆظٹط©", SRC_HARRIS)
_add_acid("C2H3ClO2", "C2H3ClO2", "ط­ظ…ط¶ ظƒظ„ظˆط±ظˆ ط§ظ„ط£ط³ظٹطھظٹظƒ", "Chloroacetic Acid", [2.87], "organic_carboxylic",
          "ط´ط¯ظٹط¯ ط§ظ„ط°ظˆط¨ط§ظ† ظپظٹ ط§ظ„ظ…ط§ط،", "ط£ظ‚ظˆظ‰ ظ…ظ† ط­ظ…ط¶ ط§ظ„ط£ط³ظٹطھظٹظƒ ط¨ط³ط¨ط¨ ط§ظ„طھط£ط«ظٹط± ط§ظ„ط­ط«ظٹ ظ„ظ€ Clط› ط³ط§ظ… ظˆط¢ظƒظ„", SRC_HARRIS)
_add_acid("C2H2Cl2O2", "C2H2Cl2O2", "ط­ظ…ط¶ ط«ظ†ط§ط¦ظٹ ظƒظ„ظˆط±ظˆ ط§ظ„ط£ط³ظٹطھظٹظƒ", "Dichloroacetic Acid", [1.35], "organic_carboxylic",
          "ظٹظ…طھط²ط¬ ظƒظ„ظٹط§ظ‹ ظ…ط¹ ط§ظ„ظ…ط§ط،", "ط­ظ…ط¶ ط¹ط¶ظˆظٹ ظ‚ظˆظٹ ظ†ط³ط¨ظٹط§ظ‹ط› ط¢ظƒظ„", SRC_CRC_ORG)
_add_acid("C2H4O3", "C2H4O3", "ط­ظ…ط¶ ط§ظ„ط¬ظ„ظٹظƒظˆظ„ظٹظƒ", "Glycolic Acid", [3.83], "organic_carboxylic",
          "ط´ط¯ظٹط¯ ط§ظ„ط°ظˆط¨ط§ظ† ظپظٹ ط§ظ„ظ…ط§ط،", "ط£ط¨ط³ط· ط­ظ…ط¶ ط£ظ„ظپط§-ظ‡ظٹط¯ط±ظˆظƒط³ظٹ", SRC_HARRIS)
_add_acid("C3H4O3", "C3H4O3", "ط­ظ…ط¶ ط§ظ„ط¨ظٹط±ظˆظپظٹظƒ", "Pyruvic Acid", [2.50], "organic_other",
          "ظٹظ…طھط²ط¬ ظƒظ„ظٹط§ظ‹ ظ…ط¹ ط§ظ„ظ…ط§ط،", "ط­ظ…ط¶ ط£ظ„ظپط§-ظƒظٹطھظˆظ†ظٹط› ظˆط³ظٹط· ط£ظٹط¶ظٹ ط£ط³ط§ط³ظٹ", SRC_HARRIS)
_add_acid("C7H6O3", "C7H6O3", "ط­ظ…ط¶ ط§ظ„ط³ط§ظ„ظٹط³ظٹظ„ظٹظƒ", "Salicylic Acid", [2.97], "organic_carboxylic",
          "ظ‚ظ„ظٹظ„ ط§ظ„ط°ظˆط¨ط§ظ† (~2 g/L ط¹ظ†ط¯ 25 آ°C)", "ط­ظ…ط¶ ظپظٹظ†ظˆظ„ظٹ ظƒط±ط¨ظˆظƒط³ظٹظ„ظٹط› ظ…ط§ط¯ط© ط£ظˆظ„ظٹط© ظ„ظ„ط£ط³ط¨ط±ظٹظ†", SRC_HARRIS)
_add_acid("C9H8O4", "C9H8O4", "ط­ظ…ط¶ ط§ظ„ط£ط³ظٹطھظٹظ„ ط³ط§ظ„ظٹط³ظٹظ„ظٹظƒ (ط§ظ„ط£ط³ط¨ط±ظٹظ†)", "Acetylsalicylic Acid (Aspirin)", [3.5], "organic_carboxylic",
          "ظ‚ظ„ظٹظ„ ط§ظ„ط°ظˆط¨ط§ظ† (~3 g/L ط¹ظ†ط¯ 20 آ°C)", "ظٹطھط­ظ„ظ„ ظ…ط§ط¦ظٹط§ظ‹ ط¨ط¨ط·ط، ط¥ظ„ظ‰ ط­ظ…ط¶ ط§ظ„ط³ط§ظ„ظٹط³ظٹظ„ظٹظƒ ظˆط§ظ„ط£ط³ظٹطھظٹظƒ", SRC_CRC_ORG)
_add_acid("C6H5OH", "C6H6O", "ط§ظ„ظپظٹظ†ظˆظ„", "Phenol", [9.99], "organic_other",
          "ط°ظˆظ‘ط§ط¨ (~83 g/L ط¹ظ†ط¯ 20 آ°C)", "ط­ظ…ط¶ ط¶ط¹ظٹظپ ط¬ط¯ط§ظ‹ط› ط³ط§ظ… ظˆط¢ظƒظ„ ظ„ظ„ط¬ظ„ط¯", SRC_HARRIS)
_add_acid("C3H4O4", "C3H4O4", "ط­ظ…ط¶ ط§ظ„ظ…ط§ظ„ظˆظ†ظٹظƒ", "Malonic Acid", [2.85, 5.70], "organic_carboxylic",
          "ط´ط¯ظٹط¯ ط§ظ„ط°ظˆط¨ط§ظ† (~763 g/L ط¹ظ†ط¯ 20 آ°C)", "ط­ظ…ط¶ ط«ظ†ط§ط¦ظٹ ط§ظ„ظƒط±ط¨ظˆظƒط³ظٹظ„", SRC_HARRIS)
_add_acid("C4H6O4", "C4H6O4", "ط­ظ…ط¶ ط§ظ„ط³ظƒط³ظٹظ†ظٹظƒ (ط§ظ„ظƒظ‡ط±ظ…ط§ظ†)", "Succinic Acid", [4.21, 5.64], "organic_carboxylic",
          "ظ…طھظˆط³ط· ط§ظ„ط°ظˆط¨ط§ظ† (~58 g/L ط¹ظ†ط¯ 20 آ°C)", "ظ‚ظٹظ…طھط§ pKa ظ…طھظ‚ط§ط±ط¨طھط§ظ† ظپظ„ط§ ظٹط¸ظ‡ط± ظ‚ظپط² ظˆط§ط¶ط­ ظ„ظ„طھظƒط§ظپط¤ ط§ظ„ط£ظˆظ„", SRC_HARRIS)
_add_acid("C4H4O4 (maleic)", "C4H4O4", "ط­ظ…ط¶ ط§ظ„ظ…ط§ظ„ط¦ظٹظƒ", "Maleic Acid", [1.91, 6.33], "organic_carboxylic",
          "ط´ط¯ظٹط¯ ط§ظ„ط°ظˆط¨ط§ظ† (~479 g/L ط¹ظ†ط¯ 20 آ°C)", "ط§ظ„ط£ظٹط²ظˆظ…ط± cisط› ظ‚ظٹظ…طھط§ pKa ظ…طھط¨ط§ط¹ط¯طھط§ظ† ط¨ط³ط¨ط¨ ط§ظ„ط±ط§ط¨ط·ط© ط§ظ„ظ‡ظٹط¯ط±ظˆط¬ظٹظ†ظٹط© ط§ظ„ط¯ط§ط®ظ„ظٹط©", SRC_HARRIS)
_add_acid("C4H4O4 (fumaric)", "C4H4O4", "ط­ظ…ط¶ ط§ظ„ظپظˆظ…ط§ط±ظٹظƒ", "Fumaric Acid", [3.05, 4.49], "organic_carboxylic",
          "ظ‚ظ„ظٹظ„ ط§ظ„ط°ظˆط¨ط§ظ† (~4.9 g/L ط¹ظ†ط¯ 25 آ°C)", "ط§ظ„ط£ظٹط²ظˆظ…ط± transط› ظ‚ظٹظ…طھط§ pKa ط£ظ‚ط±ط¨ ظ…ظ† ط§ظ„ظ…ط§ظ„ط¦ظٹظƒ", SRC_HARRIS)
_add_acid("C8H6O4", "C8H6O4", "ط­ظ…ط¶ ط§ظ„ظپط«ط§ظ„ظٹظƒ", "Phthalic Acid", [2.95, 5.41], "organic_carboxylic",
          "ظ‚ظ„ظٹظ„ ط§ظ„ط°ظˆط¨ط§ظ† (~6 g/L ط¹ظ†ط¯ 25 آ°C)", "ظ…ظ„ط­ظ‡ ط§ظ„ط­ظ…ط¶ظٹ (ظپط«ط§ظ„ط§طھ ط§ظ„ط¨ظˆطھط§ط³ظٹظˆظ… ط§ظ„ظ‡ظٹط¯ط±ظˆط¬ظٹظ†ظٹط©) ظ…ط¹ظٹط§ط± ط£ظˆظ„ظٹ ط´ط§ط¦ط¹", SRC_HARRIS)
_add_acid("C4H6O5", "C4H6O5", "ط­ظ…ط¶ ط§ظ„ظ…ط§ظ„ظٹظƒ (ط§ظ„طھظپط§ط­ظٹظƒ)", "Malic Acid", [3.46, 5.10], "organic_carboxylic",
          "ط´ط¯ظٹط¯ ط§ظ„ط°ظˆط¨ط§ظ† (~558 g/L ط¹ظ†ط¯ 20 آ°C)", "ط­ظ…ط¶ ظ‡ظٹط¯ط±ظˆظƒط³ظٹ ط«ظ†ط§ط¦ظٹ ط§ظ„ظƒط±ط¨ظˆظƒط³ظٹظ„ ظپظٹ ط§ظ„طھظپط§ط­", SRC_HARRIS)
_add_acid("C4H6O6", "C4H6O6", "ط­ظ…ط¶ ط§ظ„ط·ط±ط·ط±ظٹظƒ", "Tartaric Acid", [3.04, 4.37], "organic_carboxylic",
          "ط´ط¯ظٹط¯ ط§ظ„ط°ظˆط¨ط§ظ† (~139 g/100 mL ط¹ظ†ط¯ 20 آ°C)", "ط­ظ…ط¶ ط«ظ†ط§ط¦ظٹ ط§ظ„ظ‡ظٹط¯ط±ظˆظƒط³ظٹ ط«ظ†ط§ط¦ظٹ ط§ظ„ظƒط±ط¨ظˆظƒط³ظٹظ„ط› ط«ظ†ط§ط¦ظٹ pKa ظ…طھظ‚ط§ط±ط¨ط§ظ†", SRC_HARRIS)
_add_acid("C6H8O6", "C6H8O6", "ط­ظ…ط¶ ط§ظ„ط£ط³ظƒظˆط±ط¨ظٹظƒ (ظپظٹطھط§ظ…ظٹظ† C)", "Ascorbic Acid (Vitamin C)", [4.10, 11.6], "organic_other",
          "ط´ط¯ظٹط¯ ط§ظ„ط°ظˆط¨ط§ظ† (~330 g/L)", "pKa2 ظ…ط±طھظپط¹ ط¬ط¯ط§ظ‹ ظپظ„ط§ ظٹظڈط¹ط§ظٹط± ط¹ظ…ظ„ظٹط§ظ‹ ط¥ظ„ط§ ط§ظ„ط¨ط±ظˆطھظˆظ† ط§ظ„ط£ظˆظ„ط› ظٹطھط£ظƒط³ط¯ ط¨ط³ظ‡ظˆظ„ط©", SRC_HARRIS)


# --------------------------------------------------------------------------------------
# Extra bases
# --------------------------------------------------------------------------------------
EXTRA_BASES: dict[str, dict[str, Any]] = {}


def _add_base(*args: Any, **kwargs: Any) -> None:
    spec = _base(*args, **kwargs)
    EXTRA_BASES[spec["formula"]] = spec


# Strong hydroxides (complete dissociation in dilute solution)
_add_base("LiOH", "LiOH", "ظ‡ظٹط¯ط±ظˆظƒط³ظٹط¯ ط§ظ„ظ„ظٹط«ظٹظˆظ…", "Lithium Hydroxide", [0.0], "hydroxide_strong",
          "ط°ظˆظ‘ط§ط¨ (~128 g/L ط¹ظ†ط¯ 20 آ°C)", "ط£ط¶ط¹ظپ ظ‡ظٹط¯ط±ظˆظƒط³ظٹط¯ط§طھ ط§ظ„ظ‚ظ„ظˆظٹط§طھ ط°ظˆط¨ط§ظ†ط§ظ‹ط› ط¢ظƒظ„", SRC_CRC_INORG, strong=True)
_add_base("RbOH", "RbOH", "ظ‡ظٹط¯ط±ظˆظƒط³ظٹط¯ ط§ظ„ط±ظˆط¨ظٹط¯ظٹظˆظ…", "Rubidium Hydroxide", [0.0], "hydroxide_strong",
          "ط´ط¯ظٹط¯ ط§ظ„ط°ظˆط¨ط§ظ†", "ظ‚ط§ط¹ط¯ط© ظ‚ظ„ظˆظٹط© ظ‚ظˆظٹط© ط¬ط¯ط§ظ‹ط› ط¢ظƒظ„ ظˆظ…طھظ…ظٹظ‘ط¹", SRC_CRC_INORG, strong=True)
_add_base("CsOH", "CsOH", "ظ‡ظٹط¯ط±ظˆظƒط³ظٹط¯ ط§ظ„ط³ظٹط²ظٹظˆظ…", "Cesium Hydroxide", [0.0], "hydroxide_strong",
          "ط´ط¯ظٹط¯ ط§ظ„ط°ظˆط¨ط§ظ†", "ظ…ظ† ط£ظ‚ظˆظ‰ ط§ظ„ظ‚ظˆط§ط¹ط¯ ط§ظ„ظ‚ظ„ظˆظٹط©ط› ط¢ظƒظ„ ظˆظ…طھظ…ظٹظ‘ط¹", SRC_CRC_INORG, strong=True)
_add_base("Sr(OH)2", "Sr(OH)2", "ظ‡ظٹط¯ط±ظˆظƒط³ظٹط¯ ط§ظ„ط³طھط±ظˆظ†ط´ظٹظˆظ…", "Strontium Hydroxide", [0.0, 0.0], "hydroxide_strong",
          "ظ‚ظ„ظٹظ„ ط¥ظ„ظ‰ ظ…طھظˆط³ط· ط§ظ„ط°ظˆط¨ط§ظ† (ظٹط²ط¯ط§ط¯ ط¨ط§ظ„ط­ط±ط§ط±ط©)", "ظ‚ط§ط¹ط¯ط© ظ‚ظˆظٹط© ط«ظ†ط§ط¦ظٹط© ط§ظ„ظ‡ظٹط¯ط±ظˆظƒط³ظٹط¯ط› ط¢ظƒظ„ط©", SRC_CRC_INORG,
          hydroxides=2, strong=True)
# Ammonia-like and amines (pKb = 14 - pKa(BH+), pKa(BH+) from the named source)
_add_base("N2H4", "N2H4", "ط§ظ„ظ‡ظٹط¯ط±ط§ط²ظٹظ†", "Hydrazine", _pkb_from_conj_pka([8.07]), "ammonia_hydrazine_like",
          "ظٹظ…طھط²ط¬ ظƒظ„ظٹط§ظ‹ ظ…ط¹ ط§ظ„ظ…ط§ط،", "ظ‚ط§ط¹ط¯ط© ط¶ط¹ظٹظپط©ط› ط´ط¯ظٹط¯ ط§ظ„ط³ظ…ظٹط© ظˆظ…ط³ط±ط·ظ† ظ…ط­طھظ…ظ„ط› ط¹ط§ظ…ظ„ ظ…ط®طھط²ظ„", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("NH2OH", "NH3O", "ط§ظ„ظ‡ظٹط¯ط±ظˆظƒط³ظٹظ„ ط£ظ…ظٹظ†", "Hydroxylamine", _pkb_from_conj_pka([5.96]), "ammonia_hydrazine_like",
          "ط°ظˆظ‘ط§ط¨ ظپظٹ ط§ظ„ظ…ط§ط،", "ظ‚ط§ط¹ط¯ط© ط¶ط¹ظٹظپط©ط› ط؛ظٹط± ظ…ط³طھظ‚ط± (ظ‚ط¯ ظٹظ†ظپط¬ط± ط¹ظ†ط¯ ط§ظ„طھط³ط®ظٹظ†)", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("(CH3)2NH", "C2H7N", "ط«ظ†ط§ط¦ظٹ ظ…ظٹط«ظٹظ„ ط£ظ…ظٹظ†", "Dimethylamine", _pkb_from_conj_pka([10.77]), "amine_aliphatic",
          "ط´ط¯ظٹط¯ ط§ظ„ط°ظˆط¨ط§ظ†", "ط£ظ…ظٹظ† ط«ط§ظ†ظˆظٹ ط£ظ„ظٹظپط§طھظٹط› ط؛ط§ط² ط±ط§ط¦ط­طھظ‡ ط³ظ…ظƒظٹط©", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("(CH3)3N", "C3H9N", "ط«ظ„ط§ط«ظٹ ظ…ظٹط«ظٹظ„ ط£ظ…ظٹظ†", "Trimethylamine", _pkb_from_conj_pka([9.80]), "amine_aliphatic",
          "ط´ط¯ظٹط¯ ط§ظ„ط°ظˆط¨ط§ظ†", "ط£ظ…ظٹظ† ط«ط§ظ„ط«ظٹط› ط؛ط§ط² ط±ط§ط¦ط­طھظ‡ ط³ظ…ظƒظٹط© ظ‚ظˆظٹط©", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("C2H4(NH2)2", "C2H8N2", "ط¥ظٹط«ظٹظ„ظٹظ† ط«ظ†ط§ط¦ظٹ ط§ظ„ط£ظ…ظٹظ†", "Ethylenediamine", _pkb_from_conj_pka([9.93, 6.85]),
          "amine_aliphatic", "ظٹظ…طھط²ط¬ ظƒظ„ظٹط§ظ‹ ظ…ط¹ ط§ظ„ظ…ط§ط،", "ظ‚ط§ط¹ط¯ط© ط¶ط¹ظٹظپط© ط«ظ†ط§ط¦ظٹط©ط› ظ„ظٹط¬ط§ظ†ط¯ ط´ط§ط¦ط¹", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("C5H11N", "C5H11N", "ط§ظ„ط¨ظٹط¨ظٹط±ظٹط¯ظٹظ†", "Piperidine", _pkb_from_conj_pka([11.12]), "amine_aliphatic",
          "ظٹظ…طھط²ط¬ ظƒظ„ظٹط§ظ‹ ظ…ط¹ ط§ظ„ظ…ط§ط،", "ط£ظ…ظٹظ† ط­ظ„ظ‚ظٹ ط«ط§ظ†ظˆظٹط› ط£ظ‚ظˆظ‰ ظ‚ط§ط¹ط¯ظٹط© ظ…ظ† ط§ظ„ط£ظ…ظˆظ†ظٹط§ ط¨ظƒط«ظٹط±", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("C2H7NO", "C2H7NO", "ط¥ظٹط«ط§ظ†ظˆظ„ ط£ظ…ظٹظ†", "Ethanolamine", _pkb_from_conj_pka([9.50]), "amine_aliphatic",
          "ظٹظ…طھط²ط¬ ظƒظ„ظٹط§ظ‹ ظ…ط¹ ط§ظ„ظ…ط§ط،", "ط£ظ…ظٹظ† ط£ظˆظ„ظٹ ظƒط­ظˆظ„ظٹط› ظٹط³طھط¹ظ…ظ„ ظپظٹ ط§ظ…طھطµط§طµ CO2", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("C5H5N", "C5H5N", "ط§ظ„ط¨ظٹط±ظٹط¯ظٹظ†", "Pyridine", _pkb_from_conj_pka([5.20]), "amine_aromatic_heterocyclic",
          "ظٹظ…طھط²ط¬ ظƒظ„ظٹط§ظ‹ ظ…ط¹ ط§ظ„ظ…ط§ط،", "ظ‚ط§ط¹ط¯ط© ط¶ط¹ظٹظپط© ط­ظ„ظ‚ظٹط© ط¹ط·ط±ظٹط©ط› ظ…ط°ظٹط¨ ط³ط§ظ… ظˆط±ط§ط¦ط­طھظ‡ ظ†ظپط§ط°ط©", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("C6H5NH2", "C6H7N", "ط§ظ„ط£ظ†ظٹظ„ظٹظ†", "Aniline", _pkb_from_conj_pka([4.60]), "amine_aromatic_heterocyclic",
          "ظ‚ظ„ظٹظ„ ط§ظ„ط°ظˆط¨ط§ظ† (~36 g/L ط¹ظ†ط¯ 20 آ°C)", "ظ‚ط§ط¹ط¯ط© ط¶ط¹ظٹظپط© ط¬ط¯ط§ظ‹ ظ„ط£ظ† ط²ظˆط¬ ط§ظ„ط¥ظ„ظƒطھط±ظˆظ†ط§طھ ظ…طھظ…ط±ظƒط² ظپظٹ ط§ظ„ط­ظ„ظ‚ط©ط› ط³ط§ظ…", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("(HOCH2)3CNH2", "C4H11NO3", "طھط±ظٹط³ (Tris)", "Tris(hydroxymethyl)aminomethane", _pkb_from_conj_pka([8.07]), "buffer_base",
          "ط´ط¯ظٹط¯ ط§ظ„ط°ظˆط¨ط§ظ† (~550 g/L ط¹ظ†ط¯ 25 آ°C)", "ظ‚ط§ط¹ط¯ط© ظ…ظ†ط¸ظ…ط© ط´ط§ط¦ط¹ط© ظپظٹ ط§ظ„ظƒظٹظ…ظٹط§ط، ط§ظ„ط­ظٹظˆظٹط©ط› ظ…ط¹ظٹط§ط± ط£ظˆظ„ظٹ ظ„ظ„ط£ط­ظ…ط§ط¶", SRC_HARRIS + "; " + SRC_DERIVED_KW)
# Basic salts: pKb from the pKa chain of the conjugate acid
_add_base("K2CO3", "K2CO3", "ظƒط±ط¨ظˆظ†ط§طھ ط§ظ„ط¨ظˆطھط§ط³ظٹظˆظ…", "Potassium Carbonate", _pkb_from_conj_pka([10.33, 6.35]), "salt_basic",
          "ط´ط¯ظٹط¯ ط§ظ„ط°ظˆط¨ط§ظ† (~1120 g/L ط¹ظ†ط¯ 20 آ°C)", "ظ‚ط§ط¹ط¯ط© ط¶ط¹ظٹظپط© ط«ظ†ط§ط¦ظٹط©ط› طھظ…طھطµ ط§ظ„ط±ط·ظˆط¨ط©", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("KHCO3", "KHCO3", "ط¨ظٹظƒط±ط¨ظˆظ†ط§طھ ط§ظ„ط¨ظˆطھط§ط³ظٹظˆظ…", "Potassium Bicarbonate", _pkb_from_conj_pka([6.35]), "salt_basic",
          "ط°ظˆظ‘ط§ط¨ (~224 g/L ط¹ظ†ط¯ 20 آ°C)", "ظ‚ط§ط¹ط¯ط© ط¶ط¹ظٹظپط© ط¬ط¯ط§ظ‹ ط£ظ…ظپظˆطھظٹط±ظٹط©", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("Na3PO4", "Na3PO4", "ظپظˆط³ظپط§طھ ط§ظ„طµظˆط¯ظٹظˆظ… ط§ظ„ط«ظ„ط§ط«ظٹ", "Trisodium Phosphate", _pkb_from_conj_pka([12.35, 7.20, 2.15]),
          "salt_basic", "ط°ظˆظ‘ط§ط¨ ظپظٹ ط§ظ„ظ…ط§ط،", "ظ‚ط§ط¹ط¯ط© ط«ظ„ط§ط«ظٹط©ط› pKb1 طµط؛ظٹط± ظپطھط³ظ„ظƒ ظƒظ‚ط§ط¹ط¯ط© ظ‚ظˆظٹط© ظ†ط³ط¨ظٹط§ظ‹ ظپظٹ ط§ظ„ط®ط·ظˆط© ط§ظ„ط£ظˆظ„ظ‰", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("CH3COONa", "C2H3NaO2", "ط£ط³ظٹطھط§طھ ط§ظ„طµظˆط¯ظٹظˆظ…", "Sodium Acetate", _pkb_from_conj_pka([4.756]), "salt_basic",
          "ط´ط¯ظٹط¯ ط§ظ„ط°ظˆط¨ط§ظ† (~464 g/L ط¹ظ†ط¯ 20 آ°C)", "ظ‚ط§ط¹ط¯ط© ط¶ط¹ظٹظپط© ط¬ط¯ط§ظ‹ط› طھط³طھط¹ظ…ظ„ ظ…ط¹ ط­ظ…ط¶ ط§ظ„ط£ط³ظٹطھظٹظƒ ظƒظ…ط­ظ„ظˆظ„ ظ…ظ†ط¸ظ…", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("NaCN", "CNNa", "ط³ظٹط§ظ†ظٹط¯ ط§ظ„طµظˆط¯ظٹظˆظ…", "Sodium Cyanide", _pkb_from_conj_pka([9.21]), "salt_basic",
          "ط´ط¯ظٹط¯ ط§ظ„ط°ظˆط¨ط§ظ†", "ظ‚ط§ط¹ط¯ط© ط¶ط¹ظٹظپط©ط› ط³ط§ظ… ط¬ط¯ط§ظ‹ط› ظٹط­ط±ظ‘ط± HCN ط§ظ„ظ‚ط§طھظ„ ظ…ط¹ ط§ظ„ط£ط­ظ…ط§ط¶", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("NaF", "FNa", "ظپظ„ظˆط±ظٹط¯ ط§ظ„طµظˆط¯ظٹظˆظ…", "Sodium Fluoride", _pkb_from_conj_pka([3.17]), "salt_basic",
          "ط°ظˆظ‘ط§ط¨ (~40 g/L ط¹ظ†ط¯ 20 آ°C)", "ظ‚ط§ط¹ط¯ط© ط¶ط¹ظٹظپط© ط¬ط¯ط§ظ‹ط› ط³ط§ظ… ط¹ظ†ط¯ ط§ظ„ط¬ط±ط¹ط§طھ ط§ظ„ط¹ط§ظ„ظٹط©", SRC_HARRIS + "; " + SRC_DERIVED_KW)


# --------------------------------------------------------------------------------------
# Extra indicators (acid-form -> base-form transition ranges)
# --------------------------------------------------------------------------------------
EXTRA_INDICATORS: dict[str, dict[str, Any]] = {}


def _add_indicator(
    key: str,
    name_ar: str,
    name_en: str,
    ph_range: tuple[float, float],
    color_acid: str,
    color_base: str,
    hex_acid: str,
    hex_mid: str,
    hex_base: str,
    suitability_ar: str,
    source: str,
) -> None:
    EXTRA_INDICATORS[key] = {
        "key": key,
        "name_ar": name_ar,
        "name_en": name_en,
        "pH_range": ph_range,
        "color_acid": color_acid,
        "color_base": color_base,
        "color_acid_hex": hex_acid,
        "color_transition_hex": hex_mid,
        "color_base_hex": hex_base,
        "color_hex": hex_mid,
        "suitability_ar": suitability_ar,
        "source": source,
        "verification_status": VERIFICATION_STATUS,
        "verification_note_ar": VERIFICATION_NOTE_AR,
    }


_add_indicator("methyl_violet", "ط¨ظ†ظپط³ط¬ظٹ ط§ظ„ظ…ظٹط«ظٹظ„ (Methyl Violet)", "Methyl Violet", (0.0, 1.6),
               "ط£طµظپط± (Yellow)", "ط£ط²ط±ظ‚ ط¨ظ†ظپط³ط¬ظٹ (Blue-Violet)", "#FACC15", "#22C55E", "#4338CA",
               "ظ…ط¯ظ‰ ط´ط¯ظٹط¯ ط§ظ„ط­ظ…ظˆط¶ط©ط› ظٹظ„ط§ط¦ظ… ط§ظ„ط£ط­ظ…ط§ط¶ ط§ظ„ظ‚ظˆظٹط© ط§ظ„ظ…ط±ظƒط²ط© ط¬ط¯ط§ظ‹.", SRC_INDICATORS_CRC)
_add_indicator("thymol_blue_acid", "ط£ط²ط±ظ‚ ط§ظ„ط«ظٹظ…ظˆظ„ - ط§ظ„ظ…ط¯ظ‰ ط§ظ„ط­ظ…ط¶ظٹ (Thymol Blue, acidic range)", "Thymol Blue (acid range)",
               (1.2, 2.8), "ط£ط­ظ…ط± (Red)", "ط£طµظپط± (Yellow)", "#DC2626", "#FB923C", "#FACC15",
               "ط§ظ„ظ…ط¯ظ‰ ط§ظ„ط­ظ…ط¶ظٹ ط§ظ„ط£ظˆظ„ ظ„ط£ط²ط±ظ‚ ط§ظ„ط«ظٹظ…ظˆظ„ط› ظٹظ„ط§ط¦ظ… ظ†ظ‚ط§ط· طھظƒط§ظپط¤ ط­ظ…ط¶ظٹط© ط¬ط¯ط§ظ‹.", SRC_INDICATORS_HARRIS)
_add_indicator("methyl_yellow", "ط§ظ„ط£طµظپط± ط§ظ„ظ…ظٹط«ظٹظ„ظٹ (Methyl Yellow)", "Methyl Yellow", (2.9, 4.0),
               "ط£ط­ظ…ط± (Red)", "ط£طµظپط± (Yellow)", "#DC2626", "#F97316", "#FACC15",
               "ظ…ط¯ظ‰ ظ‚ط±ظٹط¨ ظ…ظ† ط§ظ„ظ…ظٹط«ظٹظ„ ط§ظ„ط¨ط±طھظ‚ط§ظ„ظٹط› ظٹظ„ط§ط¦ظ… ط§ظ„ظ‚ظˆط§ط¹ط¯ ط§ظ„ط¶ط¹ظٹظپط© ظ…ط¹ ط§ظ„ط£ط­ظ…ط§ط¶ ط§ظ„ظ‚ظˆظٹط©.", SRC_INDICATORS_HARRIS)
_add_indicator("dinitrophenol_24", "2,4-ط«ظ†ط§ط¦ظٹ ظ†ظٹطھط±ظˆظپظٹظ†ظˆظ„ (2,4-Dinitrophenol)", "2,4-Dinitrophenol", (2.8, 4.7),
               "ط¹ط¯ظٹظ… ط§ظ„ظ„ظˆظ† (Colorless)", "ط£طµظپط± (Yellow)", "#F8FAFC", "#FEF08A", "#FACC15",
               "ظ…ط¯ظ‰ ط­ظ…ط¶ظٹ ط¶ط¹ظٹظپط› طھط؛ظٹظ‘ط± ظ„ظˆظ†ظٹ ط¶ط¹ظٹظپ ط§ظ„ظˆط¶ظˆط­ ظ…ظ‚ط§ط±ظ†ط© ط¨ط§ظ„ظƒظˆط§ط´ظپ ط§ظ„ظ…ظ„ظˆظ‘ظ†ط© ط¨ط§ظ„ظƒط§ظ…ظ„.", SRC_INDICATORS_CRC)
_add_indicator("bromophenol_blue", "ط£ط²ط±ظ‚ ط§ظ„ط¨ط±ظˆظ…ظˆظپظٹظ†ظˆظ„ (Bromophenol Blue)", "Bromophenol Blue", (3.0, 4.6),
               "ط£طµظپط± (Yellow)", "ط£ط²ط±ظ‚ ط¨ظ†ظپط³ط¬ظٹ (Blue-Violet)", "#FACC15", "#4ADE80", "#4F46E5",
               "ظ…ظ„ط§ط¦ظ… ظ„ظ…ط¹ط§ظٹط±ط© ط§ظ„ظ‚ظˆط§ط¹ط¯ ط§ظ„ط¶ط¹ظٹظپط© ظ…ط¹ ط§ظ„ط£ط­ظ…ط§ط¶ ط§ظ„ظ‚ظˆظٹط© ط¹ظ†ط¯ ظ†ظ‚ط·ط© طھظƒط§ظپط¤ ط­ظ…ط¶ظٹط©.", SRC_INDICATORS_HARRIS)
_add_indicator("congo_red", "ط£ط­ظ…ط± ط§ظ„ظƒظˆظ†ط؛ظˆ (Congo Red)", "Congo Red", (3.0, 5.0),
               "ط£ط²ط±ظ‚ ط¨ظ†ظپط³ط¬ظٹ (Blue-Violet)", "ط£ط­ظ…ط± (Red)", "#4338CA", "#A855F7", "#DC2626",
               "ظ…ط¯ظ‰ ط­ظ…ط¶ظٹ ظˆط§ط³ط¹ط› ط§ط³طھط¹ظ…ط§ظ„ طھط§ط±ظٹط®ظٹ ظƒط§ط´ظپ ظ„ظ„ط£ط­ظ…ط§ط¶ ط§ظ„ظ…ط¹ط¯ظ†ظٹط©.", SRC_INDICATORS_CRC)
_add_indicator("chlorophenol_red", "ط£ط­ظ…ط± ط§ظ„ظƒظ„ظˆط±ظˆظپظٹظ†ظˆظ„ (Chlorophenol Red)", "Chlorophenol Red", (4.8, 6.4),
               "ط£طµظپط± (Yellow)", "ط£ط­ظ…ط± (Red)", "#FACC15", "#FB923C", "#DC2626",
               "ظ…ظ„ط§ط¦ظ… ظ„ظ†ظ‚ط§ط· طھظƒط§ظپط¤ ط­ظ…ط¶ظٹط© ط®ظپظٹظپط© (pH 5-6).", SRC_INDICATORS_HARRIS)
_add_indicator("neutral_red", "ط§ظ„ط£ط­ظ…ط± ط§ظ„ظ…طھط¹ط§ط¯ظ„ (Neutral Red)", "Neutral Red", (6.8, 8.0),
               "ط£ط­ظ…ط± (Red)", "ط£طµظپط± (Yellow)", "#DC2626", "#F97316", "#FACC15",
               "ظ…ط¯ظ‰ ط­ظˆظ„ ط§ظ„طھط¹ط§ط¯ظ„ط› ظٹظ„ط§ط¦ظ… ظ†ظ‚ط·ط© طھظƒط§ظپط¤ ظ‚ط±ظٹط¨ط© ظ…ظ† 7.", SRC_INDICATORS_HARRIS)
_add_indicator("brilliant_yellow", "ط§ظ„ط£طµظپط± ط§ظ„ظ„ط§ظ…ط¹ (Brilliant Yellow)", "Brilliant Yellow", (6.6, 7.8),
               "ط£طµظپط± (Yellow)", "ط£ط­ظ…ط± (Red)", "#FACC15", "#FB923C", "#DC2626",
               "ظ…ط¯ظ‰ ط­ظˆظ„ ط§ظ„طھط¹ط§ط¯ظ„ط› ط¨ط¯ظٹظ„ ظ„ط£ط²ط±ظ‚ ط§ظ„ط¨ط±ظˆظ…ظˆط«ظٹظ…ظˆظ„.", SRC_INDICATORS_HARRIS)
_add_indicator("cresol_red", "ط£ط­ظ…ط± ط§ظ„ظƒط±ظٹط²ظˆظ„ (Cresol Red)", "Cresol Red", (7.0, 8.8),
               "ط£طµظپط± (Yellow)", "ط£ط­ظ…ط± (Red)", "#FACC15", "#FB923C", "#DC2626",
               "ظ…ظ„ط§ط¦ظ… ظ„ظ„ط£ط­ظ…ط§ط¶ ط§ظ„ط¶ط¹ظٹظپط© ظ…ط¹ ط§ظ„ظ‚ظˆط§ط¹ط¯ ط§ظ„ظ‚ظˆظٹط© ط¹ظ†ط¯ ظ†ظ‚ط·ط© طھظƒط§ظپط¤ ظ‚ط§ط¹ط¯ظٹط© ط®ظپظٹظپط©.", SRC_INDICATORS_HARRIS)
_add_indicator("m_cresol_purple", "ط£ط±ط¬ظˆط§ظ†ظٹ ط§ظ„ظ…ظٹطھط§-ظƒط±ظٹط²ظˆظ„ (m-Cresol Purple)", "m-Cresol Purple", (7.4, 9.0),
               "ط£طµظپط± (Yellow)", "ط£ط±ط¬ظˆط§ظ†ظٹ (Purple)", "#FACC15", "#F472B6", "#7E22CE",
               "ظ…ظ„ط§ط¦ظ… ظ„ظ†ظ‚ط§ط· طھظƒط§ظپط¤ ظ‚ط§ط¹ط¯ظٹط© ط®ظپظٹظپط© ط¥ظ„ظ‰ ظ…طھظˆط³ط·ط©.", SRC_INDICATORS_CRC)
_add_indicator("thymolphthalein", "ط§ظ„ط«ظٹظ…ظˆظ„ظپط«ط§ظ„ظٹظ† (Thymolphthalein)", "Thymolphthalein", (9.3, 10.5),
               "ط¹ط¯ظٹظ… ط§ظ„ظ„ظˆظ† (Colorless)", "ط£ط²ط±ظ‚ (Blue)", "#F8FAFC", "#BFDBFE", "#2563EB",
               "ظ…ظ„ط§ط¦ظ… ظ„ظ…ط¹ط§ظٹط±ط© ط§ظ„ط£ط­ظ…ط§ط¶ ط§ظ„ط¶ط¹ظٹظپط© ط¬ط¯ط§ظ‹ ظ…ط¹ ط§ظ„ظ‚ظˆط§ط¹ط¯ ط§ظ„ظ‚ظˆظٹط© (طھظƒط§ظپط¤ pH ~9.5-10).", SRC_INDICATORS_HARRIS)
_add_indicator("alizarin_yellow_r", "ط£طµظپط± ط§ظ„ط£ظ„ظٹط²ط§ط±ظٹظ† R (Alizarin Yellow R)", "Alizarin Yellow R", (10.1, 12.0),
               "ط£طµظپط± (Yellow)", "ط£ط­ظ…ط± (Red)", "#FACC15", "#FB923C", "#DC2626",
               "ظ…ط¯ظ‰ ظ‚ط§ط¹ط¯ظٹ ظ…ط±طھظپط¹ط› ظٹظ„ط§ط¦ظ… ط£ط­ظ…ط§ط¶ط§ظ‹ ط¶ط¹ظٹظپط© ط¬ط¯ط§ظ‹ ط£ظˆ ظ…ط¹ط§ظٹط±ط§طھ ظ‚ط§ط¹ط¯ظٹط© ط´ط¯ظٹط¯ط©.", SRC_INDICATORS_HARRIS)
_add_indicator("indigo_carmine", "ظƒط§ط±ظ…ظٹظ† ط§ظ„ظ†ظٹظ„ط© (Indigo Carmine)", "Indigo Carmine", (11.4, 13.0),
               "ط£ط²ط±ظ‚ (Blue)", "ط£طµظپط± (Yellow)", "#2563EB", "#65A30D", "#FACC15",
               "ظ…ط¯ظ‰ ظ‚ط§ط¹ط¯ظٹ ط´ط¯ظٹط¯ط› ط§ط³طھط¹ظ…ط§ظ„ ظ…ط­ط¯ظˆط¯.", SRC_INDICATORS_CRC)


# --------------------------------------------------------------------------------------
# Provenance for the original (legacy) entries so that every species carries a source
# --------------------------------------------------------------------------------------
LEGACY_ACID_CATEGORY: dict[str, tuple[str, str]] = {
    "HCl": ("inorganic_strong", SRC_HARRIS),
    "H2SO4": ("inorganic_strong", SRC_HARRIS),
    "HNO3": ("inorganic_strong", SRC_HARRIS),
    "H3PO4": ("inorganic_weak", SRC_HARRIS),
    "CH3COOH": ("organic_carboxylic", SRC_HARRIS),
    "HCOOH": ("organic_carboxylic", SRC_HARRIS),
    "C6H8O7": ("organic_carboxylic", SRC_HARRIS),
    "H2C2O4": ("organic_carboxylic", SRC_HARRIS),
    "C3H6O3": ("organic_carboxylic", SRC_HARRIS),
    "C7H6O2": ("organic_carboxylic", SRC_HARRIS),
    "H2CO3": ("inorganic_weak", SRC_HARRIS),
    "H3BO3": ("inorganic_weak", SRC_HARRIS),
    "HF": ("inorganic_weak", SRC_HARRIS),
}
LEGACY_BASE_CATEGORY: dict[str, tuple[str, str]] = {
    "NaOH": ("hydroxide_strong", SRC_CRC_INORG),
    "KOH": ("hydroxide_strong", SRC_CRC_INORG),
    "Ca(OH)2": ("hydroxide_strong", SRC_CRC_INORG),
    "Ba(OH)2": ("hydroxide_strong", SRC_CRC_INORG),
    "NH3": ("ammonia_hydrazine_like", SRC_HARRIS),
    "Mg(OH)2": ("hydroxide_sparingly_soluble", SRC_CRC_INORG),
    "CH3NH2": ("amine_aliphatic", SRC_HARRIS),
    "C2H5NH2": ("amine_aliphatic", SRC_HARRIS),
    "NaHCO3": ("salt_basic", SRC_HARRIS + "; " + SRC_DERIVED_KW),
    "Na2CO3": ("salt_basic", SRC_HARRIS + "; " + SRC_DERIVED_KW),
}
LEGACY_INDICATOR_SOURCE = SRC_INDICATORS_HARRIS + "; " + SRC_INDICATORS_CRC + "; " + SRC_INDICATORS_VOGEL


def finalize_species_tables(
    acids: dict[str, dict[str, Any]],
    bases: dict[str, dict[str, Any]],
    indicators: dict[str, dict[str, Any]],
) -> None:
    """Merge the extra species into the engine tables and attach provenance to legacy entries."""
    for key, (cat, src) in LEGACY_ACID_CATEGORY.items():
        spec = acids.get(key)
        if spec is None:
            continue
        spec.setdefault("molecular_formula", spec["formula"])
        spec.setdefault("category", cat)
        spec.setdefault("category_ar", CATEGORY_AR[cat])
        spec.setdefault("source", src)
        spec.setdefault("properties_ar", spec.get("strength_ar", ""))
        spec.setdefault("verification_status", VERIFICATION_STATUS)
        spec.setdefault("verification_note_ar", VERIFICATION_NOTE_AR)
        spec.setdefault("Ka_list", [10.0 ** (-p) for p in spec["pKa_list"]])
    for key, (cat, src) in LEGACY_BASE_CATEGORY.items():
        spec = bases.get(key)
        if spec is None:
            continue
        spec.setdefault("molecular_formula", spec["formula"])
        spec.setdefault("category", cat)
        spec.setdefault("category_ar", CATEGORY_AR[cat])
        spec.setdefault("source", src)
        spec.setdefault("properties_ar", spec.get("strength_ar", ""))
        spec.setdefault("verification_status", VERIFICATION_STATUS)
        spec.setdefault("verification_note_ar", VERIFICATION_NOTE_AR)
    for key, spec in indicators.items():
        spec.setdefault("source", LEGACY_INDICATOR_SOURCE)
        spec.setdefault("verification_status", VERIFICATION_STATUS)
        spec.setdefault("verification_note_ar", VERIFICATION_NOTE_AR)

    for key, spec in EXTRA_ACIDS.items():
        acids.setdefault(key, spec)
    for key, spec in EXTRA_BASES.items():
        bases.setdefault(key, spec)
    for key, spec in EXTRA_INDICATORS.items():
        indicators.setdefault(key, spec)

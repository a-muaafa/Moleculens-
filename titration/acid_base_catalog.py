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
    "القيم منقولة من الجداول المرجعية المذكورة ومحسوبة الاتساق آلياً؛ لم تُراجَع صفحة بصفحة "
    "على النسخة الورقية، ويُنصح بمطابقتها عند الاستخدام البحثي."
)
VERIFICATION_STATUS = "transcribed_from_standard_tables_not_page_verified"

CATEGORY_AR = {
    "inorganic_strong": "حمض/قاعدة غير عضوية قوية",
    "inorganic_weak": "حمض/قاعدة غير عضوية ضعيفة",
    "organic_carboxylic": "حمض كربوكسيلي عضوي",
    "organic_other": "حمض عضوي آخر (فينول، حمض كيتوني...)",
    "hydroxide_strong": "هيدروكسيد فلزي قوي",
    "hydroxide_sparingly_soluble": "هيدروكسيد شحيح الذوبان",
    "amine_aliphatic": "أمين أليفاتي",
    "amine_aromatic_heterocyclic": "أمين عطري / حلقي غير متجانس",
    "ammonia_hydrazine_like": "أمونيا ومشتقاتها غير العضوية",
    "salt_basic": "ملح قاعدي (أنيون ضعيف الحمضية المرافقة)",
    "buffer_base": "قاعدة منظمة (Buffer)",
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
    kind_ar = "قوي" if strong else "ضعيف"
    proto_ar = {1: "أحادي البروتون", 2: "ثنائي البروتون", 3: "ثلاثي البروتون"}[n]
    return {
        "formula": key,
        "molecular_formula": molecular_formula,
        "name_ar": name_ar,
        "name_en": name_en,
        "type": "strong" if strong else "weak",
        "category": category,
        "category_ar": CATEGORY_AR[category],
        "strength_ar": f"حمض {kind_ar} {proto_ar} (pKa = {_fmt_pk(pkas)})",
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
    kind_ar = "قوية" if strong else "ضعيفة"
    proto_ar = {1: "أحادية التكافؤ", 2: "ثنائية التكافؤ", 3: "ثلاثية التكافؤ"}[n]
    return {
        "formula": key,
        "molecular_formula": molecular_formula,
        "name_ar": name_ar,
        "name_en": name_en,
        "type": "strong" if strong else "weak",
        "category": category,
        "category_ar": CATEGORY_AR[category],
        "strength_ar": f"قاعدة {kind_ar} {proto_ar} (pKb = {_fmt_pk(pkbs)})",
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
_add_acid("HBr", "HBr", "حمض الهيدروبروميك", "Hydrobromic Acid", [-9.0], "inorganic_strong",
          "يذوب كلياً في الماء (حتى ~49% w/w)", "حمض معدني قوي جداً؛ أقوى من HCl؛ أبخرة مهيّجة", SRC_HARRIS)
_add_acid("HI", "HI", "حمض الهيدرويوديك", "Hydroiodic Acid", [-10.0], "inorganic_strong",
          "يذوب كلياً في الماء", "حمض معدني قوي جداً؛ يتأكسد تدريجياً بالهواء إلى اليود (لون بني)", SRC_HARRIS)
_add_acid("HClO4", "HClO4", "حمض البيركلوريك", "Perchloric Acid", [-10.0], "inorganic_strong",
          "يمتزج كلياً مع الماء", "من أقوى الأحماض المعدنية؛ مؤكسد خطر عند التركيز العالي أو الساخن", SRC_HARRIS)
_add_acid("HClO3", "HClO3", "حمض الكلوريك", "Chloric Acid", [-1.0], "inorganic_strong",
          "يوجد في المحاليل المائية حتى ~40% فقط", "حمض قوي ومؤكسد؛ غير مستقر خارج المحلول", SRC_CRC_INORG)
_add_acid("H2SO3", "H2SO3", "حمض الكبريتوز", "Sulfurous Acid", [1.85, 7.20], "inorganic_weak",
          "لا يوجد نقياً؛ يتكوّن بإذابة SO2 في الماء", "حمض ضعيف ثنائي البروتون؛ عامل مختزل", SRC_CRC_INORG)
# Weak inorganic acids
_add_acid("HNO2", "HNO2", "حمض النيتروز", "Nitrous Acid", [3.25], "inorganic_weak",
          "يوجد في المحاليل المخففة الباردة فقط", "غير مستقر يتفكك إلى NO وHNO3؛ ملحه النتريت", SRC_HARRIS)
_add_acid("HClO2", "HClO2", "حمض الكلوروز", "Chlorous Acid", [1.94], "inorganic_weak",
          "يوجد في المحلول فقط", "حمض ضعيف غير مستقر ومؤكسد", SRC_CRC_INORG)
_add_acid("HOCl", "HOCl", "حمض الهيبوكلوروز", "Hypochlorous Acid", [7.53], "inorganic_weak",
          "يوجد في المحلول المائي فقط", "حمض ضعيف مؤكسد ومطهر؛ أساس مبيضات الكلور", SRC_HARRIS)
_add_acid("HOBr", "HOBr", "حمض الهيبوبروموز", "Hypobromous Acid", [8.63], "inorganic_weak",
          "يوجد في المحلول المائي فقط", "حمض ضعيف مؤكسد غير مستقر", SRC_CRC_INORG)
_add_acid("HCN", "HCN", "حمض الهيدروسيانيك", "Hydrocyanic Acid", [9.21], "inorganic_weak",
          "يمتزج كلياً مع الماء", "حمض ضعيف جداً؛ شديد السمية؛ يتطلب تهوية وإجراءات سلامة صارمة", SRC_HARRIS)
_add_acid("H2NSO3H", "H3NO3S", "حمض السلفاميك", "Sulfamic Acid", [0.99], "inorganic_weak",
          "ذوّاب (~213 g/L عند 20 °C)", "حمض صلب مستقر يستخدم معيارياً أولياً في المعايرة وإزالة الترسبات", SRC_CRC_INORG)
# Organic carboxylic and other
_add_acid("C3H6O2", "C3H6O2", "حمض البروبانويك (البروبيونيك)", "Propanoic Acid", [4.87], "organic_carboxylic",
          "يمتزج كلياً مع الماء", "حمض دهني قصير السلسلة؛ رائحة نفاذة", SRC_HARRIS)
_add_acid("C4H8O2", "C4H8O2", "حمض البيوتانويك (البيوتيريك)", "Butanoic Acid", [4.82], "organic_carboxylic",
          "يمتزج كلياً مع الماء", "حمض دهني؛ رائحة كريهة قوية", SRC_HARRIS)
_add_acid("C2H3ClO2", "C2H3ClO2", "حمض كلورو الأسيتيك", "Chloroacetic Acid", [2.87], "organic_carboxylic",
          "شديد الذوبان في الماء", "أقوى من حمض الأسيتيك بسبب التأثير الحثي لـ Cl؛ سام وآكل", SRC_HARRIS)
_add_acid("C2H2Cl2O2", "C2H2Cl2O2", "حمض ثنائي كلورو الأسيتيك", "Dichloroacetic Acid", [1.35], "organic_carboxylic",
          "يمتزج كلياً مع الماء", "حمض عضوي قوي نسبياً؛ آكل", SRC_CRC_ORG)
_add_acid("C2H4O3", "C2H4O3", "حمض الجليكوليك", "Glycolic Acid", [3.83], "organic_carboxylic",
          "شديد الذوبان في الماء", "أبسط حمض ألفا-هيدروكسي", SRC_HARRIS)
_add_acid("C3H4O3", "C3H4O3", "حمض البيروفيك", "Pyruvic Acid", [2.50], "organic_other",
          "يمتزج كلياً مع الماء", "حمض ألفا-كيتوني؛ وسيط أيضي أساسي", SRC_HARRIS)
_add_acid("C7H6O3", "C7H6O3", "حمض الساليسيليك", "Salicylic Acid", [2.97], "organic_carboxylic",
          "قليل الذوبان (~2 g/L عند 25 °C)", "حمض فينولي كربوكسيلي؛ مادة أولية للأسبرين", SRC_HARRIS)
_add_acid("C9H8O4", "C9H8O4", "حمض الأسيتيل ساليسيليك (الأسبرين)", "Acetylsalicylic Acid (Aspirin)", [3.5], "organic_carboxylic",
          "قليل الذوبان (~3 g/L عند 20 °C)", "يتحلل مائياً ببطء إلى حمض الساليسيليك والأسيتيك", SRC_CRC_ORG)
_add_acid("C6H5OH", "C6H6O", "الفينول", "Phenol", [9.99], "organic_other",
          "ذوّاب (~83 g/L عند 20 °C)", "حمض ضعيف جداً؛ سام وآكل للجلد", SRC_HARRIS)
_add_acid("C3H4O4", "C3H4O4", "حمض المالونيك", "Malonic Acid", [2.85, 5.70], "organic_carboxylic",
          "شديد الذوبان (~763 g/L عند 20 °C)", "حمض ثنائي الكربوكسيل", SRC_HARRIS)
_add_acid("C4H6O4", "C4H6O4", "حمض السكسينيك (الكهرمان)", "Succinic Acid", [4.21, 5.64], "organic_carboxylic",
          "متوسط الذوبان (~58 g/L عند 20 °C)", "قيمتا pKa متقاربتان فلا يظهر قفز واضح للتكافؤ الأول", SRC_HARRIS)
_add_acid("C4H4O4 (maleic)", "C4H4O4", "حمض المالئيك", "Maleic Acid", [1.91, 6.33], "organic_carboxylic",
          "شديد الذوبان (~479 g/L عند 20 °C)", "الأيزومر cis؛ قيمتا pKa متباعدتان بسبب الرابطة الهيدروجينية الداخلية", SRC_HARRIS)
_add_acid("C4H4O4 (fumaric)", "C4H4O4", "حمض الفوماريك", "Fumaric Acid", [3.05, 4.49], "organic_carboxylic",
          "قليل الذوبان (~4.9 g/L عند 25 °C)", "الأيزومر trans؛ قيمتا pKa أقرب من المالئيك", SRC_HARRIS)
_add_acid("C8H6O4", "C8H6O4", "حمض الفثاليك", "Phthalic Acid", [2.95, 5.41], "organic_carboxylic",
          "قليل الذوبان (~6 g/L عند 25 °C)", "ملحه الحمضي (فثالات البوتاسيوم الهيدروجينية) معيار أولي شائع", SRC_HARRIS)
_add_acid("C4H6O5", "C4H6O5", "حمض الماليك (التفاحيك)", "Malic Acid", [3.46, 5.10], "organic_carboxylic",
          "شديد الذوبان (~558 g/L عند 20 °C)", "حمض هيدروكسي ثنائي الكربوكسيل في التفاح", SRC_HARRIS)
_add_acid("C4H6O6", "C4H6O6", "حمض الطرطريك", "Tartaric Acid", [3.04, 4.37], "organic_carboxylic",
          "شديد الذوبان (~139 g/100 mL عند 20 °C)", "حمض ثنائي الهيدروكسي ثنائي الكربوكسيل؛ ثنائي pKa متقاربان", SRC_HARRIS)
_add_acid("C6H8O6", "C6H8O6", "حمض الأسكوربيك (فيتامين C)", "Ascorbic Acid (Vitamin C)", [4.10, 11.6], "organic_other",
          "شديد الذوبان (~330 g/L)", "pKa2 مرتفع جداً فلا يُعاير عملياً إلا البروتون الأول؛ يتأكسد بسهولة", SRC_HARRIS)


# --------------------------------------------------------------------------------------
# Extra bases
# --------------------------------------------------------------------------------------
EXTRA_BASES: dict[str, dict[str, Any]] = {}


def _add_base(*args: Any, **kwargs: Any) -> None:
    spec = _base(*args, **kwargs)
    EXTRA_BASES[spec["formula"]] = spec


# Strong hydroxides (complete dissociation in dilute solution)
_add_base("LiOH", "LiOH", "هيدروكسيد الليثيوم", "Lithium Hydroxide", [0.0], "hydroxide_strong",
          "ذوّاب (~128 g/L عند 20 °C)", "أضعف هيدروكسيدات القلويات ذوباناً؛ آكل", SRC_CRC_INORG, strong=True)
_add_base("RbOH", "RbOH", "هيدروكسيد الروبيديوم", "Rubidium Hydroxide", [0.0], "hydroxide_strong",
          "شديد الذوبان", "قاعدة قلوية قوية جداً؛ آكل ومتميّع", SRC_CRC_INORG, strong=True)
_add_base("CsOH", "CsOH", "هيدروكسيد السيزيوم", "Cesium Hydroxide", [0.0], "hydroxide_strong",
          "شديد الذوبان", "من أقوى القواعد القلوية؛ آكل ومتميّع", SRC_CRC_INORG, strong=True)
_add_base("Sr(OH)2", "Sr(OH)2", "هيدروكسيد السترونشيوم", "Strontium Hydroxide", [0.0, 0.0], "hydroxide_strong",
          "قليل إلى متوسط الذوبان (يزداد بالحرارة)", "قاعدة قوية ثنائية الهيدروكسيد؛ آكلة", SRC_CRC_INORG,
          hydroxides=2, strong=True)
# Ammonia-like and amines (pKb = 14 - pKa(BH+), pKa(BH+) from the named source)
_add_base("N2H4", "N2H4", "الهيدرازين", "Hydrazine", _pkb_from_conj_pka([8.07]), "ammonia_hydrazine_like",
          "يمتزج كلياً مع الماء", "قاعدة ضعيفة؛ شديد السمية ومسرطن محتمل؛ عامل مختزل", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("NH2OH", "NH3O", "الهيدروكسيل أمين", "Hydroxylamine", _pkb_from_conj_pka([5.96]), "ammonia_hydrazine_like",
          "ذوّاب في الماء", "قاعدة ضعيفة؛ غير مستقر (قد ينفجر عند التسخين)", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("(CH3)2NH", "C2H7N", "ثنائي ميثيل أمين", "Dimethylamine", _pkb_from_conj_pka([10.77]), "amine_aliphatic",
          "شديد الذوبان", "أمين ثانوي أليفاتي؛ غاز رائحته سمكية", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("(CH3)3N", "C3H9N", "ثلاثي ميثيل أمين", "Trimethylamine", _pkb_from_conj_pka([9.80]), "amine_aliphatic",
          "شديد الذوبان", "أمين ثالثي؛ غاز رائحته سمكية قوية", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("C2H4(NH2)2", "C2H8N2", "إيثيلين ثنائي الأمين", "Ethylenediamine", _pkb_from_conj_pka([9.93, 6.85]),
          "amine_aliphatic", "يمتزج كلياً مع الماء", "قاعدة ضعيفة ثنائية؛ ليجاند شائع", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("C5H11N", "C5H11N", "البيبيريدين", "Piperidine", _pkb_from_conj_pka([11.12]), "amine_aliphatic",
          "يمتزج كلياً مع الماء", "أمين حلقي ثانوي؛ أقوى قاعدية من الأمونيا بكثير", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("C2H7NO", "C2H7NO", "إيثانول أمين", "Ethanolamine", _pkb_from_conj_pka([9.50]), "amine_aliphatic",
          "يمتزج كلياً مع الماء", "أمين أولي كحولي؛ يستعمل في امتصاص CO2", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("C5H5N", "C5H5N", "البيريدين", "Pyridine", _pkb_from_conj_pka([5.20]), "amine_aromatic_heterocyclic",
          "يمتزج كلياً مع الماء", "قاعدة ضعيفة حلقية عطرية؛ مذيب سام ورائحته نفاذة", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("C6H5NH2", "C6H7N", "الأنيلين", "Aniline", _pkb_from_conj_pka([4.60]), "amine_aromatic_heterocyclic",
          "قليل الذوبان (~36 g/L عند 20 °C)", "قاعدة ضعيفة جداً لأن زوج الإلكترونات متمركز في الحلقة؛ سام", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("(HOCH2)3CNH2", "C4H11NO3", "تريس (Tris)", "Tris(hydroxymethyl)aminomethane", _pkb_from_conj_pka([8.07]), "buffer_base",
          "شديد الذوبان (~550 g/L عند 25 °C)", "قاعدة منظمة شائعة في الكيمياء الحيوية؛ معيار أولي للأحماض", SRC_HARRIS + "; " + SRC_DERIVED_KW)
# Basic salts: pKb from the pKa chain of the conjugate acid
_add_base("K2CO3", "K2CO3", "كربونات البوتاسيوم", "Potassium Carbonate", _pkb_from_conj_pka([10.33, 6.35]), "salt_basic",
          "شديد الذوبان (~1120 g/L عند 20 °C)", "قاعدة ضعيفة ثنائية؛ تمتص الرطوبة", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("KHCO3", "KHCO3", "بيكربونات البوتاسيوم", "Potassium Bicarbonate", _pkb_from_conj_pka([6.35]), "salt_basic",
          "ذوّاب (~224 g/L عند 20 °C)", "قاعدة ضعيفة جداً أمفوتيرية", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("Na3PO4", "Na3PO4", "فوسفات الصوديوم الثلاثي", "Trisodium Phosphate", _pkb_from_conj_pka([12.35, 7.20, 2.15]),
          "salt_basic", "ذوّاب في الماء", "قاعدة ثلاثية؛ pKb1 صغير فتسلك كقاعدة قوية نسبياً في الخطوة الأولى", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("CH3COONa", "C2H3NaO2", "أسيتات الصوديوم", "Sodium Acetate", _pkb_from_conj_pka([4.756]), "salt_basic",
          "شديد الذوبان (~464 g/L عند 20 °C)", "قاعدة ضعيفة جداً؛ تستعمل مع حمض الأسيتيك كمحلول منظم", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("NaCN", "CNNa", "سيانيد الصوديوم", "Sodium Cyanide", _pkb_from_conj_pka([9.21]), "salt_basic",
          "شديد الذوبان", "قاعدة ضعيفة؛ سام جداً؛ يحرّر HCN القاتل مع الأحماض", SRC_HARRIS + "; " + SRC_DERIVED_KW)
_add_base("NaF", "FNa", "فلوريد الصوديوم", "Sodium Fluoride", _pkb_from_conj_pka([3.17]), "salt_basic",
          "ذوّاب (~40 g/L عند 20 °C)", "قاعدة ضعيفة جداً؛ سام عند الجرعات العالية", SRC_HARRIS + "; " + SRC_DERIVED_KW)


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


_add_indicator("methyl_violet", "بنفسجي الميثيل (Methyl Violet)", "Methyl Violet", (0.0, 1.6),
               "أصفر (Yellow)", "أزرق بنفسجي (Blue-Violet)", "#FACC15", "#22C55E", "#4338CA",
               "مدى شديد الحموضة؛ يلائم الأحماض القوية المركزة جداً.", SRC_INDICATORS_CRC)
_add_indicator("thymol_blue_acid", "أزرق الثيمول - المدى الحمضي (Thymol Blue, acidic range)", "Thymol Blue (acid range)",
               (1.2, 2.8), "أحمر (Red)", "أصفر (Yellow)", "#DC2626", "#FB923C", "#FACC15",
               "المدى الحمضي الأول لأزرق الثيمول؛ يلائم نقاط تكافؤ حمضية جداً.", SRC_INDICATORS_HARRIS)
_add_indicator("methyl_yellow", "الأصفر الميثيلي (Methyl Yellow)", "Methyl Yellow", (2.9, 4.0),
               "أحمر (Red)", "أصفر (Yellow)", "#DC2626", "#F97316", "#FACC15",
               "مدى قريب من الميثيل البرتقالي؛ يلائم القواعد الضعيفة مع الأحماض القوية.", SRC_INDICATORS_HARRIS)
_add_indicator("dinitrophenol_24", "2,4-ثنائي نيتروفينول (2,4-Dinitrophenol)", "2,4-Dinitrophenol", (2.8, 4.7),
               "عديم اللون (Colorless)", "أصفر (Yellow)", "#F8FAFC", "#FEF08A", "#FACC15",
               "مدى حمضي ضعيف؛ تغيّر لوني ضعيف الوضوح مقارنة بالكواشف الملوّنة بالكامل.", SRC_INDICATORS_CRC)
_add_indicator("bromophenol_blue", "أزرق البروموفينول (Bromophenol Blue)", "Bromophenol Blue", (3.0, 4.6),
               "أصفر (Yellow)", "أزرق بنفسجي (Blue-Violet)", "#FACC15", "#4ADE80", "#4F46E5",
               "ملائم لمعايرة القواعد الضعيفة مع الأحماض القوية عند نقطة تكافؤ حمضية.", SRC_INDICATORS_HARRIS)
_add_indicator("congo_red", "أحمر الكونغو (Congo Red)", "Congo Red", (3.0, 5.0),
               "أزرق بنفسجي (Blue-Violet)", "أحمر (Red)", "#4338CA", "#A855F7", "#DC2626",
               "مدى حمضي واسع؛ استعمال تاريخي كاشف للأحماض المعدنية.", SRC_INDICATORS_CRC)
_add_indicator("chlorophenol_red", "أحمر الكلوروفينول (Chlorophenol Red)", "Chlorophenol Red", (4.8, 6.4),
               "أصفر (Yellow)", "أحمر (Red)", "#FACC15", "#FB923C", "#DC2626",
               "ملائم لنقاط تكافؤ حمضية خفيفة (pH 5-6).", SRC_INDICATORS_HARRIS)
_add_indicator("neutral_red", "الأحمر المتعادل (Neutral Red)", "Neutral Red", (6.8, 8.0),
               "أحمر (Red)", "أصفر (Yellow)", "#DC2626", "#F97316", "#FACC15",
               "مدى حول التعادل؛ يلائم نقطة تكافؤ قريبة من 7.", SRC_INDICATORS_HARRIS)
_add_indicator("brilliant_yellow", "الأصفر اللامع (Brilliant Yellow)", "Brilliant Yellow", (6.6, 7.8),
               "أصفر (Yellow)", "أحمر (Red)", "#FACC15", "#FB923C", "#DC2626",
               "مدى حول التعادل؛ بديل لأزرق البروموثيمول.", SRC_INDICATORS_HARRIS)
_add_indicator("cresol_red", "أحمر الكريزول (Cresol Red)", "Cresol Red", (7.0, 8.8),
               "أصفر (Yellow)", "أحمر (Red)", "#FACC15", "#FB923C", "#DC2626",
               "ملائم للأحماض الضعيفة مع القواعد القوية عند نقطة تكافؤ قاعدية خفيفة.", SRC_INDICATORS_HARRIS)
_add_indicator("m_cresol_purple", "أرجواني الميتا-كريزول (m-Cresol Purple)", "m-Cresol Purple", (7.4, 9.0),
               "أصفر (Yellow)", "أرجواني (Purple)", "#FACC15", "#F472B6", "#7E22CE",
               "ملائم لنقاط تكافؤ قاعدية خفيفة إلى متوسطة.", SRC_INDICATORS_CRC)
_add_indicator("thymolphthalein", "الثيمولفثالين (Thymolphthalein)", "Thymolphthalein", (9.3, 10.5),
               "عديم اللون (Colorless)", "أزرق (Blue)", "#F8FAFC", "#BFDBFE", "#2563EB",
               "ملائم لمعايرة الأحماض الضعيفة جداً مع القواعد القوية (تكافؤ pH ~9.5-10).", SRC_INDICATORS_HARRIS)
_add_indicator("alizarin_yellow_r", "أصفر الأليزارين R (Alizarin Yellow R)", "Alizarin Yellow R", (10.1, 12.0),
               "أصفر (Yellow)", "أحمر (Red)", "#FACC15", "#FB923C", "#DC2626",
               "مدى قاعدي مرتفع؛ يلائم أحماضاً ضعيفة جداً أو معايرات قاعدية شديدة.", SRC_INDICATORS_HARRIS)
_add_indicator("indigo_carmine", "كارمين النيلة (Indigo Carmine)", "Indigo Carmine", (11.4, 13.0),
               "أزرق (Blue)", "أصفر (Yellow)", "#2563EB", "#65A30D", "#FACC15",
               "مدى قاعدي شديد؛ استعمال محدود.", SRC_INDICATORS_CRC)


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

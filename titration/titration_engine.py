"""
Expanded Titration Engine for MolecuLens Virtual Titration Lab.

Supports:
1. 41 Acids (Strong, Weak, Monoprotic, Diprotic, Triprotic):
   HCl, H2SO4, HNO3, H3PO4, CH3COOH, HCOOH, C6H8O7, H2C2O4, C3H6O3, C7H6O2, H2CO3, H3BO3, HF.
2. 30 Bases (Strong, Weak, Mono/Di-functional Hydroxides, Amines, Carbonates):
   NaOH, KOH, Ca(OH)2, Ba(OH)2, NH3, Mg(OH)2, CH3NH2, C2H5NH2, NaHCO3, Na2CO3.
3. 23 Color Indicators with acidic/basic/transition hex colors and pH transition intervals:
   Phenolphthalein, Methyl Orange, Methyl Red, Bromothymol Blue, Bromocresol Green,
   Bromocresol Purple, Thymol Blue, Litmus, Phenol Red.
4. Automatic optimal indicator selection, polyprotic stepwise equivalence points,
   endpoint vs. equivalence point detection, and dynamic Erlenmeyer flask color SVG.
"""

from __future__ import annotations

import html
import math
from typing import Any
try:
    import plotly.graph_objects as go
except ImportError:
    from core._plotly_fallback import register_plotly_fallback_in_sys_modules
    register_plotly_fallback_in_sys_modules()
    import plotly.graph_objects as go

KW_25C = 1.0e-14

ACID_SPECIES_DATABASE: dict[str, dict[str, Any]] = {
    "HCl": {
        "formula": "HCl",
        "name_ar": "حمض الهيدروكلوريك",
        "name_en": "Hydrochloric Acid",
        "type": "strong",
        "strength_ar": "حمض معدني قوي أحادي البروتون",
        "pKa": -6.3,
        "pKa_list": [-6.3],
        "Ka": 1e6,
        "protons": 1,
        "molar_mass": 36.46,
        "solubility_ar": "ذوبان كامل في الماء (> 700 g/L عند 20 °C)",
    },
    "H2SO4": {
        "formula": "H2SO4",
        "name_ar": "حمض الكبريتيك",
        "name_en": "Sulfuric Acid",
        "type": "strong",
        "strength_ar": "حمض معدني قوي ثنائي البروتون (التأين الأول تام، والثاني pKa2 = 1.99)",
        "pKa": -3.0,
        "pKa_list": [-3.0, 1.99],
        "Ka": 1e3,
        "protons": 2,
        "molar_mass": 98.079,
        "solubility_ar": "يمتزج كلياً مع الماء بجميع النسب (طارد للحرارة بشدة)",
    },
    "HNO3": {
        "formula": "HNO3",
        "name_ar": "حمض النيتريك",
        "name_en": "Nitric Acid",
        "type": "strong",
        "strength_ar": "حمض معدني قوي أحادي البروتون",
        "pKa": -1.4,
        "pKa_list": [-1.4],
        "Ka": 25.0,
        "protons": 1,
        "molar_mass": 63.01,
        "solubility_ar": "يمتزج كلياً مع الماء بجميع النسب",
    },
    "H3PO4": {
        "formula": "H3PO4",
        "name_ar": "حمض الفوسفوريك",
        "name_en": "Phosphoric Acid",
        "type": "weak",
        "strength_ar": "حمض غير عضوي ضعيف ثلاثي البروتون (pKa1=2.15, pKa2=7.20, pKa3=12.35)",
        "pKa": 2.15,
        "pKa_list": [2.15, 7.20, 12.35],
        "Ka": 7.08e-3,
        "protons": 3,
        "molar_mass": 97.994,
        "solubility_ar": "شديد الذوبان في الماء (5480 g/L عند 20 °C)",
    },
    "CH3COOH": {
        "formula": "CH3COOH",
        "name_ar": "حمض الأسيتيك (حمض الخليك)",
        "name_en": "Acetic Acid",
        "type": "weak",
        "strength_ar": "حمض كربوكسيلي ضعيف أحادي البروتون",
        "pKa": 4.756,
        "pKa_list": [4.756],
        "Ka": 1.75e-5,
        "protons": 1,
        "molar_mass": 60.052,
        "solubility_ar": "يمتزج كلياً مع الماء بجميع النسب",
    },
    "HCOOH": {
        "formula": "HCOOH",
        "name_ar": "حمض الفورميك (حمض النمليك)",
        "name_en": "Formic Acid",
        "type": "weak",
        "strength_ar": "حمض كربوكسيلي ضعيف أحادي البروتون",
        "pKa": 3.75,
        "pKa_list": [3.75],
        "Ka": 1.78e-4,
        "protons": 1,
        "molar_mass": 46.025,
        "solubility_ar": "يمتزج كلياً مع الماء بجميع النسب",
    },
    "C6H8O7": {
        "formula": "C6H8O7",
        "name_ar": "حمض الستريك (حمض الليمون)",
        "name_en": "Citric Acid",
        "type": "weak",
        "strength_ar": "حمض عضوي ضعيف ثلاثي البروتون (pKa1=3.13, pKa2=4.76, pKa3=6.40)",
        "pKa": 3.13,
        "pKa_list": [3.13, 4.76, 6.40],
        "Ka": 7.41e-4,
        "protons": 3,
        "molar_mass": 192.124,
        "solubility_ar": "شديد الذوبان في الماء (592 g/L عند 20 °C)",
    },
    "H2C2O4": {
        "formula": "H2C2O4",
        "name_ar": "حمض الأكساليك",
        "name_en": "Oxalic Acid",
        "type": "weak",
        "strength_ar": "حمض عضوي ضعيف ثنائي البروتون (pKa1=1.25, pKa2=4.27)",
        "pKa": 1.25,
        "pKa_list": [1.25, 4.27],
        "Ka": 5.62e-2,
        "protons": 2,
        "molar_mass": 90.03,
        "solubility_ar": "ذائب في الماء (143 g/L عند 25 °C)",
    },
    "C3H6O3": {
        "formula": "C3H6O3",
        "name_ar": "حمض اللاكتيك (حمض اللبن)",
        "name_en": "Lactic Acid",
        "type": "weak",
        "strength_ar": "حمض هيدروكسي كربوكسيلي ضعيف أحادي البروتون",
        "pKa": 3.86,
        "pKa_list": [3.86],
        "Ka": 1.38e-4,
        "protons": 1,
        "molar_mass": 90.078,
        "solubility_ar": "شديد الذوبان ويمتزج مع الماء",
    },
    "C7H6O2": {
        "formula": "C7H6O2",
        "name_ar": "حمض البنزويك",
        "name_en": "Benzoic Acid",
        "type": "weak",
        "strength_ar": "حمض كربوكسيلي عطري ضعيف أحادي البروتون",
        "pKa": 4.20,
        "pKa_list": [4.20],
        "Ka": 6.31e-5,
        "protons": 1,
        "molar_mass": 122.12,
        "solubility_ar": "قليل الذوبان في الماء البارد (3.4 g/L عند 25 °C)، يذوب بالماء الساخن",
    },
    "H2CO3": {
        "formula": "H2CO3",
        "name_ar": "حمض الكربونيك",
        "name_en": "Carbonic Acid",
        "type": "weak",
        "strength_ar": "حمض ضعيف ثنائي البروتون (pKa1=6.35, pKa2=10.33)",
        "pKa": 6.35,
        "pKa_list": [6.35, 10.33],
        "Ka": 4.47e-7,
        "protons": 2,
        "molar_mass": 62.024,
        "solubility_ar": "يوجد في المحلول المائي باتزان مع غاز CO2 الذائب (~1.45 g/L)",
    },
    "H3BO3": {
        "formula": "H3BO3",
        "name_ar": "حمض البوريك",
        "name_en": "Boric Acid",
        "type": "weak",
        "strength_ar": "حمض لويس ضعيف جداً أحادي التكافؤ الفعال في الماء (pKa=9.24)",
        "pKa": 9.24,
        "pKa_list": [9.24],
        "Ka": 5.75e-10,
        "protons": 1,
        "molar_mass": 61.83,
        "solubility_ar": "ذائب في الماء (57 g/L عند 25 °C)",
    },
    "HF": {
        "formula": "HF",
        "name_ar": "حمض الهيدروفلوريك",
        "name_en": "Hydrofluoric Acid",
        "type": "weak",
        "strength_ar": "حمض هالوجيني ضعيف أحادي البروتون",
        "pKa": 3.17,
        "pKa_list": [3.17],
        "Ka": 6.76e-4,
        "protons": 1,
        "molar_mass": 20.006,
        "solubility_ar": "يمتزج كلياً مع الماء",
    },
}

BASE_SPECIES_DATABASE: dict[str, dict[str, Any]] = {
    "NaOH": {
        "formula": "NaOH",
        "name_ar": "هيدروكسيد الصوديوم",
        "name_en": "Sodium Hydroxide",
        "type": "strong",
        "strength_ar": "قاعدة قلوية قوية أحادية الهيدروكسيد",
        "pKb": 0.0,
        "pKb_list": [0.0],
        "Kb": 1e6,
        "hydroxides": 1,
        "molar_mass": 39.997,
        "solubility_ar": "شديد الذوبان في الماء (1090 g/L عند 20 °C)",
    },
    "KOH": {
        "formula": "KOH",
        "name_ar": "هيدروكسيد البوتاسيوم",
        "name_en": "Potassium Hydroxide",
        "type": "strong",
        "strength_ar": "قاعدة قلوية قوية أحادية الهيدروكسيد",
        "pKb": 0.0,
        "pKb_list": [0.0],
        "Kb": 1e6,
        "hydroxides": 1,
        "molar_mass": 56.105,
        "solubility_ar": "شديد الذوبان في الماء (1210 g/L عند 25 °C)",
    },
    "Ca(OH)2": {
        "formula": "Ca(OH)2",
        "name_ar": "هيدروكسيد الكالسيوم (ماء الجير)",
        "name_en": "Calcium Hydroxide",
        "type": "strong",
        "strength_ar": "قاعدة قوية ثنائية الهيدروكسيد (تامة التأين في المحاليل المخففة)",
        "pKb": 1.37,
        "pKb_list": [1.37, 2.43],
        "Kb": 1e3,
        "hydroxides": 2,
        "molar_mass": 74.093,
        "solubility_ar": "قليل الذوبان في الماء (1.73 g/L عند 20 °C — ماء الجير الرائق)",
    },
    "Ba(OH)2": {
        "formula": "Ba(OH)2",
        "name_ar": "هيدروكسيد الباريوم",
        "name_en": "Barium Hydroxide",
        "type": "strong",
        "strength_ar": "قاعدة قلوية ترابية قوية ثنائية الهيدروكسيد",
        "pKb": 0.15,
        "pKb_list": [0.15],
        "Kb": 1e4,
        "hydroxides": 2,
        "molar_mass": 171.34,
        "solubility_ar": "متوسط الذوبان في الماء (46.8 g/L عند 25 °C)",
    },
    "NH3": {
        "formula": "NH3",
        "name_ar": "الأمونيا (محلول النشادر)",
        "name_en": "Ammonia",
        "type": "weak",
        "strength_ar": "قاعدة نيتروجينية ضعيفة أحادية التكافؤ",
        "pKb": 4.75,
        "pKb_list": [4.75],
        "Kb": 1.77e-5,
        "hydroxides": 1,
        "molar_mass": 17.031,
        "solubility_ar": "شديد الذوبان في الماء (482 g/L عند 25 °C)",
    },
    "Mg(OH)2": {
        "formula": "Mg(OH)2",
        "name_ar": "هيدروكسيد المغنيسيوم (حليب المغنيسيا)",
        "name_en": "Magnesium Hydroxide",
        "type": "weak",
        "strength_ar": "قاعدة ثنائية الهيدروكسيد شحيحة الذوبان (Ksp = 5.61e-12)",
        "pKb": 2.60,
        "pKb_list": [2.60],
        "Kb": 2.5e-3,
        "hydroxides": 2,
        "molar_mass": 58.32,
        "solubility_ar": "شحيح الذوبان جداً في الماء (0.009 g/L عند 25 °C)",
    },
    "CH3NH2": {
        "formula": "CH3NH2",
        "name_ar": "ميثيل أمين",
        "name_en": "Methylamine",
        "type": "weak",
        "strength_ar": "أمين أولي أليفاتي (قاعدة عضوية ضعيفة أقوى من الأمونيا)",
        "pKb": 3.36,
        "pKb_list": [3.36],
        "Kb": 4.37e-4,
        "hydroxides": 1,
        "molar_mass": 31.057,
        "solubility_ar": "شديد الذوبان في الماء (1080 g/L عند 25 °C)",
    },
    "C2H5NH2": {
        "formula": "C2H5NH2",
        "name_ar": "إيثيل أمين",
        "name_en": "Ethylamine",
        "type": "weak",
        "strength_ar": "أمين أولي أليفاتي (قاعدة عضوية ضعيفة)",
        "pKb": 3.25,
        "pKb_list": [3.25],
        "Kb": 5.62e-4,
        "hydroxides": 1,
        "molar_mass": 45.084,
        "solubility_ar": "يمتزج كلياً مع الماء",
    },
    "NaHCO3": {
        "formula": "NaHCO3",
        "name_ar": "بيكربونات الصوديوم (صودا الخبز)",
        "name_en": "Sodium Bicarbonate",
        "type": "weak",
        "strength_ar": "ملح قاعدي ضعيف أحادي التكافؤ الحمضي",
        "pKb": 7.65,
        "pKb_list": [7.65],
        "Kb": 2.24e-8,
        "hydroxides": 1,
        "molar_mass": 84.007,
        "solubility_ar": "ذائب في الماء (96 g/L عند 20 °C)",
    },
    "Na2CO3": {
        "formula": "Na2CO3",
        "name_ar": "كربونات الصوديوم (رماد الصودا)",
        "name_en": "Sodium Carbonate",
        "type": "weak",
        "strength_ar": "قاعدة ضعيفة ثنائية التكافؤ (pKb1=3.67, pKb2=7.65)",
        "pKb": 3.67,
        "pKb_list": [3.67, 7.65],
        "Kb": 2.14e-4,
        "hydroxides": 2,
        "molar_mass": 105.988,
        "solubility_ar": "ذائب في الماء (215 g/L عند 20 °C)",
    },
}

INDICATORS_DATABASE: dict[str, dict[str, Any]] = {
    "phenolphthalein": {
        "key": "phenolphthalein",
        "name_ar": "فينولفثالين (Phenolphthalein)",
        "name_en": "Phenolphthalein",
        "pH_range": (8.2, 10.0),
        "color_acid": "عديم اللون (Colorless)",
        "color_base": "وردي قرمزي (Vivid Pink / Magenta)",
        "color_acid_hex": "#F8FAFC",
        "color_transition_hex": "#F9A8D4",
        "color_base_hex": "#EC4899",
        "color_hex": "#EC4899",
        "suitability_ar": "مثالي لمعايرة الأحماض الضعيفة مع القواعد القوية (نقطة التكافؤ pH > 8).",
    },
    "methyl_orange": {
        "key": "methyl_orange",
        "name_ar": "ميثيل برتقالي (Methyl Orange)",
        "name_en": "Methyl Orange",
        "pH_range": (3.1, 4.4),
        "color_acid": "أحمر (Red)",
        "color_base": "أصفر برتقالي (Yellow-Orange)",
        "color_acid_hex": "#EF4444",
        "color_transition_hex": "#FB923C",
        "color_base_hex": "#FACC15",
        "color_hex": "#F97316",
        "suitability_ar": "مثالي لمعايرة القواعد الضعيفة أو الكربونات مع الأحماض القوية (نقطة التكافؤ الحمضية pH 3.1 - 4.4).",
    },
    "methyl_red": {
        "key": "methyl_red",
        "name_ar": "أحمر الميثيل (Methyl Red)",
        "name_en": "Methyl Red",
        "pH_range": (4.4, 6.2),
        "color_acid": "أحمر (Red)",
        "color_base": "أصفر (Yellow)",
        "color_acid_hex": "#DC2626",
        "color_transition_hex": "#F97316",
        "color_base_hex": "#FDE047",
        "color_hex": "#EA580C",
        "suitability_ar": "مثالي لمعايرة القواعد الضعيفة (مثل الأمونيا والأمينات) مع الأحماض القوية (نقطة التكافؤ pH ~ 5.2).",
    },
    "bromothymol_blue": {
        "key": "bromothymol_blue",
        "name_ar": "أزرق بروموثيمول (Bromothymol Blue)",
        "name_en": "Bromothymol Blue",
        "pH_range": (6.0, 7.6),
        "color_acid": "أصفر (Yellow)",
        "color_base": "أزرق (Blue)",
        "color_acid_hex": "#FACC15",
        "color_transition_hex": "#10B981",
        "color_base_hex": "#2563EB",
        "color_hex": "#3B82F6",
        "suitability_ar": "مثالي للمعايرات المتعادلة بين الأحماض القوية والقواعد القوية (نقطة التكافؤ pH ≈ 7.0).",
    },
    "bromocresol_green": {
        "key": "bromocresol_green",
        "name_ar": "أخضر بروموكريسول (Bromocresol Green)",
        "name_en": "Bromocresol Green",
        "pH_range": (3.8, 5.4),
        "color_acid": "أصفر (Yellow)",
        "color_base": "أزرق مخضر (Blue-Green)",
        "color_acid_hex": "#FDE047",
        "color_transition_hex": "#34D399",
        "color_base_hex": "#0284C7",
        "color_hex": "#0D9488",
        "suitability_ar": "مثالي لمعايرة القواعد الضعيفة وأملاح الكربونات والبيكربونات (pH 3.8 - 5.4).",
    },
    "bromocresol_purple": {
        "key": "bromocresol_purple",
        "name_ar": "أرجواني بروموكريسول (Bromocresol Purple)",
        "name_en": "Bromocresol Purple",
        "pH_range": (5.2, 6.8),
        "color_acid": "أصفر (Yellow)",
        "color_base": "أرجواني (Purple)",
        "color_acid_hex": "#FACC15",
        "color_transition_hex": "#C084FC",
        "color_base_hex": "#7E22CE",
        "color_hex": "#9333EA",
        "suitability_ar": "ملائم للمعايرات ذات نقطة التكافؤ الحمضية الخفيفة إلى المتعادلة (pH 5.2 - 6.8).",
    },
    "thymol_blue": {
        "key": "thymol_blue",
        "name_ar": "أزرق الثيمول (Thymol Blue)",
        "name_en": "Thymol Blue",
        "pH_range": (8.0, 9.6),
        "color_acid": "أصفر (Yellow)",
        "color_base": "أزرق (Blue)",
        "color_acid_hex": "#FACC15",
        "color_transition_hex": "#38BDF8",
        "color_base_hex": "#1D4ED8",
        "color_hex": "#2563EB",
        "suitability_ar": "ملائم لمعايرة الأحماض الضعيفة مع القواعد القوية في الوسط القاعدي (pH 8.0 - 9.6).",
    },
    "litmus": {
        "key": "litmus",
        "name_ar": "عباد الشمس (Litmus)",
        "name_en": "Litmus",
        "pH_range": (4.5, 8.3),
        "color_acid": "أحمر (Red)",
        "color_base": "أزرق (Blue)",
        "color_acid_hex": "#EF4444",
        "color_transition_hex": "#A855F7",
        "color_base_hex": "#3B82F6",
        "color_hex": "#6366F1",
        "suitability_ar": "كاشف عام واسع النطاق للتمييز بين الوسط الحمضي والقاعدي (pH 4.5 - 8.3).",
    },
    "phenol_red": {
        "key": "phenol_red",
        "name_ar": "أحمر الفينول (Phenol Red)",
        "name_en": "Phenol Red",
        "pH_range": (6.8, 8.4),
        "color_acid": "أصفر (Yellow)",
        "color_base": "أحمر قرمزي (Red)",
        "color_acid_hex": "#FACC15",
        "color_transition_hex": "#FB7185",
        "color_base_hex": "#DC2626",
        "color_hex": "#EF4444",
        "suitability_ar": "ملائم للمعايرات القريبة من التعادل والوسط الفسيولوجي (pH 6.8 - 8.4).",
    },
}

from titration.acid_base_catalog import finalize_species_tables

finalize_species_tables(ACID_SPECIES_DATABASE, BASE_SPECIES_DATABASE, INDICATORS_DATABASE)


def solve_quadratic_positive(a: float, b: float, c: float) -> float:
    """Solves ax^2 + bx + c = 0 for the physical positive root."""
    disc = b * b - 4 * a * c
    if disc < 0:
        return 0.0
    sqrt_d = math.sqrt(disc)
    r1 = (-b + sqrt_d) / (2 * a)
    r2 = (-b - sqrt_d) / (2 * a)
    if r1 > 0:
        return r1
    if r2 > 0:
        return r2
    return 0.0




def kw_at_temperature(temperature_c: float = 25.0) -> float:
    """
    Computes water autoionization constant Kw(T) from 0 °C to 60 °C using standard
    thermodynamic Van 't Hoff relation (ΔH°_w ≈ 55.836 kJ/mol, exact 1.0e-14 at 25.0 °C).
    """
    t_c = max(0.0, min(60.0, float(temperature_c)))
    t_k = t_c + 273.15
    t_ref = 298.15
    delta_h_j = 55836.0
    r_gas = 8.314462618
    ln_kw = math.log(1.0e-14) - (delta_h_j / r_gas) * (1.0 / t_k - 1.0 / t_ref)
    return math.exp(ln_kw)


def davies_activity_coefficient(charge: int, ionic_strength: float) -> float:
    """
    Davies equation for single-ion activity coefficient γ_z at 25 °C (valid up to I ≈ 0.5 M):
      log10(γ_z) = -0.509 * z^2 * (sqrt(I) / (1 + sqrt(I)) - 0.3 * I)
    """
    if charge == 0 or ionic_strength <= 0.0:
        return 1.0
    i_eff = min(float(ionic_strength), 0.5)
    sqrt_i = math.sqrt(i_eff)
    log10_g = -0.509 * (charge ** 2) * (sqrt_i / (1.0 + sqrt_i) - 0.3 * i_eff)
    return 10.0 ** log10_g


def _acid_avg_charge(h: float, ka_list: list[float], is_strong: bool, n_protons: int) -> float:
    """Average negative charge per mole of analytical acid C_a (mass balance across all α_j)."""
    if is_strong:
        return float(n_protons)
    # Denominator D = [H]^n + Ka1*[H]^(n-1) + Ka1*Ka2*[H]^(n-2) + ...
    n = len(ka_list)
    terms = [h ** n]
    prod = 1.0
    for j, ka in enumerate(ka_list, start=1):
        prod *= ka
        terms.append(prod * (h ** (n - j)))
    denom = sum(terms)
    if denom <= 0.0:
        return 0.0
    return sum(j * terms[j] for j in range(1, n + 1)) / denom


def _base_avg_charge(oh: float, kb_list: list[float], is_strong: bool, n_oh: int) -> float:
    """Average positive charge per mole of analytical base C_b (mass balance across all β_j)."""
    if is_strong:
        return float(n_oh)
    n = len(kb_list)
    terms = [oh ** n]
    prod = 1.0
    for j, kb in enumerate(kb_list, start=1):
        prod *= kb
        terms.append(prod * (oh ** (n - j)))
    denom = sum(terms)
    if denom <= 0.0:
        return 0.0
    return sum(j * terms[j] for j in range(1, n + 1)) / denom


def solve_equilibrium_ph_exact(
    c_acid_m: float,
    c_base_m: float,
    acid_spec: dict[str, Any],
    base_spec: dict[str, Any],
    temperature_c: float = 25.0,
    use_activity_correction: bool = False,
) -> float:
    """
    General analytical acid-base equilibrium solver using exact Mass Balance,
    Charge Balance, all stepwise Ka/Kb ionization fractions, exact Kw(T),
    and optional Davies ionic-strength activity correction.
    Does NOT artificially clip pH to [0, 14] (e.g., 2.0 M HCl yields pH ≈ -0.30).
    """
    kw = kw_at_temperature(temperature_c)
    acid_strong = acid_spec.get("type") == "strong"
    base_strong = base_spec.get("type") == "strong"
    n_p = int(acid_spec.get("protons", 1))
    n_oh = int(base_spec.get("hydroxides", 1))
    ka_list = [10.0 ** (-float(p)) for p in acid_spec.get("pKa_list", [acid_spec["pKa"]])]
    kb_list = [10.0 ** (-float(p)) for p in base_spec.get("pKb_list", [base_spec["pKb"]])]

    def charge_balance(ph: float) -> float:
        h = 10.0 ** (-ph)
        oh = kw / max(h, 1e-30)
        pos_base = c_base_m * _base_avg_charge(oh, kb_list, base_strong, n_oh)
        neg_acid = c_acid_m * _acid_avg_charge(h, ka_list, acid_strong, n_p)
        # Positive charges - Negative charges (strictly decreasing in pH!)
        return (h + pos_base) - (oh + neg_acid)

    lo, hi = -2.0, 16.0
    for _ in range(90):
        mid = 0.5 * (lo + hi)
        if charge_balance(mid) > 0.0:
            lo = mid
        else:
            hi = mid
    ph_ideal = 0.5 * (lo + hi)

    if not use_activity_correction:
        return round(ph_ideal, 3)

    # First-order Davies activity correction on a(H+) = γ_1 * [H+]
    h_val = 10.0 ** (-ph_ideal)
    oh_val = kw / max(h_val, 1e-30)
    ionic_strength = 0.5 * (h_val + oh_val + c_acid_m * n_p + c_base_m * n_oh)
    gamma_h = davies_activity_coefficient(1, ionic_strength)
    ph_activity = -math.log10(max(1e-30, gamma_h * h_val))
    return round(ph_activity, 3)


def compute_ph_from_h(h_conc: float) -> float:
    h = max(h_conc, 1e-14)
    return max(0.0, min(14.0, -math.log10(h)))


def compute_ph_from_oh(oh_conc: float) -> float:
    oh = max(oh_conc, 1e-14)
    poh = -math.log10(oh)
    return max(0.0, min(14.0, 14.0 - poh))


def select_optimal_indicator_key(ph_eq: float, case_key: str) -> str:
    """Selects the scientifically optimal color indicator based on equivalence pH and titration type."""
    if case_key == "strong_acid_strong_base" and 6.2 <= ph_eq <= 7.8:
        return "bromothymol_blue"
    candidates = []
    for k, ind in INDICATORS_DATABASE.items():
        low, high = ind.get("pH_range", (7.0, 7.0))
        mid = (low + high) / 2.0
        covers = low <= ph_eq <= high
        dist = abs(mid - ph_eq)
        candidates.append((0 if covers else 1, dist, k))
    candidates.sort()
    return candidates[0][2] if candidates else "phenolphthalein"


def get_indicator_state_at_ph(indicator_key: str, current_ph: float) -> dict[str, str]:
    """Returns the active hex color and Arabic/English color description of an indicator at current_ph."""
    ind = INDICATORS_DATABASE.get(indicator_key, INDICATORS_DATABASE["bromothymol_blue"])
    low, high = ind["pH_range"]
    if current_ph < low:
        return {
            "hex": ind.get("color_acid_hex", "#FACC15"),
            "label_ar": f"لون الوسط الحمضي: {ind['color_acid']}",
            "label_en": f"Acidic form: {ind['color_acid']}",
            "phase": "acidic",
        }
    if current_ph > high:
        return {
            "hex": ind.get("color_base_hex", "#3B82F6"),
            "label_ar": f"لون الوسط القاعدي: {ind['color_base']}",
            "label_en": f"Basic form: {ind['color_base']}",
            "phase": "basic",
        }
    return {
        "hex": ind.get("color_transition_hex", ind.get("color_hex", "#10B981")),
        "label_ar": f"منطقة التحول اللوني (نقطة النهاية): بين {ind['color_acid']} و {ind['color_base']}",
        "label_en": f"Transition zone (Endpoint): between {ind['color_acid']} and {ind['color_base']}",
        "phase": "transition",
    }


class TitrationSimulator:
    """
    Rigorously computes acid-base titration thermodynamics, polyprotic steps,
    equivalence points, and indicator endpoints.
    """

    @classmethod
    def simulate_titration(
        cls,
        analyte_formula: str,
        titrant_formula: str,
        analyte_volume_ml: float = 25.0,
        analyte_conc_m: float = 0.1,
        titrant_conc_m: float = 0.1,
        max_titrant_ratio: float = 2.2,
        points_count: int = 160,
        dilution_water_ml: float = 0.0,
    ) -> dict[str, Any]:
        is_acid_analyte = analyte_formula in ACID_SPECIES_DATABASE
        if is_acid_analyte:
            acid_spec = ACID_SPECIES_DATABASE[analyte_formula]
            base_spec = BASE_SPECIES_DATABASE.get(titrant_formula, BASE_SPECIES_DATABASE["NaOH"])
            mode = "acid_analyte_titrated_by_base"
            val_analyte = int(acid_spec.get("protons", 1))
            val_titrant = int(base_spec.get("hydroxides", 1))
        else:
            base_spec = BASE_SPECIES_DATABASE.get(analyte_formula, BASE_SPECIES_DATABASE["NaOH"])
            acid_spec = ACID_SPECIES_DATABASE.get(titrant_formula, ACID_SPECIES_DATABASE["HCl"])
            mode = "base_analyte_titrated_by_acid"
            val_analyte = int(base_spec.get("hydroxides", 1))
            val_titrant = int(acid_spec.get("protons", 1))

        c_analyte = max(0.001, float(analyte_conc_m))
        v_analyte_raw = max(1.0, float(analyte_volume_ml))
        v_dilution = max(0.0, float(dilution_water_ml or 0.0))
        v_analyte = v_analyte_raw + v_dilution
        c_titrant = max(0.001, float(titrant_conc_m))

        n_analyte_mmol = c_analyte * v_analyte_raw
        c_analyte_effective = n_analyte_mmol / v_analyte
        # First equivalence volume (1st proton/hydroxide) and total stoichiometric equivalence volume
        v_eq1_ml = (n_analyte_mmol * 1.0) / (c_titrant * val_titrant)
        v_eq_total_ml = (n_analyte_mmol * val_analyte) / (c_titrant * val_titrant)
        # Primary equivalence volume reported is v_eq_total_ml (equal to v_eq1_ml for monoprotic)
        v_eq_ml = v_eq_total_ml
        v_max_ml = v_eq_ml * max_titrant_ratio

        step_equivalence_volumes_ml = [
            round(v_eq1_ml * step_idx, 2) for step_idx in range(1, val_analyte + 1)
        ]

        acid_type = acid_spec["type"]
        base_type = base_spec["type"]
        is_polyprotic = (int(acid_spec.get("protons", 1)) > 1) or (int(base_spec.get("hydroxides", 1)) > 1)

        if acid_type == "strong" and base_type == "strong":
            case_key = "strong_acid_strong_base"
            case_label_ar = "معايرة حمض قوي مع قاعدة قوية (Strong Acid - Strong Base)"
        elif acid_type == "weak" and base_type == "strong":
            case_key = "weak_acid_strong_base"
            case_label_ar = (
                "معايرة حمض ضعيف متعدد البروتونات مع قاعدة قوية (Polyprotic Weak Acid - Strong Base)"
                if int(acid_spec.get("protons", 1)) > 1
                else "معايرة حمض ضعيف مع قاعدة قوية (Weak Acid - Strong Base)"
            )
        elif acid_type == "strong" and base_type == "weak":
            case_key = "strong_acid_weak_base"
            case_label_ar = "معايرة حمض قوي مع قاعدة ضعيفة (Strong Acid - Weak Base)"
        else:
            case_key = "weak_acid_weak_base"
            case_label_ar = "معايرة حمض ضعيف مع قاعدة ضعيفة (Weak Acid - Weak Base)"

        total_vol_at_eq = v_analyte + v_eq_ml
        c_salt_at_eq = n_analyte_mmol / total_vol_at_eq

        pka_list = acid_spec.get("pKa_list", [acid_spec["pKa"]])
        pkb_list = base_spec.get("pKb_list", [base_spec["pKb"]])
        effective_last_pka = float(pka_list[-1]) if pka_list else float(acid_spec["pKa"])
        effective_last_ka = 10.0 ** (-effective_last_pka)
        effective_last_pkb = float(pkb_list[-1]) if pkb_list else float(base_spec["pKb"])
        effective_last_kb = 10.0 ** (-effective_last_pkb)

        if case_key == "strong_acid_strong_base":
            ph_eq = 7.00
            ph_half_eq = None
        elif case_key == "weak_acid_strong_base":
            kb = KW_25C / max(1e-14, effective_last_ka)
            oh_salt = solve_quadratic_positive(1.0, kb, -kb * c_salt_at_eq)
            oh_eq = math.sqrt(oh_salt * oh_salt + KW_25C)
            ph_eq = round(compute_ph_from_oh(oh_eq), 2)
            ph_half_eq = round(float(acid_spec["pKa"]), 2)
        elif case_key == "strong_acid_weak_base":
            ka_conj = KW_25C / max(1e-14, effective_last_kb)
            h_salt = solve_quadratic_positive(1.0, ka_conj, -ka_conj * c_salt_at_eq)
            h_eq = math.sqrt(h_salt * h_salt + KW_25C)
            ph_eq = round(compute_ph_from_h(h_eq), 2)
            ph_half_eq = round(14.0 - float(base_spec["pKb"]), 2)
        else:
            ph_eq = round(7.0 + 0.5 * (float(acid_spec["pKa"]) - float(base_spec["pKb"])), 2)
            ph_half_eq = round(float(acid_spec["pKa"]), 2)

        v_points: list[float] = []
        base_step = v_max_ml / max(20, (points_count - 20))
        curr_v = 0.0
        while curr_v <= v_max_ml:
            v_points.append(round(curr_v, 4))
            curr_v += base_step

        near_eq = [
            v_eq_ml * 0.5,
            v_eq_ml * 0.90,
            v_eq_ml * 0.95,
            v_eq_ml * 0.98,
            v_eq_ml * 0.99,
            v_eq_ml * 0.995,
            v_eq_ml * 1.0,
            v_eq_ml * 1.005,
            v_eq_ml * 1.01,
            v_eq_ml * 1.02,
            v_eq_ml * 1.05,
            v_eq_ml * 1.10,
        ]
        for step_v in step_equivalence_volumes_ml:
            near_eq.extend([step_v * 0.5, step_v * 0.98, step_v, step_v * 1.02])

        for ep in near_eq:
            if 0.0 <= ep <= v_max_ml:
                v_points.append(round(ep, 4))
        v_points = sorted(set(v_points))

        curve_data: list[dict[str, Any]] = []
        eq_analyte_mmol = n_analyte_mmol * val_analyte

        for v_tit in v_points:
            v_tot = v_analyte + v_tit
            eq_tit_mmol = c_titrant * v_tit * val_titrant

            if mode == "acid_analyte_titrated_by_base":
                if v_tit == 0.0:
                    if acid_type == "strong":
                        ph_val = compute_ph_from_h(c_analyte_effective * val_analyte)
                    else:
                        ka = acid_spec["Ka"]
                        h_val = solve_quadratic_positive(1.0, ka, -ka * c_analyte_effective)
                        ph_val = compute_ph_from_h(h_val)
                    stage = "نقطة البداية (Initial Point)"
                elif v_tit < v_eq_ml - 1e-5:
                    if acid_type == "strong":
                        rem_h_mmol = eq_analyte_mmol - eq_tit_mmol
                        ph_val = compute_ph_from_h(rem_h_mmol / v_tot)
                    else:
                        # Polyprotic or monoprotic Henderson-Hasselbalch buffer region
                        equiv_progress = (eq_tit_mmol / max(1e-9, eq_analyte_mmol)) * val_analyte
                        step_idx = min(len(pka_list) - 1, max(0, int(math.floor(equiv_progress))))
                        frac_in_step = max(0.01, min(0.99, equiv_progress - step_idx))
                        active_pka = float(pka_list[step_idx])
                        ph_val = active_pka + math.log10(frac_in_step / (1.0 - frac_in_step))
                    stage = "منطقة المحلول المنظم / قبل التكافؤ (Buffer Region)"
                elif abs(v_tit - v_eq_ml) <= 1e-4:
                    ph_val = ph_eq
                    stage = "نقطة التكافؤ الكلي (Equivalence Point)"
                else:
                    excess_oh_mmol = eq_tit_mmol - eq_analyte_mmol
                    if base_type == "strong":
                        ph_val = compute_ph_from_oh(excess_oh_mmol / v_tot)
                    else:
                        c_excess = excess_oh_mmol / v_tot
                        c_salt = eq_analyte_mmol / v_tot
                        ratio = max(1e-6, c_excess / max(1e-6, c_salt))
                        ph_val = 14.0 - (base_spec["pKb"] - math.log10(ratio))
                    stage = "فائض المحلول القياسي / بعد التكافؤ (Post-Equivalence)"
            else:
                if v_tit == 0.0:
                    if base_type == "strong":
                        ph_val = compute_ph_from_oh(c_analyte_effective * val_analyte)
                    else:
                        kb = base_spec["Kb"]
                        oh_val = solve_quadratic_positive(1.0, kb, -kb * c_analyte_effective)
                        ph_val = compute_ph_from_oh(oh_val)
                    stage = "نقطة البداية (Initial Point)"
                elif v_tit < v_eq_ml - 1e-5:
                    if base_type == "strong":
                        rem_oh_mmol = eq_analyte_mmol - eq_tit_mmol
                        ph_val = compute_ph_from_oh(rem_oh_mmol / v_tot)
                    else:
                        equiv_progress = (eq_tit_mmol / max(1e-9, eq_analyte_mmol)) * val_analyte
                        step_idx = min(len(pkb_list) - 1, max(0, int(math.floor(equiv_progress))))
                        frac_in_step = max(0.01, min(0.99, equiv_progress - step_idx))
                        active_pkb = float(pkb_list[step_idx])
                        ph_val = 14.0 - (active_pkb + math.log10(frac_in_step / (1.0 - frac_in_step)))
                    stage = "منطقة المحلول المنظم / قبل التكافؤ (Buffer Region)"
                elif abs(v_tit - v_eq_ml) <= 1e-4:
                    ph_val = ph_eq
                    stage = "نقطة التكافؤ الكلي (Equivalence Point)"
                else:
                    excess_h_mmol = eq_tit_mmol - eq_analyte_mmol
                    if acid_type == "strong":
                        ph_val = compute_ph_from_h(excess_h_mmol / v_tot)
                    else:
                        c_excess = excess_h_mmol / v_tot
                        c_salt = eq_analyte_mmol / v_tot
                        ratio = max(1e-6, c_excess / max(1e-6, c_salt))
                        ph_val = acid_spec["pKa"] - math.log10(ratio)
                    stage = "فائض المحلول القياسي / بعد التكافؤ (Post-Equivalence)"

            curve_data.append({
                "volume_added_ml": round(v_tit, 2),
                "ph": round(max(0.0, min(14.0, ph_val)), 2),
                "stage_ar": stage,
            })

        best_indicator_key = select_optimal_indicator_key(ph_eq, case_key)
        indicator_info = INDICATORS_DATABASE[best_indicator_key]

        # Estimate colorimetric endpoint (where pH crosses midpoint of indicator range)
        ind_mid_ph = 0.5 * (indicator_info["pH_range"][0] + indicator_info["pH_range"][1])
        endpoint_pt = min(curve_data, key=lambda pt: abs(pt["ph"] - ind_mid_ph))

        explanation_ar = (
            f"عند معايرة **{analyte_formula}** ({c_analyte} M، حجم {v_analyte} mL، التكافؤ = {val_analyte}) "
            f"بواسطة **{titrant_formula}** ({c_titrant} M، التكافؤ = {val_titrant})، "
            f"يتحقق التكافؤ الكيميائي الكلي عند إضافة **{v_eq_ml:.2f} mL** من المحلول القياسي. "
            f"قيمة الرقم الهيدروجيني عند نقطة التكافؤ هي **pH = {ph_eq:.2f}**. "
            f"الكاشف اللوني الأنسب علمياً هو **{indicator_info['name_ar']}** لأن مدى تغير لونه "
            f"({indicator_info['pH_range'][0]} - {indicator_info['pH_range'][1]}) يتقاطع مع القفزة الهيدروجينية."
        )

        step_equivalence_points: list[dict[str, Any]] = []
        for s_idx, s_vol in enumerate(step_equivalence_volumes_ml, start=1):
            if s_idx < val_analyte:
                if mode == "acid_analyte_titrated_by_base" and s_idx < len(pka_list):
                    s_ph = round(0.5 * (float(pka_list[s_idx - 1]) + float(pka_list[s_idx])), 2)
                elif mode == "base_analyte_titrated_by_acid" and s_idx < len(pkb_list):
                    s_ph = round(14.0 - 0.5 * (float(pkb_list[s_idx - 1]) + float(pkb_list[s_idx])), 2)
                else:
                    s_ph = ph_eq
            else:
                s_ph = ph_eq
            half_v_s = round(v_eq1_ml * (s_idx - 0.5), 2)
            if mode == "acid_analyte_titrated_by_base" and acid_type == "weak":
                half_ph_s = round(float(pka_list[min(s_idx - 1, len(pka_list) - 1)]), 2)
            elif mode == "base_analyte_titrated_by_acid" and base_type == "weak":
                half_ph_s = round(14.0 - float(pkb_list[min(s_idx - 1, len(pkb_list) - 1)]), 2)
            else:
                half_ph_s = None
            step_ind_key = select_optimal_indicator_key(s_ph, case_key)
            step_equivalence_points.append({
                "step": s_idx,
                "volume_ml": s_vol,
                "ph": s_ph,
                "half_equivalence_volume_ml": half_v_s,
                "half_equivalence_ph": half_ph_s,
                "recommended_indicator_key": step_ind_key,
                "recommended_indicator_name_ar": INDICATORS_DATABASE[step_ind_key]["name_ar"],
            })

        # Quantitative Buffer Analysis (Van Slyke Buffer Capacity beta_max at half-equivalence)
        if ph_half_eq is not None:
            c_buf_half = (n_analyte_mmol / (v_analyte + 0.5 * v_eq1_ml))
            beta_max = round(0.576 * c_buf_half, 4)
            buf_range = (round(max(0.0, ph_half_eq - 1.0), 2), round(min(14.0, ph_half_eq + 1.0), 2))
            buffer_analysis = {
                "has_buffer_region": True,
                "half_equivalence_volume_ml": round(v_eq1_ml * 0.5, 2),
                "optimal_buffer_ph": ph_half_eq,
                "effective_buffer_range_ph": buf_range,
                "max_buffer_capacity_mol_l_ph": beta_max,
                "summary_ar": (
                    f"تتشكل منطقة محلول منظم (Buffer Region) فعالة في النطاق الهيدروجيني (pH {buf_range[0]} – {buf_range[1]}) "
                    f"وتبلغ أقصى سعة تنظيمية (Van Slyke Buffer Capacity β_max ≈ {beta_max} mol/L·pH) عند نصف التكافؤ "
                    f"(V = {v_eq1_ml * 0.5:.2f} mL، حيث pH = {ph_half_eq:.2f})."
                ),
            }
        else:
            buffer_analysis = {
                "has_buffer_region": False,
                "half_equivalence_volume_ml": None,
                "optimal_buffer_ph": None,
                "effective_buffer_range_ph": None,
                "max_buffer_capacity_mol_l_ph": 0.0,
                "summary_ar": "لا تتشكل منطقة محلول منظم كلاسيكية في معايرة حمض قوي مع قاعدة قوية لأن كلا الطرفين تام التأين.",
            }

        vol_diff_ml = round(float(endpoint_pt["volume_added_ml"]) - round(v_eq_ml, 2), 3)
        rel_err_pct = round(100.0 * abs(vol_diff_ml) / max(1e-6, v_eq_ml), 2)
        titration_error_analysis = {
            "equivalence_volume_ml": round(v_eq_ml, 2),
            "equivalence_ph": ph_eq,
            "endpoint_volume_ml": endpoint_pt["volume_added_ml"],
            "endpoint_ph": endpoint_pt["ph"],
            "volume_difference_ml": vol_diff_ml,
            "relative_error_pct": rel_err_pct,
            "summary_ar": (
                f"الفرق بين نقطة التكافؤ النظرية (V_eq = {v_eq_ml:.2f} mL عند pH = {ph_eq:.2f}) "
                f"ونقطة النهاية اللونية للكاشف (V_end = {endpoint_pt['volume_added_ml']:.2f} mL عند pH = {endpoint_pt['ph']:.2f}) "
                f"هو ΔV = {vol_diff_ml:+g} mL (خطأ معايرة نسبي = {rel_err_pct}%)."
            ),
            "distinction_explanation_ar": (
                "نقطة التكافؤ (Equivalence Point) هي النقطة الستوكيومترية النظرية التي تتساوى عندها مكافئات الحمض والقاعدة تماماً، "
                "بينما نقطة النهاية (Endpoint) هي النقطة التجريبية التي يتغير عندها لون الكاشف المرئي فعلياً."
            ),
        }

        step_by_step_calculation = [
            f"1) حساب كمية المادة في الدورق (Analyte Moles): n = C × V = {c_analyte:g} M × {v_analyte_raw:g} mL = {n_analyte_mmol:.4f} mmol (إجمالي المكافئات = {n_analyte_mmol * val_analyte:.4f} meq).",
            f"2) حساب حجم التكافؤ الأول وحجم التكافؤ الكلي: V_eq1 = {v_eq1_ml:.2f} mL | V_eq_total = (n × z_analyte) / (C_titrant × z_titrant) = {v_eq_ml:.2f} mL.",
            f"3) الرقم الهيدروجيني الابتدائي قبل الإضافة (V = 0 mL، حجم التخفيف المضاف = {v_dilution:g} mL): pH_0 = {curve_data[0]['ph']:.2f}.",
            (
                f"4) منطقة المحلول المنظم (Henderson-Hasselbalch): عند نصف التكافؤ الأول (V = {v_eq1_ml * 0.5:.2f} mL) يتساوى [HA] و [A⁻] فيكون pH = pKa = {ph_half_eq:.2f}."
                if ph_half_eq is not None
                else "4) قبل التكافؤ: يُحسب pH مباشرةً من فائض أيونات الهيدرونيوم أو الهيدروكسيد القوية مقسوماً على الحجم الكلي (V_analyte + V_added)."
            ),
            f"5) عند نقطة التكافؤ الكلي (V_eq = {v_eq_ml:.2f} mL، الحجم الكلي = {total_vol_at_eq:.2f} mL، تركيز الملح = {c_salt_at_eq:.4f} M): pH_eq = {ph_eq:.2f}.",
        ]

        model_limitations_ar = [
            "يفترض النموذج محاليل مائية مخففة مثالية عند 25 °C (معاملات الفعالية الأيونية γ± ≈ 1.0 و Kw = 1.0×10⁻¹⁴)؛ وفي التراكيز العالية (> 0.5 M) ينبغي تصحيح القوة الأيونية بمعادلة ديباي-هوكل (Debye-Hückel).",
            "في معايرة حمض ضعيف مع قاعدة ضعيفة (Weak Acid – Weak Base) تكون قفزة pH عند التكافؤ ضحلة جداً، فلا تعطي الكواشف اللونية البصرية نقطة نهاية حادة ويلزم استخدام مقياس pH إلكتروني (Potentiometric Titration).",
            "في الأحماض متعددة البروتونات ذات ثوابت التأين المتقاربة (ΔpKa < 3) أو الضعيفة جداً (pKa > 10 مثل البروتون الثالث لحمض الفوسفوريك pKa3 = 12.35)، تندمج قفزات التكافؤ أو يحجبها أثر التسوية للماء.",
            "للقواعد شحيحة الذوبان (مثل Mg(OH)₂)، يخضع التركيز الفعلي في المحلول المائي لاتزان حاصل الإذابة (Ksp) وليس للتفكك التام في التراكيز العالية.",
        ]

        return {
            "mode": mode,
            "case_key": case_key,
            "case_label_ar": case_label_ar,
            "is_polyprotic": is_polyprotic,
            "step_equivalence_volumes_ml": step_equivalence_volumes_ml,
            "analyte": {
                "formula": analyte_formula,
                "name_ar": acid_spec["name_ar"] if is_acid_analyte else base_spec["name_ar"],
                "concentration_m": c_analyte,
                "volume_ml": v_analyte,
                "valence": val_analyte,
                "amount_mmol": round(n_analyte_mmol, 4),
            },
            "titrant": {
                "formula": titrant_formula,
                "name_ar": base_spec["name_ar"] if is_acid_analyte else acid_spec["name_ar"],
                "concentration_m": c_titrant,
                "valence": val_titrant,
            },
            "equivalence_point": {
                "volume_ml": round(v_eq_ml, 2),
                "ph": ph_eq,
            },
            "endpoint": {
                "volume_ml": endpoint_pt["volume_added_ml"],
                "ph": endpoint_pt["ph"],
                "indicator_key": best_indicator_key,
            },
            "half_equivalence_point": {
                "volume_ml": round(v_eq1_ml * 0.5, 2),
                "ph": ph_half_eq,
            } if ph_half_eq is not None else None,
            "recommended_indicator_key": best_indicator_key,
            "recommended_indicator": indicator_info,
            "step_equivalence_points": step_equivalence_points,
            "buffer_analysis": buffer_analysis,
            "titration_error_analysis": titration_error_analysis,
            "step_by_step_calculation": step_by_step_calculation,
            "model_limitations_ar": model_limitations_ar,
            "dilution_water_ml": v_dilution,
            "curve_points": curve_data,
            "explanation_ar": explanation_ar,
        }

    @classmethod
    def build_plotly_titration_chart(
        cls,
        sim_result: dict[str, Any],
        indicator_key: str | None = None,
    ) -> go.Figure:
        curve = sim_result["curve_points"]
        v_list = [p["volume_added_ml"] for p in curve]
        ph_list = [p["ph"] for p in curve]
        v_eq = sim_result["equivalence_point"]["volume_ml"]
        ph_eq = sim_result["equivalence_point"]["ph"]

        ind_key = indicator_key or sim_result.get("recommended_indicator_key", "bromothymol_blue")
        ind = INDICATORS_DATABASE.get(ind_key, sim_result["recommended_indicator"])
        ind_low, ind_high = ind["pH_range"]
        ind_mid = 0.5 * (ind_low + ind_high)
        endpoint_pt = min(curve, key=lambda pt: abs(pt["ph"] - ind_mid))

        fig = go.Figure()

        fig.add_hrect(
            y0=ind_low,
            y1=ind_high,
            fillcolor=ind["color_hex"],
            opacity=0.16,
            line_width=1,
            line_dash="dash",
            line_color=ind["color_hex"],
            annotation_text=f"مدى تحول كاشف {ind['name_ar']} (pH {ind_low} - {ind_high})",
            annotation_position="bottom right",
        )

        fig.add_hline(y=7.0, line_dash="dot", line_color="#94a3b8", line_width=1.2)

        fig.add_trace(
            go.Scatter(
                x=v_list,
                y=ph_list,
                mode="lines",
                name="منحنى pH (Titration Curve)",
                line=dict(color="#0f4c5c", width=3.5),
                hoverinfo="text",
                hovertext=[
                    f"الحجم: {v} mL<br>pH: {ph}<br>{stg}"
                    for v, ph, stg in zip(v_list, ph_list, [p["stage_ar"] for p in curve])
                ],
            )
        )

        fig.add_trace(
            go.Scatter(
                x=[v_eq],
                y=[ph_eq],
                mode="markers+text",
                name="نقطة التكافؤ (Equivalence Point)",
                text=[f"نقطة التكافؤ<br>V = {v_eq} mL, pH = {ph_eq}"],
                textposition="top right",
                marker=dict(size=14, color="#DC2626", symbol="diamond", line=dict(color="#ffffff", width=2)),
            )
        )

        fig.add_trace(
            go.Scatter(
                x=[endpoint_pt["volume_added_ml"]],
                y=[endpoint_pt["ph"]],
                mode="markers+text",
                name="نقطة نهاية المعايرة (Indicator Endpoint)",
                text=[f"نقطة النهاية اللونية<br>V = {endpoint_pt['volume_added_ml']} mL, pH = {endpoint_pt['ph']}"],
                textposition="bottom left",
                marker=dict(size=12, color=ind["color_hex"], symbol="star", line=dict(color="#1e293b", width=1.5)),
            )
        )

        half_eq = sim_result.get("half_equivalence_point")
        if half_eq and half_eq.get("ph") is not None:
            v_half = half_eq["volume_ml"]
            ph_half = half_eq["ph"]
            fig.add_trace(
                go.Scatter(
                    x=[v_half],
                    y=[ph_half],
                    mode="markers+text",
                    name="نصف التكافؤ (pH = pKa1)",
                    text=[f"نصف التكافؤ<br>pH = {ph_half}"],
                    textposition="bottom right",
                    marker=dict(size=11, color="#059669", symbol="circle", line=dict(color="#ffffff", width=2)),
                )
            )

        titrant_label = sim_result["titrant"]["name_ar"]
        titrant_form = sim_result["titrant"]["formula"]

        fig.update_layout(
            title=dict(
                text=f"منحنى المعايرة الكيميائي: {sim_result['case_label_ar']}",
                font=dict(family="Cairo, sans-serif", size=16, color="#0f4c5c"),
                x=0.5,
            ),
            xaxis=dict(
                title=f"حجم المحلول القياسي المضاف ({titrant_label} — {titrant_form}) بوحدة mL",
                gridcolor="#f1f5f9",
                zeroline=False,
            ),
            yaxis=dict(
                title="الرقم الهيدروجيني (pH)",
                range=[0, 14],
                gridcolor="#f1f5f9",
                zeroline=False,
            ),
            plot_bgcolor="#ffffff",
            paper_bgcolor="#ffffff",
            hovermode="closest",
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
            margin=dict(l=40, r=40, t=60, b=40),
            height=480,
        )

        return fig


class TitrationEngine:
    """Unified engine for Titration views and test suites."""

    ACID_DATABASE = ACID_SPECIES_DATABASE
    BASE_DATABASE = BASE_SPECIES_DATABASE
    INDICATORS = {
        k: {
            **v,
            "transition_range": v.get("pH_range", (7.0, 7.0)),
        }
        for k, v in INDICATORS_DATABASE.items()
    }

    @classmethod
    def simulate_titration(
        cls,
        acid_key: str,
        base_key: str,
        acid_conc_m: float,
        acid_volume_ml: float,
        base_conc_m: float,
        titration_mode: str = "acid_with_base",
        base_volume_ml: float | None = None,
        dilution_water_ml: float = 0.0,
    ) -> dict[str, Any]:
        if acid_key not in cls.ACID_DATABASE:
            raise ValueError(f"Unknown acid key: {acid_key}")
        if base_key not in cls.BASE_DATABASE:
            raise ValueError(f"Unknown base key: {base_key}")
        if acid_conc_m <= 0 or acid_volume_ml <= 0 or base_conc_m <= 0:
            raise ValueError("Concentrations and volumes must be strictly positive.")

        if titration_mode == "base_with_acid":
            v_base = base_volume_ml if base_volume_ml is not None and base_volume_ml > 0 else acid_volume_ml
            sim = TitrationSimulator.simulate_titration(
                analyte_formula=base_key,
                titrant_formula=acid_key,
                analyte_volume_ml=v_base,
                analyte_conc_m=base_conc_m,
                titrant_conc_m=acid_conc_m,
                dilution_water_ml=dilution_water_ml,
            )
        else:
            sim = TitrationSimulator.simulate_titration(
                analyte_formula=acid_key,
                titrant_formula=base_key,
                analyte_volume_ml=acid_volume_ml,
                analyte_conc_m=acid_conc_m,
                titrant_conc_m=base_conc_m,
                dilution_water_ml=dilution_water_ml,
            )

        eq = sim["equivalence_point"]
        endpt = sim["endpoint"]
        half_eq = sim.get("half_equivalence_point")
        curve = sim["curve_points"]
        v_list = [p["volume_added_ml"] for p in curve]
        ph_list = [p["ph"] for p in curve]

        return {
            "success": True,
            "v_eq_ml": eq["volume_ml"],
            "ph_eq": eq["ph"],
            "endpoint_volume_ml": endpt["volume_ml"],
            "endpoint_ph": endpt["ph"],
            "step_equivalence_volumes_ml": sim["step_equivalence_volumes_ml"],
            "step_equivalence_points": sim["step_equivalence_points"],
            "buffer_analysis": sim["buffer_analysis"],
            "titration_error_analysis": sim["titration_error_analysis"],
            "step_by_step_calculation": sim["step_by_step_calculation"],
            "model_limitations_ar": sim["model_limitations_ar"],
            "dilution_water_ml": sim["dilution_water_ml"],
            "is_polyprotic": sim["is_polyprotic"],
            "initial_ph": curve[0]["ph"] if curve else 7.0,
            "half_v_eq_ml": half_eq["volume_ml"] if half_eq else None,
            "half_ph": half_eq["ph"] if half_eq else None,
            "volumes_ml": v_list,
            "ph_values": ph_list,
            "sim_result": sim,
            "recommended_indicator_key": sim.get("recommended_indicator_key", "bromothymol_blue"),
            "recommended_indicator": sim.get("recommended_indicator"),
        }

    @classmethod
    def calculate_point_ph(
        cls,
        v_add: float,
        v_eq: float,
        v_acid: float,
        ca: float,
        cb: float,
        acid_info: dict[str, Any],
        base_info: dict[str, Any],
        titration_mode: str = "acid_with_base",
    ) -> dict[str, Any]:
        if titration_mode == "base_with_acid":
            sim = TitrationSimulator.simulate_titration(
                analyte_formula=base_info["formula"],
                titrant_formula=acid_info["formula"],
                analyte_volume_ml=v_acid,
                analyte_conc_m=cb,
                titrant_conc_m=ca,
                max_titrant_ratio=max(2.2, (v_add / v_eq) + 0.2 if v_eq > 0 else 2.2),
            )
        else:
            sim = TitrationSimulator.simulate_titration(
                analyte_formula=acid_info["formula"],
                titrant_formula=base_info["formula"],
                analyte_volume_ml=v_acid,
                analyte_conc_m=ca,
                titrant_conc_m=cb,
                max_titrant_ratio=max(2.2, (v_add / v_eq) + 0.2 if v_eq > 0 else 2.2),
            )
        curve = sim["curve_points"]
        closest = min(curve, key=lambda p: abs(p["volume_added_ml"] - v_add))
        return {
            "ph": closest["ph"],
            "region": closest["stage_ar"],
        }

    @classmethod
    def create_plotly_curve(
        cls,
        sim_data: dict[str, Any],
        current_v: float = None,
        indicator_key: str = None,
    ) -> go.Figure:
        sim_result = sim_data.get("sim_result")
        if not sim_result:
            return go.Figure()

        fig = TitrationSimulator.build_plotly_titration_chart(sim_result, indicator_key=indicator_key)

        if current_v is not None:
            v_list = sim_data["volumes_ml"]
            ph_list = sim_data["ph_values"]
            closest_idx = min(range(len(v_list)), key=lambda i: abs(v_list[i] - current_v))
            cur_ph = ph_list[closest_idx]
            fig.add_trace(
                go.Scatter(
                    x=[current_v],
                    y=[cur_ph],
                    mode="markers+text",
                    name="الموضع الحالي (Current)",
                    text=[f"الحالي: {current_v:.1f} mL, pH {cur_ph:.2f}"],
                    textposition="bottom right",
                    marker=dict(size=13, color="#f59e0b", symbol="circle", line=dict(color="#ffffff", width=2)),
                )
            )

        return fig

    @classmethod
    def create_flask_svg(
        cls,
        current_ph: float,
        indicator_key: str,
        current_vol_ml: float,
        v_eq_ml: float,
        analyte_name_ar: str = "",
        titrant_name_ar: str = "",
    ) -> str:
        """
        Generates an interactive SVG laboratory titration apparatus (Burette + Erlenmeyer Flask)
        showing dynamic solution color change inside the flask as titrant is added.
        """
        color_state = get_indicator_state_at_ph(indicator_key, current_ph)
        ind = INDICATORS_DATABASE.get(indicator_key, INDICATORS_DATABASE["bromothymol_blue"])
        fill_hex = color_state["hex"]
        stroke_hex = "#334155" if fill_hex.upper() == "#F8FAFC" else fill_hex
        progress_ratio = min(1.0, max(0.0, current_vol_ml / max(1.0, v_eq_ml * 2.0)))
        liquid_top_y = int(205 - progress_ratio * 22)

        return f"""
        <svg width="100%" height="295" viewBox="0 0 540 295" xmlns="http://www.w3.org/2000/svg" style="background:#f8fafc;border-radius:14px;border:1px solid #cbd5e1;">
            <!-- Header -->
            <rect x="14" y="10" width="512" height="32" rx="8" fill="#0f4c5c"/>
            <text x="270" y="31" font-family="'Cairo', sans-serif" font-size="13" font-weight="bold" fill="#ffffff" text-anchor="middle">
                محاكاة دورق المعايرة والسحاحة — الكاشف: {html.escape(ind['name_ar'])}
            </text>

            <!-- Burette at Top -->
            <rect x="138" y="48" width="24" height="72" rx="3" fill="#e0f2fe" stroke="#1e293b" stroke-width="2"/>
            <line x1="138" y1="65" x2="148" y2="65" stroke="#475569" stroke-width="1.2"/>
            <line x1="138" y1="80" x2="150" y2="80" stroke="#475569" stroke-width="1.2"/>
            <line x1="138" y1="95" x2="148" y2="95" stroke="#475569" stroke-width="1.2"/>
            <!-- Stopcock -->
            <rect x="131" y="120" width="38" height="7" rx="3" fill="#dc2626"/>
            <line x1="150" y1="127" x2="150" y2="142" stroke="#1e293b" stroke-width="2.5"/>
            <!-- Falling Titrant Drop -->
            <circle cx="150" cy="154" r="4.5" fill="#38bdf8" stroke="#0284c7" stroke-width="1"/>

            <!-- Erlenmeyer Flask Liquid (Dynamic Indicator Color) -->
            <path d="M 112 {liquid_top_y} L 188 {liquid_top_y} L 214 252 Q 150 260 86 252 Z" fill="{fill_hex}" stroke="{stroke_hex}" stroke-width="1.5" opacity="0.88"/>

            <!-- Erlenmeyer Flask Glass Outline -->
            <path d="M 135 142 L 135 168 L 82 252 Q 150 262 218 252 L 165 168 L 165 142" fill="none" stroke="#0f172a" stroke-width="3" stroke-linecap="round"/>
            <ellipse cx="150" cy="142" rx="15" ry="4" fill="none" stroke="#0f172a" stroke-width="2.2"/>

            <!-- Right Info Panel -->
            <rect x="250" y="56" width="272" height="205" rx="12" fill="#ffffff" stroke="#cbd5e1" stroke-width="1.5"/>
            <text x="386" y="82" font-family="'Cairo', sans-serif" font-size="13" font-weight="bold" fill="#0f4c5c" text-anchor="middle">
                الحالة اللونية والهيدروجينية داخل الدورق
            </text>
            <rect x="275" y="94" width="222" height="32" rx="8" fill="{fill_hex}" stroke="#475569" stroke-width="1.5"/>
            <text x="386" y="145" font-family="'Cairo', sans-serif" font-size="12" font-weight="bold" fill="#1e293b" text-anchor="middle">
                {html.escape(color_state['label_ar'])}
            </text>
            <text x="386" y="170" font-family="'Cairo', sans-serif" font-size="12" fill="#0369a1" text-anchor="middle">
                الحجم المضاف: {current_vol_ml:.2f} mL | التكافؤ: {v_eq_ml:.2f} mL
            </text>
            <text x="386" y="194" font-family="'Cairo', sans-serif" font-size="13" font-weight="bold" fill="#0f766e" text-anchor="middle">
                الرقم الهيدروجيني الحالي: pH = {current_ph:.2f}
            </text>
            <text x="386" y="218" font-family="'Cairo', sans-serif" font-size="11" fill="#64748b" text-anchor="middle">
                مدى تحول الكاشف: pH {ind['pH_range'][0]} – {ind['pH_range'][1]}
            </text>
            <text x="386" y="240" font-family="'Cairo', sans-serif" font-size="11" fill="#475569" text-anchor="middle">
                حمضي: {html.escape(ind['color_acid'])} | قاعدي: {html.escape(ind['color_base'])}
            </text>
        </svg>
        """

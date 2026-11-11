import re
from typing import Any
from chemistry_engine.formula_analyzer import (
    formulas_are_equivalent,
    get_standard_phase_code,
    normalize_chemical_formula,
    replace_subscripts,
)

VERIFIED_REACTIONS_REGISTRY: list[dict[str, Any]] = [
    {
        "reaction_id": "rxn_water_synthesis",
        "domain": "general_chemistry",
        "domain_label_ar": "الكيمياء العامة (General Chemistry)",
        "reaction_name_ar": "تفاعل اتحاد الهيدروجين والأكسجين لتكوين الماء",
        "reaction_name_en": "Synthesis of Water from Hydrogen and Oxygen",
        "reaction_type_ar": "تفاعل اتحاد مباشر وأكسدة-اختزال (Synthesis & Redox Reaction)",
        "reactants": [
            {"formula": "H2", "name_ar": "غاز الهيدروجين", "name_en": "Hydrogen", "coefficient": 2, "state": "g"},
            {"formula": "O2", "name_ar": "غاز الأكسجين", "name_en": "Oxygen", "coefficient": 1, "state": "g"},
        ],
        "products": [
            {"formula": "H2O", "name_ar": "الماء", "name_en": "Water", "coefficient": 2, "state": "l"},
        ],
        "conditions": {
            "temperature": "شرارة بدء التفاعل أو التسخين (> 500 °C) أو وجود حفاز البلاتين (Pt)",
            "pressure": "الضغط الجوي القياسي (1 atm)",
            "catalyst": "بلاتين (Pt) اختياري لخفض طاقة التنشيط",
            "thermodynamics_ar": "تفاعل طارد للحرارة بشدة (ΔH° = -571.6 kJ لكل 2 mol من الماء السائل)",
        },
        "explanation_ar": (
            "يتفاعل جزيئان من غاز الهيدروجين (H2) مع جزيء واحد من غاز الأكسجين (O2) بنسبة مولية ثابتة (2:1) "
            "لتكوين جزيئين من الماء (H2O). تنكسر الروابط (H–H) و (O=O) وتتكون روابط تساهمية قطبية جديدة (O–H) أكثر استقرارا."
        ),
        "assumptions_and_limits_ar": (
            "يفترض هذا الحساب النظري اكتمال التفاعل بنسبة 100% في نظام مغلق دون فقدان حراري أو تفاعلات جانبية. "
            "هذه نتيجة حسابية نظرية وليست إثباتا لتجربة عملية مخبرية."
        ),
        "source": "IUPAC / NIST Chemistry WebBook — Thermochemical Kinetics (H2 + O2 Synthesis)",
    },
    {
        "reaction_id": "rxn_iron_oxidation",
        "domain": "general_chemistry",
        "domain_label_ar": "الكيمياء العامة (General Chemistry)",
        "reaction_name_ar": "أكسدة الحديد لتكوين أكسيد الحديد الثلاثي",
        "reaction_name_en": "Oxidation of Iron to Iron(III) Oxide",
        "reaction_type_ar": "تفاعل اتحاد وأكسدة-اختزال (Oxidation-Reduction)",
        "reactants": [
            {"formula": "Fe", "name_ar": "الحديد", "name_en": "Iron", "coefficient": 4, "state": "s"},
            {"formula": "O2", "name_ar": "غاز الأكسجين", "name_en": "Oxygen", "coefficient": 3, "state": "g"},
        ],
        "products": [
            {"formula": "Fe2O3", "name_ar": "أكسيد الحديد الثلاثي", "name_en": "Iron(III) Oxide", "coefficient": 2, "state": "s"},
        ],
        "conditions": {
            "temperature": "درجة حرارة الغرفة بوجود الرطوبة (بطيء) أو التسخين في الهواء (سريع)",
            "pressure": "1 atm",
            "catalyst": "بخار الماء يسرع عملية التأكسد الكهروكيميائية",
            "thermodynamics_ar": "تفاعل طارد للحرارة (ΔH° = -1648.4 kJ لتكوين 2 mol من Fe2O3)",
        },
        "explanation_ar": (
            "تتأكسد 4 مولات من ذرات الحديد (Fe) بفقد إلكترونات تكافئها إلى 3 مولات من جزيئات الأكسجين (O2) "
            "لتكوين 2 مول من أكسيد الحديد الثلاثي الصلب (Fe2O3) وفق النسبة المولية 4:3:2."
        ),
        "assumptions_and_limits_ar": "حساب ستوكيوميتري نظري بفرض تحول الحديد بالكامل إلى حالة الأكسدة +3.",
        "source": "CRC Handbook of Chemistry and Physics — Inorganic Metal Oxidation",
    },
    {
        "reaction_id": "rxn_acid_base_hcl_naoh",
        "domain": "general_chemistry",
        "domain_label_ar": "الكيمياء العامة (General Chemistry)",
        "reaction_name_ar": "تفاعل تعادل حمض الهيدروكلوريك مع هيدروكسيد الصوديوم",
        "reaction_name_en": "Neutralization of Hydrochloric Acid and Sodium Hydroxide",
        "reaction_type_ar": "تفاعل تعادل حمض قوي مع قاعدة قوية (Acid-Base Neutralization)",
        "reactants": [
            {"formula": "HCl", "name_ar": "حمض الهيدروكلوريك", "name_en": "Hydrochloric Acid", "coefficient": 1, "state": "aq"},
            {"formula": "NaOH", "name_ar": "هيدروكسيد الصوديوم", "name_en": "Sodium Hydroxide", "coefficient": 1, "state": "aq"},
        ],
        "products": [
            {"formula": "NaCl", "name_ar": "كلوريد الصوديوم", "name_en": "Sodium Chloride", "coefficient": 1, "state": "aq"},
            {"formula": "H2O", "name_ar": "الماء", "name_en": "Water", "coefficient": 1, "state": "l"},
        ],
        "conditions": {
            "temperature": "الظروف القياسية (25 °C)",
            "pressure": "1 atm",
            "solvent": "محلول مائي (Aqueous Solution)",
            "thermodynamics_ar": "تفاعل طارد للحرارة (حرارة التعادل القياسية ΔH° = -57.3 kJ/mol)",
        },
        "explanation_ar": (
            "يتحد أيون الهيدروجين (H⁺) من الحمض مع أيون الهيدروكسيد (OH⁻) من القاعدة بنسبة مولية متكافئة (1:1) "
            "لتكوين جزيء الماء (H2O) وملح كلوريد الصوديوم الذائب (NaCl)."
        ),
        "assumptions_and_limits_ar": "حساب نظري لتفاعل تام في الوسط المائي المخفف.",
        "source": "USP / Ph. Eur. Volumetric Acid-Base Titration Reference & General Chemistry",
    },
    {
        "reaction_id": "rxn_haber_ammonia",
        "domain": "general_chemistry",
        "domain_label_ar": "الكيمياء العامة (General Chemistry)",
        "reaction_name_ar": "تصنيع الأمونيا من النيتروجين والهيدروجين (طريقة هابر-بوش)",
        "reaction_name_en": "Haber-Bosch Synthesis of Ammonia",
        "reaction_type_ar": "تفاعل اتحاد غازي محفز (Catalytic Synthesis)",
        "reactants": [
            {"formula": "N2", "name_ar": "غاز النيتروجين", "name_en": "Nitrogen", "coefficient": 1, "state": "g"},
            {"formula": "H2", "name_ar": "غاز الهيدروجين", "name_en": "Hydrogen", "coefficient": 3, "state": "g"},
        ],
        "products": [
            {"formula": "NH3", "name_ar": "الأمونيا", "name_en": "Ammonia", "coefficient": 2, "state": "g"},
        ],
        "conditions": {
            "temperature": "400 – 500 °C",
            "pressure": "150 – 250 atm",
            "catalyst": "حفاز الحديد المنشط (Fe with K2O/Al2O3 promoters)",
            "thermodynamics_ar": "تفاعل طارد للحرارة عكسي (ΔH° = -92.4 kJ لكل 2 mol من NH3)",
        },
        "explanation_ar": (
            "يتفاعل مول واحد من غاز النيتروجين (N2) مع 3 مولات من غاز الهيدروجين (H2) لتكوين مولين من غاز الأمونيا (NH3) "
            "بنسبة مولية 1:3:2."
        ),
        "assumptions_and_limits_ar": "يمثل الناتج المحسوب الحد الأقصى النظري (100% تحويل) قبل الوصول إلى حالة الاتزان العكسي.",
        "source": "Ullmann's Encyclopedia of Industrial Chemistry — Ammonia Haber-Bosch Process",
    },
    {
        "reaction_id": "rxn_h2o2_decomposition",
        "domain": "general_chemistry",
        "domain_label_ar": "الكيمياء العامة (General Chemistry)",
        "reaction_name_ar": "تفكك فوق أكسيد الهيدروجين إلى ماء وأكسجين",
        "reaction_name_en": "Catalytic Decomposition of Hydrogen Peroxide",
        "reaction_type_ar": "تفاعل تفكك وعدم تناسب (Decomposition / Disproportionation)",
        "reactants": [
            {"formula": "H2O2", "name_ar": "فوق أكسيد الهيدروجين", "name_en": "Hydrogen Peroxide", "coefficient": 2, "state": "aq"},
        ],
        "products": [
            {"formula": "H2O", "name_ar": "الماء", "name_en": "Water", "coefficient": 2, "state": "l"},
            {"formula": "O2", "name_ar": "غاز الأكسجين", "name_en": "Oxygen", "coefficient": 1, "state": "g"},
        ],
        "conditions": {
            "temperature": "25 °C",
            "catalyst": "ثاني أكسيد المنغنيز (MnO2) أو إنزيم الكاتالاز الحيوي",
            "thermodynamics_ar": "تفاعل طارد للحرارة (ΔH° = -196.1 kJ لكل 2 mol من H2O2)",
        },
        "explanation_ar": "يتفكك جزيئان من فوق أكسيد الهيدروجين (H2O2) لإنتاج جزيئين من الماء (H2O) وجزيء واحد من غاز الأكسجين (O2).",
        "assumptions_and_limits_ar": "حساب نظري بفرض التفكك التام للمتفاعل.",
        "source": "Rhea Database (RHEA:15889 Catalase) & Inorganic Chemistry Reference",
    },
    {
        "reaction_id": "rxn_co2_water_hydration",
        "domain": "general_chemistry",
        "domain_label_ar": "الكيمياء العامة (General Chemistry)",
        "reaction_name_ar": "تفاعل ثاني أكسيد الكربون مع الماء لتكوين حمض الكربونيك",
        "reaction_name_en": "Hydration of Carbon Dioxide to Carbonic Acid",
        "reaction_type_ar": "تفاعل إماهة أكسيد لا فلزي (Hydration Reaction)",
        "reactants": [
            {"formula": "CO2", "name_ar": "ثاني أكسيد الكربون", "name_en": "Carbon Dioxide", "coefficient": 1, "state": "g"},
            {"formula": "H2O", "name_ar": "الماء", "name_en": "Water", "coefficient": 1, "state": "l"},
        ],
        "products": [
            {"formula": "H2CO3", "name_ar": "حمض الكربونيك", "name_en": "Carbonic Acid", "coefficient": 1, "state": "aq"},
        ],
        "conditions": {
            "temperature": "25 °C",
            "pressure": "1 atm أو ضغط مرتفع لزيادة الذوبانية",
            "catalyst": "في الأنظمة الحيوية يحفزه إنزيم الأنهيدراز الكربوني (Carbonic Anhydrase)",
        },
        "explanation_ar": "يتفاعل مول واحد من ثاني أكسيد الكربون (CO2) مع مول واحد من الماء (H2O) لتكوين مول واحد من حمض الكربونيك (H2CO3).",
        "assumptions_and_limits_ar": "يمثل الحساب القيمة النظرية الكاملة لتحول CO2 الذائب إلى H2CO3.",
        "source": "Rhea Database (RHEA:16833 Carbonic Anhydrase) & Physiological Buffer Chemistry",
    },
    {
        "reaction_id": "rxn_methane_combustion",
        "domain": "organic_chemistry",
        "domain_label_ar": "الكيمياء العضوية (Organic Chemistry)",
        "reaction_name_ar": "الاحتراق الكامل لغاز الميثان",
        "reaction_name_en": "Complete Combustion of Methane",
        "reaction_type_ar": "تفاعل احتراق هيدروكربوني كامل (Complete Hydrocarbon Combustion)",
        "reactants": [
            {"formula": "CH4", "name_ar": "الميثان", "name_en": "Methane", "coefficient": 1, "state": "g"},
            {"formula": "O2", "name_ar": "غاز الأكسجين", "name_en": "Oxygen", "coefficient": 2, "state": "g"},
        ],
        "products": [
            {"formula": "CO2", "name_ar": "ثاني أكسيد الكربون", "name_en": "Carbon Dioxide", "coefficient": 1, "state": "g"},
            {"formula": "H2O", "name_ar": "الماء", "name_en": "Water", "coefficient": 2, "state": "l"},
        ],
        "conditions": {
            "temperature": "شرارة اشتعال في وفرة من الأكسجين",
            "pressure": "1 atm",
            "thermodynamics_ar": "تفاعل طارد للحرارة بشدة (ΔH°c = -890.3 kJ/mol)",
        },
        "explanation_ar": "يحترق مول واحد من الميثان (CH4) بوجود مولين من غاز الأكسجين (O2) لينتج مولا واحدا من ثاني أكسيد الكربون (CO2) ومولين من الماء (H2O).",
        "assumptions_and_limits_ar": "يفترض الحساب احتراقا كاملا دون تكوين أول أكسيد الكربون (CO) أو السخام.",
        "source": "NIST Chemistry WebBook — Standard Enthalpy of Hydrocarbon Combustion",
    },
    {
        "reaction_id": "rxn_glucose_respiration",
        "domain": "biochemistry",
        "domain_label_ar": "الكيمياء الحيوية (Biochemistry)",
        "reaction_name_ar": "الأكسدة الهوائية الكاملة للجلوكوز (التنفس الخلوي)",
        "reaction_name_en": "Aerobic Cellular Oxidation of D-Glucose",
        "reaction_type_ar": "أكسدة حيوية متعددة المراحل (Biochemical Catabolic Oxidation)",
        "reactants": [
            {"formula": "C6H12O6", "name_ar": "الجلوكوز", "name_en": "D-Glucose", "coefficient": 1, "state": "aq"},
            {"formula": "O2", "name_ar": "الأكسجين الجزيئي", "name_en": "Oxygen", "coefficient": 6, "state": "g"},
        ],
        "products": [
            {"formula": "CO2", "name_ar": "ثاني أكسيد الكربون", "name_en": "Carbon Dioxide", "coefficient": 6, "state": "g"},
            {"formula": "H2O", "name_ar": "الماء", "name_en": "Water", "coefficient": 6, "state": "l"},
        ],
        "conditions": {
            "temperature": "درجة الحرارة الفسيولوجية (37 °C)",
            "catalyst": "سلسلة إنزيمات التحلل السكري ودورة كريبس والفسفرة التأكسدية",
            "thermodynamics_ar": "تفاعل منتج للطاقة الحرة القياسية (ΔG°' = -2870 kJ/mol)",
        },
        "explanation_ar": "يتأكسد مول واحد من الجلوكوز (C6H12O6) بوجود 6 مولات من الأكسجين (O2) لينتج 6 مولات من ثاني أكسيد الكربون (CO2) و6 مولات من الماء (H2O).",
        "assumptions_and_limits_ar": "تمثل المعادلة المحصلة الستوكيومترية الكلية للتنفس الخلوي الهوائي.",
        "source": "Lehninger Principles of Biochemistry & Rhea Biochemical Pathways",
    },
]

REACTION_CATEGORY_BY_ID = {
    "rxn_water_synthesis": "redox_synthesis",
    "rxn_iron_oxidation": "redox_synthesis",
    "rxn_acid_base_hcl_naoh": "neutralization",
    "rxn_haber_ammonia": "redox_synthesis",
    "rxn_h2o2_decomposition": "decomposition",
    "rxn_co2_water_hydration": "redox_synthesis",
    "rxn_methane_combustion": "combustion",
    "rxn_glucose_respiration": "biochemical_reaction",
}

CATEGORY_LABELS_AR = {
    "neutralization": "تفاعل تعادل وأحماض وقواعد (Acid-Base Neutralization)",
    "combustion": "تفاعل احتراق (Combustion Reaction)",
    "precipitation": "تفاعل ترسيب وتبادل مزدوج (Precipitation & Metathesis)",
    "single_displacement": "تفاعل إحلال مفرد (Single Displacement)",
    "redox_synthesis": "تفاعل اتحاد وأكسدة واختزال (Synthesis & Redox)",
    "decomposition": "تفاعل تحلل وتفكك (Decomposition)",
    "organic_reaction": "تفاعل عضوي (Organic Reaction)",
    "biochemical_reaction": "تفاعل كيميائي حيوي (Biochemical Reaction)",
}

DOMAIN_BY_CATEGORY = {
    "neutralization": ("acid_base_chemistry", "الأحماض والقواعد والتعادل (Acids, Bases & Neutralization)"),
    "precipitation": ("redox_and_precipitation", "الأكسدة والاختزال والترسيب والإحلال (Redox, Precipitation & Displacement)"),
    "single_displacement": ("redox_and_precipitation", "الأكسدة والاختزال والترسيب والإحلال (Redox, Precipitation & Displacement)"),
    "redox_synthesis": ("redox_and_precipitation", "الأكسدة والاختزال والترسيب والإحلال (Redox, Precipitation & Displacement)"),
    "combustion": ("organic_chemistry", "الكيمياء العضوية (Organic Chemistry - ORD)"),
    "organic_reaction": ("organic_chemistry", "الكيمياء العضوية (Organic Chemistry - ORD)"),
    "decomposition": ("general_chemistry", "الكيمياء العامة وغير العضوية (General & Inorganic Chemistry)"),
    "biochemical_reaction": ("biochemistry", "الكيمياء الحيوية (Biochemistry - Rhea)"),
}

CATEGORY_DEFAULT_SOURCES = {
    "neutralization": "USP / Ph. Eur. Acid-Base Volumetric Standards & IUPAC Analytical Compendium",
    "precipitation": "Vogel's Textbook of Qualitative Inorganic Analysis & CRC Solubility Rules",
    "single_displacement": "CRC Handbook of Chemistry and Physics — Electrochemical Series",
    "redox_synthesis": "Greenwood & Earnshaw, Chemistry of the Elements (Inorganic Redox Reference)",
    "combustion": "NIST Chemistry WebBook — Standard Thermochemical Combustion Data",
    "decomposition": "Inorganic Thermochemical Decomposition Reference (CRC / NIST)",
    "organic_reaction": "Open Reaction Database (ORD) & March's Advanced Organic Chemistry",
    "biochemical_reaction": "Rhea Biochemical Reactions Database & Lehninger Biochemistry",
}

EXTRA_SPECIES_NAMES = {
    "NaNO3": ("نترات الصوديوم", "Sodium Nitrate"),
    "CuCl2": ("كلوريد النحاس الثنائي", "Copper(II) Chloride"),
    "Fe(OH)3": ("هيدروكسيد الحديديك", "Iron(III) Hydroxide"),
    "Cu(OH)2": ("هيدروكسيد النحاس الثنائي", "Copper(II) Hydroxide"),
    "ZnCl2": ("كلوريد الزنك", "Zinc Chloride"),
    "FeCl2": ("كلوريد الحديدوز", "Iron(II) Chloride"),
    "Cu(NO3)2": ("نترات النحاس الثنائي", "Copper(II) Nitrate"),
    "AlCl3": ("كلوريد الألومنيوم", "Aluminium Chloride"),
    "K2MnO4": ("منغنات البوتاسيوم", "Potassium Manganate"),
    "MnO2": ("ثاني أكسيد المنغنيز", "Manganese Dioxide"),
    "C2H4Br2": ("1،2-ثنائي برومو الإيثان", "1,2-Dibromoethane"),
    "CH3Cl": ("كلورو الميثان", "Chloromethane"),
    "C6H5ONa": ("فينوكسيد الصوديوم", "Sodium Phenoxide"),
    "C6H5NH3Cl": ("كلوريد الأنيلينيوم", "Anilinium Chloride"),
    "C10H15N5O10P2": ("أدينوسين ثنائي الفوسفات", "Adenosine Diphosphate"),
}

# category, reactants, products, conditions_ar
_GENERATED_REACTIONS: list[tuple[str, str, str, str]] = [
    ("neutralization", "CH3COOH + NaOH", "CH3COONa + H2O", "محلول مائي، 25 °C"),
    ("neutralization", "CH3COOH + NaHCO3", "CH3COONa + H2O + CO2", "محلول مائي، تصاعد غاز CO2"),
    ("neutralization", "H2SO4 + 2NaOH", "Na2SO4 + 2H2O", "محلول مائي، 25 °C"),
    ("neutralization", "HNO3 + KOH", "KNO3 + H2O", "محلول مائي، 25 °C"),
    ("neutralization", "2HCl + Ca(OH)2", "CaCl2 + 2H2O", "محلول مائي، 25 °C"),
    ("neutralization", "2HCl + Mg(OH)2", "MgCl2 + 2H2O", "محلول مائي، 25 °C"),
    ("neutralization", "H3PO4 + 3NaOH", "Na3PO4 + 3H2O", "محلول مائي، 25 °C"),
    ("neutralization", "HCl + NH3", "NH4Cl", "طور غازي أو مائي، 25 °C"),
    ("neutralization", "HCl + NaHCO3", "NaCl + H2O + CO2", "محلول مائي، تصاعد غاز CO2"),
    ("neutralization", "2HCl + Na2CO3", "2NaCl + H2O + CO2", "محلول مائي، تصاعد غاز CO2"),
    ("neutralization", "2HCl + CaCO3", "CaCl2 + H2O + CO2", "محلول مائي، تصاعد غاز CO2"),
    ("neutralization", "H2SO4 + Ba(OH)2", "BaSO4 + 2H2O", "محلول مائي، يتكون راسب BaSO4"),
    ("neutralization", "HCl + KOH", "KCl + H2O", "محلول مائي، 25 °C"),
    ("redox_synthesis", "CaO + H2O", "Ca(OH)2", "تفاعل طارد للحرارة، 25 °C"),
    ("redox_synthesis", "CO2 + Ca(OH)2", "CaCO3 + H2O", "محلول مائي، يعكر ماء الجير"),
    ("redox_synthesis", "SO3 + H2O", "H2SO4", "إضافة حذرة للماء"),
    ("combustion", "2C2H6 + 7O2", "4CO2 + 6H2O", "شرارة اشتعال في وفرة من الأكسجين"),
    ("combustion", "C3H8 + 5O2", "3CO2 + 4H2O", "شرارة اشتعال في وفرة من الأكسجين"),
    ("combustion", "2C4H10 + 13O2", "8CO2 + 10H2O", "شرارة اشتعال في وفرة من الأكسجين"),
    ("combustion", "C2H4 + 3O2", "2CO2 + 2H2O", "شرارة اشتعال في وفرة من الأكسجين"),
    ("combustion", "2C2H2 + 5O2", "4CO2 + 2H2O", "شرارة اشتعال في وفرة من الأكسجين"),
    ("combustion", "C2H5OH + 3O2", "2CO2 + 3H2O", "شرارة اشتعال في وفرة من الأكسجين"),
    ("combustion", "2CH3OH + 3O2", "2CO2 + 4H2O", "شرارة اشتعال في وفرة من الأكسجين"),
    ("combustion", "2C6H6 + 15O2", "12CO2 + 6H2O", "اشتعال في وفرة من الأكسجين"),
    ("combustion", "2H2S + 3O2", "2SO2 + 2H2O", "اشتعال في وفرة من الأكسجين"),
    ("redox_synthesis", "C + O2", "CO2", "تسخين في وفرة من الأكسجين"),
    ("precipitation", "AgNO3 + NaCl", "AgCl + NaNO3", "محلول مائي، راسب أبيض من AgCl"),
    ("precipitation", "AgNO3 + KCl", "AgCl + KNO3", "محلول مائي، راسب أبيض من AgCl"),
    ("precipitation", "AgNO3 + HCl", "AgCl + HNO3", "محلول مائي، راسب أبيض من AgCl"),
    ("precipitation", "BaCl2 + Na2SO4", "BaSO4 + 2NaCl", "محلول مائي، راسب أبيض من BaSO4"),
    ("precipitation", "BaCl2 + CuSO4", "BaSO4 + CuCl2", "محلول مائي، راسب أبيض من BaSO4"),
    ("precipitation", "CaCl2 + Na2CO3", "CaCO3 + 2NaCl", "محلول مائي، راسب أبيض من CaCO3"),
    ("precipitation", "FeCl3 + 3NaOH", "Fe(OH)3 + 3NaCl", "محلول مائي، راسب بني محمر"),
    ("precipitation", "CuSO4 + 2NaOH", "Cu(OH)2 + Na2SO4", "محلول مائي، راسب أزرق"),
    ("single_displacement", "Zn + CuSO4", "ZnSO4 + Cu", "محلول مائي، 25 °C"),
    ("single_displacement", "Fe + CuSO4", "FeSO4 + Cu", "محلول مائي، 25 °C"),
    ("single_displacement", "Zn + 2HCl", "ZnCl2 + H2", "محلول حمضي مخفف، تصاعد H2"),
    ("single_displacement", "Mg + 2HCl", "MgCl2 + H2", "محلول حمضي مخفف، تصاعد H2"),
    ("single_displacement", "Fe + 2HCl", "FeCl2 + H2", "محلول حمضي مخفف، تصاعد H2"),
    ("single_displacement", "2Al + 6HCl", "2AlCl3 + 3H2", "محلول حمضي، تصاعد H2"),
    ("single_displacement", "Zn + H2SO4", "ZnSO4 + H2", "حمض كبريتيك مخفف، تصاعد H2"),
    ("single_displacement", "2Na + 2H2O", "2NaOH + H2", "تفاعل عنيف طارد للحرارة"),
    ("single_displacement", "Cu + 2AgNO3", "Cu(NO3)2 + 2Ag", "محلول مائي، ترسب الفضة"),
    ("redox_synthesis", "2Mg + O2", "2MgO", "اشتعال بلهب أبيض ساطع"),
    ("redox_synthesis", "2Na + Cl2", "2NaCl", "تفاعل طارد للحرارة"),
    ("redox_synthesis", "2Ca + O2", "2CaO", "تسخين"),
    ("redox_synthesis", "S + O2", "SO2", "اشتعال الكبريت"),
    ("redox_synthesis", "2SO2 + O2", "2SO3", "حفاز V2O5، تسخين"),
    ("redox_synthesis", "2CO + O2", "2CO2", "اشتعال"),
    ("redox_synthesis", "4Al + 3O2", "2Al2O3", "تسخين"),
    ("redox_synthesis", "2Cu + O2", "2CuO", "تسخين في الهواء"),
    ("redox_synthesis", "Fe2O3 + 3CO", "2Fe + 3CO2", "الفرن العالي، درجات حرارة مرتفعة"),
    ("redox_synthesis", "H2 + Cl2", "2HCl", "ضوء أو شرارة"),
    ("redox_synthesis", "2NO + O2", "2NO2", "درجة حرارة الغرفة"),
    ("redox_synthesis", "2Zn + O2", "2ZnO", "تسخين"),
    ("redox_synthesis", "N2 + O2", "2NO", "درجات حرارة مرتفعة جدا أو برق"),
    ("decomposition", "CaCO3", "CaO + CO2", "تسخين عند نحو 850 °C"),
    ("decomposition", "2KClO3", "2KCl + 3O2", "تسخين مع حفاز MnO2"),
    ("decomposition", "2KMnO4", "K2MnO4 + MnO2 + O2", "تسخين"),
    ("decomposition", "2H2O", "2H2 + O2", "تحليل كهربائي"),
    ("decomposition", "2NaHCO3", "Na2CO3 + H2O + CO2", "تسخين"),
    ("decomposition", "H2CO3", "H2O + CO2", "درجة حرارة الغرفة"),
    ("decomposition", "2NH3", "N2 + 3H2", "درجات حرارة مرتفعة مع حفاز"),
    ("decomposition", "Cu(OH)2", "CuO + H2O", "تسخين"),
    ("organic_reaction", "CH3COOH + C2H5OH", "CH3COOCH2CH3 + H2O", "حفاز H2SO4 المركز، تسخين (أسترة فيشر)"),
    ("organic_reaction", "CH3COOH + CH3OH", "CH3COOCH3 + H2O", "حفاز H2SO4 المركز، تسخين (أسترة فيشر)"),
    ("organic_reaction", "C2H4 + H2", "C2H6", "حفاز Pt أو Pd أو Ni (هدرجة)"),
    ("organic_reaction", "C2H4 + H2O", "C2H5OH", "حفاز H3PO4، ضغط وحرارة عاليان (إماهة)"),
    ("organic_reaction", "C2H2 + H2", "C2H4", "حفاز لندلار (هدرجة جزئية)"),
    ("organic_reaction", "CH4 + Cl2", "CH3Cl + HCl", "ضوء فوق بنفسجي (استبدال جذري حر)"),
    ("organic_reaction", "C6H6 + HNO3", "C6H5NO2 + H2O", "حفاز H2SO4 المركز (نترتة)"),
    ("organic_reaction", "C2H4 + Br2", "C2H4Br2", "إضافة إلى الرابطة المزدوجة، 25 °C"),
    ("organic_reaction", "C6H5OH + NaOH", "C6H5ONa + H2O", "محلول مائي (حمضية الفينول)"),
    ("organic_reaction", "C6H5NH2 + HCl", "C6H5NH3Cl", "محلول مائي (تكوين ملح الأمين)"),
    ("biochemical_reaction", "C6H12O6", "2C2H5OH + 2CO2", "إنزيمات الخميرة، وسط لاهوائي (تخمر كحولي)"),
    ("biochemical_reaction", "6CO2 + 6H2O", "C6H12O6 + 6O2", "طاقة ضوئية وصبغة الكلوروفيل (البناء الضوئي)"),
    ("biochemical_reaction", "C10H16N5O13P3 + H2O", "C10H15N5O10P2 + H3PO4", "إنزيمات ATPase، 37 °C (تحلل ATP)"),
    ("biochemical_reaction", "CO(NH2)2 + H2O", "CO2 + 2NH3", "إنزيم اليورياز (تحلل اليوريا)"),
    ("biochemical_reaction", "C12H22O11 + H2O", "2C6H12O6", "إنزيم السكراز أو حمض مخفف (تحلل السكروز)"),
]


def _species_index() -> dict[str, tuple[str, str]]:
    from search.local_catalog import LOCAL_SCIENTIFIC_CATALOG

    index: dict[str, tuple[str, str]] = {}
    for entry in LOCAL_SCIENTIFIC_CATALOG:
        sym = entry["symbol_or_formula"]
        index.setdefault(sym, (entry["name_ar"], entry["name_en"]))
    index.update({k: v for k, v in EXTRA_SPECIES_NAMES.items() if k not in index})
    index["C6H12O6"] = ("الجلوكوز", "D-Glucose")
    return index


def _infer_species_state(formula: str, category: str, cond: str) -> str:
    phase = get_standard_phase_code(formula)
    if phase in {"g", "l"}:
        if formula in {"HCl", "HNO3", "H2SO4", "CH3COOH"} and ("محلول" in cond or category in {"neutralization", "precipitation", "single_displacement"}):
            return "aq"
        return phase
    if "محلول" in cond and formula not in {"AgCl", "BaSO4", "CaCO3", "Fe(OH)3", "Cu(OH)2", "Zn", "Fe", "Mg", "Al", "Cu", "Ag", "Na"}:
        return "aq"
    return phase


def _parse_side(side: str, index: dict[str, tuple[str, str]], category: str, cond: str) -> list[dict[str, Any]]:
    items = []
    for token in side.split("+"):
        token = token.strip()
        m = re.fullmatch(r"(\d+)?(.+)", token)
        coeff = int(m.group(1)) if m.group(1) else 1
        formula = m.group(2).strip()
        name_ar, name_en = index.get(formula, (formula, formula))
        state = _infer_species_state(formula, category, cond)
        items.append({"formula": formula, "name_ar": name_ar, "name_en": name_en, "coefficient": coeff, "state": state})
    return items


def _extract_conditions_dict(category: str, cond_str: str) -> dict[str, str]:
    cond_dict: dict[str, str] = {
        "temperature": cond_str,
        "pressure": "1 atm ما لم يذكر خلاف ذلك",
    }
    if "محلول مائي" in cond_str or "محلول حمضي" in cond_str:
        cond_dict["solvent"] = "محلول مائي (Aqueous Solution)"
    elif "طور غازي" in cond_str or category == "combustion":
        cond_dict["solvent"] = "طور غازي أو تفاعل مباشر بدون مذيب"
    if "حفاز" in cond_str or "إنزيم" in cond_str:
        cond_dict["catalyst"] = cond_str
    return cond_dict


def _build_generated_reactions() -> list[dict[str, Any]]:
    index = _species_index()
    out = []
    for idx, (category, lhs, rhs, cond) in enumerate(_GENERATED_REACTIONS, start=1):
        reactants = _parse_side(lhs, index, category, cond)
        products = _parse_side(rhs, index, category, cond)
        domain, domain_label = DOMAIN_BY_CATEGORY[category]
        label = CATEGORY_LABELS_AR[category]
        r_ar = " + ".join(r["name_ar"] for r in reactants)
        p_ar = " + ".join(p["name_ar"] for p in products)
        r_en = " + ".join(r["name_en"] for r in reactants)
        p_en = " + ".join(p["name_en"] for p in products)
        ratio = ":".join(str(r["coefficient"]) for r in reactants)
        source_ref = CATEGORY_DEFAULT_SOURCES.get(category, "Curated Reaction Registry (IUPAC / NIST Reference)")
        out.append(
            {
                "reaction_id": f"rxn_{idx:03d}_{category}",
                "domain": domain,
                "domain_label_ar": domain_label,
                "category": category,
                "reaction_name_ar": f"{label.split(' (')[0]}: {r_ar} إلى {p_ar}",
                "reaction_name_en": f"{r_en} -> {p_en}",
                "reaction_type_ar": label,
                "reactants": reactants,
                "products": products,
                "conditions": _extract_conditions_dict(category, cond),
                "explanation_ar": f"تتفاعل المواد ({r_ar}) بنسبة مولية ({ratio}) لتكوين ({p_ar}) وفق المعادلة الموزونة المسجلة في السجل المعتمد.",
                "assumptions_and_limits_ar": "حساب نظري بفرض اكتمال التفاعل بنسبة 100% دون تفاعلات جانبية، وهو ليس إثباتا لتجربة عملية.",
                "source": source_ref,
            }
        )
    return out


for _rxn in VERIFIED_REACTIONS_REGISTRY:
    _rxn.setdefault("category", REACTION_CATEGORY_BY_ID.get(_rxn["reaction_id"], "redox_synthesis"))
    _rxn.setdefault("source", CATEGORY_DEFAULT_SOURCES.get(_rxn["category"], "Curated Reaction Registry (IUPAC / NIST)"))

VERIFIED_REACTIONS_REGISTRY.extend(_build_generated_reactions())


def _clean(f: str) -> str:
    return replace_subscripts(f or "").replace(" ", "")


def _known_exact_formulas() -> set[str]:
    from search.local_catalog import LOCAL_SCIENTIFIC_CATALOG

    return {_clean(e["symbol_or_formula"]).lower() for e in LOCAL_SCIENTIFIC_CATALOG}


def _formula_matches(reaction_formula: str, user_formula: str, known: set[str]) -> bool:
    rf, uf = _clean(reaction_formula), _clean(user_formula)
    if rf == uf:
        return True
    if rf.lower() == uf.lower() and len(uf) > 2:
        return True
    if uf.lower() in known:
        return False
    return formulas_are_equivalent(rf, uf)


def _domain_matches(rxn: dict[str, Any], domain: str | None) -> bool:
    if not domain or domain in {"all", "all_domains"}:
        return True
    rxn_dom = rxn.get("domain", "")
    rxn_cat = rxn.get("category", "")
    if rxn_dom == domain:
        return True
    # Allow general_chemistry reactions that belong to sub-categories or vice-versa when domain is general_chemistry
    if domain == "general_chemistry" and rxn_dom in {"general_chemistry", "acid_base_chemistry", "redox_and_precipitation"}:
        return True
    if domain == "acid_base_chemistry" and rxn_cat == "neutralization":
        return True
    if domain == "redox_and_precipitation" and rxn_cat in {"redox_synthesis", "precipitation", "single_displacement"}:
        return True
    return False


def find_verified_reactions(
    reactant_formulas: list[str],
    reaction_type: str | None = None,
    domain: str | None = None,
) -> list[dict[str, Any]]:
    user = [f for f in reactant_formulas if f]
    if not user:
        return []
    known = _known_exact_formulas()
    matches = []
    for rxn in VERIFIED_REACTIONS_REGISTRY:
        if reaction_type and reaction_type != "all" and rxn.get("category") != reaction_type:
            continue
        if not _domain_matches(rxn, domain):
            continue
        rf = [r["formula"] for r in rxn["reactants"]]
        if len(rf) != len(user):
            continue
        remaining = list(user)
        ok = True
        for r in rf:
            hit = next((u for u in remaining if _formula_matches(r, u, known)), None)
            if hit is None:
                ok = False
                break
            remaining.remove(hit)
        if ok and not remaining:
            matches.append(dict(rxn))
    return matches


def find_verified_reaction(
    reactant_formulas: list[str],
    domain: str | None = None,
    reaction_type: str | None = None,
) -> dict[str, Any] | None:
    found = find_verified_reactions(reactant_formulas, reaction_type=reaction_type, domain=domain)
    return found[0] if found else None


def search_reactions_by_substance(
    formula_or_key: str,
    domain: str | None = None,
    reaction_type: str | None = None,
) -> list[dict[str, Any]]:
    target = formula_or_key or ""
    known = _known_exact_formulas()
    results = []
    for rxn in VERIFIED_REACTIONS_REGISTRY:
        if reaction_type and reaction_type != "all" and rxn.get("category") != reaction_type:
            continue
        if not _domain_matches(rxn, domain):
            continue
        species = [s["formula"] for s in rxn["reactants"] + rxn["products"]]
        if any(_formula_matches(s, target, known) for s in species):
            results.append(dict(rxn))
    return results

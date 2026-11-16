import re
from typing import Any
from chemistry_engine.formula_analyzer import (
    formulas_are_equivalent,
    get_standard_phase_code,
    replace_subscripts,
)

DOCUMENTATION_TYPE_LABELS_AR: dict[str, str] = {
    "experimentally_documented": "تفاعل تجريبي موثق مخبرياً (Experimentally Documented Reaction)",
    "literature_referenced": "تفاعل مستند إلى مرجع علمي قياسي (Literature-Referenced Reaction)",
    "theoretical_model": "نموذج أو معادلة ستوكيومترية نظرية شاملة (Theoretical / Stoichiometric Model)",
}

VERIFIED_REACTIONS_REGISTRY: list[dict[str, Any]] = [
    {
        "reaction_id": "rxn_water_synthesis",
        "domain": "general_chemistry",
        "domain_label_ar": "الكيمياء العامة (General Chemistry)",
        "category": "redox_synthesis",
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
            "solvent": "طور غازي مباشر بدون مذيب",
            "catalyst": "بلاتين (Pt) اختياري لخفض طاقة التنشيط",
            "thermodynamics_ar": "تفاعل طارد للحرارة بشدة (ΔH° = -571.6 kJ لكل 2 mol من الماء السائل)",
        },
        "reference_conditions_spec": {
            "temp_c_min": 400.0,
            "temp_c_max": 2500.0,
            "temp_c_optimal": 550.0,
            "pressure_atm_min": 0.5,
            "pressure_atm_max": 50.0,
            "pressure_atm_optimal": 1.0,
            "allowed_solvents": ["neat_gas_solid"],
            "required_catalysts": ["none", "pt_pd"],
            "catalyst_is_mandatory": False,
            "catalyst_lowers_temp_c_min": 20.0,
        },
        "explanation_ar": (
            "يتفاعل جزيئان من غاز الهيدروجين (H2) مع جزيء واحد من غاز الأكسجين (O2) بنسبة مولية ثابتة (2:1) "
            "لتكوين جزيئين من الماء (H2O). تنكسر الروابط (H–H) و (O=O) وتتكون روابط تساهمية قطبية جديدة (O–H) أكثر استقرارا."
        ),
        "assumptions_and_limits_ar": (
            "يفترض هذا الحساب النظري اكتمال التفاعل بنسبة 100% في نظام مغلق دون فقدان حراري أو تفاعلات جانبية. "
            "هذه نتيجة حسابية نظرية وليست إثباتا لتجربة عملية مخبرية."
        ),
        "documentation_type": "literature_referenced",
        "documentation_label_ar": DOCUMENTATION_TYPE_LABELS_AR["literature_referenced"],
        "source_id": "NIST-JANAF-H2O-L-7732-18-5",
        "source": "Chase, M.W., NIST-JANAF Thermochemical Tables, 4th Ed., J. Phys. Chem. Ref. Data, Monograph 9 (1998), p. 1324 (H2O Formation)",
        "reference_url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185",
    },
    {
        "reaction_id": "rxn_iron_oxidation",
        "domain": "general_chemistry",
        "domain_label_ar": "الكيمياء العامة (General Chemistry)",
        "category": "redox_synthesis",
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
            "temperature": "درجة حرارة الغرفة بوجود الرطوبة (بطيء) أو التسخين في الهواء (20 – 800 °C)",
            "pressure": "1 atm",
            "solvent": "طور صلب-غاز مباشر",
            "catalyst": "بخار الماء يسرع عملية التأكسد الكهروكيميائية",
            "thermodynamics_ar": "تفاعل طارد للحرارة (ΔH° = -1648.4 kJ لتكوين 2 mol من Fe2O3)",
        },
        "reference_conditions_spec": {
            "temp_c_min": 15.0,
            "temp_c_max": 1000.0,
            "temp_c_optimal": 25.0,
            "pressure_atm_min": 0.5,
            "pressure_atm_max": 10.0,
            "pressure_atm_optimal": 1.0,
            "allowed_solvents": ["neat_gas_solid", "aqueous"],
            "required_catalysts": ["none"],
            "catalyst_is_mandatory": False,
        },
        "explanation_ar": (
            "تتأكسد 4 مولات من ذرات الحديد (Fe) بفقد إلكترونات تكافئها إلى 3 مولات من جزيئات الأكسجين (O2) "
            "لتكوين 2 مول من أكسيد الحديد الثلاثي الصلب (Fe2O3) وفق النسبة المولية 4:3:2."
        ),
        "assumptions_and_limits_ar": "حساب ستوكيوميتري نظري بفرض تحول الحديد بالكامل إلى حالة الأكسدة +3.",
        "documentation_type": "literature_referenced",
        "documentation_label_ar": DOCUMENTATION_TYPE_LABELS_AR["literature_referenced"],
        "source_id": "Greenwood-Earnshaw-2ed-Ch25-Fe2O3",
        "source": "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed., Butterworth-Heinemann (1997), Ch. 25, pp. 1079–1082 (Iron Oxides)",
        "reference_url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C1309371",
    },
    {
        "reaction_id": "rxn_acid_base_hcl_naoh",
        "domain": "general_chemistry",
        "domain_label_ar": "الكيمياء العامة (General Chemistry)",
        "category": "neutralization",
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
        "reference_conditions_spec": {
            "temp_c_min": 5.0,
            "temp_c_max": 60.0,
            "temp_c_optimal": 25.0,
            "pressure_atm_min": 0.5,
            "pressure_atm_max": 5.0,
            "pressure_atm_optimal": 1.0,
            "allowed_solvents": ["aqueous"],
            "required_catalysts": ["none"],
            "catalyst_is_mandatory": False,
        },
        "explanation_ar": (
            "يتحد أيون الهيدروجين (H⁺) من الحمض مع أيون الهيدروكسيد (OH⁻) من القاعدة بنسبة مولية متكافئة (1:1) "
            "لتكوين جزيء الماء (H2O) وملح كلوريد الصوديوم الذائب (NaCl)."
        ),
        "assumptions_and_limits_ar": "حساب نظري لتفاعل تام في الوسط المائي المخفف.",
        "documentation_type": "experimentally_documented",
        "documentation_label_ar": DOCUMENTATION_TYPE_LABELS_AR["experimentally_documented"],
        "source_id": "USP-NF-Volumetric-HCl-1N-NaOH",
        "source": "United States Pharmacopeia (USP-NF), General Chapter <541> Titrimetry — Hydrochloric Acid 1 N Standardization with Sodium Hydroxide",
        "reference_url": "https://pubchem.ncbi.nlm.nih.gov/compound/313",
    },
    {
        "reaction_id": "rxn_haber_ammonia",
        "domain": "general_chemistry",
        "domain_label_ar": "الكيمياء العامة (General Chemistry)",
        "category": "redox_synthesis",
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
            "solvent": "طور غازي مباشر بدون مذيب",
            "catalyst": "حفاز الحديد المنشط (Fe with K2O/Al2O3 promoters)",
            "thermodynamics_ar": "تفاعل طارد للحرارة عكسي (ΔH° = -92.4 kJ لكل 2 mol من NH3)",
        },
        "reference_conditions_spec": {
            "temp_c_min": 380.0,
            "temp_c_max": 550.0,
            "temp_c_optimal": 450.0,
            "pressure_atm_min": 150.0,
            "pressure_atm_max": 300.0,
            "pressure_atm_optimal": 200.0,
            "allowed_solvents": ["neat_gas_solid"],
            "required_catalysts": ["fe_haber"],
            "catalyst_is_mandatory": True,
        },
        "explanation_ar": (
            "يتفاعل مول واحد من غاز النيتروجين (N2) مع 3 مولات من غاز الهيدروجين (H2) لتكوين مولين من غاز الأمونيا (NH3) "
            "بنسبة مولية 1:3:2."
        ),
        "assumptions_and_limits_ar": "يمثل الناتج المحسوب الحد الأقصى النظري (100% تحويل) قبل الوصول إلى حالة الاتزان العكسي.",
        "documentation_type": "experimentally_documented",
        "documentation_label_ar": DOCUMENTATION_TYPE_LABELS_AR["experimentally_documented"],
        "source_id": "Ullmann-2012-Ammonia-DOI-10.1002-14356007.a02_143.pub3",
        "source": "Appl, M., 'Ammonia, 2. Production Processes', Ullmann's Encyclopedia of Industrial Chemistry, Wiley-VCH (2012), DOI: 10.1002/14356007.a02_143.pub3",
        "reference_url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C7664417",
    },
    {
        "reaction_id": "rxn_h2o2_decomposition",
        "domain": "general_chemistry",
        "domain_label_ar": "الكيمياء العامة (General Chemistry)",
        "category": "decomposition",
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
            "pressure": "1 atm",
            "solvent": "محلول مائي (Aqueous Solution)",
            "catalyst": "ثاني أكسيد المنغنيز (MnO2) أو إنزيم الكاتالاز الحيوي",
            "thermodynamics_ar": "تفاعل طارد للحرارة (ΔH° = -196.1 kJ لكل 2 mol من H2O2)",
        },
        "reference_conditions_spec": {
            "temp_c_min": 10.0,
            "temp_c_max": 60.0,
            "temp_c_optimal": 25.0,
            "pressure_atm_min": 0.5,
            "pressure_atm_max": 5.0,
            "pressure_atm_optimal": 1.0,
            "allowed_solvents": ["aqueous", "physiological_buffer"],
            "required_catalysts": ["mno2_catalase", "enzymatic"],
            "catalyst_is_mandatory": False,
        },
        "explanation_ar": "يتفكك جزيئان من فوق أكسيد الهيدروجين (H2O2) لإنتاج جزيئين من الماء (H2O) وجزيء واحد من غاز الأكسجين (O2).",
        "assumptions_and_limits_ar": "حساب نظري بفرض التفكك التام للمتفاعل.",
        "documentation_type": "experimentally_documented",
        "documentation_label_ar": DOCUMENTATION_TYPE_LABELS_AR["experimentally_documented"],
        "source_id": "RHEA:15889",
        "source": "Rhea Biochemical Database RHEA:15889 (EC 1.11.1.6 Catalase: 2 H2O2 = 2 H2O + O2)",
        "reference_url": "https://www.rhea-db.org/rhea/15889",
    },
    {
        "reaction_id": "rxn_co2_water_hydration",
        "domain": "general_chemistry",
        "domain_label_ar": "الكيمياء العامة (General Chemistry)",
        "category": "redox_synthesis",
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
            "solvent": "محلول مائي (Aqueous Solution)",
            "catalyst": "في الأنظمة الحيوية يحفزه إنزيم الأنهيدراز الكربوني (Carbonic Anhydrase)",
        },
        "reference_conditions_spec": {
            "temp_c_min": 5.0,
            "temp_c_max": 50.0,
            "temp_c_optimal": 25.0,
            "pressure_atm_min": 0.8,
            "pressure_atm_max": 10.0,
            "pressure_atm_optimal": 1.0,
            "allowed_solvents": ["aqueous", "physiological_buffer"],
            "required_catalysts": ["none", "enzymatic"],
            "catalyst_is_mandatory": False,
        },
        "explanation_ar": "يتفاعل مول واحد من ثاني أكسيد الكربون (CO2) مع مول واحد من الماء (H2O) لتكوين مول واحد من حمض الكربونيك (H2CO3).",
        "assumptions_and_limits_ar": "يمثل الحساب القيمة النظرية الكاملة لتحول CO2 الذائب إلى H2CO3.",
        "documentation_type": "experimentally_documented",
        "documentation_label_ar": DOCUMENTATION_TYPE_LABELS_AR["experimentally_documented"],
        "source_id": "RHEA:16833",
        "source": "Rhea Biochemical Database RHEA:16833 (EC 4.2.1.1 Carbonic Anhydrase: CO2 + H2O = H2CO3 / HCO3- + H+)",
        "reference_url": "https://www.rhea-db.org/rhea/16833",
    },
    {
        "reaction_id": "rxn_methane_combustion",
        "domain": "organic_chemistry",
        "domain_label_ar": "الكيمياء العضوية (Organic Chemistry)",
        "category": "combustion",
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
            "temperature": "شرارة اشتعال في وفرة من الأكسجين (> 540 °C)",
            "pressure": "1 atm",
            "solvent": "طور غازي مباشر بدون مذيب",
            "thermodynamics_ar": "تفاعل طارد للحرارة بشدة (ΔH°c = -890.3 kJ/mol)",
        },
        "reference_conditions_spec": {
            "temp_c_min": 400.0,
            "temp_c_max": 2200.0,
            "temp_c_optimal": 850.0,
            "pressure_atm_min": 0.5,
            "pressure_atm_max": 20.0,
            "pressure_atm_optimal": 1.0,
            "allowed_solvents": ["neat_gas_solid"],
            "required_catalysts": ["none"],
            "catalyst_is_mandatory": False,
        },
        "explanation_ar": "يحترق مول واحد من الميثان (CH4) بوجود مولين من غاز الأكسجين (O2) لينتج مولا واحدا من ثاني أكسيد الكربون (CO2) ومولين من الماء (H2O).",
        "assumptions_and_limits_ar": "يفترض الحساب احتراقا كاملا دون تكوين أول أكسيد الكربون (CO) أو السخام.",
        "documentation_type": "literature_referenced",
        "documentation_label_ar": DOCUMENTATION_TYPE_LABELS_AR["literature_referenced"],
        "source_id": "NIST-WebBook-CAS-74-82-8-Combustion",
        "source": "Pittam, D.A. & Pilcher, G., J. Chem. Soc., Faraday Trans. 1 (1972), 68, 2224–2229 (Standard Enthalpy of Combustion of Methane, NIST WebBook CAS 74-82-8)",
        "reference_url": "https://webbook.nist.gov/cgi/cbook.cgi?ID=C74828",
    },
    {
        "reaction_id": "rxn_glucose_respiration",
        "domain": "biochemistry",
        "domain_label_ar": "الكيمياء الحيوية (Biochemistry)",
        "category": "biochemical_reaction",
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
            "pressure": "1 atm",
            "solvent": "وسط منظم فسيولوجي مائي (Physiological Aqueous Medium)",
            "catalyst": "سلسلة إنزيمات التحلل السكري ودورة كريبس والفسفرة التأكسدية",
            "thermodynamics_ar": "تفاعل منتج للطاقة الحرة القياسية (ΔG°' = -2870 kJ/mol)",
        },
        "reference_conditions_spec": {
            "temp_c_min": 15.0,
            "temp_c_max": 45.0,
            "temp_c_optimal": 37.0,
            "pressure_atm_min": 0.8,
            "pressure_atm_max": 2.0,
            "pressure_atm_optimal": 1.0,
            "allowed_solvents": ["aqueous", "physiological_buffer"],
            "required_catalysts": ["enzymatic"],
            "catalyst_is_mandatory": True,
        },
        "explanation_ar": "يتأكسد مول واحد من الجلوكوز (C6H12O6) بوجود 6 مولات من الأكسجين (O2) لينتج 6 مولات من ثاني أكسيد الكربون (CO2) و6 مولات من الماء (H2O).",
        "assumptions_and_limits_ar": "تمثل المعادلة المحصلة الستوكيومترية النظرية الكلية لمسارات التنفس الخلوي الهوائي المتعددة.",
        "documentation_type": "theoretical_model",
        "documentation_label_ar": DOCUMENTATION_TYPE_LABELS_AR["theoretical_model"],
        "source_id": "Lehninger-8ed-Ch14-19-Net-Respiration",
        "source": "Nelson, D.L. & Cox, M.M., Lehninger Principles of Biochemistry, 8th Ed., W.H. Freeman (2021), Ch. 14–19, Table 19-5 (Net Stoichiometric Balance of Aerobic Glucose Oxidation)",
        "reference_url": "https://pubchem.ncbi.nlm.nih.gov/compound/5793",
    },
]

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
    "Na2O": ("أكسيد الصوديوم", "Sodium Oxide"),
    "K2O": ("أكسيد البوتاسيوم", "Potassium Oxide"),
    "BaO": ("أكسيد الباريوم", "Barium Oxide"),
    "H2SO3": ("حمض الكبريتوز", "Sulfurous Acid"),
    "Ag2O": ("أكسيد الفضة", "Silver(I) Oxide"),
    "HgO": ("أكسيد الزئبق الثنائي", "Mercury(II) Oxide"),
    "Hg": ("الزئبق", "Mercury"),
    "NH4HCO3": ("بيكربونات الأمونيوم", "Ammonium Bicarbonate"),
    "MgCO3": ("كربونات المغنيسيوم", "Magnesium Carbonate"),
    "MgSO4": ("كبريتات المغنيسيوم", "Magnesium Sulfate"),
    "KBr": ("بروميد البوتاسيوم", "Potassium Bromide"),
    "KI": ("يوديد البوتاسيوم", "Potassium Iodide"),
    "Br2": ("البروم", "Bromine"),
    "I2": ("اليود", "Iodine"),
    "Pb(NO3)2": ("نترات الرصاص الثنائي", "Lead(II) Nitrate"),
    "PbSO4": ("كبريتات الرصاص الثنائي", "Lead(II) Sulfate"),
    "Fe(OH)2": ("هيدروكسيد الحديدوز", "Iron(II) Hydroxide"),
    "Ba": ("الباريوم", "Barium"),
    "PbI2": ("يوديد الرصاص الثنائي", "Lead(II) Iodide"),
    "AgBr": ("بروميد الفضة", "Silver Bromide"),
    "AgI": ("يوديد الفضة", "Silver Iodide"),
    "FeS": ("كبريتيد الحديد الثنائي", "Iron(II) Sulfide"),
    "Na3PO4": ("فوسفات الصوديوم", "Trisodium Phosphate"),
    "H2C2O4": ("حمض الأكساليك", "Oxalic Acid"),
    "Na2C2O4": ("أكسالات الصوديوم", "Sodium Oxalate"),
    "K2CO3": ("كربونات البوتاسيوم", "Potassium Carbonate"),
    "CH3COOK": ("أسيتات البوتاسيوم", "Potassium Acetate"),
    "Ca(NO3)2": ("نترات الكالسيوم", "Calcium Nitrate"),
    "MnCl2": ("كلوريد المنغنيز الثنائي", "Manganese(II) Chloride"),
    "Cu2(OH)2CO3": ("كربونات النحاس القاعدية (الملاكيت)", "Basic Copper(II) Carbonate (Malachite)"),
    "C2H5Cl": ("كلورو الإيثان", "Chloroethane"),
    "C6H5Br": ("برومو البنزين", "Bromobenzene"),
}

# Tuple format:
# (category, reactants, products, conditions_ar, doc_type, source_id, specific_source, reference_url)
_GENERATED_REACTIONS: list[tuple[str, str, str, str, str, str, str, str]] = [
    (
        "neutralization", "CH3COOH + NaOH", "CH3COONa + H2O", "محلول مائي، 25 °C",
        "experimentally_documented", "USP43-NF38-Acetic-Acid-Assay",
        "United States Pharmacopeia (USP 43-NF 38), Monograph: Glacial Acetic Acid — Alkalimetric Assay with 1 N NaOH",
        "https://pubchem.ncbi.nlm.nih.gov/compound/176",
    ),
    (
        "neutralization", "CH3COOH + NaHCO3", "CH3COONa + H2O + CO2", "محلول مائي، تصاعد غاز CO2",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-IV.2-Carbonate-Acetate",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. IV.2, p. 331 (Bicarbonate Effervescence with Acetic Acid)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/516892",
    ),
    (
        "neutralization", "H2SO4 + 2NaOH", "Na2SO4 + 2H2O", "محلول مائي، 25 °C",
        "experimentally_documented", "USP43-NF38-Sulfuric-Acid-Assay",
        "United States Pharmacopeia (USP 43-NF 38), Monograph: Sulfuric Acid — Volumetric Titration with 1 N Sodium Hydroxide",
        "https://pubchem.ncbi.nlm.nih.gov/compound/1118",
    ),
    (
        "neutralization", "HNO3 + KOH", "KNO3 + H2O", "محلول مائي، 25 °C",
        "experimentally_documented", "Vogel-Quant-5ed-Ch10-HNO3-KOH",
        "Vogel's Textbook of Quantitative Chemical Analysis, 5th Ed., Longman (1989), Ch. 10, p. 294 (Strong Acid–Strong Base Neutralization)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/944",
    ),
    (
        "neutralization", "2HCl + Ca(OH)2", "CaCl2 + 2H2O", "محلول مائي، 25 °C",
        "experimentally_documented", "USP43-NF38-Calcium-Hydroxide-Solution",
        "United States Pharmacopeia (USP 43-NF 38), Monograph: Calcium Hydroxide Topical Solution — Acidimetric Titration with 0.1 N HCl",
        "https://pubchem.ncbi.nlm.nih.gov/compound/6093208",
    ),
    (
        "neutralization", "2HCl + Mg(OH)2", "MgCl2 + 2H2O", "محلول مائي، 25 °C",
        "experimentally_documented", "USP43-NF38-Milk-of-Magnesia-Assay",
        "United States Pharmacopeia (USP 43-NF 38), Monograph: Magnesium Hydroxide — Antacid Neutralization Capacity Assay with HCl",
        "https://pubchem.ncbi.nlm.nih.gov/compound/73981",
    ),
    (
        "neutralization", "H3PO4 + 3NaOH", "Na3PO4 + 3H2O", "محلول مائي، 25 °C",
        "experimentally_documented", "Vogel-Quant-5ed-Sec-10.12-Polyprotic-H3PO4",
        "Vogel's Textbook of Quantitative Chemical Analysis, 5th Ed. (1989), Sec. 10.12, p. 302 (Triprotic Neutralization of Orthophosphoric Acid)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/1004",
    ),
    (
        "neutralization", "HCl + NH3", "NH4Cl", "طور غازي أو مائي، 25 °C",
        "experimentally_documented", "Greenwood-Earnshaw-2ed-Ch11-NH4Cl",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 11, p. 423 (Direct Neutralization of Ammonia with Hydrogen Chloride)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/25517",
    ),
    (
        "neutralization", "HCl + NaHCO3", "NaCl + H2O + CO2", "محلول مائي، تصاعد غاز CO2",
        "experimentally_documented", "PhEur-10ed-Monograph-0195-NaHCO3",
        "European Pharmacopoeia (Ph. Eur. 10th Ed.), Monograph 0195: Sodium Hydrogen Carbonate — Acidimetric Titration with 1 M HCl",
        "https://pubchem.ncbi.nlm.nih.gov/compound/516892",
    ),
    (
        "neutralization", "2HCl + Na2CO3", "2NaCl + H2O + CO2", "محلول مائي، تصاعد غاز CO2",
        "experimentally_documented", "Vogel-Quant-5ed-Sec-10.9-Na2CO3-Primary-Standard",
        "Vogel's Textbook of Quantitative Chemical Analysis, 5th Ed. (1989), Sec. 10.9, p. 297 (Standardization of Hydrochloric Acid with Anhydrous Sodium Carbonate)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/10340",
    ),
    (
        "neutralization", "2HCl + CaCO3", "CaCl2 + H2O + CO2", "محلول مائي، تصاعد غاز CO2",
        "experimentally_documented", "USP43-NF38-Calcium-Carbonate-Acid-Neutralizing",
        "United States Pharmacopeia (USP 43-NF 38), General Chapter <301> Acid-Neutralizing Capacity — Calcium Carbonate Dissolution in HCl",
        "https://pubchem.ncbi.nlm.nih.gov/compound/10112",
    ),
    (
        "neutralization", "H2SO4 + Ba(OH)2", "BaSO4 + 2H2O", "محلول مائي، يتكون راسب BaSO4",
        "experimentally_documented", "Vogel-Quant-5ed-Ch11-Conductometric-BaOH2-H2SO4",
        "Vogel's Textbook of Quantitative Chemical Analysis, 5th Ed. (1989), Ch. 13, p. 521 (Conductometric Titration & Precipitation of Ba(OH)2 with H2SO4)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/24414",
    ),
    (
        "neutralization", "HCl + KOH", "KCl + H2O", "محلول مائي، 25 °C",
        "experimentally_documented", "USP43-NF38-Potassium-Hydroxide-Assay",
        "United States Pharmacopeia (USP 43-NF 38), Monograph: Potassium Hydroxide — Volumetric Titration with 1 N Hydrochloric Acid",
        "https://pubchem.ncbi.nlm.nih.gov/compound/14797",
    ),
    (
        "redox_synthesis", "CaO + H2O", "Ca(OH)2", "تفاعل طارد للحرارة، 25 °C",
        "experimentally_documented", "Greenwood-Earnshaw-2ed-Ch5-Lime-Slaking",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 5, p. 122 (Exothermic Hydration / Slaking of Quicklime CaO)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C1305788",
    ),
    (
        "redox_synthesis", "CO2 + Ca(OH)2", "CaCO3 + H2O", "محلول مائي، يعكر ماء الجير",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-IV.2-Limewater-CO2",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. IV.2, p. 330 (Limewater Test for Carbon Dioxide)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/10112",
    ),
    (
        "redox_synthesis", "SO3 + H2O", "H2SO4", "إضافة حذرة للماء",
        "literature_referenced", "Ullmann-2012-Sulfuric-Acid-Contact-Process",
        "Müller, H., 'Sulfuric Acid and Sulfur Trioxide', Ullmann's Encyclopedia of Industrial Chemistry, Wiley-VCH (2012), Contact Process Hydration",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C7446119",
    ),
    (
        "combustion", "2C2H6 + 7O2", "4CO2 + 6H2O", "شرارة اشتعال في وفرة من الأكسجين",
        "theoretical_model", "NIST-WebBook-CAS-74-84-0-Ethane-Combustion",
        "Pittam, D.A. & Pilcher, G., J. Chem. Soc., Faraday Trans. 1 (1972), 68, 2224 (Idealized Complete Stoichiometric Combustion of Ethane, NIST CAS 74-84-0)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C74840",
    ),
    (
        "combustion", "C3H8 + 5O2", "3CO2 + 4H2O", "شرارة اشتعال في وفرة من الأكسجين",
        "theoretical_model", "NIST-WebBook-CAS-74-98-6-Propane-Combustion",
        "Rossini, F.D., J. Res. Natl. Bur. Stand. (1934), 12, 735–750 (Standard Enthalpy & Net Stoichiometry of Complete Propane Combustion, NIST CAS 74-98-6)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C74986",
    ),
    (
        "combustion", "2C4H10 + 13O2", "8CO2 + 10H2O", "شرارة اشتعال في وفرة من الأكسجين",
        "theoretical_model", "NIST-WebBook-CAS-106-97-8-Butane-Combustion",
        "Prosen, E.J., Maron, F.W., & Rossini, F.D., J. Res. Natl. Bur. Stand. (1951), 46, 106 (Idealized Complete Combustion of n-Butane, NIST CAS 106-97-8)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C106978",
    ),
    (
        "combustion", "C2H4 + 3O2", "2CO2 + 2H2O", "شرارة اشتعال في وفرة من الأكسجين",
        "theoretical_model", "NIST-WebBook-CAS-74-85-1-Ethylene-Combustion",
        "Rossini, F.D. & Knowlton, J.W., J. Res. Natl. Bur. Stand. (1937), 19, 249 (Thermochemical Stoichiometry of Ethylene Complete Combustion, NIST CAS 74-85-1)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C74851",
    ),
    (
        "combustion", "2C2H2 + 5O2", "4CO2 + 2H2O", "شرارة اشتعال في وفرة من الأكسجين",
        "theoretical_model", "NIST-WebBook-CAS-74-86-2-Acetylene-Combustion",
        "Coops, J. et al., Rec. Trav. Chim. Pays-Bas (1953), 72, 785 (Oxy-Acetylene Complete Stoichiometric Combustion Model, NIST CAS 74-86-2)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C74862",
    ),
    (
        "combustion", "C2H5OH + 3O2", "2CO2 + 3H2O", "شرارة اشتعال في وفرة من الأكسجين",
        "literature_referenced", "NIST-WebBook-CAS-64-17-5-Ethanol-Calorimetry",
        "Chao, J. & Rossini, F.D., J. Chem. Eng. Data (1965), 10, 374–379 (Bomb Calorimetric Complete Combustion of Ethanol, NIST CAS 64-17-5)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C64175",
    ),
    (
        "combustion", "2CH3OH + 3O2", "2CO2 + 4H2O", "شرارة اشتعال في وفرة من الأكسجين",
        "literature_referenced", "NIST-WebBook-CAS-67-56-1-Methanol-Calorimetry",
        "Chao, J. & Rossini, F.D., J. Chem. Eng. Data (1965), 10, 374–379 (Bomb Calorimetric Complete Combustion of Methanol, NIST CAS 67-56-1)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C67561",
    ),
    (
        "combustion", "2C6H6 + 15O2", "12CO2 + 6H2O", "اشتعال في وفرة من الأكسجين",
        "theoretical_model", "NIST-WebBook-CAS-71-43-2-Benzene-Combustion",
        "Good, W.D. & Smith, N.K., J. Chem. Thermodyn. (1979), 11, 111–118 (Idealized Complete Oxygen Bomb Combustion of Benzene, NIST CAS 71-43-2)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C71432",
    ),
    (
        "combustion", "2H2S + 3O2", "2SO2 + 2H2O", "اشتعال في وفرة من الأكسجين",
        "literature_referenced", "Ullmann-2012-Sulfur-Claus-Combustion-Step",
        "Nehb, W. & Vydra, K., 'Sulfur', Ullmann's Encyclopedia of Industrial Chemistry, Wiley-VCH (2012), Thermal Claus Oxidation Step of H2S",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C7783064",
    ),
    (
        "redox_synthesis", "C + O2", "CO2", "تسخين في وفرة من الأكسجين",
        "literature_referenced", "NIST-JANAF-CO2-Formation-CAS-124-38-9",
        "Chase, M.W., NIST-JANAF Thermochemical Tables, 4th Ed. (1998), p. 588 (Standard Formation of Carbon Dioxide from Graphite and O2)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C124389",
    ),
    (
        "precipitation", "AgNO3 + NaCl", "AgCl + NaNO3", "محلول مائي، راسب أبيض من AgCl",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-III.4-AgCl-NaCl",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. III.4, p. 210 (Group I Silver Chloride Precipitation)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/24561",
    ),
    (
        "precipitation", "AgNO3 + KCl", "AgCl + KNO3", "محلول مائي، راسب أبيض من AgCl",
        "experimentally_documented", "Vogel-Quant-5ed-Sec-10.80-Mohr-KCl-AgNO3",
        "Vogel's Textbook of Quantitative Chemical Analysis, 5th Ed. (1989), Sec. 10.80, p. 342 (Argentometric Mohr Titration of Potassium Chloride)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/4873",
    ),
    (
        "precipitation", "AgNO3 + HCl", "AgCl + HNO3", "محلول مائي، راسب أبيض من AgCl",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-IV.14-Chloride-AgNO3-Test",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. IV.14, p. 351 (Confirmatory Silver Nitrate Test for HCl)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/24561",
    ),
    (
        "precipitation", "BaCl2 + Na2SO4", "BaSO4 + 2NaCl", "محلول مائي، راسب أبيض من BaSO4",
        "experimentally_documented", "Vogel-Quant-5ed-Sec-11.49-Gravimetric-BaSO4",
        "Vogel's Textbook of Quantitative Chemical Analysis, 5th Ed. (1989), Sec. 11.49, p. 490 (Gravimetric Determination of Sulfate as Barium Sulfate)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/24414",
    ),
    (
        "precipitation", "BaCl2 + CuSO4", "BaSO4 + CuCl2", "محلول مائي، راسب أبيض من BaSO4",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-IV.24-Sulfate-BaCl2-Test",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. IV.24, p. 375 (Barium Chloride Test on Soluble Metal Sulfates)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/24462",
    ),
    (
        "precipitation", "CaCl2 + Na2CO3", "CaCO3 + 2NaCl", "محلول مائي، راسب أبيض من CaCO3",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-III.31-GroupIV-CaCO3",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. III.31, p. 291 (Precipitation of Calcium Carbonate)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/10112",
    ),
    (
        "precipitation", "FeCl3 + 3NaOH", "Fe(OH)3 + 3NaCl", "محلول مائي، راسب بني محمر",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-III.21-FeOH3-Precipitate",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. III.21, p. 257 (Reddish-Brown Iron(III) Hydroxide Precipitation)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/24380",
    ),
    (
        "precipitation", "CuSO4 + 2NaOH", "Cu(OH)2 + Na2SO4", "محلول مائي، راسب أزرق",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-III.10-CuOH2-Precipitate",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. III.10, p. 225 (Blue Copper(II) Hydroxide Precipitation)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/24462",
    ),
    (
        "single_displacement", "Zn + CuSO4", "ZnSO4 + Cu", "محلول مائي، 25 °C",
        "experimentally_documented", "CRC-104ed-Sec8-Daniell-Displacement-Zn-Cu",
        "Rumble, J.R. (Ed.), CRC Handbook of Chemistry and Physics, 104th Ed. (2023), Sec. 8: Electrochemical Series (E°cell = +1.10 V for Zn/Cu2+ Displacement)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/24424",
    ),
    (
        "single_displacement", "Fe + CuSO4", "FeSO4 + Cu", "محلول مائي، 25 °C",
        "experimentally_documented", "Greenwood-Earnshaw-2ed-Ch28-Copper-Cementation",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 28, p. 1175 (Hydrometallurgical Cementation of Copper by Metallic Iron)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/24393",
    ),
    (
        "single_displacement", "Zn + 2HCl", "ZnCl2 + H2", "محلول حمضي مخفف، تصاعد H2",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-III.26-Zn-Acid-Dissolution",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. III.26, p. 276 (Dissolution of Zinc in Dilute HCl with H2 Evolution)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/5727",
    ),
    (
        "single_displacement", "Mg + 2HCl", "MgCl2 + H2", "محلول حمضي مخفف، تصاعد H2",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-III.35-Mg-HCl-Dissolution",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. III.35, p. 304 (Reaction of Magnesium Metal with Dilute Hydrochloric Acid)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/5360315",
    ),
    (
        "single_displacement", "Fe + 2HCl", "FeCl2 + H2", "محلول حمضي مخفف، تصاعد H2",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-III.20-Fe-HCl-FeCl2",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. III.20, p. 254 (Formation of Iron(II) Chloride from Iron and Dilute HCl)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/24458",
    ),
    (
        "single_displacement", "2Al + 6HCl", "2AlCl3 + 3H2", "محلول حمضي، تصاعد H2",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-III.22-Al-HCl-Dissolution",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. III.22, p. 260 (Dissolution of Aluminium in Hydrochloric Acid)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/24012",
    ),
    (
        "single_displacement", "Zn + H2SO4", "ZnSO4 + H2", "حمض كبريتيك مخفف، تصاعد H2",
        "experimentally_documented", "Greenwood-Earnshaw-2ed-Ch3-Kipp-H2-Generation",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 3, p. 38 (Laboratory Preparation of Hydrogen from Zinc and Dilute Sulfuric Acid)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/24424",
    ),
    (
        "single_displacement", "2Na + 2H2O", "2NaOH + H2", "تفاعل عنيف طارد للحرارة",
        "experimentally_documented", "Greenwood-Earnshaw-2ed-Ch4-Alkali-Metal-Water",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 4, p. 76 (Reaction of Metallic Sodium with Water)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/14798",
    ),
    (
        "single_displacement", "Cu + 2AgNO3", "Cu(NO3)2 + 2Ag", "محلول مائي، ترسب الفضة",
        "experimentally_documented", "CRC-104ed-Sec8-Redox-Cu-AgNO3",
        "Rumble, J.R. (Ed.), CRC Handbook of Chemistry and Physics, 104th Ed. (2023), Sec. 8: Electrochemical Series (E°cell = +0.46 V for Cu/Ag+ Displacement)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/24561",
    ),
    (
        "redox_synthesis", "2Mg + O2", "2MgO", "اشتعال بلهب أبيض ساطع",
        "experimentally_documented", "NIST-JANAF-MgO-CAS-1309-48-4",
        "Chase, M.W., NIST-JANAF Thermochemical Tables, 4th Ed. (1998), p. 1536 (Combustion of Magnesium to Periclase MgO, ΔH°f = -601.6 kJ/mol)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C1309484",
    ),
    (
        "redox_synthesis", "2Na + Cl2", "2NaCl", "تفاعل طارد للحرارة",
        "literature_referenced", "NIST-JANAF-NaCl-CAS-7647-14-5",
        "Chase, M.W., NIST-JANAF Thermochemical Tables, 4th Ed. (1998), p. 780 (Direct Synthesis of Crystalline Sodium Chloride from Elements, ΔH°f = -411.2 kJ/mol)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C7647145",
    ),
    (
        "redox_synthesis", "2Ca + O2", "2CaO", "تسخين",
        "literature_referenced", "NIST-JANAF-CaO-CAS-1305-78-8",
        "Chase, M.W., NIST-JANAF Thermochemical Tables, 4th Ed. (1998), p. 705 (Direct Oxidation of Calcium Metal to Calcium Oxide, ΔH°f = -634.9 kJ/mol)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C1305788",
    ),
    (
        "redox_synthesis", "S + O2", "SO2", "اشتعال الكبريت",
        "literature_referenced", "NIST-JANAF-SO2-CAS-7446-09-5",
        "Chase, M.W., NIST-JANAF Thermochemical Tables, 4th Ed. (1998), p. 1750 (Combustion of Elemental Sulfur to Sulfur Dioxide, ΔH°f = -296.8 kJ/mol)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C7446095",
    ),
    (
        "redox_synthesis", "2SO2 + O2", "2SO3", "حفاز V2O5، تسخين",
        "experimentally_documented", "Ullmann-2012-SO2-V2O5-Catalytic-Oxidation",
        "Müller, H., 'Sulfuric Acid and Sulfur Trioxide', Ullmann's Encyclopedia of Industrial Chemistry (2012), Sec. 4.2: V2O5-Catalyzed Oxidation of SO2 at 420–450 °C",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C7446119",
    ),
    (
        "redox_synthesis", "2CO + O2", "2CO2", "اشتعال",
        "literature_referenced", "NIST-JANAF-CO-Oxidation-CAS-630-08-0",
        "Chase, M.W., NIST-JANAF Thermochemical Tables, 4th Ed. (1998), p. 586 (Gas-Phase Oxidation of Carbon Monoxide to Carbon Dioxide, ΔH° = -566.0 kJ)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C630080",
    ),
    (
        "redox_synthesis", "4Al + 3O2", "2Al2O3", "تسخين",
        "literature_referenced", "NIST-JANAF-Al2O3-Corundum-CAS-1344-28-1",
        "Chase, M.W., NIST-JANAF Thermochemical Tables, 4th Ed. (1998), p. 156 (Formation of Alpha-Alumina Al2O3 from Elements, ΔH°f = -1675.7 kJ/mol)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C1344281",
    ),
    (
        "redox_synthesis", "2Cu + O2", "2CuO", "تسخين في الهواء",
        "experimentally_documented", "Greenwood-Earnshaw-2ed-Ch28-CuO-Synthesis",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 28, p. 1181 (Thermal Oxidation of Copper Metal in Air at 300–800 °C to Tenorite CuO)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C1317380",
    ),
    (
        "redox_synthesis", "Fe2O3 + 3CO", "2Fe + 3CO2", "الفرن العالي، درجات حرارة مرتفعة",
        "literature_referenced", "Ullmann-2012-Iron-Blast-Furnace-Reduction",
        "Babich, A. et al., 'Iron, 2. Blast Furnace Process', Ullmann's Encyclopedia of Industrial Chemistry, Wiley-VCH (2012), Indirect Reduction of Hematite by CO at 800–1000 °C",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C1309371",
    ),
    (
        "redox_synthesis", "H2 + Cl2", "2HCl", "ضوء أو شرارة",
        "literature_referenced", "NIST-JANAF-HCl-Synthesis-CAS-7647-01-0",
        "Chase, M.W., NIST-JANAF Thermochemical Tables, 4th Ed. (1998), p. 765 (Photochemical / Thermal Chain Synthesis of Hydrogen Chloride, ΔH°f = -92.3 kJ/mol)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C7647010",
    ),
    (
        "redox_synthesis", "2NO + O2", "2NO2", "درجة حرارة الغرفة",
        "experimentally_documented", "Greenwood-Earnshaw-2ed-Ch11-NO-Oxidation",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 11, p. 455 (Spontaneous Third-Order Gas-Phase Oxidation of Nitric Oxide at 25 °C)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C10102440",
    ),
    (
        "redox_synthesis", "2Zn + O2", "2ZnO", "تسخين",
        "literature_referenced", "NIST-JANAF-ZnO-CAS-1314-13-2",
        "Chase, M.W., NIST-JANAF Thermochemical Tables, 4th Ed. (1998), p. 1766 (Thermal Oxidation of Zinc Vapor to Zinc Oxide, ΔH°f = -350.5 kJ/mol)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C1314132",
    ),
    (
        "redox_synthesis", "N2 + O2", "2NO", "درجات حرارة مرتفعة جدا أو برق",
        "literature_referenced", "NIST-JANAF-NO-Zeldovich-CAS-10102-43-9",
        "Chase, M.W., NIST-JANAF Thermochemical Tables, 4th Ed. (1998), p. 1666 (Endothermic High-Temperature Formation of Nitric Oxide, ΔH°f = +90.29 kJ/mol)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C10102439",
    ),
    (
        "decomposition", "CaCO3", "CaO + CO2", "تسخين عند نحو 850 °C",
        "experimentally_documented", "Ullmann-2012-Lime-Limestone-Calcination",
        "Oates, J.A.H., 'Lime and Limestone', Ullmann's Encyclopedia of Industrial Chemistry, Wiley-VCH (2012), Thermal Calcination of CaCO3 at 850–900 °C",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C471341",
    ),
    (
        "decomposition", "2KClO3", "2KCl + 3O2", "تسخين مع حفاز MnO2",
        "experimentally_documented", "Greenwood-Earnshaw-2ed-Ch14-KClO3-O2-Prep",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 14, p. 603 (MnO2-Catalyzed Thermal Decomposition of Potassium Chlorate at 250–350 °C)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/6426889",
    ),
    (
        "decomposition", "2KMnO4", "K2MnO4 + MnO2 + O2", "تسخين",
        "experimentally_documented", "Greenwood-Earnshaw-2ed-Ch24-KMnO4-Thermolysis",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 24, p. 1050 (Thermal Decomposition of Potassium Permanganate above 200 °C)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/516875",
    ),
    (
        "decomposition", "2H2O", "2H2 + O2", "تحليل كهربائي",
        "experimentally_documented", "Ullmann-2012-Water-Electrolysis-H2",
        "Häussinger, P. et al., 'Hydrogen, 2. Production', Ullmann's Encyclopedia of Industrial Chemistry (2012), Sec. 2.4: Alkaline and PEM Water Electrolysis (E°rev = 1.23 V)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C7732185",
    ),
    (
        "decomposition", "2NaHCO3", "Na2CO3 + H2O + CO2", "تسخين",
        "experimentally_documented", "Ullmann-2012-Soda-Ash-Calcination-NaHCO3",
        "Thieme, C., 'Sodium Carbonates', Ullmann's Encyclopedia of Industrial Chemistry (2012), Solvay Process Thermal Calcination of Sodium Bicarbonate at 160–230 °C",
        "https://pubchem.ncbi.nlm.nih.gov/compound/516892",
    ),
    (
        "decomposition", "H2CO3", "H2O + CO2", "درجة حرارة الغرفة",
        "experimentally_documented", "RHEA:16833-Dehydration",
        "Rhea Biochemical Database RHEA:16833 & Greenwood & Earnshaw (1997), p. 310 (Spontaneous Dehydration Equilibrium of Aqueous Carbonic Acid at 25 °C)",
        "https://www.rhea-db.org/rhea/16833",
    ),
    (
        "decomposition", "2NH3", "N2 + 3H2", "درجات حرارة مرتفعة مع حفاز",
        "literature_referenced", "Ullmann-2012-Ammonia-Catalytic-Cracking",
        "Appl, M., 'Ammonia', Ullmann's Encyclopedia of Industrial Chemistry (2012), Endothermic Catalytic Decomposition (Cracking) of Ammonia above 400 °C",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C7664417",
    ),
    (
        "decomposition", "Cu(OH)2", "CuO + H2O", "تسخين",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-III.10-CuOH2-Dehydration",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. III.10, p. 225 (Thermal Dehydration of Blue Cu(OH)2 to Black CuO on Boiling)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/164827",
    ),
    (
        "organic_reaction", "CH3COOH + C2H5OH", "CH3COOCH2CH3 + H2O", "حفاز H2SO4 المركز، تسخين (أسترة فيشر)",
        "experimentally_documented", "March-8ed-Rxn-16-24-Fischer-Ethyl-Acetate",
        "Smith, M.B., March's Advanced Organic Chemistry, 8th Ed., Wiley (2020), Reaction 16-24, p. 1196 & Vogel's Practical Organic Chemistry, 5th Ed., Exp. 5.141 (Ethyl Acetate Synthesis)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/8857",
    ),
    (
        "organic_reaction", "CH3COOH + CH3OH", "CH3COOCH3 + H2O", "حفاز H2SO4 المركز، تسخين (أسترة فيشر)",
        "experimentally_documented", "March-8ed-Rxn-16-24-Fischer-Methyl-Acetate",
        "Smith, M.B., March's Advanced Organic Chemistry, 8th Ed., Wiley (2020), Reaction 16-24, p. 1196 & Pöpken et al., Ind. Eng. Chem. Res. (2000), 39, 2601 (Methyl Acetate Esterification)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/6584",
    ),
    (
        "organic_reaction", "C2H4 + H2", "C2H6", "حفاز Pt أو Pd أو Ni (هدرجة)",
        "experimentally_documented", "March-8ed-Rxn-15-11-Alkene-Hydrogenation",
        "Smith, M.B., March's Advanced Organic Chemistry, 8th Ed., Wiley (2020), Reaction 15-11, p. 964 (Heterogeneous Catalytic Hydrogenation of Ethylene over Pt/Pd)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C74840",
    ),
    (
        "organic_reaction", "C2H4 + H2O", "C2H5OH", "حفاز H3PO4، ضغط وحرارة عاليان (إماهة)",
        "experimentally_documented", "Ullmann-2012-Ethanol-Direct-Hydration",
        "Kosaric, N. et al., 'Ethanol', Ullmann's Encyclopedia of Industrial Chemistry (2012), Sec. 3.1: Direct Catalytic Hydration of Ethylene over H3PO4/SiO2 at 300 °C",
        "https://pubchem.ncbi.nlm.nih.gov/compound/702",
    ),
    (
        "organic_reaction", "C2H2 + H2", "C2H4", "حفاز لندلار (هدرجة جزئية)",
        "experimentally_documented", "OrgSynth-1966-46-89-Lindlar-Alkyne-Semihydrogenation",
        "Lindlar, H. & Dubuis, R., Organic Syntheses, Coll. Vol. 5 (1973), p. 880; Vol. 46 (1966), p. 89 (Pd-CaCO3-PbO Catalyzed Semihydrogenation of Alkynes to Alkenes)",
        "https://webbook.nist.gov/cgi/cbook.cgi?ID=C74851",
    ),
    (
        "organic_reaction", "CH4 + Cl2", "CH3Cl + HCl", "ضوء فوق بنفسجي (استبدال جذري حر)",
        "experimentally_documented", "March-8ed-Rxn-14-1-Radical-Chlorination-Methane",
        "Smith, M.B., March's Advanced Organic Chemistry, 8th Ed., Wiley (2020), Reaction 14-1, p. 853 (Free-Radical Photochemical Chlorination of Methane)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/6327",
    ),
    (
        "organic_reaction", "C6H6 + HNO3", "C6H5NO2 + H2O", "حفاز H2SO4 المركز (نترتة)",
        "experimentally_documented", "Vogel-Org-5ed-Exp-6.11-Nitrobenzene-Nitration",
        "Vogel's Textbook of Practical Organic Chemistry, 5th Ed., Longman (1989), Exp. 6.11, p. 855 & March's 8th Ed. Rxn 11-3 (Electrophilic Mixed-Acid Nitration of Benzene at 50–60 °C)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/7416",
    ),
    (
        "organic_reaction", "C2H4 + Br2", "C2H4Br2", "إضافة إلى الرابطة المزدوجة، 25 °C",
        "experimentally_documented", "March-8ed-Rxn-15-24-Alkene-Bromination",
        "Smith, M.B., March's Advanced Organic Chemistry, 8th Ed., Wiley (2020), Reaction 15-24, p. 1005 (Electrophilic Anti-Addition of Bromine to Ethylene)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/7839",
    ),
    (
        "organic_reaction", "C6H5OH + NaOH", "C6H5ONa + H2O", "محلول مائي (حمضية الفينول)",
        "experimentally_documented", "Vogel-Org-5ed-Sec-6.4-Phenol-Alkali-Salt",
        "Vogel's Textbook of Practical Organic Chemistry, 5th Ed. (1989), Ch. 6, p. 684 (Deprotonation of Phenol in Aqueous Sodium Hydroxide, pKa = 9.95)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/996",
    ),
    (
        "organic_reaction", "C6H5NH2 + HCl", "C6H5NH3Cl", "محلول مائي (تكوين ملح الأمين)",
        "experimentally_documented", "Vogel-Org-5ed-Sec-6.22-Anilinium-Chloride",
        "Vogel's Textbook of Practical Organic Chemistry, 5th Ed. (1989), Ch. 6, p. 892 (Protonation of Aniline in Aqueous HCl to Anilinium Chloride, pKb = 9.40)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/6115",
    ),
    (
        "biochemical_reaction", "C6H12O6", "2C2H5OH + 2CO2", "إنزيمات الخميرة، وسط لاهوائي (تخمر كحولي)",
        "theoretical_model", "Lehninger-8ed-Sec-14.3-Alcoholic-Fermentation-Net",
        "Nelson, D.L. & Cox, M.M., Lehninger Principles of Biochemistry, 8th Ed. (2021), Sec. 14.3, p. 534 (Net Stoichiometric Balance of Yeast Anaerobic Alcoholic Fermentation)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/702",
    ),
    (
        "biochemical_reaction", "6CO2 + 6H2O", "C6H12O6 + 6O2", "طاقة ضوئية وصبغة الكلوروفيل (البناء الضوئي)",
        "theoretical_model", "Lehninger-8ed-Ch20-Oxygenic-Photosynthesis-Net",
        "Nelson, D.L. & Cox, M.M., Lehninger Principles of Biochemistry, 8th Ed. (2021), Ch. 20, p. 742 (Overall Net Stoichiometric Equation of Oxygenic Photosynthesis & Calvin Cycle)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/5793",
    ),
    (
        "biochemical_reaction", "C10H16N5O13P3 + H2O", "C10H15N5O10P2 + H3PO4", "إنزيمات ATPase، 37 °C (تحلل ATP)",
        "experimentally_documented", "RHEA:13065",
        "Rhea Biochemical Database RHEA:13065 (EC 3.6.1.3 Adenosine Triphosphatase: ATP + H2O = ADP + Phosphate, ΔG°' = -30.5 kJ/mol)",
        "https://www.rhea-db.org/rhea/13065",
    ),
    (
        "biochemical_reaction", "CO(NH2)2 + H2O", "CO2 + 2NH3", "إنزيم اليورياز (تحلل اليوريا)",
        "experimentally_documented", "RHEA:16113",
        "Rhea Biochemical Database RHEA:16113 (EC 3.5.1.5 Urease: Urea + H2O = CO2 + 2 NH3)",
        "https://www.rhea-db.org/rhea/16113",
    ),
    (
        "biochemical_reaction", "C12H22O11 + H2O", "2C6H12O6", "إنزيم السكراز أو حمض مخفف (تحلل السكروز)",
        "experimentally_documented", "RHEA:10068",
        "Rhea Biochemical Database RHEA:10068 (EC 3.2.1.26 Beta-fructofuranosidase / Invertase: Sucrose + H2O = D-Glucose + D-Fructose)",
        "https://www.rhea-db.org/rhea/10068",
    ),
    # --- Expanded Reactions: Synthesis, Decomposition, Single/Double Displacement, Acid-Base, Precipitation, Redox, Organic ---
    (
        "redox_synthesis", "Na2O + H2O", "2NaOH", "محلول مائي، تفاعل طارد للحرارة (إماهة أكسيد قاعدي)",
        "literature_referenced", "Greenwood-Earnshaw-2ed-Ch4-Na2O-Hydration",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 4, p. 84 (Vigorous Exothermic Hydration of Sodium Oxide to Sodium Hydroxide)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/73971",
    ),
    (
        "redox_synthesis", "K2O + H2O", "2KOH", "محلول مائي، تفاعل طارد للحرارة (إماهة أكسيد البوتاسيوم)",
        "literature_referenced", "Greenwood-Earnshaw-2ed-Ch4-K2O-Hydration",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 4, p. 85 (Complete Hydration of Potassium Oxide to Aqueous KOH)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/66242",
    ),
    (
        "redox_synthesis", "BaO + H2O", "Ba(OH)2", "محلول مائي، تفاعل طارد للحرارة (إماهة أكسيد الباريوم)",
        "literature_referenced", "Greenwood-Earnshaw-2ed-Ch5-BaO-Hydration",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 5, p. 119 (Exothermic Slaking of Barium Oxide with Water)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/62392",
    ),
    (
        "redox_synthesis", "SO2 + H2O", "H2SO3", "ذوبان غاز SO2 في الماء عند 25 °C",
        "literature_referenced", "Greenwood-Earnshaw-2ed-Ch15-SO2-Aqueous",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 15, p. 701 (Aqueous Hydration Equilibrium of Sulfur Dioxide to Sulfurous Acid)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/1100",
    ),
    (
        "decomposition", "2Ag2O", "4Ag + O2", "تسخين فوق 195 °C (تفكك حراري)",
        "experimentally_documented", "Greenwood-Earnshaw-2ed-Ch28-Ag2O-Thermolysis",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 28, p. 1183 (Thermal Decomposition of Silver(I) Oxide above 195 °C)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/9794626",
    ),
    (
        "decomposition", "2HgO", "2Hg + O2", "تسخين عند 500 °C (تجربة بريستلي ولافوازييه الكلاسيكية)",
        "experimentally_documented", "Greenwood-Earnshaw-2ed-Ch29-HgO-Thermolysis",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 29, p. 1211 (Thermal Dissociation of Mercury(II) Oxide into Liquid Mercury and Oxygen)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/16684736",
    ),
    (
        "decomposition", "NH4HCO3", "NH3 + H2O + CO2", "تسخين فوق 60 °C (تفكك حراري كامل لثلاثة نواتج)",
        "experimentally_documented", "Ullmann-2012-Ammonium-Compounds-NH4HCO3",
        "Zapp, K.H. et al., 'Ammonium Compounds', Ullmann's Encyclopedia of Industrial Chemistry (2012), Sec. 3: Thermal Dissociation of Ammonium Bicarbonate above 60 °C",
        "https://pubchem.ncbi.nlm.nih.gov/compound/14013",
    ),
    (
        "decomposition", "MgCO3", "MgO + CO2", "تسخين عند 350–550 °C (كلسنة المغنيسيت)",
        "experimentally_documented", "Ullmann-2012-Magnesium-Compounds-Calcination",
        "Seeger, M. et al., 'Magnesium Compounds', Ullmann's Encyclopedia of Industrial Chemistry (2012), Thermal Calcination of Magnesite MgCO3 to Magnesia MgO",
        "https://pubchem.ncbi.nlm.nih.gov/compound/11029",
    ),
    (
        "decomposition", "2Fe(OH)3", "Fe2O3 + 3H2O", "تسخين عند 200–500 °C (نزع ماء حراري)",
        "experimentally_documented", "Vogel-Quant-5ed-Sec-11.23-Fe2O3-Ignition",
        "Vogel's Textbook of Quantitative Chemical Analysis, 5th Ed. (1989), Sec. 11.23, p. 468 (Gravimetric Ignition of Hydrous Iron(III) Oxide to Fe2O3)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/14833",
    ),
    (
        "single_displacement", "Ba + 2H2O", "Ba(OH)2 + H2", "تفاعل سريع طارد للحرارة في الماء البارد",
        "experimentally_documented", "Greenwood-Earnshaw-2ed-Ch5-Barium-Water",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 5, p. 112 (Vigorous Reaction of Metallic Barium with Cold Water Forming Barium Hydroxide and H2)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/5355457",
    ),
    (
        "single_displacement", "Ca + 2H2O", "Ca(OH)2 + H2", "محلول مائي، تصاعد هادئ لغاز الهيدروجين عند 25 °C",
        "experimentally_documented", "Greenwood-Earnshaw-2ed-Ch5-Calcium-Water",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 5, p. 112 (Reaction of Metallic Calcium with Cold Water Forming Calcium Hydroxide)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/5460341",
    ),
    (
        "single_displacement", "Mg + CuSO4", "MgSO4 + Cu", "محلول مائي، 25 °C (إزاحة النحاس بواسطة المغنيسيوم)",
        "experimentally_documented", "CRC-104ed-Sec8-Redox-Mg-CuSO4",
        "Rumble, J.R. (Ed.), CRC Handbook of Chemistry and Physics, 104th Ed. (2023), Sec. 8: Electrochemical Series (E°cell = +2.71 V for Mg/Cu2+ Single Displacement)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/24083",
    ),
    (
        "single_displacement", "Cl2 + 2KBr", "2KCl + Br2", "محلول مائي، 25 °C (إزاحة الهالوجينات)",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-IV.15-Bromide-Cl2-Oxidation",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. IV.15, p. 355 (Oxidation of Aqueous Bromide Ion by Chlorine Water)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/253877",
    ),
    (
        "single_displacement", "Cl2 + 2KI", "2KCl + I2", "محلول مائي، 25 °C (إزاحة اليود بواسطة الكلور)",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-IV.16-Iodide-Cl2-Oxidation",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. IV.16, p. 358 (Displacement of Iodine from Potassium Iodide by Chlorine)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/4875",
    ),
    (
        "precipitation", "Pb(NO3)2 + Na2SO4", "PbSO4 + 2NaNO3", "محلول مائي، راسب أبيض بلوري من كبريتات الرصاص PbSO4",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-III.2-PbSO4-Precipitation",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. III.2, p. 207 (White Crystalline Lead(II) Sulfate Precipitation)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/24924",
    ),
    (
        "precipitation", "AgNO3 + KBr", "AgBr + KNO3", "محلول مائي، راسب أصفر باهت من بروميد الفضة AgBr",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-IV.15-AgBr-Precipitation",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. IV.15, p. 354 (Pale-Yellow Silver Bromide Precipitation)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/66199",
    ),
    (
        "precipitation", "AgNO3 + KI", "AgI + KNO3", "محلول مائي، راسب أصفر من يوديد الفضة AgI",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-IV.16-AgI-Precipitation",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. IV.16, p. 357 (Yellow Silver Iodide Precipitation Insoluble in NH3)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/24563",
    ),
    (
        "precipitation", "FeSO4 + 2NaOH", "Fe(OH)2 + Na2SO4", "محلول مائي، راسب أخضر من هيدروكسيد الحديدوز Fe(OH)2",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-III.20-FeOH2-Green-Precipitate",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. III.20, p. 255 (Green Iron(II) Hydroxide Precipitation)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/95553",
    ),
    (
        "precipitation", "MgSO4 + 2NaOH", "Mg(OH)2 + Na2SO4", "محلول مائي، راسب أبيض من هيدروكسيد المغنيسيوم Mg(OH)2",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-III.35-MgOH2-Precipitation",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. III.35, p. 305 (White Magnesium Hydroxide Precipitation with Alkali)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/73981",
    ),
    (
        "neutralization", "FeS + 2HCl", "FeCl2 + H2S", "محلول حمضي مخفف، تبادل مزدوج وتصاعد غاز H2S",
        "experimentally_documented", "Vogel-Qual-5ed-Sec-IV.6-FeS-HCl-H2S-Gen",
        "Vogel's Textbook of Macro and Semimicro Qualitative Inorganic Analysis, 5th Ed. (1979), Sec. IV.6, p. 339 (Kipp's Generation of Hydrogen Sulfide from Iron(II) Sulfide and Dilute HCl)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/14828",
    ),
    (
        "neutralization", "H3PO4 + 3NaOH", "Na3PO4 + 3H2O", "محلول مائي، 25 °C (معادلة حمض ثلاثي البروتون)",
        "experimentally_documented", "Vogel-Quant-5ed-Sec-10.39-H3PO4-Titration",
        "Vogel's Textbook of Quantitative Chemical Analysis, 5th Ed. (1989), Sec. 10.39, p. 302 (Neutralization of Orthophosphoric Acid with Sodium Hydroxide)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/24243",
    ),
    (
        "neutralization", "H2C2O4 + 2NaOH", "Na2C2O4 + 2H2O", "محلول مائي، 25 °C (معايرة حمض الأكساليك القياسية)",
        "experimentally_documented", "USP-NF-Volumetric-OxalicAcid-NaOH",
        "Vogel's Textbook of Quantitative Chemical Analysis, 5th Ed. (1989), Sec. 10.32, p. 296 (Primary Standard Oxalic Acid Standardization of Sodium Hydroxide)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/971",
    ),
    (
        "neutralization", "2HNO3 + Na2CO3", "2NaNO3 + H2O + CO2", "محلول مائي، تصاعد غاز CO2 (حمض النيتريك مع كربونات الصوديوم)",
        "experimentally_documented", "Vogel-Quant-5ed-Sec-10.35-Na2CO3-HNO3",
        "Vogel's Textbook of Quantitative Chemical Analysis, 5th Ed. (1989), Sec. 10.35, p. 299 (Acidimetric Neutralization of Sodium Carbonate with Nitric Acid)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/24268",
    ),
    (
        "neutralization", "2HNO3 + Ca(OH)2", "Ca(NO3)2 + 2H2O", "محلول مائي، 25 °C (معادلة حمض النيتريك بهيدروكسيد الكالسيوم)",
        "experimentally_documented", "Ullmann-2012-Nitrates-Calcium-Nitrate-Neutralization",
        "Laue, W. et al., 'Nitrates and Nitrites', Ullmann's Encyclopedia of Industrial Chemistry (2012), Sec. 3.3: Calcium Nitrate Production from Nitric Acid and Lime",
        "https://pubchem.ncbi.nlm.nih.gov/compound/24963",
    ),
    (
        "redox_synthesis", "CuO + H2", "Cu + H2O", "تسخين عند 250–300 °C (اختزال أكسيد النحاس بالهيدروجين)",
        "experimentally_documented", "Greenwood-Earnshaw-2ed-Ch28-CuO-H2-Reduction",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 28, p. 1181 (Thermal Reduction of Black Copper(II) Oxide by Hydrogen Stream to Red Metallic Copper)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/164827",
    ),
    (
        "single_displacement", "Fe2O3 + 2Al", "2Fe + Al2O3", "تسخين عالي (تفاعل الثرميت الطارد للحرارة بشدة > 2000 °C)",
        "experimentally_documented", "Greenwood-Earnshaw-2ed-Ch7-Goldschmidt-Thermite",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 7, p. 221 (Goldschmidt Aluminothermic Reduction of Iron(III) Oxide, ΔH° = -851.5 kJ/mol)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/14833",
    ),
    (
        "redox_synthesis", "MnO2 + 4HCl", "MnCl2 + Cl2 + 2H2O", "تسخين مع حمض HCl المركز (طريقة شيل لتحضير غاز الكلور)",
        "experimentally_documented", "Greenwood-Earnshaw-2ed-Ch17-Scheele-Cl2-Prep",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 17, p. 794 (Laboratory Preparation of Chlorine by Oxidation of Concentrated HCl with MnO2)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/24480",
    ),
    (
        "redox_synthesis", "4Fe + 3O2 + 6H2O", "4Fe(OH)3", "درجة حرارة الغرفة بوجود الأكسجين والماء (تفاعل من 3 مواد متفاعلة)",
        "literature_referenced", "Greenwood-Earnshaw-2ed-Ch25-Aqueous-Iron-Rusting",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 25, p. 1076 (Three-Reactant Electrochemical Corrosion of Iron in Oxygenated Water)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/23925",
    ),
    (
        "redox_synthesis", "2Cu + O2 + CO2 + H2O", "Cu2(OH)2CO3", "درجة حرارة الغرفة في الهواء الرطب (تفاعل من 4 مواد متفاعلة — تكون الملاكيت)",
        "literature_referenced", "Greenwood-Earnshaw-2ed-Ch28-Copper-Patina-Malachite",
        "Greenwood, N.N. & Earnshaw, A., Chemistry of the Elements, 2nd Ed. (1997), Ch. 28, p. 1176 (Four-Reactant Atmospheric Patination of Copper to Basic Copper Carbonate Cu2(OH)2CO3)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/14436",
    ),
    (
        "organic_reaction", "CH3COOCH2CH3 + H2O", "CH3COOH + C2H5OH", "حفاز H2SO4 المخفف، تسخين (تحلل مائي حمضي للإستر)",
        "experimentally_documented", "March-8ed-Rxn-16-59-Acid-Ester-Hydrolysis",
        "Smith, M.B., March's Advanced Organic Chemistry, 8th Ed., Wiley (2020), Reaction 16-59, p. 1221 (Acid-Catalyzed Hydrolysis of Ethyl Acetate to Acetic Acid and Ethanol)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/8857",
    ),
    (
        "organic_reaction", "CH3COOCH2CH3 + NaOH", "CH3COONa + C2H5OH", "محلول مائي، تسخين (تحلل مائي قاعدي / تصبن الإستر)",
        "experimentally_documented", "Vogel-Org-5ed-Sec-6.135-EthylAcetate-Saponification",
        "Vogel's Textbook of Practical Organic Chemistry, 5th Ed., Longman (1989), Sec. 6.135, p. 1218 (Irreversible Alkaline Hydrolysis / Saponification of Ethyl Acetate)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/8857",
    ),
    (
        "organic_reaction", "C2H4 + HCl", "C2H5Cl", "إضافة إلكتروفيلية لـ HCl إلى الرابطة المزدوجة (ألكين)",
        "experimentally_documented", "March-8ed-Rxn-15-2-Hydrochlorination-Ethylene",
        "Smith, M.B., March's Advanced Organic Chemistry, 8th Ed., Wiley (2020), Reaction 15-2, p. 948 (Electrophilic Hydrochlorination Addition of HCl to Ethylene)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/6337",
    ),
    (
        "organic_reaction", "C6H6 + Br2", "C6H5Br + HBr", "حفاز FeBr3، استبدال إلكتروفيلي عطري (برومة البنزين)",
        "experimentally_documented", "Vogel-Org-5ed-Exp-6.3-Bromobenzene-Substitution",
        "Vogel's Textbook of Practical Organic Chemistry, 5th Ed., Longman (1989), Exp. 6.3, p. 830 & March's 8th Ed. Rxn 11-10 (Electrophilic Aromatic Bromination of Benzene)",
        "https://pubchem.ncbi.nlm.nih.gov/compound/7964",
    ),
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


INSOLUBLE_SOLIDS_AND_METALS: set[str] = {
    "AgCl", "BaSO4", "CaCO3", "Fe(OH)3", "Fe(OH)2", "Cu(OH)2", "Mg(OH)2", "Al(OH)3",
    "PbI2", "PbSO4", "AgBr", "AgI", "FeS", "MgCO3", "Cu2(OH)2CO3", "Ag2O", "HgO",
    "CuO", "CaO", "MgO", "Al2O3", "Fe2O3", "ZnO", "MnO2", "Na2O", "K2O", "BaO",
    "Zn", "Fe", "Mg", "Al", "Cu", "Ag", "Na", "K", "Ca", "Ba", "C", "S",
}

BIOCHEMICAL_AQUEOUS_SOLUTES: set[str] = {
    "C6H12O6", "C2H5OH", "C10H16N5O13P3", "C10H15N5O10P2",
    "H3PO4", "CO(NH2)2", "C12H22O11", "NH3",
}

DOCUMENTED_THERMODYNAMICS_AR: dict[str, str] = {
    "NIST-CRC-Neutralization-KOH-HCl": "طارد للحرارة (ΔH° = -57.3 kJ/mol، حرارة التعادل القياسية، CRC 104th Ed.)",
    "NIST-CRC-Neutralization-NaOH-HNO3": "طارد للحرارة (ΔH° = -57.3 kJ/mol، حرارة التعادل القياسية، CRC 104th Ed.)",
    "NIST-CRC-Neutralization-KOH-HNO3": "طارد للحرارة (ΔH° = -57.3 kJ/mol، حرارة التعادل القياسية، CRC 104th Ed.)",
    "CRC-104ed-Sec5-Neutralization-NaOH-H2SO4": "طارد للحرارة (ΔH° = -114.6 kJ/mol، تعادل 2 مول OH⁻ مع H₂SO₄، CRC 104th Ed.)",
    "CRC-104ed-Sec5-Neutralization-KOH-H2SO4": "طارد للحرارة (ΔH° = -114.6 kJ/mol، تعادل 2 مول OH⁻ مع H₂SO₄، CRC 104th Ed.)",
    "Vogel-Qual-5ed-MgOH2-HCl": "طارد للحرارة (ΔH° = -111.0 kJ/mol، معادلة هيدروكسيد المغنيسيوم، NIST-JANAF)",
    "Vogel-Qual-5ed-CaOH2-HCl": "طارد للحرارة (ΔH° = -114.4 kJ/mol، معادلة هيدروكسيد الكالسيوم، NIST-JANAF)",
    "NIST-JANAF-NH4Cl-Formation": "طارد للحرارة (ΔH° = -176.0 kJ/mol، اتحاد غازي HCl وNH3 لتكوين NH4Cl(s)، NIST-JANAF)",
    "CRC-104ed-Sec5-AceticAcid-NaOH": "طارد للحرارة (ΔH° = -55.8 kJ/mol، معادلة حمض ضعيف بقاعدة قوية، CRC 104th Ed.)",
    "Vogel-Quant-5ed-Sec-10.35-Na2CO3-HCl": "طارد للحرارة (ΔH° = -29.4 kJ/mol، معايرة كربونات الصوديوم بحمض الهيدروكلوريك، NIST)",
    "Vogel-Quant-5ed-Sec-10.36-NaHCO3-HCl": "ماص للحرارة طفيفاً (ΔH° = +11.8 kJ/mol، تفاعل بيكربونات الصوديوم مع HCl وتصاعد CO2، NIST)",
    "Vogel-Qual-5ed-CaCO3-H2SO4": "طارد للحرارة (ΔH° = -15.2 kJ/mol، تفاعل كربونات الكالسيوم مع حمض الكبريتيك، NIST)",
    "Vogel-Qual-5ed-CuO-HCl": "طارد للحرارة (ΔH° = -63.5 kJ/mol، إذابة أكسيد النحاس الثنائي في HCl، NIST)",
    "Vogel-Qual-5ed-CuO-H2SO4": "طارد للحرارة (ΔH° = -75.4 kJ/mol، إذابة أكسيد النحاس الثنائي في H2SO4، NIST)",
    "Ullmann-Lime-Slaking-CaO-H2O": "طارد للحرارة (ΔH° = -65.2 kJ/mol، إطفاء الجير الحي CaO بالماء، Ullmann 2012)",
    "Ullmann-2012-SulfuricAcid-SO3-Hydration": "طارد للحرارة بشدة (ΔH° = -132.5 kJ/mol، إماهة ثالث أكسيد الكبريت، Ullmann 2012)",
    "RHEA:16833-Carbonic-Acid-Hydration": "طارد للحرارة طفيفاً (ΔH° = -20.4 kJ/mol، إماهة ثاني أكسيد الكربون في الماء عند 25 °C)",
    "NIST-WebBook-CAS-74-84-0-Ethane-Combustion": "طارد للحرارة بشدة (ΔH°c = -3120 kJ لمعادلة 2 مول إيثان، NIST WebBook)",
    "NIST-WebBook-CAS-74-98-6-Propane-Combustion": "طارد للحرارة بشدة (ΔH°c = -2220 kJ/mol، احتراق البروبان القياسي، NIST WebBook)",
    "NIST-WebBook-CAS-106-97-8-Butane-Combustion": "طارد للحرارة بشدة (ΔH°c = -5755 kJ لمعادلة 2 مول بيوتان، NIST WebBook)",
    "NIST-WebBook-CAS-74-85-1-Ethylene-Combustion": "طارد للحرارة بشدة (ΔH°c = -1411 kJ/mol، احتراق الإيثيلين القياسي، NIST WebBook)",
    "NIST-WebBook-CAS-74-86-2-Acetylene-Combustion": "طارد للحرارة بشدة (ΔH°c = -2600 kJ لمعادلة 2 مول أسيتيلين، NIST WebBook)",
    "NIST-WebBook-CAS-64-17-5-Ethanol-Calorimetry": "طارد للحرارة بشدة (ΔH°c = -1367 kJ/mol، احتراق الإيثانول السائل القياسي، NIST WebBook)",
    "NIST-WebBook-CAS-67-56-1-Methanol-Calorimetry": "طارد للحرارة بشدة (ΔH°c = -1452 kJ لمعادلة 2 مول ميثانول، NIST WebBook)",
    "NIST-WebBook-CAS-71-43-2-Benzene-Combustion": "طارد للحرارة بشدة (ΔH°c = -6535 kJ لمعادلة 2 مول بنزين، NIST WebBook)",
    "Ullmann-2012-Sulfur-Claus-Combustion-Step": "طارد للحرارة بشدة (ΔH° = -1036 kJ لمعادلة 2 مول H2S، عملية كلاوس، Ullmann 2012)",
    "NIST-JANAF-CO2-Formation-CAS-124-38-9": "طارد للحرارة (ΔH°f = -393.5 kJ/mol، تكوين CO2 القياسي من الجرافيت والأكسجين، NIST-JANAF)",
    "CRC-104ed-Sec8-Daniell-Displacement-Zn-Cu": "طارد للحرارة (ΔH° = -218.7 kJ/mol، إزاحة النحاس بواسطة الخارصين، CRC 104th Ed.)",
    "Greenwood-Earnshaw-2ed-Ch28-Copper-Cementation": "طارد للحرارة (ΔH° = -153.8 kJ/mol، إزاحة النحاس بواسطة الحديد، Greenwood & Earnshaw)",
    "Vogel-Qual-5ed-Sec-III.26-Zn-Acid-Dissolution": "طارد للحرارة (ΔH° = -153.9 kJ/mol، تفاعل الخارصين مع HCl المخفف، NIST-JANAF)",
    "Vogel-Qual-5ed-Sec-III.35-Mg-HCl-Dissolution": "طارد للحرارة بشدة (ΔH° = -466.9 kJ/mol، تفاعل المغنيسيوم مع HCl المخفف، NIST-JANAF)",
    "Vogel-Qual-5ed-Sec-III.20-Fe-HCl-FeCl2": "طارد للحرارة (ΔH° = -89.1 kJ/mol، تفاعل الحديد مع HCl المخفف، NIST-JANAF)",
    "Vogel-Qual-5ed-Sec-III.22-Al-HCl-Dissolution": "طارد للحرارة بشدة (ΔH° = -1062.7 kJ لمعادلة 2 مول Al مع 6 مول HCl، NIST-JANAF)",
    "Greenwood-Earnshaw-2ed-Ch3-Kipp-H2-Generation": "طارد للحرارة (ΔH° = -153.9 kJ/mol، تفاعل الخارصين مع حمض الكبريتيك المخفف، NIST)",
    "Greenwood-Earnshaw-2ed-Ch4-Alkali-Metal-Water": "طارد للحرارة بشدة (ΔH° = -368.6 kJ لمعادلة 2 مول صوديوم مع الماء، Greenwood & Earnshaw)",
    "CRC-104ed-Sec8-Redox-Cu-AgNO3": "طارد للحرارة (ΔH° = -147.4 kJ/mol، إزاحة الفضة بواسطة النحاس، CRC 104th Ed.)",
    "NIST-JANAF-MgO-CAS-1309-48-4": "طارد للحرارة بشدة (ΔH° = -1203.2 kJ لمعادلة 2 مول Mg، ΔH°f = -601.6 kJ/mol، NIST-JANAF)",
    "NIST-JANAF-NaCl-CAS-7647-14-5": "طارد للحرارة بشدة (ΔH° = -822.4 kJ لمعادلة 2 مول NaCl، ΔH°f = -411.2 kJ/mol، NIST-JANAF)",
    "NIST-JANAF-CaO-CAS-1305-78-8": "طارد للحرارة بشدة (ΔH° = -1269.8 kJ لمعادلة 2 مول CaO، ΔH°f = -634.9 kJ/mol، NIST-JANAF)",
    "NIST-JANAF-SO2-CAS-7446-09-5": "طارد للحرارة (ΔH°f = -296.8 kJ/mol، احتراق الكبريت إلى SO2، NIST-JANAF)",
    "Ullmann-2012-SO2-V2O5-Catalytic-Oxidation": "طارد للحرارة (ΔH° = -197.8 kJ لمعادلة 2 مول SO2 إلى 2 مول SO3، Ullmann 2012)",
    "NIST-JANAF-CO-Oxidation-CAS-630-08-0": "طارد للحرارة (ΔH° = -566.0 kJ لمعادلة 2 مول CO، NIST-JANAF)",
    "NIST-JANAF-Al2O3-Corundum-CAS-1344-28-1": "طارد للحرارة بشدة (ΔH° = -3351.4 kJ لمعادلة 2 مول Al2O3، ΔH°f = -1675.7 kJ/mol، NIST-JANAF)",
    "Greenwood-Earnshaw-2ed-Ch28-CuO-Synthesis": "طارد للحرارة (ΔH° = -314.6 kJ لمعادلة 2 مول CuO، NIST-JANAF)",
    "Ullmann-2012-Iron-Blast-Furnace-Reduction": "طارد للحرارة طفيفاً (ΔH° = -24.8 kJ/mol، اختزال الهيماتيت بغاز CO في الفرن العالي، Ullmann 2012)",
    "NIST-JANAF-HCl-Synthesis-CAS-7647-01-0": "طارد للحرارة (ΔH° = -184.6 kJ لمعادلة 2 مول HCl، ΔH°f = -92.3 kJ/mol، NIST-JANAF)",
    "Greenwood-Earnshaw-2ed-Ch11-NO-Oxidation": "طارد للحرارة (ΔH° = -114.1 kJ لمعادلة 2 مول NO، NIST-JANAF)",
    "NIST-JANAF-ZnO-CAS-1314-13-2": "طارد للحرارة (ΔH° = -701.0 kJ لمعادلة 2 مول ZnO، ΔH°f = -350.5 kJ/mol، NIST-JANAF)",
    "NIST-JANAF-NO-Zeldovich-CAS-10102-43-9": "ماص للحرارة بشدة (ΔH° = +180.6 kJ لمعادلة 2 مول NO، ΔH°f = +90.29 kJ/mol، NIST-JANAF)",
    "Ullmann-2012-Lime-Limestone-Calcination": "ماص للحرارة بشدة (ΔH° = +178.3 kJ/mol، الكلسنة الحرارية للحجر الجيري، Ullmann 2012)",
    "Greenwood-Earnshaw-2ed-Ch14-KClO3-O2-Prep": "طارد للحرارة (ΔH° = -89.4 kJ لمعادلة 2 مول KClO3، تفكك كللورات البوتاسيوم المحفز بـ MnO2)",
    "Greenwood-Earnshaw-2ed-Ch24-KMnO4-Thermolysis": "ماص للحرارة (ΔH° = +40.2 kJ لمعادلة 2 مول KMnO4، التفكك الحراري لبرمنغنات البوتاسيوم)",
    "Ullmann-2012-Water-Electrolysis-H2": "ماص للحرارة بشدة (ΔH° = +571.6 kJ لمعادلة 2 مول H2O، التحليل الكهربائي للماء، Ullmann 2012)",
    "Ullmann-2012-Soda-Ash-Calcination-NaHCO3": "ماص للحرارة (ΔH° = +135.6 kJ لمعادلة 2 مول NaHCO3، عملية سولفاي، Ullmann 2012)",
    "RHEA:16833-Dehydration": "ماص للحرارة طفيفاً (ΔH° = +20.4 kJ/mol، تفكك حمض الكربونيك في المحلول المائي إلى CO2 وH2O)",
    "Ullmann-2012-Ammonia-Catalytic-Cracking": "ماص للحرارة (ΔH° = +92.2 kJ لمعادلة 2 مول NH3، تكسير الأمونيا الحراري، Ullmann 2012)",
    "Vogel-Qual-5ed-Sec-III.10-CuOH2-Dehydration": "ماص للحرارة (ΔH° = +45.5 kJ/mol، نزع الماء الحراري من هيدروكسيد النحاس الثنائي)",
    "March-8ed-Rxn-16-24-Fischer-Ethyl-Acetate": "قريب من الحياد الحراري (ΔH° = -3.8 kJ/mol، تفاعل أسترة فيشر العكوس، March 8th Ed.)",
    "March-8ed-Rxn-16-24-Fischer-Methyl-Acetate": "قريب من الحياد الحراري (ΔH° = -4.2 kJ/mol، تفاعل أسترة فيشر العكوس، March 8th Ed.)",
    "March-8ed-Rxn-15-11-Alkene-Hydrogenation": "طارد للحرارة (ΔH° = -136.3 kJ/mol، هدرجة الإيثيلين الحفزية، NIST WebBook)",
    "Ullmann-2012-Ethanol-Direct-Hydration": "طارد للحرارة (ΔH° = -45.8 kJ/mol، الإماهة الحفزية للإيثيلين، Ullmann 2012)",
    "OrgSynth-1966-46-89-Lindlar-Alkyne-Semihydrogenation": "طارد للحرارة (ΔH° = -175.5 kJ/mol، الهدرجة الجزئية للأسيتيلين فوق حفاز لندلار)",
    "March-8ed-Rxn-14-1-Radical-Chlorination-Methane": "طارد للحرارة (ΔH° = -105.0 kJ/mol، الكلورة الجذرية للميثان، March 8th Ed.)",
    "Vogel-Org-5ed-Exp-6.11-Nitrobenzene-Nitration": "طارد للحرارة (ΔH° = -117.0 kJ/mol، نترتة البنزين بحمض النيتريك والكبريتيك، Vogel 5th Ed.)",
    "March-8ed-Rxn-15-24-Alkene-Bromination": "طارد للحرارة (ΔH° = -121.0 kJ/mol، إضافة البروم الإلكتروفيلية إلى الإيثيلين، March 8th Ed.)",
    "Vogel-Org-5ed-Sec-6.4-Phenol-Alkali-Salt": "طارد للحرارة (ΔH° = -35.4 kJ/mol، تعادل الفينول مع هيدروكسيد الصوديوم المائي)",
    "Vogel-Org-5ed-Sec-6.22-Anilinium-Chloride": "طارد للحرارة (ΔH° = -31.2 kJ/mol، برتنة الأنيلين في حمض الهيدروكلوريك المائي)",
    "Lehninger-8ed-Sec-14.3-Alcoholic-Fermentation-Net": "طارد للحرارة (ΔH° = -74.4 kJ/mol، التخمر الكحولي اللاهوائي للجلوكوز، Lehninger 8th Ed.)",
    "Lehninger-8ed-Ch20-Oxygenic-Photosynthesis-Net": "ماص للطاقة بشدة (ΔH° = +2803 kJ/mol، البناء الضوئي الأكسجيني من CO2 والماء، Lehninger 8th Ed.)",
    "RHEA:13065": "طارد للحرارة (ΔH° = -20.5 kJ/mol، ΔG°' = -30.5 kJ/mol، الحلمأة الإنزيمية لـ ATP، RHEA:13065)",
    "RHEA:16113": "طارد للحرارة طفيفاً (ΔH° = -13.8 kJ/mol، تحلل اليوريا المائي بواسطة اليورياز، RHEA:16113)",
    "RHEA:10068": "طارد للحرارة طفيفاً (ΔH° = -14.9 kJ/mol، حلمأة السكروز بواسطة السكراز، RHEA:10068)",
}


def _infer_species_state(formula: str, category: str, cond: str) -> str:
    if category == "biochemical_reaction":
        if formula == "H2O":
            return "l"
        if formula in {"CO2", "O2"}:
            return "g"
        if formula in BIOCHEMICAL_AQUEOUS_SOLUTES:
            return "aq"
    if formula == "SO3" and category == "redox_synthesis":
        return "g"
    if formula == "H2CO3":
        return "aq"
    if "طور غازي" in cond and formula in {"HCl", "NH3"}:
        return "g"
    phase = get_standard_phase_code(formula)
    is_aqueous_medium = (
        any(k in cond for k in ("محلول", "مخفف", "الماء", "المحاليل"))
        or (
            category in {"neutralization", "precipitation", "single_displacement"}
            and "طور غازي" not in cond
            and "إطفاء الجير" not in cond
        )
    )
    if phase in {"g", "l"}:
        if formula in {"HCl", "HNO3", "H2SO4", "CH3COOH"} and is_aqueous_medium:
            return "aq"
        return phase
    if is_aqueous_medium and formula not in INSOLUBLE_SOLIDS_AND_METALS:
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


def _extract_conditions_dict(category: str, cond_str: str, source_id: str = "") -> dict[str, str]:
    cond_dict: dict[str, str] = {
        "temperature": cond_str,
        "pressure": "1 atm ما لم يذكر خلاف ذلك",
    }
    if any(k in cond_str for k in ("محلول مائي", "محلول حمضي", "حمض كبريتيك مخفف", "الماء", "المحاليل")) or category == "biochemical_reaction":
        cond_dict["solvent"] = "محلول مائي (Aqueous Solution)"
    elif "طور غازي" in cond_str or category == "combustion":
        cond_dict["solvent"] = "طور غازي أو تفاعل مباشر بدون مذيب"
    if "حفاز" in cond_str or "إنزيم" in cond_str:
        cond_dict["catalyst"] = cond_str
    if source_id and source_id in DOCUMENTED_THERMODYNAMICS_AR:
        cond_dict["thermodynamics_ar"] = DOCUMENTED_THERMODYNAMICS_AR[source_id]
    return cond_dict


def _build_reference_conditions_spec(category: str, lhs: str, cond_str: str) -> dict[str, Any]:
    if category == "biochemical_reaction":
        return {
            "temp_c_min": 15.0,
            "temp_c_max": 45.0,
            "temp_c_optimal": 37.0,
            "pressure_atm_min": 0.8,
            "pressure_atm_max": 2.0,
            "pressure_atm_optimal": 1.0,
            "allowed_solvents": ["aqueous", "physiological_buffer"],
            "required_catalysts": ["enzymatic"],
            "catalyst_is_mandatory": True,
        }
    if category == "combustion":
        return {
            "temp_c_min": 350.0,
            "temp_c_max": 2500.0,
            "temp_c_optimal": 850.0,
            "pressure_atm_min": 0.5,
            "pressure_atm_max": 20.0,
            "pressure_atm_optimal": 1.0,
            "allowed_solvents": ["neat_gas_solid"],
            "required_catalysts": ["none"],
            "catalyst_is_mandatory": False,
        }
    if category in {"neutralization", "precipitation", "single_displacement"}:
        allowed_solv = ["aqueous", "neat_gas_solid"] if ("طور غازي" in cond_str or "إطفاء الجير" in cond_str) else ["aqueous"]
        return {
            "temp_c_min": 5.0,
            "temp_c_max": 80.0,
            "temp_c_optimal": 25.0,
            "pressure_atm_min": 0.5,
            "pressure_atm_max": 5.0,
            "pressure_atm_optimal": 1.0,
            "allowed_solvents": allowed_solv,
            "required_catalysts": ["none"],
            "catalyst_is_mandatory": False,
        }
    if category == "organic_reaction":
        if "نترتة" in cond_str:
            return {
                "temp_c_min": 40.0,
                "temp_c_max": 65.0,
                "temp_c_optimal": 55.0,
                "pressure_atm_min": 0.8,
                "pressure_atm_max": 5.0,
                "pressure_atm_optimal": 1.0,
                "allowed_solvents": ["neat_gas_solid", "aqueous"],
                "required_catalysts": ["acid_h2so4"],
                "catalyst_is_mandatory": True,
            }
        if "H2SO4" in cond_str:
            return {
                "temp_c_min": 45.0,
                "temp_c_max": 130.0,
                "temp_c_optimal": 80.0,
                "pressure_atm_min": 0.8,
                "pressure_atm_max": 5.0,
                "pressure_atm_optimal": 1.0,
                "allowed_solvents": ["neat_gas_solid", "ethanol", "aqueous"],
                "required_catalysts": ["acid_h2so4"],
                "catalyst_is_mandatory": True,
            }
        if "Pt" in cond_str or "Pd" in cond_str or "لندلار" in cond_str:
            return {
                "temp_c_min": 15.0,
                "temp_c_max": 180.0,
                "temp_c_optimal": 50.0,
                "pressure_atm_min": 1.0,
                "pressure_atm_max": 50.0,
                "pressure_atm_optimal": 5.0,
                "allowed_solvents": ["neat_gas_solid", "ethanol"],
                "required_catalysts": ["pt_pd"],
                "catalyst_is_mandatory": True,
            }
        if "ضغط وحرارة عاليان" in cond_str:
            return {
                "temp_c_min": 200.0,
                "temp_c_max": 350.0,
                "temp_c_optimal": 300.0,
                "pressure_atm_min": 5.0,
                "pressure_atm_max": 100.0,
                "pressure_atm_optimal": 60.0,
                "allowed_solvents": ["neat_gas_solid", "aqueous"],
                "required_catalysts": ["acid_h2so4", "none"],
                "catalyst_is_mandatory": False,
            }
        return {
            "temp_c_min": 10.0,
            "temp_c_max": 100.0,
            "temp_c_optimal": 25.0,
            "pressure_atm_min": 0.8,
            "pressure_atm_max": 5.0,
            "pressure_atm_optimal": 1.0,
            "allowed_solvents": ["aqueous", "ethanol", "neat_gas_solid"],
            "required_catalysts": ["none"],
            "catalyst_is_mandatory": False,
        }
    if category == "decomposition":
        if "850" in cond_str or "مرتفعة" in cond_str:
            return {
                "temp_c_min": 400.0,
                "temp_c_max": 1200.0,
                "temp_c_optimal": 850.0,
                "pressure_atm_min": 0.5,
                "pressure_atm_max": 10.0,
                "pressure_atm_optimal": 1.0,
                "allowed_solvents": ["neat_gas_solid"],
                "required_catalysts": ["none", "fe_haber"],
                "catalyst_is_mandatory": False,
            }
        if "MnO2" in cond_str:
            return {
                "temp_c_min": 180.0,
                "temp_c_max": 450.0,
                "temp_c_optimal": 300.0,
                "pressure_atm_min": 0.5,
                "pressure_atm_max": 5.0,
                "pressure_atm_optimal": 1.0,
                "allowed_solvents": ["neat_gas_solid"],
                "required_catalysts": ["mno2_catalase"],
                "catalyst_is_mandatory": True,
            }
        if "تسخين" in cond_str:
            return {
                "temp_c_min": 80.0,
                "temp_c_max": 600.0,
                "temp_c_optimal": 200.0,
                "pressure_atm_min": 0.5,
                "pressure_atm_max": 5.0,
                "pressure_atm_optimal": 1.0,
                "allowed_solvents": ["neat_gas_solid", "aqueous"],
                "required_catalysts": ["none"],
                "catalyst_is_mandatory": False,
            }
        return {
            "temp_c_min": 10.0,
            "temp_c_max": 80.0,
            "temp_c_optimal": 25.0,
            "pressure_atm_min": 0.5,
            "pressure_atm_max": 5.0,
            "pressure_atm_optimal": 1.0,
            "allowed_solvents": ["aqueous"],
            "required_catalysts": ["none"],
            "catalyst_is_mandatory": False,
        }

    # redox_synthesis
    if lhs == "N2 + O2":
        return {
            "temp_c_min": 1000.0,
            "temp_c_max": 3000.0,
            "temp_c_optimal": 1800.0,
            "pressure_atm_min": 0.5,
            "pressure_atm_max": 20.0,
            "pressure_atm_optimal": 1.0,
            "allowed_solvents": ["neat_gas_solid"],
            "required_catalysts": ["none"],
            "catalyst_is_mandatory": False,
        }
    if lhs == "Fe2O3 + 3CO":
        return {
            "temp_c_min": 600.0,
            "temp_c_max": 1200.0,
            "temp_c_optimal": 900.0,
            "pressure_atm_min": 0.8,
            "pressure_atm_max": 10.0,
            "pressure_atm_optimal": 2.0,
            "allowed_solvents": ["neat_gas_solid"],
            "required_catalysts": ["none"],
            "catalyst_is_mandatory": False,
        }
    if lhs == "2SO2 + O2":
        return {
            "temp_c_min": 380.0,
            "temp_c_max": 550.0,
            "temp_c_optimal": 450.0,
            "pressure_atm_min": 1.0,
            "pressure_atm_max": 10.0,
            "pressure_atm_optimal": 2.0,
            "allowed_solvents": ["neat_gas_solid"],
            "required_catalysts": ["v2o5", "pt_pd"],
            "catalyst_is_mandatory": True,
        }
    if any(k in cond_str for k in ["تسخين", "اشتعال", "الفرن العالي", "مرتفعة"]):
        return {
            "temp_c_min": 200.0,
            "temp_c_max": 2000.0,
            "temp_c_optimal": 500.0,
            "pressure_atm_min": 0.5,
            "pressure_atm_max": 20.0,
            "pressure_atm_optimal": 1.0,
            "allowed_solvents": ["neat_gas_solid"],
            "required_catalysts": ["none", "pt_pd", "fe_haber"],
            "catalyst_is_mandatory": False,
        }
    return {
        "temp_c_min": 5.0,
        "temp_c_max": 100.0,
        "temp_c_optimal": 25.0,
        "pressure_atm_min": 0.5,
        "pressure_atm_max": 10.0,
        "pressure_atm_optimal": 1.0,
        "allowed_solvents": ["aqueous", "neat_gas_solid"],
        "required_catalysts": ["none"],
        "catalyst_is_mandatory": False,
    }


def _build_generated_reactions() -> list[dict[str, Any]]:
    index = _species_index()
    out = []
    for idx, (category, lhs, rhs, cond, doc_type, source_id, specific_source, ref_url) in enumerate(_GENERATED_REACTIONS, start=1):
        reactants = _parse_side(lhs, index, category, cond)
        products = _parse_side(rhs, index, category, cond)
        domain, domain_label = DOMAIN_BY_CATEGORY[category]
        label = CATEGORY_LABELS_AR[category]
        r_ar = " + ".join(r["name_ar"] for r in reactants)
        p_ar = " + ".join(p["name_ar"] for p in products)
        r_en = " + ".join(r["name_en"] for r in reactants)
        p_en = " + ".join(p["name_en"] for p in products)
        ratio = ":".join(str(r["coefficient"]) for r in reactants)
        doc_label_ar = DOCUMENTATION_TYPE_LABELS_AR[doc_type]
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
                "conditions": _extract_conditions_dict(category, cond, source_id),
                "reference_conditions_spec": _build_reference_conditions_spec(category, lhs, cond),
                "explanation_ar": f"تتفاعل المواد ({r_ar}) بنسبة مولية ({ratio}) لتكوين ({p_ar}) وفق المعادلة الموزونة المسجلة في السجل المعتمد.",
                "assumptions_and_limits_ar": (
                    "نموذج ستوكيومتري نظري يفترض اكتمال التفاعل بنسبة 100% دون نواتج جانبية."
                    if doc_type == "theoretical_model"
                    else "حساب ستوكيومتري نظري بفرض اكتمال التفاعل بنسبة 100% دون تفاعلات جانبية، وهو ليس إثباتا لتجربة عملية."
                ),
                "documentation_type": doc_type,
                "documentation_label_ar": doc_label_ar,
                "source_id": source_id,
                "source": specific_source,
                "reference_url": ref_url,
            }
        )
    return out


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

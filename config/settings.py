import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
LOCAL_SQLITE_DB_PATH = os.path.join(BASE_DIR, "data", "moleculens_local.db")

APP_NAME = "MolecuLens"
APP_SUBTITLE = "Scientific Learning & Exploration Platform"
APP_SLOGAN = "See Molecules. Understand Science."
TEAM_NAME = "BioCode"

ENTITY_CHEMICAL_ELEMENT = "chemical_element"
ENTITY_CHEMICAL_COMPOUND = "chemical_compound"
ENTITY_DRUG_SUBSTANCE = "drug_substance"
ENTITY_BIOMOLECULE_PROTEIN = "biomolecule_protein"

ENTITY_TYPE_LABELS_AR = {
    ENTITY_CHEMICAL_ELEMENT: "عنصر كيميائي",
    ENTITY_CHEMICAL_COMPOUND: "مركب كيميائي",
    ENTITY_DRUG_SUBSTANCE: "مادة دوائية",
    ENTITY_BIOMOLECULE_PROTEIN: "بروتين أو جزيء حيوي",
}

ENTITY_TYPE_LABELS_EN = {
    ENTITY_CHEMICAL_ELEMENT: "Chemical Element",
    ENTITY_CHEMICAL_COMPOUND: "Chemical Compound",
    ENTITY_DRUG_SUBSTANCE: "Drug Substance",
    ENTITY_BIOMOLECULE_PROTEIN: "Protein / Biomolecule",
}

CHEMICAL_SUBCLASS_LABELS_AR = {
    "element": "عنصر كيميائي",
    "acid": "حمض كيميائي",
    "base": "قاعدة كيميائية",
    "salt": "ملح كيميائي",
    "ion": "أيون كيميائي",
    "inorganic_compound": "مركب غير عضوي",
    "organic_compound": "مركب عضوي",
    "aromatic_compound": "مركب عطري",
    "alcohol": "كحول عضوي",
    "aldehyde": "ألدهيد عضوي",
    "ketone": "كيتون عضوي",
    "ester": "إستر عضوي",
    "carboxylic_acid": "حمض كربوكسيلي عضوي",
    "drug": "مادة دوائية صيدلانية",
    "biomolecule": "جزيء أو بروتين حيوي",
}

OP_SCIENTIFIC_RETRIEVAL = "scientific_retrieval"
OP_SCIENTIFIC_CALCULATION = "scientific_calculation"
OP_SCIENTIFIC_PREDICTION = "scientific_prediction"

OPERATION_LABELS_AR = {
    OP_SCIENTIFIC_RETRIEVAL: "تفاعل أو سجل علمي موثق مسترجع من قاعدة بيانات معتمدة",
    OP_SCIENTIFIC_CALCULATION: "نتيجة محسوبة رياضيا وكيميائيا وفق القوانين الحتمية",
    OP_SCIENTIFIC_PREDICTION: "توقع نظري محكوم بحدود وافتراضات علمية",
}

REACTION_DOMAINS = [
    ("all_domains", "جميع المجالات الكيميائية (All Domains)"),
    ("general_chemistry", "الكيمياء العامة وغير العضوية (General & Inorganic Chemistry)"),
    ("acid_base_chemistry", "تفاعلات الأحماض والقواعد والتعادل (Acids, Bases & Neutralization)"),
    ("redox_and_precipitation", "تفاعلات الأكسدة والاختزال والترسيب والإحلال (Redox, Precipitation & Displacement)"),
    ("organic_chemistry", "الكيمياء العضوية (Organic Chemistry - ORD)"),
    ("biochemistry", "الكيمياء الحيوية (Biochemistry - Rhea)"),
]

SUPPORTED_AMOUNT_UNITS = ["mol", "mmol", "g", "mg", "L_STP"]

CLOSED_DROPDOWN_OPTIONS = {
    "formula_types": [
        ("molecular", "الصيغة الجزيئية والعامة (Molecular Formula)"),
        ("condensed", "الصيغة البنائية المكثفة (Condensed Formula مثل CH3COOH)"),
        ("hill", "صيغة نظام هيل القياسي (Hill System Formula)"),
        ("smiles", "التمثيل الخطي الجزيئي (Canonical SMILES)"),
    ],
    "substance_categories": [
        ("all", "جميع فئات المواد العلمية (All Categories)"),
        ("element", "العناصر الكيميائية (Elements)"),
        ("acid", "الأحماض الكيميائية (Acids)"),
        ("base", "القواعد الكيميائية (Bases)"),
        ("salt", "الأملاح والأيونات (Salts & Ions)"),
        ("inorganic_compound", "المركبات غير العضوية (Inorganic Compounds)"),
        ("organic_compound", "المركبات العضوية والعطرية والكحولات والكيتونات والإسترات (Organic)"),
        ("drug", "الأدوية والمركبات الصيدلانية (Drugs)"),
        ("biomolecule", "البروتينات والجزيئات الحيوية (Biomolecules)"),
    ],
    "quantity_types": [
        ("molar_amount", "كمية المادة بالمولات (Molar Amount - n)"),
        ("mass_amount", "الكتلة الوزنية (Mass - m)"),
        ("gas_volume_stp", "حجم الغاز في الظروف القياسية (Gas Volume at STP - V)"),
    ],
    "amount_units": [
        ("mol", "مول (mol)"),
        ("mmol", "مللي مول (mmol)"),
        ("g", "جرام (g)"),
        ("mg", "مللي جرام (mg)"),
        ("L_STP", "لتر عند الظروف القياسية (L_STP)"),
    ],
    "reaction_types": [
        ("all", "تلقائي بحسب المواد المتفاعلة (Auto-detect Reaction Type)"),
        ("neutralization", "تفاعل تعادل حمض وقاعدة (Acid-Base Neutralization)"),
        ("combustion", "تفاعل احتراق (Combustion)"),
        ("precipitation", "تفاعل ترسيب وتبادل مزدوج (Precipitation)"),
        ("single_displacement", "تفاعل إحلال مفرد (Single Displacement)"),
        ("redox_synthesis", "تفاعل اتحاد وأكسدة واختزال (Synthesis & Redox)"),
        ("decomposition", "تفاعل تحلل وتفكك (Decomposition)"),
        ("organic_reaction", "تفاعل عضوي وأسترة (Organic Reaction)"),
        ("biochemical_reaction", "تفاعل كيميائي حيوي (Biochemical Reaction)"),
    ],
    "temperature_conditions": [
        ("25 °C (الظروف القياسية STP/SATP)", "25 °C — درجة حرارة الغرفة القياسية (298.15 K)"),
        ("37 °C (درجة الحرارة الفسيولوجية)", "37 °C — درجة الحرارة الفسيولوجية الحيوية (310.15 K)"),
        ("100 °C (تسخين مائي / ارتداد حراري)", "100 °C — تسخين مائي أو ارتداد حراري (Reflux)"),
        ("450 °C (تسخين صناعي محفز)", "450 °C — تسخين صناعي عالي (مثل عملية هابر-بوش)"),
        ("850 °C (تحلل حراري / كلسنة)", "850 °C — تحلل حراري عالي أو شرارة اشتعال"),
    ],
    "pressure_conditions": [
        ("1 atm (الضغط الجوي القياسي)", "1 atm — الضغط الجوي القياسي (101.325 kPa)"),
        ("5 atm (ضغط متوسط)", "5 atm — ضغط غازي متوسط"),
        ("200 atm (ضغط صناعي عالي)", "200 atm — ضغط صناعي مرتفع"),
    ],
    "solvent_conditions": [
        ("aqueous", "محلول مائي (Aqueous Solution - H2O)"),
        ("neat_gas_solid", "بدون مذيب / طور غازي أو صلب مباشر (Solvent-free / Gas-Solid)"),
        ("ethanol", "مذيب عضوي قطبي — إيثانول (Ethanol)"),
        ("physiological_buffer", "وسط منظم فسيولوجي (Physiological Buffer pH 7.4)"),
    ],
    "catalyst_conditions": [
        ("none", "بدون عامل حفاز إضافي (No Catalyst)"),
        ("acid_h2so4", "حفاز حمضي — حمض الكبريتيك المركز (H2SO4 Catalyst)"),
        ("pt_pd", "حفاز فلزي انتقالي — بلاتين أو بالاديوم (Pt / Pd)"),
        ("fe_haber", "حفاز الحديد المنشط (Promoted Iron Catalyst - Fe)"),
        ("mno2_catalase", "ثاني أكسيد المنغنيز أو إنزيم الكاتالاز (MnO2 / Catalase)"),
        ("enzymatic", "إنزيم حيوي متخصص (Enzymatic Catalysis - Rhea)"),
    ],
}

MOLAR_VOLUME_STP_LITERS = 22.414

PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
CHEBI_OLS_SEARCH_URL = "https://www.ebi.ac.uk/ols4/api/search"
CHEMBL_BASE_URL = "https://www.ebi.ac.uk/chembl/api/data"
UNIPROT_BASE_URL = "https://rest.uniprot.org/uniprotkb"
PDB_BASE_URL = "https://data.rcsb.org/rest/v1/core/entry"
RHEA_BASE_URL = "https://www.rhea-db.org/rhea"
ORD_BASE_URL = "https://open-reaction-database.org"

EXTERNAL_API_TIMEOUT_SECONDS = 5
CACHE_TTL_SECONDS = 3600

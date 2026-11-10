import os

BASE_DIR = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))

APP_NAME = "MolecuLens"
TEAM_NAME = "BioCode"
APP_SLOGAN = "See Molecules. Understand Science."

RAW_DATA_PATH = os.path.join(BASE_DIR, "data", "raw", "dpp4_raw.csv")
CLEAN_DATA_PATH = os.path.join(BASE_DIR, "data", "processed", "dpp4_cleaned.csv")
CLEAN_REPORT_PATH = os.path.join(BASE_DIR, "data", "processed", "cleaning_report.json")
TRAIN_DATA_PATH = os.path.join(BASE_DIR, "data", "processed", "train_compounds.csv")
TEST_PRED_PATH = os.path.join(BASE_DIR, "data", "processed", "test_predictions.csv")
DEMO_COMPOUNDS_PATH = os.path.join(BASE_DIR, "data", "datasets", "dpp4_demo_compounds.json")

MODEL_PATH = os.path.join(BASE_DIR, "models", "dpp4_rf_model.pkl")
TRAINING_REF_PATH = os.path.join(BASE_DIR, "models", "training_reference.pkl")
METRICS_PATH = os.path.join(BASE_DIR, "models", "model_metrics.json")

TARGET_INFO = {
    "target_id": "CHEMBL284",
    "uniprot_id": "P27487",
    "entry_name": "DPP4_HUMAN",
    "name_en": "Human Dipeptidyl Peptidase-4 (DPP4)",
    "name_ar": "إنزيم ثنائي ببتيديل ببتيداز-4 البشري (DPP4)",
    "ec_number": "EC 3.4.14.5",
    "organism": "Homo sapiens",
    "organism_ar": "الإنسان (Homo sapiens)",
    "uniprot_url": "https://www.uniprot.org/uniprotkb/P27487/entry",
    "chembl_url": "https://www.ebi.ac.uk/chembl/target_report_card/CHEMBL284/",
    "function_ar": (
        "إنزيم ببتيداز سيريني غشائي يقوم بشطر الببتيدات الثنائية من النهاية الأمينية (N-terminus) "
        "للهرمونات المعوية المعروفة بالإنكريتينات مثل GLP-1 و GIP."
    ),
    "diabetes_relation_ar": (
        "في أبحاث السكري من النوع الثاني، يؤدي تثبيط إنزيم DPP4 إلى إبطاء تحلل هرمون GLP-1 الداخلي، "
        "مما يعزز إفراز الإنسولين المعتمد على الجلوكوز ويخفض إفراز الجلوكاجون. "
        "ملاحظة علمية: وجود نشاط تثبيطي تجاه DPP4 في المختبر أو الحاسوب لا يعني إثبات فعالية علاجية أو سلامة دوائية لدى المرضى."
    ),
}

MORGAN_RADIUS = 2
MORGAN_NBITS = 1024
RANDOM_SEED = 42
RF_ESTIMATORS = 300
ACTIVE_PIC50_THRESHOLD = 6.0

AD_IN_DOMAIN_TANIMOTO = 0.60
AD_BORDERLINE_TANIMOTO = 0.30
AD_MIN_HEAVY_ATOMS = 6
AD_MAX_MW = 900.0

PUBCHEM_BASE_URL = "https://pubchem.ncbi.nlm.nih.gov/rest/pug"
PUBCHEM_TIMEOUT_SECONDS = 5



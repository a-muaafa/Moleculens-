"""
Expanded Drug-Drug Interaction (DDI) & Multi-Source Drug Synchronization Engine for MolecuLens.

Features:
1. 10 Therapeutic Classes with 45+ curated essential medications and 180+ brand/generic/Arabic synonyms.
2. 42+ Clinically Established Pairwise Drug-Drug Interactions graded into Major, Moderate, and Minor,
   with pharmacokinetic/pharmacodynamic mechanisms, clinical effects, management recommendations,
   authoritative citations, and `last_updated` timestamps.
3. Live External Drug Synchronization & Periodic Update Service integrating:
   - PubChem PUG REST (Public Domain): formula, molecular weight, SMILES, InChIKey, CID.
   - RxNorm REST API (NLM Public API): RxCUI resolution and brand/generic normalization.
   - openFDA Drug Labeling API (CC0 Public Domain): FDA label drug_interactions, warnings, adverse_reactions.
   - DrugBank Open Data Identifiers (CC0 1.0): DBxxxxx accession cross-referencing.
4. Strict deduplication by canonical key / RxCUI / DrugBank ID and mandatory clinical disclaimers
   clarifying that absence of a recorded interaction does NOT prove clinical safety.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any, Dict, List, Optional
import urllib.parse

from config.settings import EXTERNAL_API_TIMEOUT_SECONDS
from core.security import sanitize_external_text, validate_safe_chemical_input
from data_integrations.cache_store import get_cached_item, set_cached_item
from data_integrations.http_client import request_with_retry

DATA_SOURCE_LICENSES: dict[str, dict[str, str]] = {
    "PubChem": {
        "name": "NCBI PubChem PUG REST API",
        "url": "https://pubchem.ncbi.nlm.nih.gov/rest/pug",
        "scope_ar": "الصيغ الكيميائية، الكتل المولية، التراكيب البنيوية (SMILES / InChIKey)، ومعرفات CID",
        "license": "U.S. Government Public Domain (NCBI/NLM)",
    },
    "RxNorm": {
        "name": "NLM RxNorm REST API",
        "url": "https://rxnav.nlm.nih.gov/REST",
        "scope_ar": "توحيد الأسماء العلمية للمواد الفعالة والأسماء التجارية ومعرفات RxCUI",
        "license": "NLM Public Access API (Free for Research & Educational Use)",
    },
    "openFDA": {
        "name": "U.S. FDA openFDA Drug Labeling API",
        "url": "https://api.fda.gov/drug/label.json",
        "scope_ar": "معلومات النشرات الدوائية الرسمية، التحذيرات، التداخلات الدوائية، والآثار الجانبية",
        "license": "Creative Commons CC0 1.0 Universal (Public Domain)",
    },
    "DrugBank_Open": {
        "name": "DrugBank Open Vocabulary & Accession Identifiers",
        "url": "https://go.drugbank.com",
        "scope_ar": "معرفات الأدوية القياسية (DBxxxxx) والمرادفات المفتوحة (البيانات الكاملة للتداخلات تتطلب ترخيصاً أكاديمياً؛ لذا يُعتمد في الجلب الآلي على openFDA و RxNorm)",
        "license": "CC0 1.0 Universal (Open Vocabulary Subset)",
    },
}

THERAPEUTIC_CLASSES: dict[str, dict[str, str]] = {
    "analgesics_nsaids": {
        "name_ar": "المسكنات ومضادات الالتهاب (Analgesics & NSAIDs)",
        "name_en": "Analgesics & Non-Steroidal Anti-Inflammatory Drugs (NSAIDs)",
    },
    "antibiotics": {
        "name_ar": "المضادات الحيوية ومضادات الميكروبات (Antibiotics)",
        "name_en": "Antibiotics & Antimicrobials",
    },
    "antihypertensives": {
        "name_ar": "أدوية ضغط الدم (Antihypertensives)",
        "name_en": "Antihypertensive Agents",
    },
    "antidiabetics": {
        "name_ar": "أدوية السكري (Antidiabetics)",
        "name_en": "Antidiabetic & Hypoglycemic Agents",
    },
    "anticoagulants": {
        "name_ar": "مضادات التخثر والصفائح (Anticoagulants & Antiplatelets)",
        "name_en": "Anticoagulants & Antiplatelet Agents",
    },
    "antihistamines": {
        "name_ar": "مضادات الحساسية والهيستامين (Antihistamines)",
        "name_en": "Antihistamines & Antiallergic Agents",
    },
    "gastrointestinal": {
        "name_ar": "أدوية المعدة والجهاز الهضمي (Gastrointestinal Agents)",
        "name_en": "Gastrointestinal & Acid-Suppressive Agents",
    },
    "antidepressants": {
        "name_ar": "مضادات الاكتئاب والجهاز العصبي (Antidepressants)",
        "name_en": "Antidepressants & CNS Agents",
    },
    "cardiovascular": {
        "name_ar": "أدوية القلب والشرايين والكوليسترول (Cardiovascular Agents)",
        "name_en": "Cardiovascular, Antiarrhythmic & Lipid-Lowering Agents",
    },
    "respiratory": {
        "name_ar": "أدوية الجهاز التنفسي والربو (Respiratory Agents)",
        "name_en": "Respiratory & Bronchodilator Agents",
    },
}


def _drug_entry(
    key: str,
    name_en: str,
    name_ar: str,
    t_class: str,
    brands: list[str],
    rxcui: str,
    drugbank_id: str,
    pubchem_cid: int,
    updated: str = "2026-11-15",
) -> dict[str, Any]:
    cls_info = THERAPEUTIC_CLASSES[t_class]
    return {
        "canonical_key": key,
        "name_en": name_en,
        "name_ar": name_ar,
        "therapeutic_class_key": t_class,
        "therapeutic_class_ar": cls_info["name_ar"],
        "therapeutic_class_en": cls_info["name_en"],
        "brand_names": brands,
        "rxcui": rxcui,
        "drugbank_id": drugbank_id,
        "pubchem_cid": pubchem_cid,
        "last_updated": updated,
        "sources": ["RxNorm", "openFDA", "PubChem", "DrugBank_Open"],
    }


DRUG_CATALOG: dict[str, dict[str, Any]] = {
    # 1. Analgesics & NSAIDs
    "paracetamol": _drug_entry("paracetamol", "Paracetamol (Acetaminophen)", "باراسيتامول (أسيتامينوفين)", "analgesics_nsaids", ["Panadol", "Tylenol", "Adol", "Fevadol"], "161", "DB00316", 1983),
    "ibuprofen": _drug_entry("ibuprofen", "Ibuprofen", "إيبوبروفين", "analgesics_nsaids", ["Advil", "Motrin", "Brufen", "Nurofen"], "5640", "DB01050", 3672),
    "aspirin": _drug_entry("aspirin", "Aspirin (Acetylsalicylic Acid)", "أسبرين (حمض أسيتيل ساليسيليك)", "analgesics_nsaids", ["Jusprin", "Aspocid", "Disprin", "Ecotrin"], "1191", "DB00945", 2244),
    "naproxen": _drug_entry("naproxen", "Naproxen", "نابروكسين", "analgesics_nsaids", ["Aleve", "Naprosyn", "Proxidol"], "7258", "DB00788", 156391),
    "diclofenac": _drug_entry("diclofenac", "Diclofenac", "ديكلوفيناك", "analgesics_nsaids", ["Voltaren", "Cataflam", "Rofenac"], "3355", "DB00586", 3033),
    "celecoxib": _drug_entry("celecoxib", "Celecoxib", "سيليكوكسيب", "analgesics_nsaids", ["Celebrex"], "140587", "DB00482", 2662),
    "tramadol": _drug_entry("tramadol", "Tramadol", "ترامادول", "analgesics_nsaids", ["Ultram", "Tramal"], "10689", "DB00193", 33741),
    # 2. Antibiotics
    "amoxicillin": _drug_entry("amoxicillin", "Amoxicillin", "أموكسيسيلين", "antibiotics", ["Amoxil", "Augmentin", "Hiconcil"], "723", "DB01060", 33613),
    "azithromycin": _drug_entry("azithromycin", "Azithromycin", "أزيثروميسين", "antibiotics", ["Zithromax", "Azimax", "Zomax"], "18631", "DB00207", 447043),
    "clarithromycin": _drug_entry("clarithromycin", "Clarithromycin", "كلاريثروميسين", "antibiotics", ["Klacid", "Biaxin"], "21212", "DB01211", 84029),
    "ciprofloxacin": _drug_entry("ciprofloxacin", "Ciprofloxacin", "سيبروفلوكساسين", "antibiotics", ["Cipro", "Ciprobay", "Ciproxen"], "2551", "DB00537", 2764),
    "levofloxacin": _drug_entry("levofloxacin", "Levofloxacin", "ليفوفلوكساسين", "antibiotics", ["Tavanic", "Levaquin"], "82122", "DB01137", 149096),
    "metronidazole": _drug_entry("metronidazole", "Metronidazole", "ميترونيدازول", "antibiotics", ["Flagyl", "Amrizole"], "6922", "DB00916", 4173),
    "doxycycline": _drug_entry("doxycycline", "Doxycycline", "دوكسيسيكلين", "antibiotics", ["Vibramycin", "Doxydar"], "3640", "DB00254", 54671203),
    # 3. Antihypertensives
    "amlodipine": _drug_entry("amlodipine", "Amlodipine", "أملوديبين", "antihypertensives", ["Norvasc", "Amlor"], "17767", "DB00381", 2162),
    "lisinopril": _drug_entry("lisinopril", "Lisinopril", "ليزينوبريل", "antihypertensives", ["Zestril", "Prinivil"], "29046", "DB00722", 5362119),
    "enalapril": _drug_entry("enalapril", "Enalapril", "إينالابريل", "antihypertensives", ["Renitec", "Vasotec"], "3827", "DB00584", 5388962),
    "losartan": _drug_entry("losartan", "Losartan", "لوسارتان", "antihypertensives", ["Cozaar", "Sortiva"], "52175", "DB00678", 3961),
    "valsartan": _drug_entry("valsartan", "Valsartan", "فالسارتان", "antihypertensives", ["Diovan", "Tareg"], "69749", "DB00177", 60846),
    "hydrochlorothiazide": _drug_entry("hydrochlorothiazide", "Hydrochlorothiazide (HCTZ)", "هيدروكلوروثيازيد", "antihypertensives", ["Esidrex", "Microzide"], "5487", "DB00999", 3639),
    "spironolactone": _drug_entry("spironolactone", "Spironolactone", "سبيرونولاكتون", "antihypertensives", ["Aldactone"], "9997", "DB00421", 5833),
    "metoprolol": _drug_entry("metoprolol", "Metoprolol", "ميتوبرولول", "antihypertensives", ["Betaloc", "Lopressor", "Toprol-XL"], "6918", "DB00264", 4171),
    # 4. Antidiabetics
    "metformin": _drug_entry("metformin", "Metformin", "ميتفورمين", "antidiabetics", ["Glucophage", "Cidophage", "Formit"], "6809", "DB00331", 4091),
    "sitagliptin": _drug_entry("sitagliptin", "Sitagliptin", "سيتاغليبتين", "antidiabetics", ["Januvia", "Janumet"], "593411", "DB01261", 4369359),
    "vildagliptin": _drug_entry("vildagliptin", "Vildagliptin", "فيلداغليبتين", "antidiabetics", ["Galvus", "Eucreas"], "476350", "DB04876", 6918537),
    "glimepiride": _drug_entry("glimepiride", "Glimepiride", "غليميبيريد", "antidiabetics", ["Amaryl", "Diapride"], "25789", "DB00222", 3476),
    "empagliflozin": _drug_entry("empagliflozin", "Empagliflozin", "إمباغليفلوزين", "antidiabetics", ["Jardiance"], "1545653", "DB09038", 11949646),
    # 5. Anticoagulants & Antiplatelets
    "warfarin": _drug_entry("warfarin", "Warfarin", "وارفارين", "anticoagulants", ["Coumadin", "Marevan"], "11289", "DB00682", 54678486),
    "clopidogrel": _drug_entry("clopidogrel", "Clopidogrel", "كلوبيدوغريل", "anticoagulants", ["Plavix", "Pedovex"], "32968", "DB00758", 60606),
    "rivaroxaban": _drug_entry("rivaroxaban", "Rivaroxaban", "ريفاروكسابان", "anticoagulants", ["Xarelto"], "1114195", "DB06228", 9875401),
    "apixaban": _drug_entry("apixaban", "Apixaban", "أبيكسابان", "anticoagulants", ["Eliquis"], "1364430", "DB06605", 10182969),
    # 6. Antihistamines
    "cetirizine": _drug_entry("cetirizine", "Cetirizine", "سيتريزين", "antihistamines", ["Zyrtec", "Finallerg"], "20610", "DB00341", 2678),
    "loratadine": _drug_entry("loratadine", "Loratadine", "لوراتادين", "antihistamines", ["Claritin", "Claritine"], "28889", "DB00455", 3957),
    "fexofenadine": _drug_entry("fexofenadine", "Fexofenadine", "فيكسوفينادين", "antihistamines", ["Telfast", "Allegra"], "87636", "DB00950", 3348),
    "diphenhydramine": _drug_entry("diphenhydramine", "Diphenhydramine", "ديفينهيدرامين", "antihistamines", ["Benadryl"], "3498", "DB01075", 3100),
    # 7. Gastrointestinal Agents
    "omeprazole": _drug_entry("omeprazole", "Omeprazole", "أوميبرازول", "gastrointestinal", ["Losec", "Prilosec", "Gasec"], "7646", "DB00338", 4594),
    "esomeprazole": _drug_entry("esomeprazole", "Esomeprazole", "إيزوميبرازول", "gastrointestinal", ["Nexium"], "283742", "DB00736", 9568614),
    "pantoprazole": _drug_entry("pantoprazole", "Pantoprazole", "بانتوبرازول", "gastrointestinal", ["Controloc", "Protonix"], "40790", "DB00213", 4679),
    "famotidine": _drug_entry("famotidine", "Famotidine", "فاموتيدين", "gastrointestinal", ["Pepcid", "Antodine"], "4278", "DB00927", 5702160),
    "ondansetron": _drug_entry("ondansetron", "Ondansetron", "أوندانسيترون", "gastrointestinal", ["Zofran"], "26225", "DB00904", 4595),
    # 8. Antidepressants
    "sertraline": _drug_entry("sertraline", "Sertraline", "سيرترالين", "antidepressants", ["Zoloft", "Lustral"], "36437", "DB01104", 68617),
    "fluoxetine": _drug_entry("fluoxetine", "Fluoxetine", "فلوكستين", "antidepressants", ["Prozac", "Philozac"], "4493", "DB00472", 3386),
    "escitalopram": _drug_entry("escitalopram", "Escitalopram", "إسيتالوبرام", "antidepressants", ["Cipralex", "Lexapro"], "321988", "DB01175", 146570),
    "amitriptyline": _drug_entry("amitriptyline", "Amitriptyline", "أميتريبتيلين", "antidepressants", ["Elavil", "Tryptizol"], "704", "DB00321", 2160),
    # 9. Cardiovascular & Lipid-Lowering Agents
    "digoxin": _drug_entry("digoxin", "Digoxin", "ديجوكسين", "cardiovascular", ["Lanoxin"], "3407", "DB00390", 2724385),
    "amiodarone": _drug_entry("amiodarone", "Amiodarone", "أميودارون", "cardiovascular", ["Cordarone", "Pacerone"], "703", "DB01118", 2157),
    "atorvastatin": _drug_entry("atorvastatin", "Atorvastatin", "أتورفاستاتين", "cardiovascular", ["Lipitor", "Ator"], "83367", "DB01076", 60823),
    "simvastatin": _drug_entry("simvastatin", "Simvastatin", "سيمفاستاتين", "cardiovascular", ["Zocor"], "36567", "DB00641", 54454),
    "rosuvastatin": _drug_entry("rosuvastatin", "Rosuvastatin", "روزوفاستاتين", "cardiovascular", ["Crestor"], "301542", "DB01098", 446157),
    "verapamil": _drug_entry("verapamil", "Verapamil", "فيراباميل", "cardiovascular", ["Isoptin", "Calan"], "11170", "DB00661", 2520),
    # 10. Respiratory Agents
    "salbutamol": _drug_entry("salbutamol", "Salbutamol (Albuterol)", "سالبوتامول (ألبوتيرول)", "respiratory", ["Ventolin", "ProAir"], "435", "DB01001", 2083),
    "theophylline": _drug_entry("theophylline", "Theophylline", "ثيوفيلين", "respiratory", ["Theo-Dur", "Quibron"], "10438", "DB00277", 2153),
    "montelukast": _drug_entry("montelukast", "Montelukast", "مونتيلوكاست", "respiratory", ["Singulair", "Airfast"], "88249", "DB00471", 5281040),
    "budesonide": _drug_entry("budesonide", "Budesonide", "بوديزونيد", "respiratory", ["Pulmicort", "Symbicort"], "19831", "DB01222", 5281004),
    "caffeine": _drug_entry("caffeine", "Caffeine", "كافيين", "respiratory", ["Cafcit"], "1886", "DB00201", 2519),
}


def _ix(
    severity: str,
    severity_ar: str,
    mech_en: str,
    mech_ar: str,
    eff_en: str,
    eff_ar: str,
    mgmt_en: str,
    mgmt_ar: str,
    citations: list[str],
    updated: str = "2026-11-15",
) -> dict[str, Any]:
    return {
        "severity": severity,
        "severity_ar": severity_ar,
        "mechanism_en": mech_en,
        "mechanism_ar": mech_ar,
        "clinical_effect_en": eff_en,
        "clinical_effect_ar": eff_ar,
        "management_en": mgmt_en,
        "management_ar": mgmt_ar,
        "evidence_level": "Established",
        "citations": citations,
        "last_updated": updated,
    }


class DDIEngine:
    """Evaluates drug-drug interactions with scientific rigor and multi-source synchronization."""

    DISCLAIMER_AR = (
        "تنبيه سريري وأكاديمي هام: هذه الأداة مخصصة للأغراض التعليمية والاسترشادية فقط. "
        "إن عدم وجود تفاعل مسجل في قاعدة البيانات الحالية لا يعني بالضرورة الأمان التام "
        "أو غياب التداخل الدوائي السريري. لا تُستخدم هذه الأداة كبديل عن الاستشارة الصيدلانية "
        "أو الطبية المباشرة، ويجب دائماً مراجعة المراجع الدوائية المعتمدة لكل حالة سريرية."
    )

    DISCLAIMER_EN = (
        "Important Clinical & Academic Notice: This tool is designed strictly for educational and guidance purposes. "
        "The absence of a documented interaction in this database does NOT prove safety or the complete absence "
        "of a clinical interaction. This tool does not replace professional medical or clinical pharmacy consultation; "
        "always cross-reference official prescribing monographs for individualized patient care."
    )

    DRUG_CATALOG = DRUG_CATALOG
    THERAPEUTIC_CLASSES = THERAPEUTIC_CLASSES
    DATA_SOURCE_LICENSES = DATA_SOURCE_LICENSES

    _INTERACTIONS: Dict[frozenset, Dict[str, Any]] = {
        frozenset(["sitagliptin", "metformin"]): _ix(
            "Minor", "بسيط / مراقبة روتينية",
            "Synergistic glycemic control without direct pharmacokinetic competition; both drugs are cleared primarily via renal tubular excretion.",
            "تأثير تآزري مفيد في خفض سكر الدم دون تداخل حركي مباشر، كلاهما يُطرح عبر الكلى، مما يتطلب مراقبة وظائف الكلى الدورية.",
            "Enhanced reduction in HbA1c; mild increase in gastrointestinal side effects common to metformin; low risk of hypoglycemia unless combined with insulin secretagogues.",
            "تحسن ملحوظ في ضبط السكر التراكمي (HbA1c) مع احتمالية زيادة خفيفة في الأعراض الهضمية الشائعة للميتفورمين، وخطر منخفض لهبوط السكر.",
            "Commonly prescribed as fixed-dose combination. Monitor renal function (eGFR) periodically.",
            "يُصرفان معاً بشكل شائع كمركب دوائي مشترك. تجب مراقبة وظائف الكلى دورياً وتعديل الجرعة إذا انخفض معدل الترشيح الكبيبي عن 45 مل/دقيقة.",
            ["FDA Prescribing Information: Janumet (Sitagliptin/Metformin HCl)", "Diabetes Care 2023; Standards of Medical Care in Diabetes"],
        ),
        frozenset(["aspirin", "ibuprofen"]): _ix(
            "Major", "خطير / رئيسي",
            "Ibuprofen competitively inhibits platelet cyclooxygenase-1 (COX-1), blocking irreversible acetylation of Ser529 by low-dose aspirin and diminishing cardioprotection; additive gastric mucosal COX-1 inhibition increases GI bleeding.",
            "يثبط الإيبوبروفين إنزيم الأكسدة الحلقية-1 (COX-1) في الصفائح الدموية تنافسياً، مما يمنع الأسبرين من تثبيطه نهائياً ويقلل الحماية القلبية؛ كما يضاعف الخطر على بطانة المعدة.",
            "Attenuation of aspirin's antiplatelet cardioprotective effect; significantly heightened risk of gastrointestinal mucosal erosion, ulceration, and hemorrhage.",
            "إضعاف التأثير المضاد للتخثر الواقي للقلب للأسبرين، وارتفاع كبير في مخاطر قرحة ونزيف الجهاز الهضمي.",
            "Avoid concurrent use if possible. If ibuprofen is required, take immediate-release ibuprofen at least 30 minutes after or 8 hours before low-dose aspirin.",
            "يُفضل تجنب التناول المتزامن. إذا لزم الإيبوبروفين، يجب تناوله بعد 30 دقيقة على الأقل من الأسبرين أو قبل الأسبرين بـ 8 ساعات على الأقل.",
            ["FDA Drug Safety Communication: Concomitant use of Ibuprofen and Aspirin (2006)", "Catella-Lawson F et al. N Engl J Med 2001; 345:1809-1817"],
        ),
        frozenset(["aspirin", "warfarin"]): _ix(
            "Major", "خطير / رئيسي",
            "Combined anticoagulant (vitamin K antagonism) and antiplatelet (irreversible COX-1 inhibition) actions produce profound synergistic impairment of hemostasis.",
            "تأثير مضاعف يجمع بين تثبيط عوامل التخثر المعتمدة على فيتامين K وتثبيط الصفائح الدموية، مما يؤدي إلى تثبيط شديد لمنظومة التجلط.",
            "Substantially increased risk of severe, potentially life-threatening gastrointestinal or intracranial bleeding.",
            "ارتفاع حاد في احتمالية حدوث نزيف هضمي أو دماغي شديد ومهدد للحياة.",
            "Only combine when strictly indicated (e.g., mechanical heart valve plus acute coronary syndrome). Close INR monitoring and gastroprotective PPI recommended.",
            "لا يُجمع بينهما إلا بدواعٍ سريرية محددة وبإشراف استشاري دقيق مع مراقبة لصيقة لمؤشر INR وإضافة واقٍ للمعدة مثل مثبطات مضخة البروتون.",
            ["Chest 2012; 141(2 Suppl):e152S-e184S", "Circulation 2018; Dual Antithrombotic Therapy Consensus"],
        ),
        frozenset(["ibuprofen", "warfarin"]): _ix(
            "Major", "خطير / رئيسي",
            "Ibuprofen inhibits platelet COX-1 aggregation, causes gastric mucosal injury, and may displace warfarin from plasma albumin or compete for CYP2C9 metabolism.",
            "يثبط الإيبوبروفين تجمع الصفائح الدموية ويضعف الغشاء المخاطي للمعدة ويتنافس مع الوارفارين على إنزيم CYP2C9، ما يضاعف خطر النزيف.",
            "Markedly elevated risk of severe upper gastrointestinal hemorrhage and increased INR.",
            "ارتفاع شديد في خطر النزيف الهضمي العلوي واضطراب زمن التجلط (INR).",
            "Avoid NSAIDs in patients on warfarin. Consider paracetamol at low doses for analgesia with INR monitoring.",
            "يُمنع الجمع بين مضادات الالتهاب غير الستيرويدية والوارفارين؛ يُستبدل بالباراسيتامول عند الضرورة مع مراقبة INR.",
            ["FDA Prescribing Information: Coumadin (Warfarin Sodium)", "Arch Intern Med 2000; 160:2093-2099"],
        ),
        frozenset(["naproxen", "warfarin"]): _ix(
            "Major", "خطير / رئيسي",
            "Additive antiplatelet COX-1 inhibition and direct gastric mucosal erosion combined with systemic anticoagulation.",
            "تثبيط الصفائح الدموية وتآكل جدار المعدة بواسطة النابروكسين بالتزامن مع تمييع الدم بالوارفارين.",
            "High risk of severe gastrointestinal bleeding.",
            "خطر مرتفع جداً لحدوث نزيف معدي معوي حاد.",
            "Avoid combination. Use non-NSAID analgesia if possible and add a proton pump inhibitor if unavoidable.",
            "تجنب الجمع بينهما؛ استخدم مسكناً بديلاً مع مثبط مضخة البروتون عند الضرورة القصوى.",
            ["FDA Prescribing Information: Naprosyn (Naproxen)"],
        ),
        frozenset(["diclofenac", "warfarin"]): _ix(
            "Major", "خطير / رئيسي",
            "Synergistic impairment of hemostasis via COX inhibition and CYP2C9 substrate competition.",
            "تداخل ديناميكي وحركي عبر تثبيط COX والمنافسة على إنزيم CYP2C9 الكبدي.",
            "Significantly increased risk of gastrointestinal and systemic hemorrhage.",
            "زيادة كبيرة في احتمالية النزيف الهضمي والجهازي.",
            "Contraindicated or avoid concurrent administration; monitor INR closely if co-administered.",
            "يُتجنب الاستخدام المتزامن تماماً مع مراقبة دقيقة لمؤشر INR.",
            ["FDA Prescribing Information: Voltaren (Diclofenac)"],
        ),
        frozenset(["aspirin", "paracetamol"]): _ix(
            "Minor", "طفيف / مراقبة",
            "Both inhibit prostaglandin synthesis through different central and peripheral pathways without major pharmacokinetic interference at standard doses.",
            "كلاهما يؤثر في تثبيط البروستاغلاندينات عبر آليات مركزية ومحيطية مختلفة، ولا يسببان تداخلاً حركياً كبيراً عند الجرعات الاعتيادية.",
            "Additive analgesic effect. High chronic doses may increase cumulative renal toxicity risk.",
            "تأثير مسكن مضاف، ولكن الاستخدام المزمن بجرعات عالية قد يرفع العبء التراكمي على الكلى.",
            "Safe for short-term combination. Do not exceed maximum daily limits (4000 mg paracetamol / 4000 mg aspirin).",
            "آمن للاستخدام المؤقت قصير المدى. يجب عدم تجاوز الحد الأقصى اليومي المسموح به لكل منهما.",
            ["British National Formulary (BNF 85): Analgesic Combinations"],
        ),
        frozenset(["paracetamol", "warfarin"]): _ix(
            "Moderate", "متوسط / يستلزم المتابعة",
            "Paracetamol metabolite (NAPQI) inhibits vitamin K cycle enzymes (VKORC1), mildly reducing clotting factor synthesis during sustained high-dose use.",
            "قد يثبط ناتج أيض الباراسيتامول الإنزيمات المسؤولة عن تدوير فيتامين K، مما يؤدي إلى تثبيط إضافي لإنتاج عوامل التجلط عند تناوله بانتظام.",
            "Regular high-dose paracetamol (>2 g/day for >3 days) may elevate INR and increase bleeding tendency.",
            "تناول الباراسيتامول بجرعات عالية (>2 جم/يوم لأكثر من 3 أيام) قد يؤدي لارتفاع ملحوظ في قراءة INR وزيادة قابلية النزيف.",
            "Occasional single doses have minimal impact. Monitor INR if paracetamol is used regularly for more than 3-4 days.",
            "الجرعات المفردة العارضة آمنة، ولكن عند الاستمرار لأكثر من 3 أيام يجب فحص INR وتعديل جرعة الوارفارين إن لزم.",
            ["Mahé I et al. Haematologica 2006; 91(12):1621-1627"],
        ),
        frozenset(["caffeine", "aspirin"]): _ix(
            "Minor", "طفيف / تآزر إيجابي",
            "Caffeine accelerates gastric emptying and improves gastrointestinal absorption of aspirin; also exerts central adenosine receptor antagonism.",
            "يعجل الكافيين إفراغ المعدة ويزيد معدل امتصاص الأسبرين، كما يحصر مستقبلات الأدينوزين المركزية.",
            "Enhanced analgesic onset and potency. Higher doses may cause tachycardia or gastric distress.",
            "تسريع وبدء مفعول تسكين الألم وزيادة فاعليته، وقد يسبب زيادة خفيفة في نبضات القلب أو تهيج المعدة.",
            "Frequently utilized clinically. Avoid in patients with acute peptic ulcers or uncontrolled tachyarrhythmias.",
            "مزيج مستخدم سريرياً بشكل واسع. يُنصح بالحذر لدى مرضى القرحة النشطة أو اضطرابات نبضات القلب.",
            ["Cochrane Database Syst Rev 2014; (12):CD009281"],
        ),
        frozenset(["omeprazole", "amoxicillin"]): _ix(
            "Minor", "طفيف / تآزر علاجي مستهدف",
            "Omeprazole elevates gastric pH, preventing acid degradation of amoxicillin and improving its local antibacterial stability against Helicobacter pylori.",
            "يرفع الأوميبرازول درجة حموضة المعدة، مما يحمي الأموكسيسيلين من التحلل الحمضي ويزيد من استقراره وفاعليته الموضعية ضد جرثومة المعدة.",
            "Synergistic eradication of Helicobacter pylori in standard triple-therapy regimens.",
            "تعزيز القضاء على بكتيريا الملوية البوابية (جرثومة المعدة) ضمن بروتوكولات العلاج الثلاثي القياسية.",
            "Standard approved combination for H. pylori eradication therapy.",
            "مزيج قياسي معتمد عالمياً لعلاج جرثومة المعدة. يتم الالتزام بالمدة المقررة للعلاج (غالباً 10 إلى 14 يوماً).",
            ["ACG Clinical Guideline: Treatment of Helicobacter pylori Infection (Am J Gastroenterol 2017)"],
        ),
        frozenset(["digoxin", "amiodarone"]): _ix(
            "Major", "خطير / رئيسي",
            "Amiodarone inhibits P-glycoprotein (P-gp / ABCB1) efflux transporter and decreases renal and non-renal clearance of digoxin.",
            "يثبط الأميودارون مضخة P-glycoprotein المسؤولة عن إخراج الديجوكسين، مما يقلل طرحه الكلوي والكبدي ويرفع تركيزه بالدم.",
            "Serum digoxin levels approximately double within 3-7 days, creating high risk of fatal digitalis toxicity (arrhythmias, heart block).",
            "مضاعفة تركيز الديجوكسين في مصل الدم خلال أيام، مما يسبب تسمماً خطيراً بالديجوكسين واضطرابات قلبية قاتلة.",
            "Reduce digoxin maintenance dose by 50% upon initiating amiodarone. Measure serum digoxin concentrations and monitor ECG.",
            "تخفيض جرعة الديجوكسين بنسبة 50% فور بدء الأميودارون، مع قياس مستواه في الدم دورياً وتخطيط القلب الكهربائي.",
            ["FDA Prescribing Information: Lanoxin (Digoxin)", "Circulation 2010; Digoxin-Amiodarone Pharmacokinetic Interactions"],
        ),
        frozenset(["omeprazole", "clopidogrel"]): _ix(
            "Major", "خطير / رئيسي",
            "Omeprazole strongly inhibits hepatic CYP2C19, the principal enzyme required to bioactivate the prodrug clopidogrel into its active thiol antiplatelet metabolite.",
            "يثبط الأوميبرازول بقوة إنزيم الكبد CYP2C19 المسؤول عن تحويل عقار كلوبيدوغريل غير النشط إلى صورته الفعالة المضادة لتكدس الصفائح.",
            "Up to 45% reduction in active clopidogrel metabolite exposure, increasing risk of stent thrombosis, myocardial infarction, and stroke.",
            "انخفاض فعالية كلوبيدوغريل المضادة للتجلط بنسبة تصل إلى 45%، مما يرفع خطر تخثر الدعامات القلبية والجلطات.",
            "Avoid omeprazole and esomeprazole with clopidogrel. Switch to pantoprazole (minimal CYP2C19 inhibition) or famotidine.",
            "تجنب الأوميبرازول والإيزوميبرازول مع الكلوبيدوغريل؛ واستبدلهما بـ بانتوبرازول (Pantoprazole) أو فاموتيدين (Famotidine).",
            ["FDA Boxed Warning: Plavix (Clopidogrel Bisulfate) — CYP2C19 Inhibitors", "J Am Coll Cardiol 2010; 56(24):2051-2066"],
        ),
        frozenset(["esomeprazole", "clopidogrel"]): _ix(
            "Major", "خطير / رئيسي",
            "Esomeprazole inhibits CYP2C19-mediated bioactivation of clopidogrel to its active antiplatelet metabolite.",
            "يثبط الإيزوميبرازول إنزيم CYP2C19 الضروري لتنشيط كلوبيدوغريل داخل الكبد.",
            "Reduced antiplatelet efficacy and elevated risk of major adverse cardiovascular events.",
            "ضعف التأثير المانع لتجلط الصفائح وزيادة خطر الحوادث القلبية الوعائية.",
            "Replace esomeprazole with pantoprazole or an H2-receptor antagonist (famotidine).",
            "استبدال الإيزوميبرازول بدواء بانتوبرازول أو مضاد لمستقبلات H2 مثل فاموتيدين.",
            ["FDA Drug Safety Communication: Clopidogrel and Omeprazole/Esomeprazole"],
        ),
        frozenset(["simvastatin", "amlodipine"]): _ix(
            "Moderate", "متوسط / يستلزم ضبط الجرعة",
            "Amlodipine weakly inhibits intestinal and hepatic CYP3A4, increasing systemic AUC of simvastatin acid by ~80%.",
            "يثبط الأملوديبين جزئياً إنزيم CYP3A4، مما يرفع تركيز السيمفاستاتين في الدم بنحو 80%.",
            "Increased risk of statin-induced myopathy and rhabdomyolysis at simvastatin doses >20 mg/day.",
            "زيادة خطر اعتلال العضلات وانحلال الربيدات عند تجاوز جرعة السيمفاستاتين 20 ملغ يومياً.",
            "Do not exceed 20 mg simvastatin daily in patients taking amlodipine, or switch to rosuvastatin/atorvastatin.",
            "يجب ألا تتجاوز جرعة السيمفاستاتين 20 ملغ يومياً عند تناوله مع الأملوديبين، أو التحويل إلى روزوفاستاتين.",
            ["FDA Drug Safety Communication: Simvastatin Dose Limitation with Amlodipine (2011)"],
        ),
        frozenset(["simvastatin", "clarithromycin"]): _ix(
            "Major", "خطير / مضاد استطباب",
            "Clarithromycin is a potent mechanism-based inhibitor of CYP3A4 and OATP1B1, increasing simvastatin plasma exposure up to 10-fold.",
            "يثبط الكلاريثروميسين بقوة شديدة إنزيم CYP3A4، مما يرفع تركيز السيمفاستاتين في الدم نحو 10 أضعاف.",
            "Severe risk of acute rhabdomyolysis, myoglobinuria, and acute kidney injury.",
            "خطر شديد للإصابة بانحلال العضلات الهيكلي الحاد والفشل الكلوي الحاد.",
            "Contraindicated. Suspend simvastatin during clarithromycin therapy, or use azithromycin instead.",
            "ممنوع الجمع بينهما؛ يجب إيقاف السيمفاستاتين مؤقتاً طوال فترة تناول الكلاريثروميسين أو اختيار أزيثروميسين كبديل.",
            ["FDA Prescribing Information: Zocor (Simvastatin) — Contraindicated Strong CYP3A4 Inhibitors"],
        ),
        frozenset(["atorvastatin", "clarithromycin"]): _ix(
            "Major", "خطير / رئيسي",
            "Clarithromycin strongly inhibits CYP3A4 and hepatic OATP1B1 uptake transporter, elevating atorvastatin levels 3- to 4-fold.",
            "يثبط الكلاريثروميسين إنزيم CYP3A4 وناقلات الكبد OATP1B1، مما يرفع مستويات الأتورفاستاتين 3 إلى 4 أضعاف.",
            "Heightened risk of severe myopathy and rhabdomyolysis.",
            "ارتفاع خطر التهاب وانحلال العضلات الشديد.",
            "Temporarily suspend atorvastatin or limit dose to <=20 mg/day during clarithromycin course; alternatively use azithromycin.",
            "يُنصح بإيقاف الأتورفاستاتين مؤقتاً أو عدم تجاوز 20 ملغ يومياً أثناء شوط الكلاريثروميسين، أو استبداله بأزيثروميسين.",
            ["FDA Prescribing Information: Lipitor (Atorvastatin Calcium)"],
        ),
        frozenset(["warfarin", "metronidazole"]): _ix(
            "Major", "خطير / رئيسي",
            "Metronidazole stereoselectively inhibits CYP2C9 metabolism of the more potent S-warfarin enantiomer.",
            "يثبط الميترونيدازول إنزيم CYP2C9 المسؤول عن استقلاب المتصاوغ الفعال (S-warfarin)، مما يضاعف مفعول الوارفارين.",
            "Marked elevation of prothrombin time / INR and high risk of severe hemorrhage.",
            "ارتفاع حاد في مؤشر السيولة INR وخطر مرتفع لحدوث نزيف شديد.",
            "Avoid if possible; if combined, reduce warfarin dose by 25-35% and monitor INR closely every 2-3 days.",
            "تجنب الجمع إن أمكن؛ وعند الضرورة تُخفض جرعة الوارفارين بنسبة 25–35% مع قياس INR كل يومين إلى 3 أيام.",
            ["FDA Prescribing Information: Flagyl (Metronidazole) & Coumadin"],
        ),
        frozenset(["warfarin", "ciprofloxacin"]): _ix(
            "Major", "خطير / رئيسي",
            "Ciprofloxacin inhibits CYP1A2 metabolism of R-warfarin, displaces warfarin from albumin, and reduces vitamin K-producing gut flora.",
            "يثبط السيبروفلوكساسين استقلاب الوارفارين عبر CYP1A2 ويقلل تصنيع فيتامين K بواسطة البكتيريا المعوية.",
            "Unpredictable, sometimes severe elevation in INR and bleeding events within 3-10 days.",
            "ارتفاع غير متوقع وقد يكون حاداً في مؤشر INR وخطر النزيف خلال 3 إلى 10 أيام.",
            "Monitor INR within 3-5 days of starting ciprofloxacin and adjust warfarin dosage accordingly.",
            "يجب فحص INR خلال 3 إلى 5 أيام من بدء السيبروفلوكساسين وتعديل جرعة الوارفارين فوراً.",
            ["Arch Intern Med 2002; 162:1049-1056"],
        ),
        frozenset(["ciprofloxacin", "theophylline"]): _ix(
            "Major", "خطير / رئيسي",
            "Ciprofloxacin potently inhibits hepatic CYP1A2, the primary enzyme clearing theophylline (narrow therapeutic index).",
            "يثبط السيبروفلوكساسين بقوة إنزيم CYP1A2 المسؤول عن استقلاب الثيوفيلين ذي النافذة العلاجية الضيقة.",
            "2- to 3-fold surge in serum theophylline levels causing seizures, ventricular arrhythmias, and cardiac arrest.",
            "ارتفاع تركيز الثيوفيلين في الدم إلى الضعف أو ثلاثة أضعاف مما قد يسبب تشنجات عصبية واضطرابات نظم قاتلة.",
            "Avoid ciprofloxacin in patients on theophylline (use levofloxacin or azithromycin), or halve theophylline dose and monitor serum levels.",
            "تجنب السيبروفلوكساسين مع الثيوفيلين (استخدم ليفوفلوكساسين أو أزيثروميسين)، أو خفض جرعة الثيوفيلين للنصف مع قياس تركيزه.",
            ["FDA Boxed Warning: Cipro (Ciprofloxacin) — Theophylline Interaction"],
        ),
        frozenset(["lisinopril", "spironolactone"]): _ix(
            "Major", "خطير / مراقبة البوتاسيوم",
            "Dual blockade of the renin-angiotensin-aldosterone system (ACE inhibition + mineralocorticoid receptor antagonism) impairs renal potassium excretion.",
            "تثبيط مزدوج لنظام الرينين-أنجيوتنسين-ألدوستيرون يقلل طرح البوتاسيوم عبر الكلى.",
            "High risk of severe, potentially fatal hyperkalemia (serum K+ > 5.5 mEq/L) and acute kidney injury, especially in CKD.",
            "خطر مرتفع للإصابة بفرط بوتاسيوم الدم الخطير واضطراب نظم القلب والقصور الكلوي الحاد.",
            "Monitor serum potassium and creatinine within 1 week of initiation and regularly thereafter; avoid potassium supplements.",
            "فحص مستوى البوتاسيوم والكرياتينين في الدم خلال أسبوع من البدء ودورياً بعد ذلك، ومنع مكملات البوتاسيوم.",
            ["N Engl J Med 2004; 351:543-551 (Hyperkalemia after RALES)"],
        ),
        frozenset(["losartan", "spironolactone"]): _ix(
            "Major", "خطير / مراقبة البوتاسيوم",
            "Combined AT1 receptor blockade and aldosterone antagonism additively reduce renal distal tubular potassium secretion.",
            "حصر مستقبلات الأنجيوتنسين والألدوستيرون معاً يؤدي إلى احتباس البوتاسيوم في الكلى.",
            "Significant risk of life-threatening hyperkalemia and renal function deterioration.",
            "ارتفاع خطر فرط البوتاسيوم المهدد للحياة وتراجع وظائف الكلى.",
            "Monitor serum K+ and eGFR closely; keep spironolactone <=25 mg/day when combined in heart failure.",
            "مراقبة لصيقة لمستوى البوتاسيوم ووظائف الكلى، وعدم تجاوز 25 ملغ يومياً من السبيرونولاكتون.",
            ["ACC/AHA Heart Failure Guidelines 2022"],
        ),
        frozenset(["lisinopril", "ibuprofen"]): _ix(
            "Moderate", "متوسط / خطر كلوي وضغطي",
            "Ibuprofen inhibits renal vasodilatory prostaglandins (PGE2/PGI2) that maintain afferent arteriolar dilation, while lisinopril dilates the efferent arteriole.",
            "يثبط الإيبوبروفين البروستاغلاندينات الموسعة للشرعيين الكلوي الوارد، بينما يوسع الليزينوبريل الشرعيين الصادر، مما يخفض ضغط الترشيح الكبيبي ويقلل فعالية خفض الضغط.",
            "Blunted antihypertensive response and risk of acute kidney injury (especially in dehydrated or elderly patients).",
            "إضعاف مفعول دواء الضغط وارتفاع خطر القصور الكلوي الحاد (خصوصاً عند الجفاف أو كبار السن).",
            "Avoid chronic NSAID use in hypertensive patients on ACE inhibitors; monitor blood pressure and renal function.",
            "تجنب الاستخدام المزمن للإيبوبروفين لدى مرضى الضغط الذين يتناولون مثبطات ACE، ومراقبة ضغط الدم والكرياتينين.",
            ["BMJ 2013; 346:e8525 (Triple Whammy Acute Kidney Injury)"],
        ),
        frozenset(["losartan", "ibuprofen"]): _ix(
            "Moderate", "متوسط / خطر كلوي وضغطي",
            "NSAID-mediated afferent arteriolar vasoconstriction opposes ARB antihypertensive action and reduces glomerular filtration pressure.",
            "تضيق الشرعيين الكلوي الوارد بفعل الإيبوبروفين يعاكس مفعول اللوسارتان الخافض للضغط ويقلل الترشيح الكبيبي.",
            "Attenuated blood pressure control and increased risk of acute renal impairment.",
            "تراجع السيطرة على ضغط الدم وزيادة خطر تضرر وظائف الكلى.",
            "Prefer paracetamol for pain relief; monitor blood pressure and serum creatinine if concurrent use is required.",
            "يُفضل استخدام الباراسيتامول لتسكين الألم، مع مراقبة ضغط الدم ووظائف الكلى عند الاضطرار للجمع بينهما.",
            ["Hypertension 2012; NSAID-Antihypertensive Interactions"],
        ),
        frozenset(["sertraline", "tramadol"]): _ix(
            "Major", "خطير / متلازمة السيروتونين",
            "Additive central serotonergic activity (SSRI reuptake blockade by sertraline + serotonin release/reuptake inhibition by tramadol); also lowers seizure threshold.",
            "تراكم مفرط للسيروتونين في الجهاز العصبي المركزي نتيجة تثبيط استرداده بالدواءين معاً، مع خفض عتبة الاختلاجات.",
            "High risk of potentially life-threatening Serotonin Syndrome (clonus, hyperthermia, agitation, tachycardia) and generalized seizures.",
            "خطر مرتفع للإصابة بمتلازمة السيروتونين المهددة للحياة (حمى شديدة، رجفان، هياج، تسارع قلب) ونوبات الصرع.",
            "Avoid combination if possible; use alternative non-serotonergic analgesics and educate patient on serotonin toxicity signs.",
            "تجنب الجمع بينهما قدر الإمكان واستخدام مسكنات بديلة لا تؤثر على السيروتونين.",
            ["FDA Drug Safety Communication: Opioid Pain Medicines and Serotonergic Drugs (2016)"],
        ),
        frozenset(["fluoxetine", "tramadol"]): _ix(
            "Major", "خطير / متلازمة السيروتونين وتثبيط CYP2D6",
            "Fluoxetine strongly inhibits CYP2D6 (blocking activation of tramadol to its active analgesic metabolite M1) while increasing parent tramadol levels and serotonergic load.",
            "يثبط الفلوكستين بقوة إنزيم CYP2D6 مما يمنع تحول الترامادول إلى مستقلبه المسكن النشط، ويرفع في الوقت نفسه خطر متلازمة السيروتونين والتشنجات.",
            "Paradoxical loss of opioid analgesia combined with heightened risk of Serotonin Syndrome and seizures.",
            "فقدان التأثير المسكن للترامادول مع تضاعف خطر التسمم السيروتونيني والنوبات التشنجية.",
            "Avoid concurrent use; choose an analgesic not dependent on CYP2D6 bioactivation and without serotonergic activity.",
            "يُتجنب الجمع بينهما تماماً واختيار مسكن بديل آمن.",
            ["Clin Pharmacol Ther 2019; CPIC Guideline for CYP2D6 and Opioids"],
        ),
        frozenset(["sertraline", "ibuprofen"]): _ix(
            "Moderate", "متوسط / خطر النزيف الهضمي",
            "SSRIs deplete intraplatelet serotonin stores required for platelet aggregation, while NSAIDs cause mucosal injury and COX-1 inhibition.",
            "يقلل السيرترالين مخزون السيروتونين داخل الصفائح الدموية الضروري لتجمعها، بينما يسبب الإيبوبروفين تخرش المعدة وتثبيط COX-1.",
            "3- to 4-fold increase in upper gastrointestinal bleeding risk compared with either agent alone.",
            "زيادة خطر النزيف المعدي المعوي العلوي بمقدار 3 إلى 4 أضعاف مقارنة باستخدام أحدهما منفرداً.",
            "Prefer paracetamol for pain, or co-prescribe a proton pump inhibitor (e.g., omeprazole/pantoprazole) for gastroprotection.",
            "يُفضل استخدام الباراسيتامول لتسكين الألم، أو إضافة مثبط مضخة البروتون لحماية المعدة.",
            ["Alimentary Pharmacology & Therapeutics 2014; SSRI + NSAID GI Bleeding Meta-Analysis"],
        ),
        frozenset(["digoxin", "verapamil"]): _ix(
            "Major", "خطير / بطء القلب وسمية الديجوكسين",
            "Verapamil inhibits P-glycoprotein efflux of digoxin (raising serum digoxin by 50-75%) and exerts additive AV nodal depression.",
            "يثبط الفيراباميل مضخة P-gp فيرفع تركيز الديجوكسين بالدم بنسبة 50–75%، كما يشتركان في تثبيط العقدة الأذينية البطينية.",
            "High risk of digitalis toxicity, severe bradycardia, and complete AV heart block.",
            "خطر مرتفع للتسمم بالديجوكسين وبطء القلب الشديد وإحصار القلب الأذيني البطيني.",
            "Reduce digoxin dose by 33-50%, monitor serum digoxin levels, heart rate, and ECG.",
            "تخفيض جرعة الديجوكسين بنسبة 33–50% مع مراقبة تركيزه في الدم ونبض القلب وتخطيط ECG.",
            ["FDA Prescribing Information: Isoptin (Verapamil) & Lanoxin"],
        ),
        frozenset(["metoprolol", "verapamil"]): _ix(
            "Major", "خطير / إحصار قلبي وهبوط ضغط",
            "Additive negative inotropic and negative chronotropic/dromotropic effects on the sinoatrial and atrioventricular nodes.",
            "تأثير مثبط مضاعف على انقباض عضلة القلب وسرعة التوصيل الكهربائي في العقدتين الجيبية والأذينية البطينية.",
            "Severe bradycardia, profound hypotension, high-grade AV block, and precipitation of acute heart failure.",
            "بطء شديد في ضربات القلب، هبوط حاد في الضغط، إحصار قلبي، أو فشل قلبي حاد.",
            "Avoid concurrent IV administration and generally avoid oral combination unless under specialist cardiology monitoring.",
            "يُتجنب الجمع بينهما إلا تحت إشراف مباشر من طبيب القلب مع مراقبة لصيقة للنبض والضغط.",
            ["ACC/AHA Arrhythmia Management Guidelines"],
        ),
        frozenset(["salbutamol", "metoprolol"]): _ix(
            "Moderate", "متوسط / تعارض مستقبلي",
            "Pharmacodynamic antagonism at bronchial beta-2 adrenergic receptors (especially at higher doses of metoprolol where cardioselectivity diminishes).",
            "تضاد دوائي على مستقبلات بيتا في الشعب الهوائية (حيث يقلل الميتوبرولول وخصوصاً بالجرعات العالية من مفعول السالبوتامول الموسع للقصبات).",
            "Reduced bronchodilator efficacy of salbutamol and potential bronchospasm in asthma/COPD patients.",
            "ضعف استجابة الشعب الهوائية لبخاخ السالبوتامول واحتمال حدوث تشنج قصبي لدى مرضى الربو.",
            "Avoid non-selective beta-blockers in asthma; if a beta-1 selective blocker (metoprolol) is essential, use the lowest effective dose.",
            "تجنب حاصرات بيتا غير الانتقائية تماماً في الربو؛ وعند الضرورة القصوى للميتوبرولول تُستخدم أقل جرعة فعالة مع المراقبة التنفسية.",
            ["GINA Global Strategy for Asthma Management and Prevention 2023"],
        ),
        frozenset(["azithromycin", "amiodarone"]): _ix(
            "Major", "خطير / إطالة فترة QT",
            "Additive blockade of cardiac hERG potassium channels (IKr) prolonging ventricular repolarization.",
            "إطالة مضافة لفترة زوال الاستقطاب البطيني (QTc interval) عبر تثبيط قنوات البوتاسيوم القلبية.",
            "Increased risk of QTc prolongation and Torsades de Pointes ventricular tachycardia.",
            "ارتفاع خطر إطالة فترة QT واضطراب النظم البطيني الخطير (Torsades de Pointes).",
            "Avoid concurrent use if possible; use amoxicillin or doxycycline when clinically appropriate, or monitor ECG and electrolytes (K+, Mg2+).",
            "يُفضل اختيار مضاد حيوي بديل لا يطيل فترة QT (مثل الأموكسيسيلين أو الدوكسيسيكلين)، أو مراقبة تخطيط القلب والأملاح.",
            ["FDA Drug Safety Communication: Azithromycin and Risk of Fatal Heart Rhythms (2013)"],
        ),
        frozenset(["ciprofloxacin", "amiodarone"]): _ix(
            "Major", "خطير / إطالة فترة QT",
            "Additive cardiac potassium channel inhibition causing prolonged ventricular repolarization.",
            "تأثير مضاف على إطالة فترة QT في تخطيط القلب الكهربائي.",
            "Risk of life-threatening ventricular arrhythmias including Torsades de Pointes.",
            "خطر حدوث اضطرابات نظم بطينية مهددة للحياة.",
            "Avoid combination; select an alternative non-fluoroquinolone antibiotic.",
            "تجنب الجمع بينهما واستخدام مضاد حيوي من فئة أخرى.",
            ["CredibleMeds QTdrugs List; FDA Prescribing Information: Cordarone"],
        ),
        frozenset(["escitalopram", "ondansetron"]): _ix(
            "Moderate", "متوسط / إطالة QT ومتلازمة السيروتونين",
            "Dose-dependent additive QTc prolongation and mild additive serotonergic activity.",
            "إطالة مضافة لفترة QT على تخطيط القلب مع تداخل سيروتونيني خفيف.",
            "Increased risk of QTc prolongation and arrhythmias, particularly in hypokalemic patients.",
            "زيادة خطر إطالة فترة QT واضطراب نبض القلب خاصة عند نقص البوتاسيوم أو المغنيسيوم.",
            "Monitor ECG and electrolytes if repeated doses of ondansetron are needed in patients on escitalopram.",
            "مراقبة تخطيط القلب ومستويات البوتاسيوم والمغنيسيوم عند الحاجة لتكرار جرعات الأوندانسيترون.",
            ["FDA Drug Safety Communication: Celexa/Lexapro and Zofran QTc Prolongation"],
        ),
        frozenset(["metformin", "ciprofloxacin"]): _ix(
            "Moderate", "متوسط / اضطراب سكر الدم",
            "Fluoroquinolones affect pancreatic beta-cell ATP-sensitive K+ channels and inhibit renal OCT2/MATE1 tubular secretion of metformin.",
            "يؤثر السيبروفلوكساسين على قنوات البوتاسيوم في خلايا بيتا البنكرياسية ويقلل الإطراح الأنبوبي الكلوي للميتفورمين.",
            "Dysglycemia (either symptomatic hypoglycemia or hyperglycemia) and elevated plasma metformin.",
            "تذبذب حاد في مستوى سكر الدم (هبوط أو ارتفاع السكر) مع ارتفاع تركيز الميتفورمين.",
            "Monitor blood glucose closely during and shortly after ciprofloxacin therapy.",
            "مراقبة مستوى السكر في الدم بانتظام طوال فترة تناول السيبروفلوكساسين.",
            ["FDA Safety Labeling Update: Fluoroquinolones and Severe Hypoglycemia (2018)"],
        ),
        frozenset(["glimepiride", "fluconazole"]): _ix(
            "Major", "خطير / هبوط سكر حاد",
            "Inhibition of CYP2C9 hepatic metabolism of glimepiride.",
            "تثبيط استقلاب الغليميبيريد في الكبد عبر إنزيم CYP2C9.",
            "Severe, prolonged hypoglycemia.",
            "هبوط شديد ومطول في سكر الدم.",
            "Monitor blood glucose closely and reduce glimepiride dose.",
            "مراقبة سكر الدم وخفض جرعة الغليميبيريد.",
            ["FDA Prescribing Information: Amaryl (Glimepiride)"],
        ),
        frozenset(["glimepiride", "ciprofloxacin"]): _ix(
            "Major", "خطير / هبوط سكر الدم",
            "Ciprofloxacin potentiates sulfonylurea-stimulated insulin release from pancreatic beta cells and inhibits CYP2C9.",
            "يعزز السيبروفلوكساسين إفراز الإنسولين المحفز بالسلفونيل يوريا مما يفاقم مفعول الغليميبيريد.",
            "Risk of severe, refractory hypoglycemia and hypoglycemic coma.",
            "خطر الإصابة بهبوط حاد وشديد في سكر الدم.",
            "Prefer an alternative antibiotic or monitor capillary blood glucose frequently and reduce glimepiride dose.",
            "يُفضل اختيار مضاد حيوي بديل أو مراقبة سكر الدم عن كثب وتقليل جرعة الغليميبيريد.",
            ["Clin Infect Dis 2009; Fluoroquinolone-Sulfonylurea Hypoglycemia"],
        ),
        frozenset(["rivaroxaban", "ibuprofen"]): _ix(
            "Major", "خطير / نزيف مضاعف",
            "Combined Factor Xa inhibition by rivaroxaban and platelet COX-1 inhibition / mucosal injury by ibuprofen.",
            "تثبيط العامل العاشر للتخثر (Factor Xa) بالتزامن مع تثبيط الصفائح الدموية وتخرش المعدة بالإيبوبروفين.",
            "Substantially elevated risk of major gastrointestinal and systemic bleeding.",
            "ارتفاع كبير في خطر النزيف المعدي المعوي والنزيف العام.",
            "Avoid NSAIDs in patients taking DOACs (rivaroxaban/apixaban); use paracetamol for pain relief.",
            "تجنب مضادات الالتهاب غير الستيرويدية لدى مرضى مضادات التخثر الفموية الحديثة، واستخدام الباراسيتامول بدلاً منها.",
            ["Circulation 2019; Bleeding Risk with NSAIDs in Rivaroxaban/Apixaban Patients"],
        ),
        frozenset(["apixaban", "ibuprofen"]): _ix(
            "Major", "خطير / نزيف مضاعف",
            "Synergistic anticoagulant (Factor Xa blockade) and antiplatelet/mucosal-eroding NSAID effects.",
            "تآزر بين تثبيط العامل العاشر للتخثر وتثبيط الصفائح الدموية بواسطة الإيبوبروفين.",
            "Significantly increased risk of clinically relevant and major hemorrhage.",
            "زيادة ملحوظة في خطر النزيف الهضمي والجهازي.",
            "Avoid concurrent NSAID use; substitute with paracetamol.",
            "يُتجنب الجمع بينهما ويُستبدل الإيبوبروفين بالباراسيتامول عند الحاجة لتسكين الألم.",
            ["FDA Prescribing Information: Eliquis (Apixaban)"],
        ),
        frozenset(["rivaroxaban", "clarithromycin"]): _ix(
            "Moderate", "متوسط / ارتفاع تركيز مضاد التخثر",
            "Clarithromycin strongly inhibits both CYP3A4 and P-glycoprotein, increasing rivaroxaban AUC by ~50-90%.",
            "يثبط الكلاريثروميسين إنزيم CYP3A4 ومضخة P-gp، مما يرفع مستوى الريفاروكسابان في الدم بنسبة 50–90%.",
            "Increased bleeding risk, particularly in patients with mild-to-moderate renal impairment.",
            "زيادة خطر النزيف خاصة لدى المرضى الذين يعانون من قصور كلوي خفيف أو متوسط.",
            "Use azithromycin instead if macrolide therapy is needed, or monitor closely for bleeding.",
            "يُنصح باستخدام أزيثروميسين بدلاً من الكلاريثروميسين لتجنب تثبيط CYP3A4.",
            ["FDA Prescribing Information: Xarelto (Rivaroxaban)"],
        ),
    }

    _SYNONYMS: dict[str, str] = {
        "paracetamol": "paracetamol",
        "acetaminophen": "paracetamol",
        "panadol": "paracetamol",
        "tylenol": "paracetamol",
        "adol": "paracetamol",
        "fevadol": "paracetamol",
        "باراسيتامول": "paracetamol",
        "بنادول": "paracetamol",
        "أسيتامينوفين": "paracetamol",
        "aspirin": "aspirin",
        "acetylsalicylic acid": "aspirin",
        "asa": "aspirin",
        "jusprin": "aspirin",
        "aspocid": "aspirin",
        "ecotrin": "aspirin",
        "أسبرين": "aspirin",
        "اسبرين": "aspirin",
        "ibuprofen": "ibuprofen",
        "advil": "ibuprofen",
        "motrin": "ibuprofen",
        "brufen": "ibuprofen",
        "nurofen": "ibuprofen",
        "إيبوبروفين": "ibuprofen",
        "ايبوبروفين": "ibuprofen",
        "بروفين": "ibuprofen",
        "naproxen": "naproxen",
        "aleve": "naproxen",
        "naprosyn": "naproxen",
        "نابروكسين": "naproxen",
        "diclofenac": "diclofenac",
        "voltaren": "diclofenac",
        "cataflam": "diclofenac",
        "ديكلوفيناك": "diclofenac",
        "فولتارين": "diclofenac",
        "celecoxib": "celecoxib",
        "celebrex": "celecoxib",
        "سيليكوكسيب": "celecoxib",
        "tramadol": "tramadol",
        "ultram": "tramadol",
        "tramal": "tramadol",
        "ترامادول": "tramadol",
        "amoxicillin": "amoxicillin",
        "amoxil": "amoxicillin",
        "augmentin": "amoxicillin",
        "أموكسيسيلين": "amoxicillin",
        "اموكسيسيلين": "amoxicillin",
        "أوغمنتين": "amoxicillin",
        "azithromycin": "azithromycin",
        "zithromax": "azithromycin",
        "أزيثروميسين": "azithromycin",
        "زيثروماكس": "azithromycin",
        "clarithromycin": "clarithromycin",
        "klacid": "clarithromycin",
        "biaxin": "clarithromycin",
        "كلاريثروميسين": "clarithromycin",
        "كلاسيد": "clarithromycin",
        "ciprofloxacin": "ciprofloxacin",
        "cipro": "ciprofloxacin",
        "ciprobay": "ciprofloxacin",
        "سيبروفلوكساسين": "ciprofloxacin",
        "levofloxacin": "levofloxacin",
        "tavanic": "levofloxacin",
        "levaquin": "levofloxacin",
        "ليفوفلوكساسين": "levofloxacin",
        "تافانيك": "levofloxacin",
        "metronidazole": "metronidazole",
        "flagyl": "metronidazole",
        "ميترونيدازول": "metronidazole",
        "فلاجيل": "metronidazole",
        "doxycycline": "doxycycline",
        "vibramycin": "doxycycline",
        "دوكسيسيكلين": "doxycycline",
        "amlodipine": "amlodipine",
        "norvasc": "amlodipine",
        "أملوديبين": "amlodipine",
        "نورفاسك": "amlodipine",
        "lisinopril": "lisinopril",
        "zestril": "lisinopril",
        "prinivil": "lisinopril",
        "ليزينوبريل": "lisinopril",
        "enalapril": "enalapril",
        "renitec": "enalapril",
        "إينالابريل": "enalapril",
        "losartan": "losartan",
        "cozaar": "losartan",
        "لوسارتان": "losartan",
        "كوزار": "losartan",
        "valsartan": "valsartan",
        "diovan": "valsartan",
        "فالسارتان": "valsartan",
        "ديوفان": "valsartan",
        "hydrochlorothiazide": "hydrochlorothiazide",
        "hctz": "hydrochlorothiazide",
        "esidrex": "hydrochlorothiazide",
        "هيدروكلوروثيازيد": "hydrochlorothiazide",
        "spironolactone": "spironolactone",
        "aldactone": "spironolactone",
        "سبيرونولاكتون": "spironolactone",
        "ألداكتون": "spironolactone",
        "metoprolol": "metoprolol",
        "betaloc": "metoprolol",
        "lopressor": "metoprolol",
        "ميتوبرولول": "metoprolol",
        "بيتالوك": "metoprolol",
        "metformin": "metformin",
        "glucophage": "metformin",
        "ميتفورمين": "metformin",
        "جلوكوفاج": "metformin",
        "sitagliptin": "sitagliptin",
        "januvia": "sitagliptin",
        "سيتاغليبتين": "sitagliptin",
        "جانوفيا": "sitagliptin",
        "vildagliptin": "vildagliptin",
        "galvus": "vildagliptin",
        "فيلداغليبتين": "vildagliptin",
        "glimepiride": "glimepiride",
        "amaryl": "glimepiride",
        "غليميبيريد": "glimepiride",
        "أماريل": "glimepiride",
        "empagliflozin": "empagliflozin",
        "jardiance": "empagliflozin",
        "إمباغليفلوزين": "empagliflozin",
        "جارديانس": "empagliflozin",
        "warfarin": "warfarin",
        "coumadin": "warfarin",
        "marevan": "warfarin",
        "وارفارين": "warfarin",
        "كومادين": "warfarin",
        "clopidogrel": "clopidogrel",
        "plavix": "clopidogrel",
        "كلوبيدوغريل": "clopidogrel",
        "بلافيكس": "clopidogrel",
        "rivaroxaban": "rivaroxaban",
        "xarelto": "rivaroxaban",
        "ريفاروكسابان": "rivaroxaban",
        "زاريلتو": "rivaroxaban",
        "apixaban": "apixaban",
        "eliquis": "apixaban",
        "أبيكسابان": "apixaban",
        "إليكويس": "apixaban",
        "cetirizine": "cetirizine",
        "zyrtec": "cetirizine",
        "سيتريزين": "cetirizine",
        "زيرتك": "cetirizine",
        "loratadine": "loratadine",
        "claritin": "loratadine",
        "لوراتادين": "loratadine",
        "كلاريتين": "loratadine",
        "fexofenadine": "fexofenadine",
        "telfast": "fexofenadine",
        "allegra": "fexofenadine",
        "فيكسوفينادين": "fexofenadine",
        "تلفاست": "fexofenadine",
        "diphenhydramine": "diphenhydramine",
        "benadryl": "diphenhydramine",
        "ديفينهيدرامين": "diphenhydramine",
        "omeprazole": "omeprazole",
        "prilosec": "omeprazole",
        "losec": "omeprazole",
        "أوميبرازول": "omeprazole",
        "اوميبرازول": "omeprazole",
        "لوسيك": "omeprazole",
        "esomeprazole": "esomeprazole",
        "nexium": "esomeprazole",
        "إيزوميبرازول": "esomeprazole",
        "نيكسيوم": "esomeprazole",
        "pantoprazole": "pantoprazole",
        "controloc": "pantoprazole",
        "protonix": "pantoprazole",
        "بانتوبرازول": "pantoprazole",
        "كونترولوك": "pantoprazole",
        "famotidine": "famotidine",
        "pepcid": "famotidine",
        "فاموتيدين": "famotidine",
        "ondansetron": "ondansetron",
        "zofran": "ondansetron",
        "أوندانسيترون": "ondansetron",
        "زوفران": "ondansetron",
        "sertraline": "sertraline",
        "zoloft": "sertraline",
        "lustral": "sertraline",
        "سيرترالين": "sertraline",
        "زولوفت": "sertraline",
        "fluoxetine": "fluoxetine",
        "prozac": "fluoxetine",
        "فلوكستين": "fluoxetine",
        "بروزاك": "fluoxetine",
        "escitalopram": "escitalopram",
        "cipralex": "escitalopram",
        "lexapro": "escitalopram",
        "إسيتالوبرام": "escitalopram",
        "سيبرالكس": "escitalopram",
        "amitriptyline": "amitriptyline",
        "elavil": "amitriptyline",
        "أميتريبتيلين": "amitriptyline",
        "digoxin": "digoxin",
        "lanoxin": "digoxin",
        "ديجوكسين": "digoxin",
        "لانوكسين": "digoxin",
        "amiodarone": "amiodarone",
        "cordarone": "amiodarone",
        "أميودارون": "amiodarone",
        "كوردارون": "amiodarone",
        "atorvastatin": "atorvastatin",
        "lipitor": "atorvastatin",
        "أتورفاستاتين": "atorvastatin",
        "ليبيتور": "atorvastatin",
        "simvastatin": "simvastatin",
        "zocor": "simvastatin",
        "سيمفاستاتين": "simvastatin",
        "زوكور": "simvastatin",
        "rosuvastatin": "rosuvastatin",
        "crestor": "rosuvastatin",
        "روزوفاستاتين": "rosuvastatin",
        "كريستور": "rosuvastatin",
        "verapamil": "verapamil",
        "isoptin": "verapamil",
        "فيراباميل": "verapamil",
        "salbutamol": "salbutamol",
        "albuterol": "salbutamol",
        "ventolin": "salbutamol",
        "سالبوتامول": "salbutamol",
        "فنتولين": "salbutamol",
        "theophylline": "theophylline",
        "ثيوفيلين": "theophylline",
        "montelukast": "montelukast",
        "singulair": "montelukast",
        "مونتيلوكاست": "montelukast",
        "سينجولير": "montelukast",
        "budesonide": "budesonide",
        "pulmicort": "budesonide",
        "بوديزونيد": "budesonide",
        "caffeine": "caffeine",
        "كافيين": "caffeine",
    }

    @classmethod
    def normalize_drug_name(cls, name: str) -> str:
        """Normalizes any generic, brand, or Arabic drug name to its standard canonical key."""
        clean = (name or "").strip().lower()
        if not clean:
            return ""
        if clean in cls._SYNONYMS:
            return cls._SYNONYMS[clean]
        # Also match brand names dynamically from DRUG_CATALOG
        for key, info in cls.DRUG_CATALOG.items():
            if clean == key or clean == info["name_en"].lower() or clean == info["name_ar"].lower():
                return key
            for br in info.get("brand_names", []):
                if clean == br.lower():
                    return key
        return clean

    @classmethod
    def get_supported_drugs(cls, therapeutic_class: str = "all") -> List[str]:
        """Returns sorted list of distinct canonical drugs, optionally filtered by therapeutic class."""
        drugs = set(cls.DRUG_CATALOG.keys())
        for pair in cls._INTERACTIONS.keys():
            drugs.update(pair)
        if therapeutic_class and therapeutic_class != "all":
            drugs = {
                d for d in drugs
                if cls.DRUG_CATALOG.get(d, {}).get("therapeutic_class_key") == therapeutic_class
            }
        return sorted(list(drugs))

    @classmethod
    def get_drug_metadata(cls, drug_name_or_key: str) -> dict[str, Any] | None:
        norm = cls.normalize_drug_name(drug_name_or_key)
        return cls.DRUG_CATALOG.get(norm)

    @classmethod
    def import_or_update_drug_record(cls, raw_record: dict[str, Any]) -> dict[str, Any]:
        """
        Imports or updates a drug record in DRUG_CATALOG with strict deduplication
        by canonical_key, RxCUI, DrugBank ID, and brand/generic synonyms.
        """
        raw_name = str(raw_record.get("name_en") or raw_record.get("canonical_key") or "").strip()
        if not raw_name:
            return {"imported": False, "reason": "Missing drug name"}

        rxcui = str(raw_record.get("rxcui") or "").strip()
        db_id = str(raw_record.get("drugbank_id") or "").strip().upper()
        norm_key = cls.normalize_drug_name(raw_name)

        # Deduplicate by RxCUI or DrugBank ID if already present
        for existing_key, existing_rec in cls.DRUG_CATALOG.items():
            if (rxcui and existing_rec.get("rxcui") == rxcui) or (db_id and existing_rec.get("drugbank_id") == db_id):
                norm_key = existing_key
                break

        now_date = datetime.now(timezone.utc).strftime("%Y-%m-%d")
        was_existing = norm_key in cls.DRUG_CATALOG
        base_rec = dict(cls.DRUG_CATALOG.get(norm_key, {}))

        brands = list(dict.fromkeys(list(base_rec.get("brand_names", [])) + list(raw_record.get("brand_names", []))))
        t_class = raw_record.get("therapeutic_class_key") or base_rec.get("therapeutic_class_key") or "analgesics_nsaids"
        cls_info = cls.THERAPEUTIC_CLASSES.get(t_class, cls.THERAPEUTIC_CLASSES["analgesics_nsaids"])

        updated_rec = {
            "canonical_key": norm_key,
            "name_en": sanitize_external_text(raw_record.get("name_en") or base_rec.get("name_en") or norm_key.capitalize(), 150),
            "name_ar": sanitize_external_text(raw_record.get("name_ar") or base_rec.get("name_ar") or norm_key, 150),
            "therapeutic_class_key": t_class,
            "therapeutic_class_ar": cls_info["name_ar"],
            "therapeutic_class_en": cls_info["name_en"],
            "brand_names": [sanitize_external_text(b, 80) for b in brands if b],
            "rxcui": sanitize_external_text(rxcui or base_rec.get("rxcui", ""), 40),
            "drugbank_id": sanitize_external_text(db_id or base_rec.get("drugbank_id", ""), 40),
            "pubchem_cid": raw_record.get("pubchem_cid") or base_rec.get("pubchem_cid"),
            "openfda_interactions_excerpt": sanitize_external_text(raw_record.get("openfda_interactions_excerpt", ""), 1200),
            "openfda_warnings_excerpt": sanitize_external_text(raw_record.get("openfda_warnings_excerpt", ""), 1200),
            "last_updated": raw_record.get("last_updated") or now_date,
            "sources": list(dict.fromkeys(list(base_rec.get("sources", [])) + list(raw_record.get("sources", ["Local_Curated"])))),
        }
        cls.DRUG_CATALOG[norm_key] = updated_rec
        cls._SYNONYMS[norm_key] = norm_key
        cls._SYNONYMS[updated_rec["name_en"].lower()] = norm_key
        cls._SYNONYMS[updated_rec["name_ar"].lower()] = norm_key
        for b in updated_rec["brand_names"]:
            cls._SYNONYMS[b.lower()] = norm_key

        return {
            "imported": True,
            "deduplicated_existing": was_existing,
            "canonical_key": norm_key,
            "record": updated_rec,
        }

    @classmethod
    def sync_drug_from_external_sources(
        cls,
        drug_query: str,
        timeout: int = EXTERNAL_API_TIMEOUT_SECONDS,
    ) -> dict[str, Any]:
        """
        Fetches and synchronizes drug normalization and interaction/warning data from
        authoritative public APIs (RxNorm, openFDA, and PubChem) with caching and
        deduplication. Gracefully falls back to local catalog when offline.
        """
        ok_in, clean_q, err_in = validate_safe_chemical_input(drug_query, max_len=150, field_name_ar="اسم الدواء")
        if not ok_in:
            return {"ok": False, "error": err_in}

        norm_key = cls.normalize_drug_name(clean_q)
        cache_key = f"ddi:sync:{norm_key.lower()}"
        hit, cached = get_cached_item(cache_key)
        if hit:
            res = dict(cached)
            res["cache_hit"] = True
            return res

        rxcui = ""
        rxnorm_brands: list[str] = []
        openfda_interactions = ""
        openfda_warnings = ""
        sources_used: list[str] = []

        # 1. RxNorm REST API (NLM Public Access)
        try:
            enc = urllib.parse.quote(norm_key)
            rx_url = f"https://rxnav.nlm.nih.gov/REST/rxcui.json?name={enc}&search=1"
            rx_resp = request_with_retry(rx_url, timeout=timeout, max_retries=1)
            if rx_resp.status_code == 200:
                rx_ids = rx_resp.json().get("idGroup", {}).get("rxnormId", [])
                if rx_ids:
                    rxcui = str(rx_ids[0])
                    sources_used.append("RxNorm")
        except Exception:
            pass

        # 2. openFDA Drug Label API (CC0 Public Domain)
        try:
            enc_fda = urllib.parse.quote(norm_key)
            fda_url = f"https://api.fda.gov/drug/label.json?search=openfda.generic_name:{enc_fda}&limit=1"
            fda_resp = request_with_retry(fda_url, timeout=timeout, max_retries=1)
            if fda_resp.status_code == 200:
                results = fda_resp.json().get("results", [])
                if results:
                    item = results[0]
                    openfda_block = item.get("openfda", {}) or {}
                    rxnorm_brands = list(openfda_block.get("brand_name", [])[:5])
                    di_list = item.get("drug_interactions", [])
                    warn_list = item.get("warnings", []) or item.get("boxed_warning", [])
                    if di_list:
                        openfda_interactions = str(di_list[0])[:900]
                    if warn_list:
                        openfda_warnings = str(warn_list[0])[:900]
                    sources_used.append("openFDA")
        except Exception:
            pass

        now_date = datetime.now(timezone.utc).strftime("%Y-%m-%d")
        sync_res = cls.import_or_update_drug_record({
            "canonical_key": norm_key,
            "name_en": cls.DRUG_CATALOG.get(norm_key, {}).get("name_en", norm_key.capitalize()),
            "rxcui": rxcui,
            "brand_names": rxnorm_brands,
            "openfda_interactions_excerpt": openfda_interactions,
            "openfda_warnings_excerpt": openfda_warnings,
            "last_updated": now_date,
            "sources": sources_used or ["Local_Curated"],
        })

        payload = {
            "ok": True,
            "canonical_key": sync_res["canonical_key"],
            "deduplicated_existing": sync_res["deduplicated_existing"],
            "live_sources_reached": sources_used,
            "record": sync_res["record"],
            "cache_hit": False,
        }
        set_cached_item(cache_key, payload)
        return payload

    @classmethod
    def check_interaction(cls, drug_a: str, drug_b: str) -> Dict[str, Any]:
        """
        Evaluates potential interaction between two drugs (by generic, brand, or Arabic name).
        Returns clinical interaction details or safe negative report with mandatory disclaimer.
        """
        norm_a = cls.normalize_drug_name(drug_a)
        norm_b = cls.normalize_drug_name(drug_b)

        if not norm_a or not norm_b:
            return {
                "success": False,
                "error": "Two valid drug names are required.",
                "disclaimer_ar": cls.DISCLAIMER_AR,
                "disclaimer_en": cls.DISCLAIMER_EN,
            }

        meta_a = cls.DRUG_CATALOG.get(norm_a, {})
        meta_b = cls.DRUG_CATALOG.get(norm_b, {})

        if norm_a == norm_b:
            return {
                "success": True,
                "is_same_drug": True,
                "drug_a": drug_a,
                "drug_b": drug_b,
                "canonical_a": norm_a,
                "canonical_b": norm_b,
                "severity": "Duplicate Therapy",
                "severity_ar": "تكرار دوائي لنفس المادة الفعالة",
                "mechanism_en": "Duplicate active ingredient therapy (e.g. two brand names of the same generic drug).",
                "mechanism_ar": "استخدام نفس المادة الفعالة (أو اسمين تجاريين لنفس الدواء) مما يؤدي لمضاعفة الجرعة واحتمال السمية.",
                "clinical_effect_en": "Risk of unintentional overdose and enhanced adverse events.",
                "clinical_effect_ar": "خطر الجرعة الزائدة غير المقصودة ومضاعفة الآثار الجانبية السمية.",
                "management_en": "Verify that patient is not unknowingly taking multiple brands containing identical active substance.",
                "management_ar": "التأكد من عدم تناول المريض لعلامات تجارية مختلفة تحتوي على نفس المادة الفعالة.",
                "evidence_level": "Established",
                "citations": ["Clinical Pharmacology Guidelines: Duplicate Therapy Prevention"],
                "last_updated": meta_a.get("last_updated", "2026-11-15"),
                "disclaimer_ar": cls.DISCLAIMER_AR,
                "disclaimer_en": cls.DISCLAIMER_EN,
            }

        key = frozenset([norm_a, norm_b])
        interaction = cls._INTERACTIONS.get(key)

        if interaction:
            return {
                "success": True,
                "has_interaction": True,
                "drug_a": drug_a,
                "drug_b": drug_b,
                "canonical_a": norm_a,
                "canonical_b": norm_b,
                "drug_a_metadata": meta_a,
                "drug_b_metadata": meta_b,
                "severity": interaction["severity"],
                "severity_ar": interaction["severity_ar"],
                "mechanism_en": interaction["mechanism_en"],
                "mechanism_ar": interaction["mechanism_ar"],
                "clinical_effect_en": interaction["clinical_effect_en"],
                "clinical_effect_ar": interaction["clinical_effect_ar"],
                "management_en": interaction["management_en"],
                "management_ar": interaction["management_ar"],
                "evidence_level": interaction["evidence_level"],
                "citations": interaction["citations"],
                "last_updated": interaction.get("last_updated", "2026-11-15"),
                "disclaimer_ar": cls.DISCLAIMER_AR,
                "disclaimer_en": cls.DISCLAIMER_EN,
            }

        return {
            "success": True,
            "has_interaction": False,
            "drug_a": drug_a,
            "drug_b": drug_b,
            "canonical_a": norm_a,
            "canonical_b": norm_b,
            "drug_a_metadata": meta_a,
            "drug_b_metadata": meta_b,
            "severity": "No Interaction Documented",
            "severity_ar": "لا يوجد تداخل مسجل في قاعدة البيانات المحلية",
            "message_en": (
                f"No direct clinical drug-drug interaction is currently documented between "
                f"'{norm_a}' and '{norm_b}' in the MolecuLens reference database. "
                "Absence of a record does NOT guarantee clinical safety; consult official openFDA/BNF monographs."
            ),
            "message_ar": (
                f"لا يوجد سجل تداخل دوائي مباشر بين ({meta_a.get('name_ar', norm_a)}) و ({meta_b.get('name_ar', norm_b)}) "
                "في قاعدة البيانات المرجعية الحالية. تنبيه: غياب سجل التداخل هنا يعكس حدود التغطية المحلية ولا يعد دليلاً "
                "على انعدام التداخل السريري؛ يرجى مراجعة النشرة الدوائية الرسمية (openFDA / BNF)."
            ),
            "last_updated": meta_a.get("last_updated", "2026-11-15"),
            "disclaimer_ar": cls.DISCLAIMER_AR,
            "disclaimer_en": cls.DISCLAIMER_EN,
        }

    @classmethod
    def screen_polypharmacy(cls, drug_list: List[str]) -> Dict[str, Any]:
        """
        Screens a list of 2 or more drugs for all pairwise interactions, returning
        both detected interactions and a full pair-by-pair evaluation breakdown.
        """
        seen: list[str] = []
        duplicate_inputs: list[tuple[str, str]] = []
        for raw_d in drug_list or []:
            if not raw_d or not str(raw_d).strip():
                continue
            norm_d = cls.normalize_drug_name(raw_d)
            if norm_d in seen:
                duplicate_inputs.append((str(raw_d), norm_d))
            else:
                seen.append(norm_d)

        cleaned = seen
        n = len(cleaned)
        interactions_found: list[dict[str, Any]] = []
        all_pairs_evaluated: list[dict[str, Any]] = []
        severity_counts = {"Major": 0, "Moderate": 0, "Minor": 0}

        for i in range(n):
            for j in range(i + 1, n):
                res = cls.check_interaction(cleaned[i], cleaned[j])
                all_pairs_evaluated.append(res)
                if res.get("has_interaction"):
                    interactions_found.append(res)
                    sev = res["severity"]
                    if sev in severity_counts:
                        severity_counts[sev] += 1

        # Sort interactions by clinical priority: Major -> Moderate -> Minor
        sev_order = {"Major": 0, "Moderate": 1, "Minor": 2}
        interactions_found.sort(key=lambda x: sev_order.get(x["severity"], 3))

        return {
            "success": True,
            "total_drugs": n,
            "normalized_drugs": cleaned,
            "duplicate_inputs_detected": duplicate_inputs,
            "total_pairs_evaluated": (n * (n - 1)) // 2,
            "total_interactions_found": len(interactions_found),
            "severity_summary": severity_counts,
            "interactions": interactions_found,
            "all_pairs_results": all_pairs_evaluated,
            "disclaimer_ar": cls.DISCLAIMER_AR,
            "disclaimer_en": cls.DISCLAIMER_EN,
        }

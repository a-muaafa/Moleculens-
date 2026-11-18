from typing import Any
from chemistry_engine.structure_analyzer import compute_rdkit_molecular_metrics
from config.settings import OP_SCIENTIFIC_CALCULATION, OP_SCIENTIFIC_RETRIEVAL
from core.schemas import build_provenance_tag

# Laboratory safety & pre-lab preparation profiles for pharmaceutical and chemical substances
PHARMACY_PRELAB_NOTES: dict[str, dict[str, str]] = {
    "aspirin": {
        "hazard_ar": "مهيج خفيف للعينين والجهاز التنفسي عند استنشاق المسحوق الجاف (GHS07 — Irritant).",
        "ppe_ar": "معطف المختبر (Lab Coat)، قفازات النتريل (Nitrile Gloves)، ونظارات الوقاية (Safety Goggles).",
        "quality_test_ar": "اختبار النقاوة بكاشف كلوريد الحديديك (FeCl₃ Test): عدم ظهور لون بنفسجي يدل على خلو الأسبرين من حمض الساليسيليك الحر غير المتفاعل.",
        "handling_ar": "يحفظ في عبوات محكمة بعيداً عن الرطوبة والحرارة لمنع تحلله المائي إلى حمض الساليسيليك وحمض الخليك.",
    },
    "paracetamol": {
        "hazard_ar": "آمن نسبياً عند التعامل المخبري المعتاد؛ ضار عند الابتلاع بجرعات كبيرة (سمية كبدية Hepatotoxicity).",
        "ppe_ar": "معطف المختبر، قفازات واقية، وتجنب تطاير الغبار أثناء الوزن.",
        "quality_test_ar": "التحقق من درجة الانصهار القياسية (169–170.5 °C) وامتصاص الأشعة فوق البنفسجية (UV Spectrophotometry at λmax ≈ 243 nm).",
        "handling_ar": "يُوزن بدقة على ميزان تحليلي حساس ويحمى من الضوء المباشر.",
    },
    "ibuprofen": {
        "hazard_ar": "مهيج للعينين والأغشية المخاطية عند تطاير البودرة الفعالة (API Dust Hazard).",
        "ppe_ar": "معطف المختبر، قفازات النتريل، وكمامة واقية من الغبار عند وزن المسحوق الخام.",
        "quality_test_ar": "معايرة حمض-قاعدة في وسط كحولي باستخدام هيدروكسيد الصوديوم القياسي (NaOH Titration) ودليل الفينولفثالين.",
        "handling_ar": "قليل الذوبان جداً في الماء؛ يذاب في المذيبات العضوية القطبية مثل الإيثانول أو المحاليل القلوية.",
    },
    "sitagliptin": {
        "hazard_ar": "مركب دوائي نشط عالي الانتقائية؛ يجب تجنب ملامسة الجلد أو استنشاق غباره في المختبر.",
        "ppe_ar": "معطف المختبر، قفازات النتريل، ونظارات السلامة.",
        "quality_test_ar": "التحقق من النقاوة الفراغية (Chiral HPLC للمتصاوغ النشط 3R) ومطابقة طيف الأشعة تحت الحمراء (FT-IR).",
        "handling_ar": "يتعامل معه غالباً بصورة ملح فوسفات أحادي الإماهة (Sitagliptin Phosphate Monohydrate) لتحسين الذوبانية والثبات.",
    },
    "vildagliptin": {
        "hazard_ar": "مادة دوائية فعالة تحتوي على مجموعة نتريل؛ تجنب ملامسة المسحوق للجلد والعينين.",
        "ppe_ar": "معطف المختبر، قفازات النتريل، ونظارات الوقاية.",
        "quality_test_ar": "الفحص الكروماتوغرافي السائل عالي الأداء (RP-HPLC) للتحقق من عدم تحلل مجموعة السيانو (CN).",
        "handling_ar": "يحفظ في مكان جاف وبارد بعيداً عن العوامل المؤكسدة والقلويات القوية.",
    },
    "metformin": {
        "hazard_ar": "مهيج خفيف للجلد والعينين في صورته الخام.",
        "ppe_ar": "معطف المختبر وقفازات واقية.",
        "quality_test_ar": "المعايرة في وسط لا مائي (Non-aqueous Titration) أو القياس الطيفي عند 232 nm.",
        "handling_ar": "شديد الذوبان في الماء في صورة هيدروكلوريد الميتفورمين (Metformin HCl).",
    },
    "caffeine": {
        "hazard_ar": "ضار عند الابتلاع بكميات مركزة (GHS07 — Acute Oral Toxicity Category 4).",
        "ppe_ar": "معطف المختبر، قفازات النتريل، ونظارات السلامة.",
        "quality_test_ar": "اختبار الموريكسيد اللوني للبيورينات (Murexide Test) وقياس درجة الانصهار (235–238 °C).",
        "handling_ar": "قابل للتسامي عند التسخين العالي؛ يستخلص مخبرياً بالمذيبات العضوية ثنائية الطور.",
    },
    "amoxicillin": {
        "hazard_ar": "محسس تنفسي وجلدي قوي (GHS08 — Respiratory & Skin Sensitizer) لمن لديهم حساسية للبنسلينات.",
        "ppe_ar": "معطف المختبر، قفازات مزدوجة، كمامة واقية، والعمل داخل خزانة سحب الغازات (Fume Hood) عند تحضير المعلقات.",
        "quality_test_ar": "المعايرة اليودية أو الكروماتوغرافيا السائلة (HPLC) للتحقق من سلامة حلقة بيتا-لاكتام (β-Lactam Ring).",
        "handling_ar": "حساس للرطوبة والحرارة والأحماض التي تسبب فتح حلقة بيتا-لاكتام وفقدان الفعالية.",
    },
    "omeprazole": {
        "hazard_ar": "مهيج خفيف للعينين والجلد؛ حساس جداً للأوساط الحمضية.",
        "ppe_ar": "معطف المختبر، قفازات النتريل، ونظارات الوقاية.",
        "quality_test_ar": "اختبار الثبات في الأوساط المنظمة (Buffer pH 7.4 vs pH 2.0) والتحليل بـ HPLC.",
        "handling_ar": "يتحلل سريعاً في الأوساط الحمضية؛ لذلك يُحضر صيدلانياً في كبسولات مغلفة معوياً (Enteric-coated).",
    },
}


def build_pharmacology_section(entity: dict[str, Any]) -> dict[str, Any]:
    pharm = dict(entity.get("pharmacology_details", {}))
    smiles = entity.get("canonical_smiles", "")
    rdkit_metrics = compute_rdkit_molecular_metrics(smiles) if smiles else None
    key = entity.get("key", "")
    prelab = PHARMACY_PRELAB_NOTES.get(
        key,
        {
            "hazard_ar": "مادة دوائية فعالة حيوياً (Active Pharmaceutical Ingredient — API): يلزم تجنب الاستنشاق أو التلامس المباشر.",
            "ppe_ar": "معطف المختبر (Lab Coat)، قفازات النتريل (Nitrile Gloves)، ونظارات الوقاية (Safety Goggles).",
            "quality_test_ar": "التحقق من الهوية والنقاوة عبر مطابقة درجة الانصهار والتحليل الطيفي والكروماتوغرافي (USP / BP Pharmacopoeia).",
            "handling_ar": "تُوزن بدقة في بيئة جافة ونظيفة مع الالتزام ببروتوكول السلامة المخبرية لطلاب الصيدلة.",
        },
    )

    return {
        "therapeutic_class_ar": pharm.get("therapeutic_class_ar", "مادة دوائية نشطة حيوياً (Active Drug Substance)"),
        "indication_ar": pharm.get("indication_ar", "يستخدم وفق الاستطبابات الدوائية المعتمدة"),
        "target_protein_ar": pharm.get("target_protein_ar", "هدف بروتيني جزيئي (Molecular Target)"),
        "mechanism_of_action_ar": pharm.get("mechanism_of_action_ar", "يرتبط بالهدف الجزيئي لتعديل المسار الفسيولوجي."),
        "documented_bioactivity_ar": pharm.get("documented_bioactivity_ar", "موثق في قاعدة بيانات ChEMBL"),
        "approval_status_ar": pharm.get("approval_status_ar", "مسجل في قواعد البيانات الدوائية المعتمدة"),
        "pdb_complex_id": pharm.get("pdb_complex_id"),
        "calculated_drug_descriptors": rdkit_metrics,
        "descriptors_available": rdkit_metrics is not None,
        "descriptors_status_ar": (
            "تم حساب الخصائص الجزيئية الدوائية حتمياً من البنية الكيميائية عبر RDKit."
            if rdkit_metrics is not None
            else "الخصائص الجزيئية الدوائية المحسوبة (MW, LogP, HBD, HBA, TPSA) غير متاحة حالياً لعدم توفر محرك التحليل البنيوي (RDKit)؛ يلتزم النظام بعدم عرض قيم افتراضية."
        ),
        "prelab_safety_and_qc": prelab,
        "retrieval_provenance": build_provenance_tag(
            OP_SCIENTIFIC_RETRIEVAL,
            sources=["ChEMBL", "PubChem"],
            notes_ar="البيانات الدوائية والنشاط الحيوي مسترجعة من السجلات المرجعية الموثقة وليست تنبؤاً حاسوبياً.",
        ),
        "calculation_provenance": build_provenance_tag(
            OP_SCIENTIFIC_CALCULATION,
            sources=["RDKit Chemistry Engine"] if rdkit_metrics is not None else [],
            notes_ar=(
                "الخصائص الجزيئية الدوائية (MW, LogP, HBD, HBA, TPSA, Lipinski Rule of Five) محسوبة حتمياً من البنية الكيميائية."
                if rdkit_metrics is not None
                else "تعذر حساب الخصائص الدوائية لعدم توفر محرك RDKit؛ لم يتم توليد أي قيم افتراضية."
            ),
        ),
    }


PHARMACOLOGY_EDUCATIONAL_DISCLAIMER_AR = (
    "فصل صارم بين المعلومات التعليمية والإرشادات الطبية: جميع المحتويات الدوائية في منصة MolecuLens "
    "مخصصة للتعليم الأكاديمي لطلاب الصيدلة والعلوم الصحية فقط، ولا تُعد استشارة طبية أو وصفة علاجية أو بديلاً عن المراجع السريرية الرسمية."
)

PHARMACOLOGY_PILLARS: list[dict[str, str]] = [
    {
        "code": "classification",
        "title_ar": "1. تصنيف الأدوية والمجموعات العلاجية",
        "title_en": "1. Drug Classification & Therapeutic Classes",
        "desc_ar": (
            "تصنيف منهجي للأدوية وفق نظام ATC العالمي والمجموعات العلاجية والدوائية "
            "(المسكنات، المضادات الحيوية، أدوية القلب والضغط، السكري، الجهاز الهضمي، والتنفسي)."
        ),
        "desc_en": (
            "Systematic organization of medications by WHO ATC classification and therapeutic families "
            "(analgesics, antibiotics, cardiovascular, antidiabetics, GI, and respiratory agents)."
        ),
        "planned_sources": "WHO ATC Index, RxNorm, ChEMBL",
    },
    {
        "code": "moa",
        "title_ar": "2. آليات عمل الأدوية (Mechanism of Action)",
        "title_en": "2. Mechanism of Action (MoA)",
        "desc_ar": (
            "تفسير جزيئي وخلوي لكيفية إحداث الدواء لتأثيره العلاجي: تثبيط الإنزيمات، "
            "تنشيط أو حصر المستقبلات، تعديل القنوات الأيونية، أو التداخل مع المسارات الأيضية."
        ),
        "desc_en": (
            "Molecular and cellular explanation of therapeutic action: enzyme inhibition, "
            "receptor agonism/antagonism, ion channel modulation, and metabolic pathway targeting."
        ),
        "planned_sources": "ChEMBL Mechanism API, IUPHAR/BPS Guide to Pharmacology",
    },
    {
        "code": "targets",
        "title_ar": "3. المستقبلات والأهداف الحيوية (Receptors & Biological Targets)",
        "title_en": "3. Receptors & Biological Targets",
        "desc_ar": (
            "ربط المادة الدوائية بالبروتينات والمستقبلات المستهدفة (GPCRs، القنوات الأيونية، "
            "الإنزيمات، والمستقبلات النووية) مع معرفات UniProt و PDB."
        ),
        "desc_en": (
            "Linking active pharmaceutical ingredients to biological macromolecular targets "
            "(GPCRs, ion channels, enzymes, nuclear receptors) with UniProt and PDB cross-references."
        ),
        "planned_sources": "UniProtKB, RCSB PDB, Open Targets",
    },
    {
        "code": "adme",
        "title_ar": "4. الحرائك الدوائية (Pharmacokinetics — ADME)",
        "title_en": "4. Pharmacokinetics (ADME Profile)",
        "desc_ar": (
            "دراسة رحلة الدواء داخل الجسم: الامتصاص (Absorption) والتوافر الحيوي، التوزيع (Distribution) "
            "والارتباط ببروتينات البلازما، الاستقلاب الكبدي (Metabolism / CYP450)، والإطراح (Excretion) وعمر النصف."
        ),
        "desc_en": (
            "Quantitative profiling of Absorption (bioavailability), Distribution (volume of distribution, "
            "plasma protein binding), Metabolism (CYP450 isoenzymes), and Excretion (clearance, half-life)."
        ),
        "planned_sources": "FDA Clinical Pharmacology Monographs, ChEMBL ADMET",
    },
    {
        "code": "dose_response",
        "title_ar": "5. العلاقة بين الجرعة والاستجابة (Dose–Response Relationship)",
        "title_en": "5. Dose–Response Relationship & Therapeutic Index",
        "desc_ar": (
            "نمذجة منحنيات التركيز والاستجابة (EC50, IC50, Emax)، الفعالية القصوى (Efficacy) مقابل القوة الدوائية (Potency)، "
            "والنافذة العلاجية (Therapeutic Window) ومؤشر الأمان الدوائي."
        ),
        "desc_en": (
            "Modeling graded and quantal dose-response curves (EC50, IC50, Emax), potency vs. efficacy, "
            "and therapeutic index safety margins."
        ),
        "planned_sources": "Goodman & Gilman's Pharmacological Basis of Therapeutics",
    },
    {
        "code": "adverse_effects",
        "title_ar": "6. الآثار الجانبية وموانع الاستخدام (Adverse Effects & Contraindications)",
        "title_en": "6. Adverse Effects & Contraindications",
        "desc_ar": (
            "توثيق الآثار الجانبية الشائعة والخطيرة، التحذيرات المؤطرة (Boxed Warnings)، "
            "موانع الاستخدام، والاحتياطات في القصور الكلوي أو الكبدي والحمل والرضاعة."
        ),
        "desc_en": (
            "Structured catalog of adverse drug reactions, boxed warnings, contraindications, "
            "and special population precautions (renal/hepatic impairment, pregnancy)."
        ),
        "planned_sources": "openFDA Drug Labeling API, DailyMed",
    },
    {
        "code": "references",
        "title_ar": "7. المصادر الدوائية الموثوقة وفصل التعليم عن الإرشاد الطبي (Trusted Sources & Educational Scope)",
        "title_en": "7. Drug Monographs, Trusted References & Educational Scope",
        "desc_ar": (
            "ربط كل معلومة دوائية بمراجعها الدستورية والسريرية المعتمدة (FDA، EMA، BNF، USP-NF) "
            "مع فصل صارم بين المحتوى التعليمي الأكاديمي والإرشادات الطبية السريرية."
        ),
        "desc_en": (
            "Direct linkage of each drug monograph to authoritative regulatory and pharmacopeial "
            "references (FDA, EMA, BNF, USP-NF) with strict separation between education and clinical advice."
        ),
        "planned_sources": "openFDA (CC0), RxNorm (NLM), PubChem (Public Domain)",
    },
]

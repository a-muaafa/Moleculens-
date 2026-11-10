from functools import lru_cache
from rdkit import Chem
from rdkit.Chem import Descriptors, Lipinski, QED, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from chemistry.molecular_structure import get_molecular_formula, get_murcko_scaffold, validate_smiles


@lru_cache(maxsize=1)
def get_pains_catalog():
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    return FilterCatalog(params)


PROPERTY_EDUCATIONAL_GUIDE = {
    "molecular_weight": {
        "title_ar": "الوزن الجزيئي (Molecular Weight)",
        "unit": "Da (g/mol)",
        "reference_range": "<= 500 Da",
        "meaning": "مجموع الكتل الذرية لجميع الذرات المكونة للجزيء مقاسا بوحدة الدالتون (Da).",
        "calculation": "يحسب عبر جمع الأوزان الذرية القياسية لكل عنصر في الصيغة الجزيئية باستخدام مكتبة RDKit.",
        "importance": "تتميز الجزيئات ذات الوزن الأقل من 500 دالتون بقدرة أعلى على الانتشار السلبي عبر الأغشية الخلوية وجدار الأمعاء.",
        "limitations": "الوزن الجزيئي المناسب لا يضمن ارتباط المركب بالهدف البروتيني، كما أن بعض الأدوية الفعالة تتجاوز هذا الحد إذا كانت تنتقل عبر نواقل حيوية.",
    },
    "logp": {
        "title_ar": "معامل التوزع الدهني المائي (LogP)",
        "unit": "بدون وحدة (لوغاريتمي)",
        "reference_range": "<= 5.0",
        "meaning": "مقياس لدرجة محبة الجزيء للدهون (Lipophilicity) مقارنة بالماء.",
        "calculation": "يحسب حاسوبيا بطريقة Wildman-Crippen (cLogP) بجمع مساهمات المجموعات الذرية في الجزيء.",
        "importance": "يتحكم في توازن ذوبانية المركب في السوائل الحيوية المائية وقدرته على اختراق الطبقات الدهنية الثنائية للأغشية الخلوية.",
        "limitations": "القيمة المحسوبة تقديرية وقد تختلف عن القياس المخبري (LogD عند درجة حموضة فسيولوجية) خصوصا للمركبات المتأينة.",
    },
    "hbd": {
        "title_ar": "مانحات روابط الهيدروجين (Hydrogen Bond Donors - HBD)",
        "unit": "عدد المجموعات",
        "reference_range": "<= 5",
        "meaning": "عدد ذرات الهيدروجين المرتبطة بذرات ذات كهرسلبية عالية مثل الأكسجين (OH) أو النيتروجين (NH).",
        "calculation": "يحصى عبر عد مجموعات NH و OH القادرة على منح ذرة هيدروجين في البنية الجزيئية.",
        "importance": "تسهم في تثبيت المركب داخل موقع الارتباط في البروتين، لكن زيادتها المفرطة تعيق عبور الغشاء الدهني بسبب طاقة التميؤ العالية.",
        "limitations": "العدد الإجمالي لا يوضح ما إذا كانت هذه المجموعات تشكل روابط هيدروجينية داخلية (Intramolecular) تحسن النفاذية.",
    },
    "hba": {
        "title_ar": "مستقبلات روابط الهيدروجين (Hydrogen Bond Acceptors - HBA)",
        "unit": "عدد الذرات",
        "reference_range": "<= 10",
        "meaning": "عدد الذرات الغنية بالإلكترونات (خصوصا النيتروجين والأكسجين) القادرة على استقبال رابطة هيدروجينية.",
        "calculation": "يحسب في قاعدة Lipinski بعد ذرات النيتروجين والأكسجين في الجزيء.",
        "importance": "ضرورية للذوبان المائي والتفاعل الانتقائي مع الأحماض الأمينية في الهدف البروتيني.",
        "limitations": "لا تتساوى جميع ذرات النيتروجين والأكسجين في قوة استقبالها للروابط الهيدروجينية.",
    },
    "rotatable_bonds": {
        "title_ar": "الروابط القابلة للدوران (Rotatable Bonds)",
        "unit": "عدد الروابط",
        "reference_range": "<= 10",
        "meaning": "عدد الروابط الأحادية غير الحلقية وغير الطرفية التي تسمح بدوران حر لأجزاء الجزيء.",
        "calculation": "تحصى عبر تحديد الروابط الأحادية بين الذرات الثقيلة غير الموجودة في حلقات وغير المرتبطة بذرات طرفية (مع استثناء روابط الأميد).",
        "importance": "تعبر عن المرونة الفراغية للجزيء؛ انخفاضها يقلل خسارة الإنتروبيا الحرارية عند ارتباط المركب بموقع الإنزيم.",
        "limitations": "بعض المركبات المرنة ترتبط بفعالية عالية إذا كانت تشكل تفاعلات أيونية قوية تعوض فقدان الإنتروبيا.",
    },
    "tpsa": {
        "title_ar": "المساحة السطحية القطبية (TPSA)",
        "unit": "Å²",
        "reference_range": "< 140 Å²",
        "meaning": "مجموع المساحات السطحية للذرات القطبية (الأكسجين والنيتروجين والهيدروجين المرتبط بهما) في الجزيء.",
        "calculation": "تحسب بطريقة Ertl السريعة بالاعتماد على مساهمات الأجزاء القطبية الطوبولوجية دون الحاجة لتوليد شكل ثلاثي الأبعاد.",
        "importance": "مؤشر رئيسي لتقدير امتصاص الأمعاء (< 140 Å²) والقدرة على عبور الحاجز الدموي الدماغي (< 90 Å²).",
        "limitations": "تعتمد على البنية ثنائية الأبعاد ولا تراعي انطواء الجزيء في الفراغ ثلاثي الأبعاد.",
    },
}


def calculate_descriptors(mol_or_smiles) -> dict:
    if isinstance(mol_or_smiles, str):
        mol, canonical_smi, err = validate_smiles(mol_or_smiles)
        if mol is None:
            raise ValueError(err)
    else:
        mol = mol_or_smiles
        if mol is None:
            raise ValueError("كائن الجزيء غير صالح.")
        canonical_smi = Chem.MolToSmiles(mol, canonical=True)

    mw = float(Descriptors.MolWt(mol))
    logp = float(Descriptors.MolLogP(mol))
    hbd = int(Lipinski.NumHDonors(mol))
    hba = int(Lipinski.NumHAcceptors(mol))
    rot_bonds = int(Lipinski.NumRotatableBonds(mol))
    tpsa = float(Descriptors.TPSA(mol))
    heavy_atoms = int(mol.GetNumHeavyAtoms())
    ring_count = int(rdMolDescriptors.CalcNumRings(mol))
    qed_val = float(QED.qed(mol))
    formula = get_molecular_formula(mol)
    _, scaffold_smi = get_murcko_scaffold(mol)

    lipinski_checks = {
        "mw_ok": mw <= 500.0,
        "logp_ok": logp <= 5.0,
        "hbd_ok": hbd <= 5,
        "hba_ok": hba <= 10,
    }
    violations = int(sum(not ok for ok in lipinski_checks.values()))
    passes_lipinski = violations <= 1

    catalog = get_pains_catalog()
    match = catalog.GetFirstMatch(mol)
    has_pains = match is not None
    pains_desc = match.GetDescription() if has_pains else None
    pains_status = (
        f"تم رصد تنبيه هيكلي محتمل للتداخل المخبري ({pains_desc})"
        if has_pains
        else "لم ترصد تنبيهات PAINS المعروفة"
    )

    return {
        "canonical_smiles": canonical_smi,
        "molecular_formula": formula,
        "murcko_scaffold": scaffold_smi,
        "molecular_weight": round(mw, 2),
        "logp": round(logp, 2),
        "hbd": hbd,
        "hba": hba,
        "rotatable_bonds": rot_bonds,
        "tpsa": round(tpsa, 2),
        "heavy_atoms": heavy_atoms,
        "ring_count": ring_count,
        "qed": round(qed_val, 3),
        "lipinski": {
            "checks": lipinski_checks,
            "violations_count": violations,
            "passes_rule_of_five": passes_lipinski,
            "summary_ar": (
                f"عدد مخالفات قاعدة Lipinski: {violations} من 4 "
                f"({'ضمن النطاق المعتاد للامتصاص الفموي' if passes_lipinski else 'خارج النطاق المعتاد للامتصاص الفموي'})."
            ),
            "disclaimer_ar": (
                "مطابقة قاعدة Lipinski تشير فقط إلى ملاءمة الخصائص الفيزيائية للامتصاص الفموي، "
                "ولا تعد دليلا على فعالية المركب تجاه الهدف البروتيني أو سلامته الدوائية."
            ),
        },
        "pains": {
            "has_pains_alert": has_pains,
            "alert_description": pains_desc,
            "status_ar": pains_status,
            "disclaimer_ar": (
                "عدم رصد تنبيهات PAINS يعني خلو البنية من الأنماط المسجلة كمسببات شائعة للنتائج الإيجابية الكاذبة في الفحوصات المخبرية، "
                "ولا يعني أن المركب آمن حيويا أو فعال علاجيا."
            ),
        },
    }


from chemistry.descriptors import calculate_descriptors
from chemistry.molecular_structure import validate_smiles
from config.settings import (
    AD_BORDERLINE_TANIMOTO,
    AD_IN_DOMAIN_TANIMOTO,
    AD_MAX_MW,
    AD_MIN_HEAVY_ATOMS,
)


def classify_applicability_domain(
    smiles_or_mol,
    max_tanimoto: float,
    mean_top5_tanimoto: float = 0.0,
    train_bounds: dict | None = None,
) -> dict:
    if isinstance(smiles_or_mol, str):
        mol, canon_smi, err = validate_smiles(smiles_or_mol)
        if mol is None:
            return {
                "status_code": "Out-of-Domain",
                "status_ar": "خارج النطاق – Out-of-Domain",
                "pill_class": "pill-low",
                "can_show_reliable_prediction": False,
                "max_tanimoto": 0.0,
                "mean_top5_tanimoto": 0.0,
                "reasons_ar": [err or "بنية كيميائية غير صالحة."],
                "thresholds_justification_ar": (
                    "يتطلب دخول نطاق التطبيق بنية جزيئية صحيحة وتشابها هيكليا كافيا مع مركبات التدريب."
                ),
                "scientific_disclaimer_ar": (
                    "هذا التصنيف يعبر عن مدى انتماء المركب إلى نطاق تطبيق النموذج، وليس احتمال صحة التنبؤ أو مجال ثقة إحصائيا."
                ),
            }
    else:
        mol = smiles_or_mol

    desc = calculate_descriptors(mol)
    heavy_atoms = desc["heavy_atoms"]
    mw = desc["molecular_weight"]
    logp = desc["logp"]
    tpsa = desc["tpsa"]

    bounds = train_bounds or {
        "mw_min": 140.0,
        "mw_max": 750.0,
        "logp_min": -2.5,
        "logp_max": 6.5,
        "tpsa_min": 15.0,
        "tpsa_max": 190.0,
    }

    structural_pass = (heavy_atoms >= AD_MIN_HEAVY_ATOMS) and (mw <= AD_MAX_MW)
    props_within_train = (
        bounds["mw_min"] <= mw <= bounds["mw_max"]
        and bounds["logp_min"] <= logp <= bounds["logp_max"]
        and bounds["tpsa_min"] <= tpsa <= bounds["tpsa_max"]
    )

    reasons = []
    if not structural_pass:
        if heavy_atoms < AD_MIN_HEAVY_ATOMS:
            reasons.append(
                f"عدد الذرات الثقيلة ({heavy_atoms}) أقل من الحد الأدنى لجزيئات مثبطات DPP4 العضوية ({AD_MIN_HEAVY_ATOMS} ذرات)."
            )
        if mw > AD_MAX_MW:
            reasons.append(f"الوزن الجزيئي ({mw:.1f} Da) يتجاوز الحد الأعلى لنطاق مركبات التدريب ({AD_MAX_MW:.0f} Da).")

    if max_tanimoto < AD_BORDERLINE_TANIMOTO:
        reasons.append(
            f"أقصى تشابه Tanimoto مع مركبات التدريب ({max_tanimoto:.3f}) أقل من عتبة الحد الأدنى لنطاق التطبيق ({AD_BORDERLINE_TANIMOTO:.2f})."
        )
    elif max_tanimoto < AD_IN_DOMAIN_TANIMOTO:
        reasons.append(
            f"أقصى تشابه Tanimoto مع مركبات التدريب ({max_tanimoto:.3f}) يقع في المنطقة الحدية بين {AD_BORDERLINE_TANIMOTO:.2f} و {AD_IN_DOMAIN_TANIMOTO:.2f}."
        )
    else:
        reasons.append(
            f"أقصى تشابه Tanimoto مع مركبات التدريب ({max_tanimoto:.3f}) يحقق عتبة التمثيل الهيكلي الجيد (>= {AD_IN_DOMAIN_TANIMOTO:.2f})."
        )

    if not props_within_train and structural_pass:
        reasons.append("بعض الخصائص الفيزيائية الكيميائية (الوزن الجزيئي أو LogP أو TPSA) تقع على أطراف توزيع بيانات التدريب.")

    if (not structural_pass) or (max_tanimoto < AD_BORDERLINE_TANIMOTO):
        status_code = "Out-of-Domain"
        status_ar = "خارج النطاق – Out-of-Domain"
        pill_class = "pill-low"
        can_show_reliable = False
        summary_ar = (
            "المركب مختلف بدرجة كبيرة عن المركبات التي تدرب عليها النموذج أو لا يستوفي الشروط البنائية لنطاق التطبيق. "
            "لا يجوز اعتماد تنبؤ النموذج لهذا المركب كنتيجة موثوقة."
        )
    elif (max_tanimoto >= AD_IN_DOMAIN_TANIMOTO) and props_within_train:
        status_code = "In-Domain"
        status_ar = "داخل النطاق – In-Domain"
        pill_class = "pill-high"
        can_show_reliable = True
        summary_ar = (
            "المركب قريب بنيويا وفيزيائيا من نطاق المركبات التي تدرب عليها النموذج. "
            "يمكن الاسترشاد بالتنبؤ مع التأكيد على أنه توقع حاسوبي أولي وليس قياسا مخبريا."
        )
    else:
        status_code = "Borderline"
        status_ar = "حدّي – Borderline"
        pill_class = "pill-mid"
        can_show_reliable = True
        summary_ar = (
            "المركب قريب جزئيا من نطاق التدريب، لكن توجد مؤشرات على ضعف التمثيل الهيكلي المباشر أو ارتفاع عدم اليقين. "
            "يعرض التنبؤ مع تحذير واضح بضرورة التحقق المخبري."
        )

    return {
        "status_code": status_code,
        "status_ar": status_ar,
        "pill_class": pill_class,
        "can_show_reliable_prediction": can_show_reliable,
        "max_tanimoto": round(float(max_tanimoto), 3),
        "mean_top5_tanimoto": round(float(mean_top5_tanimoto), 3),
        "structural_pass": structural_pass,
        "props_within_train": props_within_train,
        "summary_ar": summary_ar,
        "reasons_ar": reasons,
        "thresholds_justification_ar": (
            "تستند عتبات بصمات Morgan (نصف القطر 2، 1024 بت) في الكيمياء الحاسوبية إلى أن Tanimoto >= 0.60 يعكس تقاربا "
            "في النواة الكيميائية والمجموعات الوظيفية، بينما تشير القيم بين 0.30 و 0.59 إلى تشابه جزئي، وتدل القيم الأقل من 0.30 "
            "على غياب نظائر هيكلية قريبة في بيانات التدريب."
        ),
        "scientific_disclaimer_ar": (
            "هذا التصنيف يعبر عن مدى انتماء المركب إلى نطاق تطبيق النموذج (Applicability Domain)، "
            "وليس احتمال صحة التنبؤ أو مجال ثقة إحصائيا."
        ),
    }


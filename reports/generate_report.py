from config.settings import APP_NAME, APP_SLOGAN, TARGET_INFO, TEAM_NAME


def generate_educational_report(
    compound_info: dict,
    descriptors: dict,
    published_records: list[dict],
    prediction: dict,
    metrics: dict,
    student_guess: dict | None = None,
) -> str:
    ad = prediction["applicability_domain"]
    nearest_df = prediction["nearest_train_compounds"]

    lines = [
        f"# التقرير التعليمي الشامل - منصة {APP_NAME}",
        f"**الفريق:** {TEAM_NAME} | **الشعار:** {APP_SLOGAN}",
        "**حوكمة المشروع:** تقرير تعليمي للفرز الحاسوبي الأولي؛ لا يعد أداة لاكتشاف أدوية جديدة أو إثبات فعالية علاجية أو تقديم استشارة طبية.",
        "",
        "---",
        "## أولا: هوية المركب والبنية الجزيئية والهدف البروتيني",
        f"- **الاسم العربي:** {compound_info.get('name_ar', 'غير محدد')}",
        f"- **الاسم الإنجليزي:** {compound_info.get('name_en', 'غير محدد')}",
        f"- **الاسم العلمي (IUPAC):** `{compound_info.get('iupac_name', 'غير محدد')}`",
        f"- **الصيغة الجزيئية:** `{descriptors['molecular_formula']}`",
        f"- **معرف PubChem CID:** `{compound_info.get('pubchem_cid') or 'غير متوفر'}`",
        f"- **معرف ChEMBL ID:** `{compound_info.get('chembl_id') or 'غير متوفر'}`",
        f"- **صيغة SMILES القياسية:** `{descriptors['canonical_smiles']}`",
        f"- **الهيكل الأساسي (Bemis-Murcko Scaffold):** `{descriptors['murcko_scaffold'] or 'لا يوجد نظام حلقي'}`",
        f"- **مصدر بيانات الهوية:** {compound_info.get('source', 'القاموس المحلي / RDKit')}",
        f"- **الهدف البروتيني المدروس:** {TARGET_INFO['name_en']} ({TARGET_INFO['name_ar']}) | UniProt: `{TARGET_INFO['uniprot_id']}` | ChEMBL: `{TARGET_INFO['target_id']}` | الكائن الحي: `{TARGET_INFO['organism']}`",
        "",
        "---",
        "## ثانيا: الخصائص الكيميائية المحسوبة (RDKit - منفصلة عن القياسات والتنبؤات)",
        f"- **الوزن الجزيئي (Molecular Weight):** `{descriptors['molecular_weight']} Da`",
        f"- **معامل التوزع الدهني (LogP):** `{descriptors['logp']}`",
        f"- **مانحات روابط الهيدروجين (HBD):** `{descriptors['hbd']}`",
        f"- **مستقبلات روابط الهيدروجين (HBA):** `{descriptors['hba']}`",
        f"- **الروابط القابلة للدوران (Rotatable Bonds):** `{descriptors['rotatable_bonds']}`",
        f"- **المساحة السطحية القطبية (TPSA):** `{descriptors['tpsa']} Å²`",
        f"- **تقييم قاعدة Lipinski:** {descriptors['lipinski']['summary_ar']}",
        f"  - *ملاحظة علمية:* {descriptors['lipinski']['disclaimer_ar']}",
        f"- **فحص تنبيهات PAINS:** {descriptors['pains']['status_ar']}",
        f"  - *ملاحظة علمية:* {descriptors['pains']['disclaimer_ar']}",
        "",
        "---",
        "## ثالثا: الأدلة العلمية المنشورة من التجارب المخبرية (ChEMBL Experimental Evidence)",
    ]

    if published_records:
        for idx, rec in enumerate(published_records, start=1):
            lines.extend(
                [
                    f"### سجل تجريبي منشور ({idx}):",
                    f"- **معرف المركب:** `{rec.get('compound_id') or rec.get('molecule_chembl_id')}`",
                    f"- **الهدف البروتيني والكائن الحي:** `{rec.get('target_id', 'CHEMBL284')}` (`{rec.get('target_organism', 'Homo sapiens')}`)",
                    f"- **نوع القياس:** `{rec.get('standard_type', 'IC50')}` (العلاقة: `{rec.get('standard_relation', '=')}`)",
                    f"- **القيمة الأصلية والوحدة الأصلية:** `{rec.get('original_value', rec.get('standard_value'))} {rec.get('original_units', rec.get('standard_units', 'nM'))}`",
                    f"- **عدد القياسات المتوافقة المجمعة:** `{rec.get('measurement_count', 1)}`",
                    f"- **وسيط النشاط الموحد (Median Activity):** `{float(rec.get('median_activity', rec.get('standard_value'))):.2f} nM` (`pIC50 = {float(rec.get('pIC50', 0.0)):.2f}`)",
                    f"- **مراجع التجارب والوثائق الأصلية (Assay & Document IDs):** `{rec.get('source_ids') or (str(rec.get('assay_id')) + ' / ' + str(rec.get('document_id')))}`",
                ]
            )
    else:
        lines.append("- لا يوجد سجل تجريبي مطابق تماما لهذا المركب في مجموعة بيانات ChEMBL المحلية المنظفة لهدف DPP4 البشري.")

    lines.extend(
        [
            "",
            "---",
            "## رابعا: نتائج التنبؤ الحاسوبي لنموذج الذكاء الاصطناعي وتصنيف نطاق التطبيق",
            "- **تنبيه الفصل المنهجي:** القيم التالية هي مخرجات حاسوبية لنموذج Random Forest وليست نتائج تجارب مخبرية.",
            f"- **نوع القياس المتنبأ به:** {prediction['measurement_type']}",
            f"- **وحدة القياس:** {prediction['units']}",
            f"- **تصنيف نطاق تطبيق النموذج (Applicability Domain):** **{ad['status_ar']}**",
            f"- **أقصى تشابه Tanimoto مع مركبات التدريب:** `{prediction['max_tanimoto']:.3f}` (متوسط أقرب 5 مركبات: `{prediction['mean_top5_tanimoto']:.3f}`)",
        ]
    )

    if ad["can_show_reliable_prediction"]:
        lines.extend(
            [
                f"- **قيمة التنبؤ اللوغاريتمي (Predicted pIC50):** `{prediction['pred_pic50']:.2f}`",
                f"- **قيمة تركيز التثبيط النصفي المكافئة (Predicted IC50):** `{prediction['pred_ic50_nm']:,.2f} nM`",
                f"- **التفسير الاسترشادي للنشاط:** {prediction['category_ar']} — {prediction['explanation_ar']}",
            ]
        )
    else:
        lines.extend(
            [
                "- **حالة عرض التنبؤ العددي:** محجوب كقيمة معتمدة لأن المركب يقع **خارج نطاق تطبيق النموذج (Out-of-Domain)**.",
                f"  - *(القيمة الحاسوبية الخام لأغراض المقارنة التعليمية فقط مع التحذير: pIC50 = `{prediction['pred_pic50']:.2f}` / IC50 = `{prediction['pred_ic50_nm']:,.2f} nM` — غير موثوقة بسبب غياب التغطية الكيميائية في بيانات التدريب)*.",
            ]
        )

    lines.append("- **أسباب ومبررات تصنيف النطاق:**")
    for r in ad["reasons_ar"]:
        lines.append(f"  - {r}")
    lines.append(f"- **إقرار منهجي:** {ad['scientific_disclaimer_ar']}")

    lines.extend(
        [
            "",
            "### مؤشرات أداء النموذج على مجموعة الاختبار (Bemis-Murcko Scaffold Split):",
            f"- **MAE (pIC50):** `{metrics['mae_pic50']}`",
            f"- **RMSE (pIC50):** `{metrics['rmse_pic50']}`",
            f"- **R² (التنبؤ العددي):** `{metrics['r2_pic50']}`",
            f"- **ROC-AUC (للتمييز التصنيفي عند عتبة 1000 nM):** `{metrics['roc_auc']}`",
            f"  - *توضيح:* {metrics.get('roc_auc_note_ar', '')}",
            "",
            "### أقرب 5 مركبات مرجعية في بيانات التدريب فقط:",
        ]
    )

    for _, r_near in nearest_df.iterrows():
        lines.append(
            f"- `{r_near['molecule_chembl_id']}` | Tanimoto: `{r_near['tanimoto']:.3f}` | "
            f"IC50: `{r_near['standard_value']:.2f} nM` (`pIC50 = {r_near['pIC50']:.2f}`) | "
            f"Assay: `{r_near['assay_chembl_id']}` | Document: `{r_near['document_chembl_id']}`"
        )

    if student_guess and student_guess.get("submitted"):
        lines.extend(
            [
                "",
                "---",
                "## خامسا: مقارنة تخمين الطالب (Guess Before You See) بنتيجة النموذج",
                f"- **تخمين الطالب للنشاط الحيوي:** {student_guess.get('activity_guess')}",
                f"- **تخمين الطالب لنطاق التطبيق:** {student_guess.get('domain_guess')}",
                f"- **نتيجة النموذج الفعلية:** {prediction['category_ar']} | النطاق: {ad['status_ar']}",
                f"- **التفسير التعليمي للمقارنة:** {student_guess.get('comparison_explanation', '')}",
            ]
        )

    lines.extend(
        [
            "",
            "---",
            "## سادسا: القيود العلمية ومصادر البيانات",
            "- يقتصر النموذج على قياسات الربط المخبري المباشر (IC50) لإنزيم DPP4 البشري ولا يجوز خلطها بقياسات Ki دون معادلة تحويل حركية موثقة.",
            "- التنبؤ الحاسوبي لا يثبت الفعالية الدوائية أو الانتقائية مقابل إنزيمي DPP8 و DPP9 أو السلامة السريرية.",
            "- **المصادر المرجعية:** UniProt (`P27487`), ChEMBL (`CHEMBL284`), PubChem PUG REST API, RDKit Open-Source Cheminformatics.",
        ]
    )

    return "\n".join(lines)


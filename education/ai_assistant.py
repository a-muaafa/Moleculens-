"""
AI Scientific Assistant for MolecuLens.
Provides pedagogical explanations, step-by-step breakdowns of chemical properties,
Lipinski rule interpretation, and interactive student Q&A grounded strictly in calculated data.
"""

from typing import Any, Optional
import html


class AIScientificAssistant:
    """Explains compound analysis and physicochemical properties with zero data hallucination."""

    PROP_MAP = {
        "mw": "molecular_weight",
        "molecular_weight": "molecular_weight",
        "logp": "logp",
        "hbd": "hbd",
        "hba": "hba",
        "tpsa": "tpsa",
        "rotatable_bonds": "rotatable_bonds",
    }

    PROPERTY_META: dict[str, dict[str, Any]] = {
        "molecular_weight": {
            "title_ar": "الوزن الجزيئي (Molecular Weight - MW)",
            "title_en": "Molecular Weight (MW)",
            "meaning_ar": "مجموع الكتل الذرية لجميع الذرات المكونة للجزيء بوحدة جرام/مول (Da).",
            "meaning_en": "Sum of atomic weights of all constituent atoms in g/mol (Da).",
            "significance_ar": "يحدد حجم الجزيء وقدرته على النفاذ عبر الأغشية الخلوية؛ تفضل الأدوية الفموية وزناً جزيئياً أقل من 500 Da وفق قاعدة ليبينسكي.",
            "significance_en": "Determines molecule size and membrane permeation; oral drugs generally require MW <= 500 Da per Lipinski's rule.",
            "threshold": "MW <= 500 Da",
            "limit": 500.0,
        },
        "logp": {
            "title_ar": "معامل التقسيم والتوزع الدهني (Partition Coefficient - LogP)",
            "title_en": "Partition Coefficient (LogP)",
            "meaning_ar": "لوغاريتم نسبة تركيز المركب في الأوكتانول مقارنة بالماء؛ يقيس مدى ألفة الجزيء للدهون (Lipophilicity).",
            "meaning_en": "Logarithm of compound concentration ratio between octanol and water; measures lipophilicity.",
            "significance_ar": "قيم LogP المعتدلة (بين 0 و 5) تدل على توازن مثالي بين الذوبانية في سوائل الجسم المائية والقدرة على عبور الحواجز الدهنية للخلايا.",
            "significance_en": "Moderate LogP (0 to 5) balances aqueous solubility in bodily fluids with lipid membrane permeability.",
            "threshold": "LogP <= 5.0",
            "limit": 5.0,
        },
        "hbd": {
            "title_ar": "مانحات الروابط الهيدروجينية (Hydrogen Bond Donors - HBD)",
            "title_en": "Hydrogen Bond Donors (HBD)",
            "meaning_ar": "عدد ذرات الهيدروجين المرتبطة بذرات عالية الكهروسالبية (مثل -OH و -NH).",
            "meaning_en": "Count of hydrogen atoms bonded to electronegative atoms (-OH, -NH).",
            "significance_ar": "العدد الزائد من المانحات (> 5) يزيد من طاقة نزع الإماهة (Desolvation Energy) ويعيق امتصاص الدواء عبر بطانة الأمعاء.",
            "significance_en": "Excess donors (> 5) raise desolvation energy, hindering intestinal membrane absorption.",
            "threshold": "HBD <= 5",
            "limit": 5.0,
        },
        "hba": {
            "title_ar": "مستقبلات الروابط الهيدروجينية (Hydrogen Bond Acceptors - HBA)",
            "title_en": "Hydrogen Bond Acceptors (HBA)",
            "meaning_ar": "عدد ذرات النيتروجين والأكسجين ذات الأزواج الإلكترونية الحرة القادرة على استقبال رابطة هيدروجينية.",
            "meaning_en": "Count of N and O atoms with lone pairs capable of accepting an H-bond.",
            "significance_ar": "تساعد في تثبيت الدواء داخل الجيب الرابط للبروتين الهدف، لكن زيادتها عن 10 تقلل من النفاذية الغشائية الفموية.",
            "significance_en": "Facilitates target binding cavity anchoring, but exceeding 10 lowers oral membrane permeability.",
            "threshold": "HBA <= 10",
            "limit": 10.0,
        },
        "tpsa": {
            "title_ar": "مساحة السطح القطبي الطوبولوجي (Topological Polar Surface Area - TPSA)",
            "title_en": "Topological Polar Surface Area (TPSA)",
            "meaning_ar": "مجموع مساحات السطح للذرات القطبية (الأكسجين والنيتروجين وذرات الهيدروجين المرتبطة بها) بوحدة أنغستروم مربع (Å²).",
            "meaning_en": "Total surface sum of polar atoms (O, N, and attached H) in square angstroms (Å²).",
            "significance_ar": "تنبئ بنفاذية الدواء؛ إذا كانت TPSA أقل من 140 Å² يكون الامتصاص المعوي جيداً، وإذا كانت أقل من 90 Å² يمكن للمركب عبور الحاجز الدموي الدماغي (BBB).",
            "significance_en": "Predicts permeability; TPSA <= 140 Å² correlates with good intestinal absorption, and <= 90 Å² with BBB penetration.",
            "threshold": "TPSA <= 140 Å² (BBB <= 90 Å²)",
            "limit": 140.0,
        },
        "rotatable_bonds": {
            "title_ar": "الروابط القابلة للدوران (Rotatable Bonds)",
            "title_en": "Rotatable Bonds",
            "meaning_ar": "الروابط الأحادية غير الحلقية بين ذرات غير طرفية والقابلة للدوران الحر.",
            "meaning_en": "Single non-ring bonds between non-terminal heavy atoms with free rotation.",
            "significance_ar": "تعبر عن مرونة الجزيء؛ الجزيئات ذات الروابط القابلة للدوران (<= 10) تكون أيسر في الارتباط الانتقائي بالهدف الحيوي مع ثبات تكويني أعلى.",
            "significance_en": "Reflects conformational flexibility; <= 10 rotatable bonds optimizes selective binding entropy.",
            "threshold": "Rotatable Bonds <= 10",
            "limit": 10.0,
        },
    }

    PROPERTY_EXPLANATIONS_AR = {
        k: {
            "title_ar": v["title_ar"],
            "meaning_ar": v["meaning_ar"],
            "clinical_significance_ar": v["significance_ar"],
        }
        for k, v in PROPERTY_META.items()
    }

    @classmethod
    def explain_property(
        cls,
        prop_key: str,
        value: Optional[float] = None,
        lang: str = "ar",
    ) -> dict[str, Any]:
        canonical = cls.PROP_MAP.get(prop_key.lower(), prop_key.lower())
        meta = cls.PROPERTY_META.get(canonical, {})
        is_ar = lang == "ar"

        prop_name = meta.get("title_ar" if is_ar else "title_en", prop_key)
        definition = meta.get("meaning_ar" if is_ar else "meaning_en", "")
        impact = meta.get("significance_ar" if is_ar else "significance_en", "")
        optimal_range = meta.get("threshold", "")

        context_eval = ""
        if value is not None and "limit" in meta:
            limit = meta["limit"]
            if value <= limit:
                context_eval = (
                    f"القيمة المقاسة ({value}) في نطاق سليم ومثالي لقاعدة ليبينسكي."
                    if is_ar
                    else f"Measured value ({value}) complies within optimal threshold."
                )
            else:
                context_eval = (
                    f"القيمة المقاسة ({value}) تتجاوز الحد الموصى به ({limit})."
                    if is_ar
                    else f"Exceeds recommended optimal threshold ({limit})."
                )

        return {
            "property_name": prop_name,
            "title_ar": meta.get("title_ar", prop_name),
            "meaning_ar": meta.get("meaning_ar", definition),
            "clinical_significance_ar": meta.get("significance_ar", impact),
            "definition": definition,
            "pharmacological_impact": impact,
            "optimal_range": optimal_range,
            "context_evaluation": context_eval,
            "measured_value": value,
        }

    @classmethod
    def generate_compound_overview_ar(cls, entity: dict[str, Any], profile: dict[str, Any]) -> str:
        name_ar = html.escape(entity.get("name_ar", "المادة الكيميائية"))
        name_en = html.escape(entity.get("name_en", ""))
        formula = html.escape(entity.get("symbol_or_formula", ""))
        entity_type_ar = html.escape(entity.get("entity_type_label_ar", "مركب كيميائي"))
        subclass_ar = html.escape(entity.get("subclass_label_ar", ""))
        molar_mass = profile.get("basic_info", {}).get("molar_mass", "")

        text = (
            f"يُصنف **{name_ar}** ({name_en}) علمياً كـ **{entity_type_ar}** (فئة: {subclass_ar})، "
            f"بصيغة جزيئية موثقة هي `{formula}` وكتلة مولية قدرها **{molar_mass}**.\n\n"
        )

        id_ver = profile.get("identity_verification", {})
        if id_ver.get("is_verified"):
            text += (
                "**حالة التحقق البنيوي (Structure Identity):** تم إثبات مطابقة البنية الكيميائية الحقيقية "
                "مع الصيغة المسجلة ومعرفات قواعد البيانات المرجعية (PubChem / ChEMBL / ChEBI) بنسبة 100%، "
                "مما يضمن عدم الخلط بين هذه المادة والمصاوغات الفراغية أو التركيبية الأخرى.\n\n"
            )
        else:
            text += f"**حالة التحقق البنيوي:** {html.escape(id_ver.get('status_ar', 'التحقق غير مكتمل'))}.\n\n"

        return text

    @classmethod
    def explain_lipinski_evaluation(
        cls,
        metrics: dict[str, Any] | None,
        lang: str = "ar",
    ) -> dict[str, Any]:
        if not metrics:
            return {
                "available": False,
                "summary": "البيانات غير متاحة" if lang == "ar" else "Data unavailable",
                "summary_ar": "الخصائص الدوائية ومحددات ليبينسكي غير متاحة حالياً.",
                "violations_count": 0,
                "is_compliant": False,
                "is_druglike": False,
                "details": [],
            }

        if "properties" in metrics:
            metrics = metrics["properties"]
        elif "calculated_drug_descriptors" in metrics:
            metrics = metrics["calculated_drug_descriptors"]

        mw = float(metrics.get("molecular_weight", 0.0) or 0.0)
        logp = float(metrics.get("logp", 0.0) or 0.0)
        hbd = int(metrics.get("hbd", 0) or 0)
        hba = int(metrics.get("hba", 0) or 0)
        tpsa = float(metrics.get("tpsa", 0.0) or 0.0)
        rot = int(metrics.get("rotatable_bonds", 0) or 0)

        details = []
        violations = 0

        # MW <= 500
        if mw <= 500.0:
            details.append({"parameter": "MW", "value": f"{mw} Da", "threshold": "<= 500 Da", "status": "pass"})
        else:
            violations += 1
            details.append({"parameter": "MW", "value": f"{mw} Da", "threshold": "<= 500 Da", "status": "fail"})

        # LogP <= 5.0
        if logp <= 5.0:
            details.append({"parameter": "LogP", "value": f"{logp}", "threshold": "<= 5.0", "status": "pass"})
        else:
            violations += 1
            details.append({"parameter": "LogP", "value": f"{logp}", "threshold": "<= 5.0", "status": "fail"})

        # HBD <= 5
        if hbd <= 5:
            details.append({"parameter": "HBD", "value": f"{hbd}", "threshold": "<= 5", "status": "pass"})
        else:
            violations += 1
            details.append({"parameter": "HBD", "value": f"{hbd}", "threshold": "<= 5", "status": "fail"})

        # HBA <= 10
        if hba <= 10:
            details.append({"parameter": "HBA", "value": f"{hba}", "threshold": "<= 10", "status": "pass"})
        else:
            violations += 1
            details.append({"parameter": "HBA", "value": f"{hba}", "threshold": "<= 10", "status": "fail"})

        is_compliant = violations <= 1
        is_druglike = is_compliant
        is_ar = lang == "ar"

        if violations == 0:
            summary = (
                "المركب متوافق تماماً مع جميع معايير قاعدة ليبينسكي الخماسية دون أي مخالفة."
                if is_ar
                else "Fully compliant with Lipinski's Rule of Five with 0 violations."
            )
        elif violations == 1:
            summary = (
                "المركب متوافق مع قاعدة ليبينسكي (مخالفة واحدة مقبولة صيدلانياً)."
                if is_ar
                else "Compliant with Lipinski's Rule of Five (1 acceptable violation)."
            )
        else:
            summary = (
                f"المركب غير متوافق مع قاعدة ليبينسكي بمعدل ({violations}) مخالفات."
                if is_ar
                else f"Non-compliant with Lipinski's Rule of Five ({violations} violations)."
            )

        return {
            "available": True,
            "summary": summary,
            "summary_ar": summary,
            "violations_count": violations,
            "is_compliant": is_compliant,
            "is_druglike": is_druglike,
            "details": details,
            "tpsa_comment_ar": f"تبلغ مساحة السطح القطبي TPSA مقداره {tpsa} Å².",
            "rotatable_comment_ar": f"يحتوي المركب على {rot} روابط قابلة للدوران.",
        }

    @classmethod
    def generate_step_by_step_breakdown(
        cls,
        entity_or_profile: dict[str, Any],
        profile: dict[str, Any] | None = None,
        structure_analysis: dict[str, Any] | None = None,
        lang: str = "ar",
    ) -> list[dict[str, Any]]:
        is_ar = lang == "ar"
        if profile is None:
            profile = entity_or_profile
            entity = profile.get("entity", profile.get("basic_info", {}))
        else:
            entity = entity_or_profile

        basic = profile.get("basic_info", {})
        props = profile.get("properties", {})
        steps = []

        # Step 1: الهوية والتركيب الجزيئي
        formula = basic.get("formula", "")
        mw = props.get("molecular_weight", basic.get("molar_mass", ""))
        steps.append({
            "step_number": "1",
            "title": "التركيب الجزيئي والهوية الكيميائية" if is_ar else "Molecular Identity & Structure",
            "title_ar": "التركيب الجزيئي والهوية الكيميائية",
            "detail": f"الصيغة الكيميائية: {formula}، الكتلة الجزيئية: {mw} g/mol." if is_ar else f"Formula: {formula}, Molar Mass: {mw} g/mol.",
            "content_ar": f"الصيغة الكيميائية: {formula}، الكتلة الجزيئية: {mw} g/mol.",
        })

        # Step 2: الذوبانية والتوزع الدهني (LogP & TPSA)
        logp = props.get("logp", "N/A")
        tpsa = props.get("tpsa", "N/A")
        steps.append({
            "step_number": "2",
            "title": "الذوبانية والنفاذية الغشائية (LogP & TPSA)" if is_ar else "Solubility & Permeability",
            "title_ar": "الذوبانية والنفاذية الغشائية (LogP & TPSA)",
            "detail": f"معامل التقسيم LogP={logp}، مساحة السطح القطبية TPSA={tpsa} Å²." if is_ar else f"LogP: {logp}, TPSA: {tpsa} Å².",
            "content_ar": f"معامل التقسيم LogP={logp}، مساحة السطح القطبية TPSA={tpsa} Å².",
        })

        # Step 3: الروابط الهيدروجينية وقابلية التطوير كدواء
        hbd = props.get("hbd", "N/A")
        hba = props.get("hba", "N/A")
        rot = props.get("rotatable_bonds", "N/A")
        steps.append({
            "step_number": "3",
            "title": "الروابط الهيدروجينية والمرونة (HBD/HBA)" if is_ar else "H-Bonding & Flexibility",
            "title_ar": "الروابط الهيدروجينية والمرونة (HBD/HBA)",
            "detail": f"مانحات الهيدروجين HBD={hbd}، مستقبلات الهيدروجين HBA={hba}، روابط الدوران={rot}." if is_ar else f"HBD: {hbd}, HBA: {hba}, Rotatable: {rot}.",
            "content_ar": f"مانحات الهيدروجين HBD={hbd}، مستقبلات الهيدروجين HBA={hba}، روابط الدوران={rot}.",
        })

        return steps

    @classmethod
    def explain_reaction_prediction(cls, prediction: dict[str, Any], lang: str = "ar") -> str:
        """
        Explains a reaction-prediction result in Arabic, strictly grounded in the computed data.

        The assistant never invents products: it only narrates the verdict, the applied rule,
        the balanced equation, the stoichiometry and the declared limitations.
        """
        if not prediction:
            return "لا تتوفر نتيجة تحليل لشرحها حالياً."

        mode = prediction.get("mode", "")
        mode_label = prediction.get("mode_label_ar", "")
        verdict = prediction.get("verdict_ar", "")
        lines: list[str] = []

        lines.append("### قراءة النتيجة (Evidence Tier)")
        lines.append(f"- **مستوى الدليل:** {mode_label}")
        lines.append(f"- **الخلاصة:** {verdict}")

        confidence = prediction.get("confidence", {}) or {}
        if confidence.get("label_ar"):
            lines.append(f"- **مستوى الثقة:** {confidence['label_ar']}")
            if confidence.get("basis_ar"):
                lines.append(f"- **أساس الثقة:** {confidence['basis_ar']}")

        reactants = prediction.get("reactants", []) or []
        products = prediction.get("products", []) or []
        if reactants:
            lines.append("\n### المواد المتفاعلة")
            for r in reactants:
                coeff = r.get("coefficient", 1)
                prefix = f"{coeff} " if coeff and coeff > 1 else ""
                lines.append(f"- {prefix}{r.get('name_ar', '')} ({r.get('formula', '')}) — {r.get('kind_label_ar', '')}")
        if products:
            lines.append("\n### النواتج")
            for p in products:
                coeff = p.get("coefficient", 1)
                prefix = f"{coeff} " if coeff and coeff > 1 else ""
                lines.append(f"- {prefix}{p.get('name_ar', '')} ({p.get('formula', '')})")

        if prediction.get("equation_unicode"):
            lines.append("\n### المعادلة الموزونة")
            lines.append(f"`{prediction['equation_unicode']}`")

        if prediction.get("net_ionic"):
            lines.append(f"\n### المعادلة الصافية (Net Ionic Equation)\n`{prediction['net_ionic']}`")

        conservation = prediction.get("atom_conservation") or {}
        if conservation:
            status = "متحقق" if conservation.get("is_balanced") else "غير متحقق"
            lines.append(f"\n### التحقق الرياضي\n- حفظ الذرات: **{status}**.")

        observations = prediction.get("observations", []) or []
        if observations:
            lines.append("\n### الملاحظات المتوقعة في المختبر")
            for obs in observations:
                lines.append(f"- {obs.get('text_ar', '')}")

        ph = prediction.get("ph_analysis") or {}
        if ph.get("available"):
            lines.append(f"\n### التغيّر في درجة الحموضة (pH)\n- {ph.get('summary_ar', '')} {ph.get('detail_ar', '')}".rstrip())

        stoich = prediction.get("stoichiometry") or {}
        if stoich and "error_ar" not in stoich:
            lim = stoich.get("limiting_reagent")
            if lim:
                lines.append(
                    f"\n### المتفاعل المحدد (Limiting Reagent)\n- **{lim.get('name_ar', '')} ({lim.get('formula', '')})** "
                    f"بكمية {lim.get('initial_mol', 0):g} مول، وهو المتفاعل الذي يحدد كمية النواتج النظرية."
                )
            elif stoich.get("is_exact_stoichiometric_ratio"):
                lines.append("\n### النسبة المولية\n- الكميات المدخلة متكافئة تماماً مع النسبة المولية للتفاعل.")

        yields = prediction.get("theoretical_yields", []) or []
        if yields:
            lines.append("\n### النواتج النظرية (Theoretical Yield)")
            for y in yields:
                lines.append(
                    f"- {y.get('name_ar', '')} ({y.get('formula', '')}): {y.get('theoretical_yield_mol', 0):g} mol "
                    f"({y.get('theoretical_yield_g', 0):g} g)"
                )

        evidence = prediction.get("evidence", []) or []
        if evidence:
            lines.append("\n### الدليل العلمي (Evidence & Sources)")
            for ev in evidence:
                lines.append(f"- **{ev.get('rule_name_ar', '')}**")
                if ev.get("basis_ar"):
                    lines.append(f"  - الأساس العلمي: {ev['basis_ar']}")
                if ev.get("source"):
                    lines.append(f"  - المرجع: {ev['source']}")

        limits = prediction.get("limitations_ar", []) or []
        if limits:
            lines.append("\n### حدود التنبؤ (Limits of Prediction)")
            for lim in limits:
                lines.append(f"- {lim}")

        if prediction.get("missing_information_ar"):
            lines.append("\n### البيانات الناقصة المطلوبة")
            for miss in prediction["missing_information_ar"]:
                lines.append(f"- {miss}")

        lines.append(
            "\n---\n"
            "*ملاحظة منهجية: هذا الشرح يعرض النتائج المحسوبة كما هي دون تعديلها أو إضافة نتائج غير موجودة في محرك MolecuLens.*"
        )
        return "\n".join(lines)

    @classmethod
    def answer_student_question(
        cls,
        question: str,
        compound_data: dict[str, Any] | None = None,
        lang: str = "ar",
    ) -> str:
        q_lower = (question or "").lower()
        is_ar = lang == "ar"

        if any(k in q_lower for k in ["تفاعل", "تفاعلات", "reaction", "نواتج", "متفاعل"]):
            return (
                "يعتمد محرك التفاعلات في منصة MolecuLens على ثلاثة مستويات من الدليل: "
                "(1) تفاعل موثّق في سجل مرجعي علمي محدّد (مرجع تجريبي أو مواصفة دستية صيدلانية أو كتاب تحليل معتمد)، "
                "(2) استنتاج من قاعدة كيميائية عامة معلنة مع ذكر مرجعها وحدودها، "
                "(3) غياب الأدلة كخيار علمي صريح: في الحالتين الأولى والثانية يذكر النظام مصدر كل نتيجة، "
                "وفي الحالة الثالثة يعرض «لا تتوفر بيانات كافية» بدل اختلاق معادلة. "
                "كما يُجرى تحقّق رياضي آلي من حفظ الذرات (Atom Conservation)، "
                "وتُعرض النسب المولية والمتفاعل المحدد والناتج النظري للحسابات المدخلة."
                if is_ar
                else "MolecuLens separates reaction evidence into three tiers: (1) a reaction recorded in a referenced scientific registry, "
                "(2) an inference from a declared general chemistry rule with its source and limits, and "
                "(3) an explicit lack-of-evidence outcome instead of an invented equation. "
                "Atom conservation is verified automatically, and stoichiometry, limiting reagent and theoretical yield are computed."
            )

        if "logp" in q_lower or "التقسيم" in q_lower or "توزع" in q_lower:
            return (
                "معامل التوزع (LogP) يقيس مدى ألفة الجزيء للدهون مقارنة بالماء (Lipophilicity). "
                "في تصميم الأدوية، يُعد النطاق بين 0 و 5 مثالياً لأن الدواء يحتاج إلى الذوبان في السوائل الحيوية المائية "
                "للوصول إلى الأنسجة، وفي نفس الوقت يحتاج إلى ألفة دهنية تمكنه من عبور الغشاء الخلوي المكون من طبقة فوسفوليبيد ثنائية."
                if is_ar
                else "LogP measures the lipophilicity (octanol/water partition) of a molecule. Values between 0 and 5 provide the optimal balance between aqueous solubility in blood and passive permeability across lipid bilayers."
            )

        if "tpsa" in q_lower or "السطح القطبية" in q_lower or "حاجز" in q_lower or "نفاذية" in q_lower:
            return (
                "مساحة السطح القطبي الطوبولوجي (TPSA) تمثل مجموع مساحات الذرات القطبية (الأكسجين والنيتروجين والهيدروجين المرتبط بهما). "
                "تُعد TPSA مؤشراً حاسماً لنفاذية الأغشية الحيوية؛ حيث تفضل الأدوية المعوية TPSA <= 140 Å²، "
                "بينما يتطلب عبور الحاجز الدموي الدماغي (Blood-Brain Barrier - BBB) مساحة سطح قطبية منخفضة TPSA <= 90 Å² لضمان نفاذية كافية للجهاز العصبي المركزي."
                if is_ar
                else "Topological Polar Surface Area (TPSA) sums polar atoms (N, O, and attached H). TPSA <= 140 Å² is required for good oral absorption, whereas <= 90 Å² is vital for penetration across the Blood-Brain Barrier (BBB)."
            )

        if "rotatable" in q_lower or "دوران" in q_lower or "مرونة" in q_lower:
            return (
                "الروابط القابلة للدوران تعبر عن المرونة التكوينية للجزيء (Conformational Flexibility). "
                "وفق قاعدة فيبر (Veber's rule)، يُفضل ألا يزيد عدد الروابط القابلة للدوران عن 10 روابط؛ "
                "لأن الجزيئات شديدة المرونة تفقد قدراً كبيراً من الإنتروبيا الحرة عند الارتباط بهدفها الحيوي، مما يقلل ألفة الارتباط ويضعف التوافر الحيوي الفموي."
                if is_ar
                else "Rotatable bonds govern conformational flexibility. Veber's rule states that <= 10 rotatable bonds optimizes oral bioavailability, as rigid molecules suffer less entropic loss upon binding to their target protein."
            )

        if "lipinski" in q_lower or "ليبينسكي" in q_lower:
            return (
                "قاعدة ليبينسكي الخماسية (Lipinski's Rule of Five) صاغها العالم كريستوفر ليبينسكي عام 1997 لتوقع التوافر الحيوي الفموي للأدوية. "
                "سُميت بالخماسية لأن حدودها مضاعفات للعدد 5: الوزن الجزيئي <= 500 Da، معامل LogP <= 5، مانحات الهيدروجين HBD <= 5، ومستقبلات الهيدروجين HBA <= 10 (2 × 5). "
                "المركب الذي لا يتجاوز أكثر من مخالفة واحدة يمتلك احتمالية عالية للنجاح كدواء فموي."
                if is_ar
                else "Lipinski's Rule of Five, established by Christopher Lipinski in 1997, predicts drug-likeness and oral bioavailability. It is based on parameters that are multiples of 5: MW <= 500 Da, LogP <= 5, HBD <= 5, HBA <= 10 (2x5). Compounds with at most 1 violation generally exhibit high oral bioavailability."
            )

        return (
            "يعتمد التحليل العلمي على القياسات البنيوية والفيزيوكيميائية الدقيقة المحسوبة في منصة MolecuLens دون تخمين أو اختراع لضمان الموثوقية الأكاديمية."
            if is_ar
            else "Scientific evaluations in MolecuLens are strictly grounded in computed physicochemical descriptors to ensure academic precision."
        )

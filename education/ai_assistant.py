"""
Grounded Scientific Educational Assistant for MolecuLens.
Explicitly structured as a 3-layer Rule-Based & Retrieval-Grounded Educational Assistant
(NOT a generative LLM):
  Layer 1 — Scientific Data: computed physicochemical descriptors, equations, and verified records.
  Layer 2 — Retrieval & Context: retrieves grounded context from compound_data and local verified catalogs with citations.
  Layer 3 — Explanation: explains computed results and concepts in Arabic/English without inventing numbers,
            reactions, or unverified drug claims, and explicitly distinguishes reference facts from inferences.
"""

from typing import Any, Optional
import html


class AIScientificAssistant:
    """
    Rule-based and retrieval-grounded educational assistant.
    Never invents numbers, reaction products, or unverified medical/pharmacological claims.
    """

    ASSISTANT_ARCHITECTURE = "rule_based_retrieval_grounded"
    ASSISTANT_LABEL_AR = "مساعد تعليمي قواعدي ومبني على الاسترجاع الموثق (Rule-Based & Retrieval-Grounded Assistant)"
    ASSISTANT_LABEL_EN = "Rule-Based & Retrieval-Grounded Educational Assistant"

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
            "meaning_ar": "مجموع الكتل الذرية لجميع الذرات المكونة للجزيء بوحدة دالتون / جرام لكل مول (Da).",
            "meaning_en": "Sum of atomic weights of all constituent atoms in Da (g/mol).",
            "significance_ar": (
                "مؤشر إرشادي لحجم الجزيء؛ تشير قاعدة ليبينسكي الاستكشافية إلى حد MW <= 500 Da كمرشّح أولي "
                "لكثير من الجزيئات الصغيرة الفموية، لكن توجد أدوية فموية معتمدة تتجاوز هذا الحد ولا يضمن وحده الامتصاص."
            ),
            "significance_en": (
                "An empirical size heuristic; Lipinski's Rule of Five uses MW <= 500 Da as an initial filter "
                "for small oral molecules, though approved oral exceptions exist and MW alone does not guarantee absorption."
            ),
            "threshold": "MW <= 500 Da (إرشادي)",
            "limit": 500.0,
            "source": "Lipinski et al. (1997) Adv. Drug Deliv. Rev. 23(1-3), 3-25",
        },
        "logp": {
            "title_ar": "معامل التقسيم الدهني/المائي (Partition Coefficient - LogP)",
            "title_en": "Partition Coefficient (LogP)",
            "meaning_ar": "لوغاريتم نسبة تركيز الشكل المتعادل للمركب بين الأوكتانول والماء؛ يقيس الألفة للدهون (Lipophilicity).",
            "meaning_en": "Logarithm of the neutral species concentration ratio between octanol and water; measures lipophilicity.",
            "significance_ar": (
                "يُستخدم نطاق LogP المعتدل (غالباً نحو 0 إلى 3 أو <= 5 وفق مؤشر ليبينسكي) كقاعدة إرشادية أولية "
                "لموازنة الذوبانية المائية والنفاذية السلبية عبر الأغشية، لكنه ليس نطاقاً مثالياً عاماً لكل الأدوية "
                "ولا يغني عن قياس الذوبانية والأيض وLogD عند pH فسيولوجي."
            ),
            "significance_en": (
                "Moderate LogP values (heuristic threshold LogP <= 5 in Lipinski's filter) serve as an empirical guideline "
                "balancing aqueous solubility and passive membrane diffusion, not a universal ideal range or guarantee of bioavailability."
            ),
            "threshold": "LogP <= 5.0 (إرشادي)",
            "limit": 5.0,
            "source": "Lipinski et al. (1997); Sangster (1989) J. Phys. Chem. Ref. Data",
        },
        "hbd": {
            "title_ar": "مانحات الروابط الهيدروجينية (Hydrogen Bond Donors - HBD)",
            "title_en": "Hydrogen Bond Donors (HBD)",
            "meaning_ar": "عدد ذرات الهيدروجين المرتبطة بذرات عالية الكهروسالبية (مثل -OH و -NH).",
            "meaning_en": "Count of hydrogen atoms bonded to electronegative atoms (-OH, -NH).",
            "significance_ar": (
                "مؤشر إرشادي: ارتفاع عدد المانحات عن 5 يرتبط إحصائياً بزيادة طاقة نزع الإماهة وصعوبة النفاذ السلبي، "
                "مع وجود استثناءات تعتمد على النقل النشط أو طبيعة الهدف الدوائي."
            ),
            "significance_en": (
                "An empirical heuristic: HBD > 5 statistically correlates with higher desolvation penalty and lower passive permeability, "
                "though active transport and clinical exceptions exist."
            ),
            "threshold": "HBD <= 5 (إرشادي)",
            "limit": 5.0,
            "source": "Lipinski et al. (1997) Adv. Drug Deliv. Rev.",
        },
        "hba": {
            "title_ar": "مستقبلات الروابط الهيدروجينية (Hydrogen Bond Acceptors - HBA)",
            "title_en": "Hydrogen Bond Acceptors (HBA)",
            "meaning_ar": "عدد ذرات النيتروجين والأكسجين القادرة على استقبال رابطة هيدروجينية.",
            "meaning_en": "Count of N and O atoms capable of accepting a hydrogen bond.",
            "significance_ar": (
                "مؤشر إرشادي ضمن قاعدة ليبينسكي (<= 10)؛ يساعد في التفاعل القطبي مع الأهداف الحيوية، "
                "وزيادته المفرطة قد تقلل النفاذية السلبية دون أن تمنع التطوير الدوائي بصورة مطلقة."
            ),
            "significance_en": (
                "An empirical heuristic within Lipinski's filter (<= 10); excess acceptors often correlate with lower passive "
                "permeability but do not rigidly preclude drug development."
            ),
            "threshold": "HBA <= 10 (إرشادي)",
            "limit": 10.0,
            "source": "Lipinski et al. (1997) Adv. Drug Deliv. Rev.",
        },
        "tpsa": {
            "title_ar": "مساحة السطح القطبي الطوبولوجي (Topological Polar Surface Area - TPSA)",
            "title_en": "Topological Polar Surface Area (TPSA)",
            "meaning_ar": "مجموع مساحات السطح التقريبية للذرات القطبية (الأكسجين والنيتروجين والهيدروجين المرتبط بهما) بوحدة (Å²).",
            "meaning_en": "Approximate surface sum over polar atoms (O, N, and attached H) in Å².",
            "significance_ar": (
                "مؤشر إرشادي يرتبط إحصائياً بالنفاذية السلبية: القيم <= 140 Å² (قاعدة فيبر) ترتبط غالباً بفرص أفضل "
                "للامتصاص المعوي السلبي، والقيم <= 90 Å² تُستخدم كمؤشر مبدئي لاحتمال عبور الحاجز الدموي الدماغي (BBB)، "
                "لكنها ليست شروطاً لازمة أو حاسمة ولا تضمن النفاذية الفعلية (التي تتأثر بمضخات الطرد مثل P-gp والتأين)."
            ),
            "significance_en": (
                "An empirical heuristic statistically associated with passive permeability: TPSA <= 140 Å² (Veber) is a guideline "
                "for oral passive absorption and <= 90 Å² is a rough guideline for CNS/BBB penetration, NOT a strict requirement or guarantee."
            ),
            "threshold": "TPSA <= 140 Å² إرشادي للامتصاص (<= 90 Å² إرشادي للـ BBB)",
            "limit": 140.0,
            "source": "Ertl et al. (2000) J. Med. Chem.; Veber et al. (2002) J. Med. Chem. 45(12), 2615-2623",
        },
        "rotatable_bonds": {
            "title_ar": "الروابط القابلة للدوران (Rotatable Bonds)",
            "title_en": "Rotatable Bonds",
            "meaning_ar": "الروابط الأحادية غير الحلقية بين ذرات ثقيلة غير طرفية والقابلة للدوران.",
            "meaning_en": "Single non-ring bonds between non-terminal heavy atoms.",
            "significance_ar": (
                "مؤشر إرشادي للمرونة التكوينية؛ تقترح قاعدة فيبر <= 10 روابط كدليل عام لتقليل الخسارة الإنتروبية "
                "عند الارتباط، لكنها قاعدة إرشادية وليست شرطاً قاطعاً."
            ),
            "significance_en": (
                "An empirical flexibility heuristic; Veber's guideline suggests <= 10 rotatable bonds to limit entropic loss, "
                "not a strict requirement."
            ),
            "threshold": "Rotatable Bonds <= 10 (إرشادي)",
            "limit": 10.0,
            "source": "Veber et al. (2002) J. Med. Chem.",
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
                    f"القيمة المحسوبة ({value}) تقع ضمن الحد الإرشادي الشائع (<= {limit})؛ وهذا مؤشر أولي وليس ضماناً للأداء الحيوي."
                    if is_ar
                    else f"Computed value ({value}) falls within the common heuristic guideline (<= {limit}); indicative only, not a guarantee."
                )
            else:
                context_eval = (
                    f"القيمة المحسوبة ({value}) تتجاوز الحد الإرشادي الشائع ({limit})، مما يستدعي مراجعة أوسع للخصائص الحركية."
                    if is_ar
                    else f"Computed value ({value}) exceeds the common heuristic guideline ({limit})."
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
            "source": meta.get("source", ""),
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
                "**حالة التحقق البنيوي (Structure Identity):** تم التحقق من تطابق الصيغة المستخرجة من البنية "
                "مع الصيغة المسجلة ومعرفات قواعد البيانات المرجعية.\n\n"
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
        """
        Evaluates Lipinski's Rule of Five strictly without coercing missing values to 0.0.
        - Missing required descriptors produce status='unavailable' (NEVER 'pass').
        - Explicitly separates Lipinski rule compliance from broader Drug-likeness.
        """
        is_ar = lang == "ar"
        caveat_ar = (
            "تنبيه منهجي: قاعدة ليبينسكي مرشّح إرشادي أولي للتوافق البنيوي فقط (Lipinski Compliance) "
            "وليست مكافئة لمفهوم الصلاحية الدوائية الكاملة (Drug-likeness)، ولا تتنبأ وحدها بالفعالية أو السلامة أو التوافر الحيوي الفموي الفعلي."
        )
        caveat_en = (
            "Methodological caveat: Lipinski compliance is an empirical heuristic filter only, distinct from overall "
            "drug-likeness, and does not by itself predict efficacy, safety, or actual oral bioavailability."
        )

        if not metrics:
            return {
                "available": False,
                "status": "unavailable",
                "summary": "البيانات غير متاحة لتقييم قاعدة ليبينسكي." if is_ar else "Data unavailable for Lipinski evaluation.",
                "summary_ar": "الخصائص الفيزيوكيميائية غير متاحة لتقييم قاعدة ليبينسكي.",
                "violations_count": None,
                "missing_parameters": ["MW", "LogP", "HBD", "HBA"],
                "is_compliant": None,
                "lipinski_compliant": None,
                "is_druglike": None,
                "drug_likeness_note_ar": caveat_ar,
                "drug_likeness_note_en": caveat_en,
                "details": [],
            }

        if "properties" in metrics:
            metrics = metrics["properties"]
        elif "calculated_drug_descriptors" in metrics:
            metrics = metrics["calculated_drug_descriptors"]

        def _parse_num(val: Any) -> Optional[float]:
            if val is None or isinstance(val, bool):
                return None
            try:
                f = float(val)
                return f if f == f else None  # reject NaN
            except (TypeError, ValueError):
                return None

        mw = _parse_num(metrics.get("molecular_weight"))
        logp = _parse_num(metrics.get("logp"))
        hbd = _parse_num(metrics.get("hbd"))
        hba = _parse_num(metrics.get("hba"))
        tpsa = _parse_num(metrics.get("tpsa"))
        rot = _parse_num(metrics.get("rotatable_bonds"))

        checks = [
            ("MW", mw, 500.0, "<= 500 Da", " Da"),
            ("LogP", logp, 5.0, "<= 5.0", ""),
            ("HBD", hbd, 5.0, "<= 5", ""),
            ("HBA", hba, 10.0, "<= 10", ""),
        ]

        details: list[dict[str, Any]] = []
        missing_params: list[str] = []
        violations = 0

        for label, val, limit, threshold_str, unit in checks:
            if val is None:
                missing_params.append(label)
                details.append({
                    "parameter": label,
                    "value": "غير متاح (N/A)" if is_ar else "Unavailable",
                    "threshold": threshold_str,
                    "status": "unavailable",
                })
            elif val <= limit:
                disp = f"{int(val)}" if label in {"HBD", "HBA"} else f"{val}{unit}"
                details.append({"parameter": label, "value": disp, "threshold": threshold_str, "status": "pass"})
            else:
                violations += 1
                disp = f"{int(val)}" if label in {"HBD", "HBA"} else f"{val}{unit}"
                details.append({"parameter": label, "value": disp, "threshold": threshold_str, "status": "fail"})

        # Missing required descriptors MUST produce unavailable, never a false pass!
        if missing_params:
            miss_str = ", ".join(missing_params)
            summary_ar = (
                f"تقييم قاعدة ليبينسكي غير مكتمل (unavailable) بسبب نقص الخصائص التالية: {miss_str}. "
                "لا تُعتبر البيانات الناقصة متوافقة مع القاعدة."
            )
            summary_en = (
                f"Lipinski evaluation unavailable due to missing descriptors: {miss_str}. "
                "Missing data is never treated as a pass."
            )
            return {
                "available": False,
                "status": "unavailable",
                "summary": summary_ar if is_ar else summary_en,
                "summary_ar": summary_ar,
                "violations_count": None,
                "missing_parameters": missing_params,
                "is_compliant": None,
                "lipinski_compliant": None,
                "is_druglike": None,
                "drug_likeness_note_ar": caveat_ar,
                "drug_likeness_note_en": caveat_en,
                "details": details,
                "tpsa_comment_ar": f"تبلغ مساحة السطح القطبي TPSA مقداره {tpsa} Å² (مؤشر إرشادي فقط)." if tpsa is not None else "قيمة TPSA غير متاحة.",
                "rotatable_comment_ar": f"يحتوي المركب على {int(rot)} روابط قابلة للدوران." if rot is not None else "عدد الروابط القابلة للدوران غير متاح.",
            }

        lipinski_compliant = violations <= 1

        if violations == 0:
            summary = (
                "المركب يحقق جميع حدود قاعدة ليبينسكي الإرشادية الأربعة (0 مخالفات). "
                "يُعد هذا مؤشراً أولياً فقط ولا يثبت وحده النجاح الدوائي أو التوافر الحيوي."
                if is_ar
                else "Meets all 4 empirical Lipinski Rule of Five thresholds (0 violations); an initial heuristic only, not proof of oral bioavailability or drug efficacy."
            )
        elif violations == 1:
            summary = (
                "المركب يسجل مخالفة واحدة لحدود قاعدة ليبينسكي الإرشادية (ضمن هامش التسامح التقليدي <= 1 مخالفة)، "
                "وهذا لا يضمن وحده النجاح كدواء فموي ويتطلب تقييماً حركياً وسمّياً شاملاً."
                if is_ar
                else "Has 1 Lipinski threshold violation (within the traditional <= 1 heuristic allowance); does not guarantee oral drug success."
            )
        else:
            summary = (
                f"المركب يتجاوز حدود قاعدة ليبينسكي الإرشادية بعدد ({violations}) مخالفات؛ "
                "مما يشير إحصائياً إلى تحديات محتملة في الامتصاص السلبي (مع وجود استثناءات دوائية معروفة)."
                if is_ar
                else f"Exceeds empirical Lipinski thresholds ({violations} violations), indicating potential passive permeability challenges."
            )

        return {
            "available": True,
            "status": "evaluated",
            "summary": summary,
            "summary_ar": summary,
            "violations_count": violations,
            "missing_parameters": [],
            "is_compliant": lipinski_compliant,
            "lipinski_compliant": lipinski_compliant,
            # Explicitly decoupled: Lipinski compliance alone is NOT overall drug-likeness
            "is_druglike": None,
            "drug_likeness_note_ar": caveat_ar,
            "drug_likeness_note_en": caveat_en,
            "details": details,
            "tpsa_comment_ar": f"تبلغ مساحة السطح القطبي TPSA مقداره {tpsa} Å² (مؤشر إرشادي للنفاذية السلبية)." if tpsa is not None else "قيمة TPSA غير متاحة.",
            "rotatable_comment_ar": f"يحتوي المركب على {int(rot)} روابط قابلة للدوران." if rot is not None else "عدد الروابط القابلة للدوران غير متاح.",
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
        basic = profile.get("basic_info", {})
        props = profile.get("properties", {})
        steps = []

        formula = basic.get("formula", "")
        mw = props.get("molecular_weight", basic.get("molar_mass", "غير متاح"))
        steps.append({
            "step_number": "1",
            "title": "التركيب الجزيئي والهوية الكيميائية" if is_ar else "Molecular Identity & Structure",
            "title_ar": "التركيب الجزيئي والهوية الكيميائية",
            "detail": f"الصيغة الكيميائية: {formula}، الكتلة الجزيئية: {mw} g/mol." if is_ar else f"Formula: {formula}, Molar Mass: {mw} g/mol.",
            "content_ar": f"الصيغة الكيميائية: {formula}، الكتلة الجزيئية: {mw} g/mol.",
        })

        logp = props.get("logp", "N/A")
        tpsa = props.get("tpsa", "N/A")
        steps.append({
            "step_number": "2",
            "title": "المؤشرات الإرشادية للذوبانية والنفاذية (LogP & TPSA)" if is_ar else "Heuristic Solubility & Permeability Indicators",
            "title_ar": "المؤشرات الإرشادية للذوبانية والنفاذية (LogP & TPSA)",
            "detail": f"معامل التقسيم LogP={logp}، مساحة السطح القطبية TPSA={tpsa} Å² (مؤشرات إرشادية أولية)." if is_ar else f"LogP: {logp}, TPSA: {tpsa} Å² (empirical heuristics).",
            "content_ar": f"معامل التقسيم LogP={logp}، مساحة السطح القطبية TPSA={tpsa} Å² (مؤشرات إرشادية أولية).",
        })

        hbd = props.get("hbd", "N/A")
        hba = props.get("hba", "N/A")
        rot = props.get("rotatable_bonds", "N/A")
        steps.append({
            "step_number": "3",
            "title": "الروابط الهيدروجينية والمرونة التكوينية (HBD/HBA)" if is_ar else "H-Bonding & Flexibility",
            "title_ar": "الروابط الهيدروجينية والمرونة التكوينية (HBD/HBA)",
            "detail": f"مانحات الهيدروجين HBD={hbd}، مستقبلات الهيدروجين HBA={hba}، روابط الدوران={rot}." if is_ar else f"HBD: {hbd}, HBA: {hba}, Rotatable: {rot}.",
            "content_ar": f"مانحات الهيدروجين HBD={hbd}، مستقبلات الهيدروجين HBA={hba}، روابط الدوران={rot}.",
        })

        return steps

    @classmethod
    def explain_reaction_prediction(cls, prediction: dict[str, Any], lang: str = "ar") -> str:
        """
        Explains a reaction-prediction result in Arabic, strictly grounded in the computed data.
        """
        if not prediction:
            return "لا تتوفر نتيجة تحليل لشرحها حالياً."

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
            f"*ملاحظة منهجية: صادر عن {cls.ASSISTANT_LABEL_AR} — يعرض النتائج المحسوبة كما هي دون تعديل أو اختراع.*"
        )
        return "\n".join(lines)

    @classmethod
    def _extract_compound_context(cls, compound_data: dict[str, Any] | None) -> dict[str, Any]:
        """Layer 1 & 2 helper: extracts verified scientific data & citations from compound_data."""
        if not compound_data or not isinstance(compound_data, dict):
            return {"has_compound": False, "props": {}, "sources": []}

        entity = compound_data.get("entity", compound_data)
        profile = compound_data.get("profile", compound_data)
        props = (
            profile.get("properties")
            or compound_data.get("properties")
            or compound_data.get("calculated_drug_descriptors")
            or {}
        )
        basic = profile.get("basic_info", {})
        sources = list(entity.get("sources", []) or profile.get("sources", []) or [])
        name_ar = entity.get("name_ar") or basic.get("name_ar") or entity.get("name_en") or ""
        name_en = entity.get("name_en") or basic.get("name_en") or ""
        formula = entity.get("symbol_or_formula") or basic.get("formula") or ""
        return {
            "has_compound": bool(name_ar or name_en or formula or props),
            "name_ar": name_ar,
            "name_en": name_en,
            "formula": formula,
            "props": dict(props) if isinstance(props, dict) else {},
            "sources": sources,
        }

    @classmethod
    def answer_student_question_structured(
        cls,
        question: str,
        compound_data: dict[str, Any] | None = None,
        lang: str = "ar",
    ) -> dict[str, Any]:
        """
        3-Layer Grounded Scientific Q&A:
          - Layer 1 (scientific_data): exact computed numbers from compound_data (never invented).
          - Layer 2 (retrieval_context): retrieved facts and explicit source citations.
          - Layer 3 (explanation): clear Arabic/English explanation separating reference facts from inferences.
        """
        q_lower = (question or "").strip().lower()
        is_ar = lang == "ar"
        ctx = cls._extract_compound_context(compound_data)
        props = ctx["props"]
        sources_used: list[str] = list(ctx["sources"])
        reference_facts: list[str] = []
        inferences: list[str] = []
        grounded_in_compound = False

        comp_label = (
            f"{ctx['name_ar']} ({ctx['formula']})" if ctx.get("name_ar") and ctx.get("formula")
            else (ctx.get("name_ar") or ctx.get("name_en") or ctx.get("formula") or "")
        )

        if any(k in q_lower for k in ["تفاعل", "تفاعلات", "reaction", "نواتج", "متفاعل"]):
            reference_facts.append(
                "يميز محرك MolecuLens بين: (1) التفاعلات الموثّقة بمراجع علمية، (2) حالات عدم التفاعل المثبتة، "
                "(3) الاستنتاج القواعدي المعلن المصدر والحدود، و(4) غياب الأدلة الكافية دون اختلاق نواتج."
                if is_ar
                else "MolecuLens distinguishes documented reactions, explicit no-reaction cases, declared general-rule inferences, and insufficient evidence without inventing products."
            )
            sources_used.append("MolecuLens VERIFIED_REACTIONS_REGISTRY & Prediction Engine")
            if comp_label:
                grounded_in_compound = True
                inferences.append(
                    f"للمادة الحالية ({comp_label}): يمكن تمريرها مباشرة إلى مختبر التفاعلات (Reaction Lab) لاختبار تفاعلاتها الموثقة أو شبكة التفاعلات متعددة المواد."
                    if is_ar
                    else f"For {comp_label}: pass it into Reaction Lab to inspect verified reactions or multi-reactant networks."
                )

        elif "logp" in q_lower or "التقسيم" in q_lower or "توزع" in q_lower:
            sources_used.append(cls.PROPERTY_META["logp"]["source"])
            reference_facts.append(cls.PROPERTY_META["logp"]["significance_ar" if is_ar else "significance_en"])
            if props.get("logp") is not None:
                grounded_in_compound = True
                val = props["logp"]
                inferences.append(
                    f"القيمة المحسوبة لـ {comp_label or 'المركب الحالي'}: LogP = {val} "
                    f"({'ضمن الحد الإرشادي <= 5' if float(val) <= 5.0 else 'يتجاوز الحد الإرشادي 5.0'}؛ مؤشر أولي وليس ضماناً للذوبانية أو الامتصاص)."
                    if is_ar
                    else f"Computed value for {comp_label or 'this compound'}: LogP = {val} (heuristic indicator only)."
                )

        elif "tpsa" in q_lower or "السطح القطب" in q_lower or "حاجز" in q_lower or "نفاذية" in q_lower:
            sources_used.append(cls.PROPERTY_META["tpsa"]["source"])
            reference_facts.append(cls.PROPERTY_META["tpsa"]["significance_ar" if is_ar else "significance_en"])
            if props.get("tpsa") is not None:
                grounded_in_compound = True
                val = props["tpsa"]
                inferences.append(
                    f"القيمة المحسوبة لـ {comp_label or 'المركب الحالي'}: TPSA = {val} Å² "
                    "(مؤشر إرشادي للنفاذية السلبية وليس شرطاً حاسماً أو إثباتاً لعبور الحواجز الحيوية)."
                    if is_ar
                    else f"Computed value for {comp_label or 'this compound'}: TPSA = {val} Å² (empirical guideline only)."
                )

        elif "rotatable" in q_lower or "دوران" in q_lower or "مرونة" in q_lower:
            sources_used.append(cls.PROPERTY_META["rotatable_bonds"]["source"])
            reference_facts.append(cls.PROPERTY_META["rotatable_bonds"]["significance_ar" if is_ar else "significance_en"])
            if props.get("rotatable_bonds") is not None:
                grounded_in_compound = True
                inferences.append(
                    f"عدد الروابط القابلة للدوران المحسوب لـ {comp_label or 'المركب الحالي'}: {props['rotatable_bonds']}."
                    if is_ar
                    else f"Computed rotatable bonds for {comp_label or 'this compound'}: {props['rotatable_bonds']}."
                )

        elif "lipinski" in q_lower or "ليبينسكي" in q_lower:
            sources_used.append(cls.PROPERTY_META["molecular_weight"]["source"])
            reference_facts.append(
                "قاعدة ليبينسكي (1997) هي مرشّح إرشادي أولي (MW <= 500، LogP <= 5، HBD <= 5، HBA <= 10) "
                "ولا تعني وحدها الصلاحية الدوائية الكاملة أو ضمان التوافر الحيوي أو الفعالية أو السلامة."
                if is_ar
                else "Lipinski's Rule of Five (1997) is an empirical initial filter (MW <= 500, LogP <= 5, HBD <= 5, HBA <= 10); it does not by itself predict efficacy, safety, or actual oral bioavailability."
            )
            if props:
                lip = cls.explain_lipinski_evaluation(props, lang=lang)
                grounded_in_compound = True
                inferences.append(
                    f"تقييم {comp_label or 'المركب الحالي'}: {lip['summary']}"
                    if is_ar
                    else f"Evaluation for {comp_label or 'this compound'}: {lip['summary']}"
                )

        elif any(k in q_lower for k in ["كتلة", "وزن", "صيغة", "formula", "mw", "mass"]) and ctx["has_compound"]:
            grounded_in_compound = True
            mw_val = props.get("molecular_weight")
            reference_facts.append(
                f"البيانات المسجلة للمادة ({comp_label}): الصيغة = {ctx.get('formula') or 'غير متاحة'}، "
                f"الوزن الجزيئي المحسوب = {mw_val if mw_val is not None else 'غير متاح'}."
                if is_ar
                else f"Recorded data for {comp_label}: formula={ctx.get('formula') or 'N/A'}, MW={mw_val if mw_val is not None else 'N/A'}."
            )

        else:
            # Honest unavailability response — never hallucinate drug/medical facts!
            if ctx["has_compound"] and comp_label:
                grounded_in_compound = True
                reference_facts.append(
                    f"لا تتوفر إجابة مرجعية مباشرة لهذا السؤال بعينه داخل القاعدة الموثقة الحالية. "
                    f"البيانات المتاحة حالياً للمادة ({comp_label}) تقتصر على الصيغة والخصائص الفيزيوكيميائية المحسوبة."
                    if is_ar
                    else f"No verified direct answer is available for this specific question in the current grounded knowledge base for {comp_label}."
                )
            else:
                reference_facts.append(
                    "لا تتوفر إجابة علمية موثقة لهذا السؤال في قاعدة البيانات المحلية الحالية؛ "
                    "يمتنع المساعد عن التخمين أو إضافة معلومات دوائية دون مصدر."
                    if is_ar
                    else "No verified scientific answer is available for this question in the local knowledge base; the assistant does not fabricate facts."
                )

        dedup_sources = list(dict.fromkeys(s for s in sources_used if s))
        parts = []
        if reference_facts:
            parts.append(("**حقائق مرجعية:** " if is_ar else "**Reference facts:** ") + " ".join(reference_facts))
        if inferences:
            parts.append(("**قراءة سياقية لبيانات المركب:** " if is_ar else "**Context interpretation:** ") + " ".join(inferences))
        if dedup_sources:
            parts.append(("**المصادر المستخدمة:** " if is_ar else "**Sources:** ") + " | ".join(dedup_sources))
        else:
            parts.append(
                "**ملاحظة:** لا توجد مصادر دوائية مرفقة لهذا الاستعلام؛ يرجى الرجوع إلى المراجع المعتمدة."
                if is_ar
                else "**Note:** No pharmacological sources matched this query."
            )

        return {
            "assistant_architecture": cls.ASSISTANT_ARCHITECTURE,
            "assistant_label_ar": cls.ASSISTANT_LABEL_AR,
            "grounded_in_compound_data": grounded_in_compound,
            "layer1_scientific_data": ctx,
            "layer2_sources_used": dedup_sources,
            "layer3_reference_facts": reference_facts,
            "layer3_inferences": inferences,
            "answer_text": "\n\n".join(parts),
        }

    @classmethod
    def answer_student_question(
        cls,
        question: str,
        compound_data: dict[str, Any] | None = None,
        lang: str = "ar",
    ) -> str:
        """Returns the formatted grounded educational answer string (backwards-compatible API)."""
        return cls.answer_student_question_structured(question, compound_data=compound_data, lang=lang)["answer_text"]

from typing import Any
from chemistry_engine.formula_analyzer import format_chemical_equation_display, format_chemical_formula_display
from chemistry_engine.periodic_table import ELEMENT_NAMES_AR
from config.settings import OP_SCIENTIFIC_CALCULATION, OP_SCIENTIFIC_PREDICTION, OP_SCIENTIFIC_RETRIEVAL, REACTION_DOMAINS
from core.schemas import build_provenance_tag
from core.validators import validate_reaction_inputs
from reaction_engine.equation_balancer import (
    balance_reaction_algebraically,
    format_balanced_equation,
    verify_atom_conservation,
)
from reaction_engine.reaction_registry import CATEGORY_LABELS_AR, find_verified_reactions
from reaction_engine.stoichiometry import compute_reaction_stoichiometry
from reaction_engine.yield_calculator import compute_theoretical_yields

NOBLE_GASES = {"He", "Ne", "Ar", "Kr", "Xe", "Rn"}
DOMAIN_LABELS_DICT = dict(REACTION_DOMAINS)


def evaluate_reaction_conditions_compatibility(
    documented_conditions: dict[str, Any],
    user_conditions: dict[str, Any] | None,
    category: str = "",
) -> dict[str, Any]:
    """Compares user-specified reaction conditions (temperature, pressure, solvent, catalyst) against documented reaction requirements."""
    if not user_conditions:
        return {
            "is_compatible": True,
            "status_ar": "متوافق مع الظروف المرجعية القياسية (Compatible with Reference Conditions)",
            "warnings_ar": [],
            "details_ar": "تُطبق الظروف المرجعية المسجلة للتفاعل دون تعارض.",
        }

    doc_temp = str(documented_conditions.get("temperature", "")).lower()
    doc_press = str(documented_conditions.get("pressure", "")).lower()
    doc_solv = str(documented_conditions.get("solvent", "")).lower()
    doc_cat = str(documented_conditions.get("catalyst", "")).lower()
    combined_doc = f"{doc_temp} {doc_press} {doc_solv} {doc_cat}"

    u_temp = str(user_conditions.get("temperature", "")).strip()
    u_press = str(user_conditions.get("pressure", "")).strip()
    u_solv = str(user_conditions.get("solvent", "")).strip()
    u_cat = str(user_conditions.get("catalyst", "")).strip()

    warnings: list[str] = []

    # 1. Temperature compatibility check
    if u_temp:
        requires_high_heat = any(k in combined_doc for k in ["850", "500", "400", "فرن", "شرارة", "اشتعال", "تسخين", "حرارة مرتفعة"])
        requires_physiological = "37" in combined_doc or category == "biochemical_reaction"
        if requires_high_heat and ("25 °C" in u_temp or "37 °C" in u_temp):
            warnings.append(
                f"تنبيه درجة الحرارة (Temperature Mismatch): التفاعل الموثق يتطلب تسخيناً حرارياً أو شرارة اشتعال ({documented_conditions.get('temperature')})، بينما تم اختيار ({u_temp})."
            )
        elif requires_physiological and ("450 °C" in u_temp or "850 °C" in u_temp or "100 °C" in u_temp):
            warnings.append(
                f"تنبيه حيوي (Thermal Denaturation): التفاعل الكيميائي الحيوي يعمل عند ({documented_conditions.get('temperature', '37 °C')})، واختيار ({u_temp}) يؤدي إلى تمسخ الإنزيمات وتوقف التفاعل."
            )

    # 2. Pressure compatibility check
    if u_press:
        requires_high_press = any(k in combined_doc for k in ["150", "200", "250", "ضغط وحرارة عاليان", "ضغط مرتفع"])
        if requires_high_press and "1 atm" in u_press and "أو" not in doc_press:
            warnings.append(
                f"تنبيه الضغط (Pressure Mismatch): يتطلب هذا التفاعل ضغطاً مرتفعاً ({documented_conditions.get('pressure')}) لتحقيق مردود عملي، بينما تم اختيار ({u_press})."
            )

    # 3. Solvent / medium compatibility check
    if u_solv:
        requires_aqueous = "محلول مائي" in combined_doc or category in {"neutralization", "precipitation"}
        if requires_aqueous and ("بدون مذيب" in u_solv or "Gas-Solid" in u_solv):
            warnings.append(
                "تنبيه الوسط/المذيب (Solvent Mismatch): يتطلب هذا التفاعل وسطاً مائياً (Aqueous Solution) لتفكك الأيونات وتفاعلها، بينما تم اختيار طور جاف بدون مذيب."
            )
        elif category == "combustion" and ("محلول مائي" in u_solv or "Buffer" in u_solv):
            warnings.append(
                "تنبيه الوسط/المذيب (Solvent Mismatch): تفاعلات الاحتراق تحدث في الطور الغازي/المباشر بوجود الأكسجين، ولا تتم داخل المحاليل المائية."
            )

    # 4. Catalyst compatibility check
    if u_cat:
        u_cat_none = "بدون عامل حفاز" in u_cat or u_cat == "none"
        if ("h2so4" in combined_doc and ("أسترة" in combined_doc or "نترتة" in combined_doc)) and not ("H2SO4" in u_cat or "حمضي" in u_cat):
            warnings.append(
                f"تنبيه العامل الحفاز (Catalyst Mismatch): يتطلب هذا التفاعل العضوي حفاز حمض الكبريتيك المركز (H₂SO₄) لتنشيط التفاعل، بينما تم اختيار ({u_cat})."
            )
        elif ("إنزيم" in combined_doc or category == "biochemical_reaction") and u_cat_none:
            warnings.append(
                f"تنبيه التحفيز الحيوي (Enzymatic Requirement): يعتمد هذا التفاعل الحيوي على تحفيز إنزيمي متخصص ({documented_conditions.get('catalyst', 'إنزيم حيوي')})، واختيار (بدون عامل حفاز) يمنع سيره بالمعدل الفسيولوجي."
            )
        elif ("fe" in combined_doc and "هابر" in combined_doc) and not ("Fe" in u_cat or "الحديد" in u_cat):
            warnings.append(
                "تنبيه العامل الحفاز (Catalyst Mismatch): تتطلب طريقة هابر-بوش حفاز الحديد المنشط (Promoted Iron Catalyst - Fe) لكسر الرابطة الثلاثية في النيتروجين."
            )
        elif ("pt" in combined_doc or "pd" in combined_doc or "لندلار" in combined_doc) and "هدرجة" in combined_doc and not ("Pt" in u_cat or "Pd" in u_cat):
            warnings.append(
                "تنبيه العامل الحفاز (Catalyst Mismatch): تتطلب تفاعلات الهدرجة حفازاً فلزياً انتقالياً مثل البلاتين أو البالاديوم (Pt / Pd)."
            )

    is_compat = len(warnings) == 0
    return {
        "is_compatible": is_compat,
        "status_ar": (
            "متوافق مع ظروف التفاعل الموثقة (Compatible Conditions)"
            if is_compat
            else "يوجد اختلاف بين الظروف المختارة وظروف التفاعل المرجعية (Condition Mismatch Warning)"
        ),
        "warnings_ar": warnings,
        "details_ar": (
            "جميع الظروف المدخلة (الحرارة، الضغط، المذيب، الحفاز) متوافقة مع شروط التفاعل الموثق."
            if is_compat
            else " | ".join(warnings)
        ),
    }


def analyze_chemical_reaction(
    reactant_inputs: list[dict[str, Any]],
    domain: str = "all_domains",
    custom_conditions: dict[str, Any] | None = None,
) -> dict[str, Any]:
    is_valid, err_msg = validate_reaction_inputs(reactant_inputs)
    if not is_valid:
        return {
            "supported": False,
            "error_type": "validation_error",
            "rejection_reason_ar": err_msg,
        }

    formulas: list[str] = []
    names_ar: list[str] = []
    for item in reactant_inputs:
        sub = item.get("substance", {})
        f = sub.get("symbol_or_formula") or item.get("formula", "")
        n_ar = sub.get("name_ar") or f
        if f:
            formulas.append(f)
            names_ar.append(n_ar)

    for f in formulas:
        if f in NOBLE_GASES:
            el_ar = ELEMENT_NAMES_AR.get(f, f)
            return {
                "supported": False,
                "error_type": "chemically_inert",
                "reactant_formulas": formulas,
                "rejection_reason_ar": (
                    f"رفض علمي (Chemically Inert): المادة المدخلة ({el_ar} — {f}) تنتمي إلى الغازات النبيلة الخاملة ذات الغلاف الإلكتروني المكتمل، "
                    "ولا تدخل في تفاعلات كيميائية عادية تحت الظروف القياسية. يلتزم النظام بعدم اختلاق نواتج غير حقيقية."
                ),
            }

    rtype = (custom_conditions or {}).get("reaction_type")
    candidates = find_verified_reactions(formulas, reaction_type=rtype, domain=domain)
    rxn = candidates[0] if candidates else None

    if rxn is None:
        # Check if the reaction exists in another reaction_type or another domain
        any_domain_or_type = find_verified_reactions(formulas, reaction_type=None, domain=None)
        if any_domain_or_type:
            if rtype and rtype != "all" and not find_verified_reactions(formulas, reaction_type=rtype, domain=None):
                valid_types = "، ".join(
                    dict.fromkeys(CATEGORY_LABELS_AR.get(c.get("category", ""), c.get("category", "")) for c in any_domain_or_type)
                )
                return {
                    "supported": False,
                    "error_type": "reaction_type_mismatch",
                    "reactant_formulas": formulas,
                    "rejection_reason_ar": (
                        "عدم تطابق نوع التفاعل (Reaction Type Mismatch): يوجد تفاعل موثق لهذه المواد لكنه من نوع آخر غير النوع المختار. "
                        f"الحل المقترح: اختر خيار (تلقائي بحسب المواد المتفاعلة) أو اختر النوع المناسب: {valid_types}."
                    ),
                }
            if domain and domain not in {"all", "all_domains"}:
                valid_domains = "، ".join(dict.fromkeys(c.get("domain_label_ar", c.get("domain", "")) for c in any_domain_or_type))
                chosen_dom_ar = DOMAIN_LABELS_DICT.get(domain, domain)
                return {
                    "supported": False,
                    "error_type": "reaction_domain_mismatch",
                    "reactant_formulas": formulas,
                    "rejection_reason_ar": (
                        f"عدم تطابق مجال التفاعل (Domain Filter Mismatch): التفاعل الموثق لهذه المواد لا يندرج ضمن المجال المختار حالياً ({chosen_dom_ar}). "
                        f"الحل المقترح: اختر (جميع المجالات الكيميائية - All Domains) أو غير المجال إلى: {valid_domains}."
                    ),
                }

        joined = " + ".join(f"{n} ({format_chemical_formula_display(f)})" for n, f in zip(names_ar, formulas))
        return {
            "supported": False,
            "error_type": "unsupported_reaction",
            "reactant_formulas": formulas,
            "rejection_reason_ar": (
                f"تفاعل غير موثق في السجل المعتمد: لا يتوفر تفاعل كيميائي موثق لخلط ({joined}). "
                "تلتزم منصة MolecuLens بالدقة العلمية الصارمة وترفض تخمين أو اختلاق نواتج أو معادلات غير موثقة. "
                "الحل المقترح: اختر متفاعلات موثقة من التجارب المرجعية أو استعرض التفاعلات المتاحة لكل مادة."
            ),
        }

    r_specs = [dict(r) for r in rxn["reactants"]]
    p_specs = [dict(p) for p in rxn["products"]]

    alg_balance = balance_reaction_algebraically(
        [r["formula"] for r in r_specs],
        [p["formula"] for p in p_specs],
    )
    if alg_balance is not None:
        r_coeffs, p_coeffs = alg_balance
        for idx, c in enumerate(r_coeffs):
            r_specs[idx]["coefficient"] = c
        for idx, c in enumerate(p_coeffs):
            p_specs[idx]["coefficient"] = c

    conservation = verify_atom_conservation(r_specs, p_specs)
    eq_formatted = format_balanced_equation(r_specs, p_specs)
    eq_unicode = format_chemical_equation_display(eq_formatted["text"])
    stoich = compute_reaction_stoichiometry(r_specs, reactant_inputs)
    yields = compute_theoretical_yields(p_specs, stoich["reaction_extent_mol"])

    cond_compat = evaluate_reaction_conditions_compatibility(
        rxn.get("conditions", {}),
        custom_conditions,
        category=rxn.get("category", ""),
    )

    ratio_str = ":".join(str(r["coefficient"]) for r in r_specs)
    r_summary = " مع ".join(f"{r['name_ar']} ({format_chemical_formula_display(r['formula'])})" for r in r_specs)
    p_summary = " و ".join(
        f"{y['theoretical_yield_mol']:g} mol ({y['theoretical_yield_g']:g} g) من {y['name_ar']} ({format_chemical_formula_display(y['formula'])})"
        for y in yields
    )

    if stoich["is_exact_stoichiometric_ratio"]:
        stoich_status_ar = "حساب نظري بنسبة مولية متكافئة تماماً (Exact Stoichiometric Ratio)"
        dynamic_explanation_ar = (
            f"تفاعلت كمية {r_summary} وفق النسبة المولية ({ratio_str})، وتكوّن نظرياً (Theoretical Yield) مقدار {p_summary}. "
            f"{rxn['explanation_ar']} هذه نتيجة حسابية نظرية وليست إثباتاً لتجربة عملية."
        )
    else:
        lim = stoich["limiting_reagent"]
        lim_name = f"{lim['name_ar']} ({format_chemical_formula_display(lim['formula'])})" if lim else ""
        stoich_status_ar = f"حساب نظري محكوم بالمتفاعل المحدد (Limiting Reagent): {lim_name}"
        dynamic_explanation_ar = (
            f"تفاعلت المواد وفق النسبة المولية ({ratio_str})، وكان المتفاعل المحدد (Limiting Reagent) للتفاعل هو {lim_name} حيث استُهلك بالكامل "
            f"وتوقف التفاعل عند انتهائه، مما أدى إلى تكوّن ناتج نظري (Theoretical Yield) قدره {p_summary}. "
            "هذه نتيجة حسابية نظرية وليست إثباتاً لتجربة عملية."
        )

    merged_conditions = dict(rxn.get("conditions", {}))
    if custom_conditions:
        merged_conditions["user_specified_conditions"] = custom_conditions
    merged_conditions["compatibility_check"] = cond_compat

    # Find all documented reactions for these reactants across all types/domains to show alternatives
    all_for_reactants = find_verified_reactions(formulas, reaction_type=None, domain=None)
    alternatives = [
        {
            "reaction_id": c["reaction_id"],
            "reaction_name_ar": c["reaction_name_ar"],
            "category": c.get("category", ""),
            "domain": c.get("domain", ""),
        }
        for c in all_for_reactants
        if c["reaction_id"] != rxn["reaction_id"]
    ]
    return {
        "supported": True,
        "reaction_id": rxn["reaction_id"],
        "category": rxn.get("category", ""),
        "alternative_documented_reactions": alternatives,
        "conditions_compatibility": cond_compat,
        "documented_reaction_info": {
            "label_ar": "التفاعل الموثق",
            "label_bilingual": "التفاعل الموثق (Documented Reaction)",
            "reaction_name_ar": rxn["reaction_name_ar"],
            "reaction_type_ar": rxn["reaction_type_ar"],
            "source": rxn.get("source", "Verified Chemical Reaction Registry"),
            "conditions": dict(rxn.get("conditions", {})),
            "equation_text": eq_formatted["text"],
            "equation_unicode": eq_unicode,
        },
        "calculated_stoichiometry_result": {
            "label_ar": "النتيجة المحسوبة",
            "label_bilingual": "النتيجة المحسوبة (Calculated Result)",
            "limiting_reagent": stoich["limiting_reagent"],
            "theoretical_yields": yields,
            "is_exact_stoichiometric_ratio": stoich["is_exact_stoichiometric_ratio"],
            "atom_conservation_verified": conservation["is_balanced"],
        },
        "theoretical_prediction_notes": {
            "label_ar": "التوقع النظري",
            "label_bilingual": "التوقع النظري والافتراضات (Theoretical Prediction)",
            "notes_ar": rxn.get("assumptions_and_limits_ar", ""),
            "conditions_compatibility_ar": cond_compat["details_ar"],
        },
        "domain": rxn["domain"],
        "domain_label_ar": rxn["domain_label_ar"],
        "reaction_name_ar": rxn["reaction_name_ar"],
        "reaction_name_en": rxn["reaction_name_en"],
        "reaction_type_ar": rxn["reaction_type_ar"],
        "equation_text": eq_formatted["text"],
        "equation_unicode": eq_unicode,
        "equation_latex": eq_formatted["latex"],
        "atom_conservation": conservation,
        "stoichiometric_status_ar": stoich_status_ar,
        "is_exact_stoichiometric_ratio": stoich["is_exact_stoichiometric_ratio"],
        "limiting_reagent": stoich["limiting_reagent"],
        "excess_reagents": stoich["excess_reagents"],
        "reactants_table": stoich["reactants_table"] if "reactants_table" in stoich else stoich["reactants_breakdown"],
        "theoretical_yields": yields,
        "conditions": merged_conditions,
        "explanation_ar": dynamic_explanation_ar,
        "assumptions_and_limits_ar": rxn.get("assumptions_and_limits_ar", "هذه نتيجة حسابية نظرية وليست إثباتا لتجربة عملية."),
        "retrieval_provenance": build_provenance_tag(
            OP_SCIENTIFIC_RETRIEVAL,
            sources=[rxn.get("source", "Verified Chemical Reaction Registry")],
            notes_ar="نوع التفاعل والنواتج وظروف التفاعل مسترجعة من السجل العلمي الموثق.",
        ),
        "calculation_provenance": build_provenance_tag(
            OP_SCIENTIFIC_CALCULATION,
            sources=["Algebraic Equation Balancer", "Stoichiometry & Yield Engine"],
            notes_ar="معاملات الموازنة وحفظ الذرات والمتفاعل المحدد (Limiting Reagent) والناتج النظري (Theoretical Yield) محسوبة رياضياً بدقة.",
        ),
        "prediction_provenance": build_provenance_tag(
            OP_SCIENTIFIC_PREDICTION,
            sources=["Thermodynamic & Kinetic Condition Evaluator"],
            notes_ar=cond_compat["details_ar"],
        ),
    }

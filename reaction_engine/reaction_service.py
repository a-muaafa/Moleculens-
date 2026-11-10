from typing import Any
from chemistry_engine.periodic_table import ELEMENT_NAMES_AR
from config.settings import OP_SCIENTIFIC_CALCULATION, OP_SCIENTIFIC_RETRIEVAL
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


def analyze_chemical_reaction(
    reactant_inputs: list[dict[str, Any]],
    domain: str = "general_chemistry",
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
                    f"رفض علمي: المادة المدخلة ({el_ar} - {f}) تنتمي إلى الغازات النبيلة الخاملة ذات الغلاف الإلكتروني المكتمل، "
                    "ولا تدخل في تفاعلات كيميائية عادية تحت الظروف القياسية. يلتزم النظام بعدم اختلاق نواتج غير حقيقية."
                ),
            }

    rtype = (custom_conditions or {}).get("reaction_type")
    candidates = find_verified_reactions(formulas, reaction_type=rtype)
    rxn = candidates[0] if candidates else None
    if rxn is None and rtype and rtype != "all":
        any_type = find_verified_reactions(formulas)
        if any_type:
            return {
                "supported": False,
                "error_type": "reaction_type_mismatch",
                "reactant_formulas": formulas,
                "rejection_reason_ar": (
                    "يوجد تفاعل موثق لهذه المواد لكنه من نوع آخر غير النوع المختار. "
                    "اختر خيار (تلقائي بحسب المواد المتفاعلة) أو النوع المناسب: "
                    + "، ".join(CATEGORY_LABELS_AR.get(c.get("category", ""), c.get("category", "")) for c in any_type)
                ),
            }
    if rxn is None:
        joined = " + ".join(f"{n} ({f})" for n, f in zip(names_ar, formulas))
        return {
            "supported": False,
            "error_type": "unsupported_reaction",
            "reactant_formulas": formulas,
            "rejection_reason_ar": (
                f"لا يتوفر تفاعل كيميائي موثق في قاعدة التفاعلات المعتمدة لخلط ({joined}). "
                "تلتزم منصة MolecuLens بالدقة العلمية وترفض تخمين أو اختلاق نواتج أو معادلات غير موثقة."
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
    stoich = compute_reaction_stoichiometry(r_specs, reactant_inputs)
    yields = compute_theoretical_yields(p_specs, stoich["reaction_extent_mol"])

    ratio_str = ":".join(str(r["coefficient"]) for r in r_specs)
    r_summary = " مع ".join(f"{r['name_ar']} ({r['formula']})" for r in r_specs)
    p_summary = " و ".join(
        f"{y['theoretical_yield_mol']:g} mol ({y['theoretical_yield_g']:g} g) من {y['name_ar']} ({y['formula']})"
        for y in yields
    )

    if stoich["is_exact_stoichiometric_ratio"]:
        stoich_status_ar = "حساب نظري بنسبة مولية متكافئة تماما (تُستهلك جميع المتفاعلات بالكامل دون فائض)"
        dynamic_explanation_ar = (
            f"تفاعلت كمية {r_summary} وفق النسبة المولية ({ratio_str})، وتكوّن نظريا {p_summary}. "
            f"{rxn['explanation_ar']} هذه نتيجة حسابية نظرية وليست إثباتا لتجربة عملية."
        )
    else:
        lim = stoich["limiting_reagent"]
        lim_name = f"{lim['name_ar']} ({lim['formula']})" if lim else ""
        stoich_status_ar = f"حساب نظري محكوم بالمتفاعل المحدد: {lim_name}"
        dynamic_explanation_ar = (
            f"تفاعلت المواد وفق النسبة المولية ({ratio_str})، وكان المتفاعل المحدد للتفاعل هو {lim_name} حيث استُهلك بالكامل "
            f"وتوقف التفاعل عند انتهائه، مما أدى إلى تكوّن نظري قدره {p_summary}. "
            "هذه نتيجة حسابية نظرية وليست إثباتا لتجربة عملية."
        )

    merged_conditions = dict(rxn.get("conditions", {}))
    if custom_conditions:
        merged_conditions["user_specified_conditions"] = custom_conditions

    alternatives = [
        {"reaction_id": c["reaction_id"], "reaction_name_ar": c["reaction_name_ar"], "category": c.get("category", "")}
        for c in candidates[1:]
    ]
    return {
        "supported": True,
        "reaction_id": rxn["reaction_id"],
        "category": rxn.get("category", ""),
        "alternative_documented_reactions": alternatives,
        "documented_reaction_info": {
            "label_ar": "التفاعل الموثق",
            "reaction_name_ar": rxn["reaction_name_ar"],
            "reaction_type_ar": rxn["reaction_type_ar"],
            "source": rxn.get("source", "Verified Chemical Reaction Registry"),
            "conditions": dict(rxn.get("conditions", {})),
            "equation_text": eq_formatted["text"],
        },
        "calculated_stoichiometry_result": {
            "label_ar": "النتيجة المحسوبة",
            "limiting_reagent": stoich["limiting_reagent"],
            "theoretical_yields": yields,
            "is_exact_stoichiometric_ratio": stoich["is_exact_stoichiometric_ratio"],
            "atom_conservation_verified": conservation["is_balanced"],
        },
        "theoretical_prediction_notes": {
            "label_ar": "التوقع النظري",
            "notes_ar": rxn.get("assumptions_and_limits_ar", ""),
        },
        "domain": rxn["domain"],
        "domain_label_ar": rxn["domain_label_ar"],
        "reaction_name_ar": rxn["reaction_name_ar"],
        "reaction_name_en": rxn["reaction_name_en"],
        "reaction_type_ar": rxn["reaction_type_ar"],
        "equation_text": eq_formatted["text"],
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
            sources=["Verified Chemical Reaction Registry"],
            notes_ar="نوع التفاعل والنواتج وظروف التفاعل مسترجعة من السجل العلمي الموثق.",
        ),
        "calculation_provenance": build_provenance_tag(
            OP_SCIENTIFIC_CALCULATION,
            sources=["Algebraic Equation Balancer", "Stoichiometry & Yield Engine"],
            notes_ar="معاملات الموازنة وحفظ الذرات والمتفاعل المحدد والناتج النظري محسوبة رياضيا بدقة.",
        ),
    }


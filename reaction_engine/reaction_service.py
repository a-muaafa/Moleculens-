import re
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
from reaction_engine.reaction_registry import (
    CATEGORY_LABELS_AR,
    DOCUMENTATION_TYPE_LABELS_AR,
    find_verified_reactions,
)
from reaction_engine.stoichiometry import compute_reaction_stoichiometry
from reaction_engine.yield_calculator import compute_theoretical_yields

NOBLE_GASES = {"He", "Ne", "Ar", "Kr", "Xe", "Rn"}
DOMAIN_LABELS_DICT = dict(REACTION_DOMAINS)

SOLVENT_LABELS_AR: dict[str, str] = {
    "aqueous": "محلول مائي (Aqueous Solution)",
    "neat_gas_solid": "طور غازي/صلب مباشر بدون مذيب (Solvent-Free / Gas-Solid)",
    "ethanol": "مذيب عضوي قطبي — إيثانول (Ethanol)",
    "physiological_buffer": "وسط منظم فسيولوجي (Physiological Buffer pH 7.4)",
}

CATALYST_LABELS_AR: dict[str, str] = {
    "none": "بدون عامل حفاز إضافي (No Catalyst)",
    "acid_h2so4": "حفاز حمض الكبريتيك المركز (H₂SO₄ Catalyst)",
    "pt_pd": "حفاز فلزي انتقالي — بلاتين أو بالاديوم (Pt / Pd)",
    "fe_haber": "حفاز الحديد المنشط (Promoted Iron Catalyst - Fe)",
    "mno2_catalase": "ثاني أكسيد المنغنيز أو إنزيم الكاتالاز (MnO₂ / Catalase)",
    "enzymatic": "إنزيم حيوي متخصص (Enzymatic Catalysis)",
}


def parse_temperature_to_celsius(temp_input: Any) -> float | None:
    """
    Parses a temperature specification (float, int, dict with value/unit, or string in °C, K, °F)
    and converts it quantitatively into degrees Celsius (°C).
    """
    if temp_input is None:
        return None
    if isinstance(temp_input, (int, float)):
        return float(temp_input)
    if isinstance(temp_input, dict):
        val = temp_input.get("value")
        unit = str(temp_input.get("unit", "C")).strip().upper()
        if val is None:
            return None
        v = float(val)
        if unit in {"K", "KELVIN"}:
            return round(v - 273.15, 3)
        if unit in {"F", "°F", "FAHRENHEIT"}:
            return round((v - 32.0) * 5.0 / 9.0, 3)
        return round(v, 3)

    s = str(temp_input).strip()
    if not s:
        return None
    # Match explicit number + unit (°C, C, K, °F, F)
    m_c = re.search(r"([+\-]?\d+(?:\.\d+)?)\s*°?\s*[Cc]\b", s)
    if m_c:
        return round(float(m_c.group(1)), 3)
    m_k = re.search(r"([+\-]?\d+(?:\.\d+)?)\s*[Kk]\b", s)
    if m_k:
        return round(float(m_k.group(1)) - 273.15, 3)
    m_f = re.search(r"([+\-]?\d+(?:\.\d+)?)\s*°?\s*[Ff]\b", s)
    if m_f:
        return round((float(m_f.group(1)) - 32.0) * 5.0 / 9.0, 3)
    m_num = re.fullmatch(r"([+\-]?\d+(?:\.\d+)?)", s)
    if m_num:
        return round(float(m_num.group(1)), 3)
    return None


def parse_pressure_to_atm(press_input: Any) -> float | None:
    """
    Parses a pressure specification (float, int, dict with value/unit, or string in atm, kPa, bar, Pa, mmHg, psi)
    and converts it quantitatively into standard atmospheres (atm).
    """
    if press_input is None:
        return None
    if isinstance(press_input, (int, float)):
        return float(press_input)
    if isinstance(press_input, dict):
        val = press_input.get("value")
        unit = str(press_input.get("unit", "atm")).strip().lower()
        if val is None:
            return None
        v = float(val)
        if unit == "kpa":
            return round(v / 101.325, 4)
        if unit == "pa":
            return round(v / 101325.0, 4)
        if unit == "bar":
            return round(v / 1.01325, 4)
        if unit in {"mmhg", "torr"}:
            return round(v / 760.0, 4)
        if unit == "psi":
            return round(v / 14.6959, 4)
        return round(v, 4)

    s = str(press_input).strip()
    if not s:
        return None
    m_atm = re.search(r"([+\-]?\d+(?:\.\d+)?)\s*atm\b", s, flags=re.IGNORECASE)
    if m_atm:
        return round(float(m_atm.group(1)), 4)
    m_kpa = re.search(r"([+\-]?\d+(?:\.\d+)?)\s*kpa\b", s, flags=re.IGNORECASE)
    if m_kpa:
        return round(float(m_kpa.group(1)) / 101.325, 4)
    m_bar = re.search(r"([+\-]?\d+(?:\.\d+)?)\s*bar\b", s, flags=re.IGNORECASE)
    if m_bar:
        return round(float(m_bar.group(1)) / 1.01325, 4)
    m_mmhg = re.search(r"([+\-]?\d+(?:\.\d+)?)\s*(?:mmhg|torr)\b", s, flags=re.IGNORECASE)
    if m_mmhg:
        return round(float(m_mmhg.group(1)) / 760.0, 4)
    m_psi = re.search(r"([+\-]?\d+(?:\.\d+)?)\s*psi\b", s, flags=re.IGNORECASE)
    if m_psi:
        return round(float(m_psi.group(1)) / 14.6959, 4)
    m_num = re.fullmatch(r"([+\-]?\d+(?:\.\d+)?)", s)
    if m_num:
        return round(float(m_num.group(1)), 4)
    return None


def normalize_solvent_code(solv_input: Any) -> str | None:
    if not solv_input:
        return None
    s = str(solv_input).strip().lower()
    if not s:
        return None
    if s in {"aqueous", "neat_gas_solid", "ethanol", "physiological_buffer"}:
        return s
    if "buffer" in s or "فسيولوجي" in s or "7.4" in s:
        return "physiological_buffer"
    if "محلول مائي" in s or "aqueous" in s or "h2o" in s:
        return "aqueous"
    if "إيثانول" in s or "ethanol" in s:
        return "ethanol"
    if "بدون مذيب" in s or "gas-solid" in s or "solvent-free" in s or "طور غازي" in s:
        return "neat_gas_solid"
    return s


def normalize_catalyst_code(cat_input: Any) -> str | None:
    if not cat_input:
        return None
    s = str(cat_input).strip().lower()
    if not s:
        return None
    if s in {"none", "acid_h2so4", "pt_pd", "fe_haber", "mno2_catalase", "enzymatic"}:
        return s
    if "بدون عامل حفاز" in s or "no catalyst" in s:
        return "none"
    if "h2so4" in s or "حمض الكبريتيك" in s or "حمضي" in s:
        return "acid_h2so4"
    if "pt" in s or "pd" in s or "بلاتين" in s or "بالاديوم" in s or "لندلار" in s:
        return "pt_pd"
    if "fe" in s or "الحديد" in s:
        return "fe_haber"
    if "mno2" in s or "الكاتالاز" in s or "catalase" in s:
        return "mno2_catalase"
    if "إنزيم" in s or "enzymatic" in s or "rhea" in s:
        return "enzymatic"
    return s


def _infer_fallback_reference_spec(documented_conditions: dict[str, Any], category: str) -> dict[str, Any]:
    """Derives a structured reference conditions spec if evaluate_reaction_conditions_compatibility is called without an explicit spec."""
    doc_temp = str(documented_conditions.get("temperature", ""))
    doc_press = str(documented_conditions.get("pressure", ""))
    doc_solv = str(documented_conditions.get("solvent", ""))
    doc_cat = str(documented_conditions.get("catalyst", ""))
    combined = f"{doc_temp} {doc_press} {doc_solv} {doc_cat}".lower()

    if category == "biochemical_reaction" or "37" in doc_temp:
        return {
            "temp_c_min": 15.0,
            "temp_c_max": 45.0,
            "temp_c_optimal": 37.0,
            "pressure_atm_min": 0.8,
            "pressure_atm_max": 2.0,
            "pressure_atm_optimal": 1.0,
            "allowed_solvents": ["aqueous", "physiological_buffer"],
            "required_catalysts": ["enzymatic"],
            "catalyst_is_mandatory": True,
        }
    if any(k in combined for k in ["850", "500", "400", "فرن", "شرارة", "اشتعال"]) or category == "combustion":
        p_min = 150.0 if any(k in doc_press for k in ["150", "200", "250"]) else 0.5
        p_max = 300.0 if p_min >= 150.0 else 25.0
        req_cat = ["fe_haber"] if ("fe" in combined and "هابر" in combined) else (["none", "pt_pd"] if "pt" in combined else ["none"])
        return {
            "temp_c_min": 380.0,
            "temp_c_max": 2500.0,
            "temp_c_optimal": 550.0,
            "pressure_atm_min": p_min,
            "pressure_atm_max": p_max,
            "pressure_atm_optimal": 200.0 if p_min >= 150.0 else 1.0,
            "allowed_solvents": ["neat_gas_solid"],
            "required_catalysts": req_cat,
            "catalyst_is_mandatory": "fe_haber" in req_cat,
            "catalyst_lowers_temp_c_min": 20.0 if "pt" in combined else None,
        }
    if "h2so4" in combined and ("أسترة" in combined or "نترتة" in combined):
        return {
            "temp_c_min": 45.0,
            "temp_c_max": 130.0,
            "temp_c_optimal": 80.0,
            "pressure_atm_min": 0.8,
            "pressure_atm_max": 5.0,
            "pressure_atm_optimal": 1.0,
            "allowed_solvents": ["neat_gas_solid", "ethanol", "aqueous"],
            "required_catalysts": ["acid_h2so4"],
            "catalyst_is_mandatory": True,
        }
    return {
        "temp_c_min": 5.0,
        "temp_c_max": 80.0,
        "temp_c_optimal": 25.0,
        "pressure_atm_min": 0.5,
        "pressure_atm_max": 5.0,
        "pressure_atm_optimal": 1.0,
        "allowed_solvents": ["aqueous"] if ("محلول مائي" in combined or category in {"neutralization", "precipitation"}) else ["aqueous", "neat_gas_solid", "ethanol"],
        "required_catalysts": ["none"],
        "catalyst_is_mandatory": False,
    }


def evaluate_reaction_conditions_compatibility(
    documented_conditions: dict[str, Any],
    user_conditions: dict[str, Any] | None,
    category: str = "",
    reference_conditions_spec: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Quantitatively compares user-specified reaction conditions (temperature in °C/K/°F, pressure in atm/kPa/bar,
    solvent medium code, and catalyst code) against structured reference bounds.
    """
    spec = dict(reference_conditions_spec or _infer_fallback_reference_spec(documented_conditions, category))

    if not user_conditions:
        return {
            "is_compatible": True,
            "status_ar": "متوافق مع الظروف المرجعية القياسية (Compatible with Reference Conditions)",
            "warnings_ar": [],
            "details_ar": "تُطبق الظروف المرجعية المسجلة للتفاعل دون تعارض.",
            "reference_bounds": spec,
            "parsed_user_conditions": {},
            "parameter_checks": {},
        }

    u_temp_raw = user_conditions.get("temperature")
    u_press_raw = user_conditions.get("pressure")
    u_solv_raw = user_conditions.get("solvent")
    u_cat_raw = user_conditions.get("catalyst")

    u_temp_c = parse_temperature_to_celsius(u_temp_raw)
    u_press_atm = parse_pressure_to_atm(u_press_raw)
    u_solv_code = normalize_solvent_code(u_solv_raw)
    u_cat_code = normalize_catalyst_code(u_cat_raw)

    warnings: list[str] = []
    parameter_checks: dict[str, dict[str, Any]] = {}

    # 1. Quantitative Temperature Bounds Comparison
    t_min = float(spec.get("temp_c_min", 0.0))
    t_max = float(spec.get("temp_c_max", 2000.0))
    if u_cat_code == "pt_pd" and spec.get("catalyst_lowers_temp_c_min") is not None:
        t_min = float(spec["catalyst_lowers_temp_c_min"])

    if u_temp_c is not None:
        temp_ok = t_min <= u_temp_c <= t_max
        parameter_checks["temperature"] = {
            "ok": temp_ok,
            "user_value_c": u_temp_c,
            "min_c": t_min,
            "max_c": t_max,
            "optimal_c": spec.get("temp_c_optimal"),
        }
        if u_temp_c < t_min:
            warnings.append(
                f"تنبيه درجة الحرارة (Temperature Below Minimum Bound): القيمة المدخلة ({u_temp_c:g} °C) أقل من الحد الأدنى المرجعي "
                f"لبدء التفاعل ({t_min:g} °C؛ النطاق المرجعي المسموح: [{t_min:g} – {t_max:g} °C])."
            )
        elif u_temp_c > t_max:
            if category == "biochemical_reaction" or t_max <= 50.0:
                warnings.append(
                    f"تنبيه حيوي (Thermal Denaturation Above Maximum Bound): درجة الحرارة المدخلة ({u_temp_c:g} °C) تتجاوز الحد الأعلى الفسيولوجي "
                    f"({t_max:g} °C؛ النطاق المسموح: [{t_min:g} – {t_max:g} °C]) مما يؤدي إلى تمسخ الإنزيمات وتوقف التفاعل."
                )
            else:
                warnings.append(
                    f"تنبيه درجة الحرارة (Temperature Exceeds Maximum Bound): القيمة المدخلة ({u_temp_c:g} °C) تتجاوز الحد الأعلى المرجعي "
                    f"({t_max:g} °C؛ النطاق المسموح: [{t_min:g} – {t_max:g} °C]) مما قد يسبب تفكك النواتج أو تفاعلات جانبية."
                )

    # 2. Quantitative Pressure Bounds Comparison
    p_min = float(spec.get("pressure_atm_min", 0.1))
    p_max = float(spec.get("pressure_atm_max", 500.0))
    if u_press_atm is not None:
        press_ok = p_min <= u_press_atm <= p_max
        parameter_checks["pressure"] = {
            "ok": press_ok,
            "user_value_atm": u_press_atm,
            "min_atm": p_min,
            "max_atm": p_max,
            "optimal_atm": spec.get("pressure_atm_optimal"),
        }
        if u_press_atm < p_min:
            warnings.append(
                f"تنبيه الضغط (Pressure Below Minimum Bound): الضغط المدخل ({u_press_atm:g} atm) أقل من الحد الأدنى المطلوب "
                f"({p_min:g} atm؛ النطاق المرجعي المسموح: [{p_min:g} – {p_max:g} atm])."
            )
        elif u_press_atm > p_max:
            warnings.append(
                f"تنبيه الضغط (Pressure Exceeds Maximum Bound): الضغط المدخل ({u_press_atm:g} atm) يتجاوز الحد الأعلى المرجعي "
                f"({p_max:g} atm؛ النطاق المسموح: [{p_min:g} – {p_max:g} atm])."
            )

    # 3. Structured Solvent / Medium Compatibility Check
    allowed_solvents: list[str] = list(spec.get("allowed_solvents") or [])
    if u_solv_code and allowed_solvents:
        solv_ok = u_solv_code in allowed_solvents
        parameter_checks["solvent"] = {
            "ok": solv_ok,
            "user_solvent_code": u_solv_code,
            "allowed_solvents": allowed_solvents,
        }
        if not solv_ok:
            allowed_ar = " أو ".join(SOLVENT_LABELS_AR.get(c, c) for c in allowed_solvents)
            chosen_ar = SOLVENT_LABELS_AR.get(u_solv_code, str(u_solv_raw))
            warnings.append(
                f"تنبيه الوسط/المذيب (Solvent Medium Mismatch): الوسط المختار ({chosen_ar}) غير متوافق مع طور التفاعل المرجعي؛ "
                f"الأوساط المسموحة: ({allowed_ar})."
            )

    # 4. Structured Catalyst Requirement Check
    req_catalysts: list[str] = list(spec.get("required_catalysts") or ["none"])
    cat_mandatory = bool(spec.get("catalyst_is_mandatory", False))
    if u_cat_code:
        if cat_mandatory:
            cat_ok = u_cat_code in req_catalysts
        else:
            cat_ok = True
        parameter_checks["catalyst"] = {
            "ok": cat_ok,
            "user_catalyst_code": u_cat_code,
            "required_catalysts": req_catalysts,
            "catalyst_is_mandatory": cat_mandatory,
        }
        if not cat_ok:
            req_ar = " أو ".join(CATALYST_LABELS_AR.get(c, c) for c in req_catalysts)
            chosen_ar = CATALYST_LABELS_AR.get(u_cat_code, str(u_cat_raw))
            warnings.append(
                f"تنبيه العامل الحفاز (Catalyst Requirement Mismatch): يتطلب هذا التفاعل حفازاً إلزامياً ({req_ar})، "
                f"بينما تم اختيار ({chosen_ar})."
            )

    is_compat = len(warnings) == 0
    return {
        "is_compatible": is_compat,
        "status_ar": (
            "متوافق مع النطاقات والظروف المرجعية للتفاعل (Compatible Conditions)"
            if is_compat
            else "يوجد تعارض كمي أو نوعي بين الظروف المختارة والنطاق المرجعي للتفاعل (Condition Mismatch Warning)"
        ),
        "warnings_ar": warnings,
        "details_ar": (
            "جميع الظروف المدخلة (درجة الحرارة، الضغط، المذيب، الحفاز) تقع ضمن الحدود المرجعية المسموحة للتفاعل."
            if is_compat
            else " | ".join(warnings)
        ),
        "reference_bounds": spec,
        "parsed_user_conditions": {
            "temperature_c": u_temp_c,
            "pressure_atm": u_press_atm,
            "solvent_code": u_solv_code,
            "catalyst_code": u_cat_code,
        },
        "parameter_checks": parameter_checks,
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
        reference_conditions_spec=rxn.get("reference_conditions_spec"),
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
    merged_conditions["reference_conditions_spec"] = rxn.get("reference_conditions_spec", {})
    if custom_conditions:
        merged_conditions["user_specified_conditions"] = custom_conditions
    merged_conditions["compatibility_check"] = cond_compat

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

    doc_type = rxn.get("documentation_type", "experimentally_documented")
    doc_label_ar = rxn.get("documentation_label_ar") or DOCUMENTATION_TYPE_LABELS_AR.get(doc_type, "")

    return {
        "supported": True,
        "reaction_id": rxn["reaction_id"],
        "category": rxn.get("category", ""),
        "documentation_type": doc_type,
        "documentation_label_ar": doc_label_ar,
        "source_id": rxn.get("source_id", ""),
        "reference_url": rxn.get("reference_url", ""),
        "alternative_documented_reactions": alternatives,
        "conditions_compatibility": cond_compat,
        "documented_reaction_info": {
            "label_ar": "التفاعل الموثق",
            "label_bilingual": "التفاعل الموثق (Documented Reaction)",
            "reaction_name_ar": rxn["reaction_name_ar"],
            "reaction_type_ar": rxn["reaction_type_ar"],
            "documentation_type": doc_type,
            "documentation_label_ar": doc_label_ar,
            "source_id": rxn.get("source_id", ""),
            "source": rxn.get("source", ""),
            "reference_url": rxn.get("reference_url", ""),
            "conditions": dict(rxn.get("conditions", {})),
            "reference_conditions_spec": dict(rxn.get("reference_conditions_spec", {})),
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
            sources=[rxn.get("source", "")],
            notes_ar=f"المعادلة والظروف مصنفة كـ ({doc_label_ar}) ومسترجعة بمعرف مرجعي ({rxn.get('source_id', '')}).",
        ),
        "calculation_provenance": build_provenance_tag(
            OP_SCIENTIFIC_CALCULATION,
            sources=["Algebraic Equation Balancer", "Stoichiometry & Yield Engine"],
            notes_ar="معاملات الموازنة وحفظ الذرات والمتفاعل المحدد (Limiting Reagent) والناتج النظري (Theoretical Yield) محسوبة رياضياً بدقة.",
        ),
        "prediction_provenance": build_provenance_tag(
            OP_SCIENTIFIC_PREDICTION,
            sources=["Quantitative Thermodynamic & Kinetic Condition Evaluator"],
            notes_ar=cond_compat["details_ar"],
        ),
    }

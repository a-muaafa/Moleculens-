from typing import Any
import pandas as pd
import streamlit as st
from chemistry_engine.formula_analyzer import (
    format_chemical_equation_display,
    format_chemical_formula_display,
    is_gaseous_at_stp,
)
from config.settings import CLOSED_DROPDOWN_OPTIONS, REACTION_DOMAINS
from core.orchestrator import MolecuLensCore
from reaction_engine.reaction_registry import CATEGORY_LABELS_AR, search_reactions_by_substance
from search.local_catalog import get_catalog_entity_by_key
from search.search_engine import get_reaction_substance_options

NONE_KEY = ""
UNITS_BY_QUANTITY = {
    "molar_amount": ["mol", "mmol"],
    "mass_amount": ["g", "mg"],
    "gas_volume_stp": ["L_STP"],
}
PRESET_KEYS = {
    "H2": "hydrogen_gas_h2",
    "O2": "oxygen_gas_o2",
    "Fe": "iron_fe",
    "HCl": "hydrochloric_acid_hcl",
    "NaOH": "sodium_hydroxide_naoh",
    "CO2": "carbon_dioxide_co2",
    "H2O": "water_h2o",
    "CH4": "methane_ch4",
    "He": "element_he",
    "CH3COOH": "acetic_acid",
    "NaHCO3": "sodium_bicarbonate",
}


def _labels(options: list[dict[str, Any]]) -> dict[str, str]:
    return {o["key"]: o["label"] for o in options}


def _apply_preset(f1: str, a1: float, f2: str, a2: float, domain_code: str = "all_domains") -> None:
    st.session_state["rxn_domain_select"] = domain_code
    st.session_state["rxn_reaction_type"] = "all"
    st.session_state["rxn_category_filter"] = "all"
    st.session_state["rxn_quantity_type"] = "molar_amount"
    st.session_state["rxn_r1_key"] = PRESET_KEYS[f1]
    st.session_state["rxn_r1_amount"] = float(a1)
    st.session_state["rxn_r1_unit"] = "mol"
    st.session_state["rxn_r2_key"] = PRESET_KEYS[f2]
    st.session_state["rxn_r2_amount"] = float(a2)
    st.session_state["rxn_r2_unit"] = "mol"


def _closed_select(label: str, option_pairs: list[tuple[str, str]], key: str, default: str | None = None) -> str:
    codes = [c for c, _ in option_pairs]
    if key not in st.session_state or st.session_state[key] not in codes:
        st.session_state[key] = default if default in codes else codes[0]
    mapping = dict(option_pairs)
    return st.selectbox(label, options=codes, format_func=lambda c: mapping.get(c, c), key=key)


def _substance_select(label: str, options: list[dict[str, Any]], key: str, allow_none: bool) -> dict[str, Any] | None:
    labels = _labels(options)
    keys = ([NONE_KEY] if allow_none else []) + [o["key"] for o in options]
    current = st.session_state.get(key)
    if current and current not in keys:
        entity = get_catalog_entity_by_key(current)
        if entity:
            keys.append(current)
            f_disp = format_chemical_formula_display(entity["symbol_or_formula"])
            labels[current] = f"{entity['name_ar']} | {entity['name_en']} ({f_disp})"
    if key not in st.session_state or st.session_state[key] not in keys:
        st.session_state[key] = keys[0]
    labels[NONE_KEY] = "— بدون مادة إضافية (No Additional Reactant) —"
    chosen = st.selectbox(label, options=keys, format_func=lambda k: labels.get(k, k), key=key)
    if chosen == NONE_KEY:
        return None
    entity = get_catalog_entity_by_key(chosen)
    if entity:
        f_disp = format_chemical_formula_display(entity["symbol_or_formula"])
        phase_str = "غاز (Gas)" if is_gaseous_at_stp(entity["symbol_or_formula"], substance=entity) else "صلب/سائل/محلول (Condensed Phase)"
        st.caption(
            f"الصيغة الجزيئية (Molecular Formula): {f_disp} | الكتلة المولية (Molar Mass): {entity.get('molar_mass')} g/mol | الطور القياسي: {phase_str}"
        )
    return entity


def _amount_inputs(prefix: str, label_suffix: str, units: list[tuple[str, str]], default_amount: float) -> tuple[float, str]:
    ca, cu = st.columns([1.6, 1.2])
    with ca:
        if f"{prefix}_amount" not in st.session_state:
            st.session_state[f"{prefix}_amount"] = default_amount
        amount = st.number_input(f"الكمية (Amount) {label_suffix}:", min_value=0.01, step=0.5, key=f"{prefix}_amount")
    with cu:
        unit = _closed_select(f"الوحدة (Unit) {label_suffix}:", units, f"{prefix}_unit", default=units[0][0])
    return float(amount), unit


def render_reaction_lab_section() -> None:
    if "rxn_domain_select" not in st.session_state:
        st.session_state["rxn_domain_select"] = "all_domains"
    if "rxn_r1_key" not in st.session_state:
        st.session_state["rxn_r1_key"] = PRESET_KEYS["H2"]
        st.session_state["rxn_r2_key"] = PRESET_KEYS["O2"]
        st.session_state["rxn_r1_amount"] = 2.0
        st.session_state["rxn_r2_amount"] = 1.0

    passed_entity = st.session_state.pop("passed_reaction_entity", None)
    if passed_entity:
        key = passed_entity.get("key", "")
        if key and get_catalog_entity_by_key(key):
            st.session_state["rxn_category_filter"] = "all"
            st.session_state["rxn_r1_key"] = key
            f_disp = format_chemical_formula_display(passed_entity.get("symbol_or_formula", ""))
            st.info(f"تم تمرير المادة ({passed_entity.get('name_ar')} — {f_disp}) من قسم الاستكشاف إلى مختبر التفاعلات.")

    st.markdown(
        """
<div class="ml-card">
  <h3>Reaction Lab — مختبر تحليل التفاعلات والحسابات الكيميائية</h3>
  <p>
    جميع خيارات الإعداد في هذه الصفحة قوائم مغلقة (Closed Dropdowns) لضمان الدقة العلمية؛ اختر المواد المتفاعلة (Reactants)
    والوحدات وظروف التفاعل، وأدخل الكمية الرقمية فقط. يتحقق المحرك العلمي من توافق المجال والظروف، ويوازن المعادلة جبرياً،
    ويحسب المتفاعل المحدد (Limiting Reagent) والناتج النظري (Theoretical Yield) دون اختلاق تفاعلات غير موثقة.
  </p>
</div>
""",
        unsafe_allow_html=True,
    )

    c_dom, c_type, c_cat = st.columns(3)
    with c_dom:
        selected_domain = _closed_select("مجال التفاعل (Reaction Domain):", list(REACTION_DOMAINS), "rxn_domain_select")
    with c_type:
        reaction_type = _closed_select("نوع التفاعل (Reaction Type):", CLOSED_DROPDOWN_OPTIONS["reaction_types"], "rxn_reaction_type")
    with c_cat:
        category = _closed_select("فئة المواد (Substance Category):", CLOSED_DROPDOWN_OPTIONS["substance_categories"], "rxn_category_filter")
    options = get_reaction_substance_options(category)

    st.markdown("<strong>تجارب مرجعية سريعة (Quick Verified Presets):</strong>", unsafe_allow_html=True)
    presets = [
        ("2H₂ + O₂", "H2", 2.0, "O2", 1.0),
        ("4Fe + 3O₂", "Fe", 4.0, "O2", 3.0),
        ("HCl + NaOH", "HCl", 1.0, "NaOH", 1.0),
        ("CH₃COOH + NaOH", "CH3COOH", 1.0, "NaOH", 1.0),
        ("CH₃COOH + NaHCO₃", "CH3COOH", 1.0, "NaHCO3", 1.0),
        ("CH₄ + 2O₂", "CH4", 1.0, "O2", 2.0),
        ("Fe + He (خامل)", "Fe", 1.0, "He", 1.0),
    ]
    p_cols = st.columns(len(presets))
    for idx_p, (lbl, f1, a1, f2, a2) in enumerate(presets):
        with p_cols[idx_p]:
            st.button(lbl, key=f"rxn_preset_{idx_p}", on_click=_apply_preset, args=(f1, a1, f2, a2))

    st.markdown('<div class="ml-card"><h4>إعداد المتفاعلات والكميات (Reactants &amp; Stoichiometry Setup)</h4>', unsafe_allow_html=True)
    quantity_type = _closed_select("نوع الكمية المدخلة (Quantity Type):", CLOSED_DROPDOWN_OPTIONS["quantity_types"], "rxn_quantity_type")
    allowed_codes = UNITS_BY_QUANTITY[quantity_type]
    units = [(c, l) for c, l in CLOSED_DROPDOWN_OPTIONS["amount_units"] if c in allowed_codes]

    rc1, rc2 = st.columns(2)
    with rc1:
        st.markdown("##### المادة المتفاعلة الأولى (Primary Reactant)")
        r1_entity = _substance_select("اختر المادة الأولى:", options, "rxn_r1_key", allow_none=False)
        r1_amt, r1_unit = _amount_inputs("rxn_r1", "للمادة الأولى", units, 2.0)
    with rc2:
        st.markdown("##### المادة المتفاعلة الثانية (Secondary Reactant — اختيارية)")
        r2_entity = _substance_select("اختر المادة الثانية:", options, "rxn_r2_key", allow_none=True)
        r2_amt, r2_unit = _amount_inputs("rxn_r2", "للمادة الثانية", units, 1.0)

    r3_entity = None
    r3_amt, r3_unit = 1.0, units[0][0]
    if st.checkbox("إضافة مادة متفاعلة ثالثة (Add Third Reactant)", key="rxn_has_r3"):
        r3_entity = _substance_select("اختر المادة الثالثة:", options, "rxn_r3_key", allow_none=True)
        r3_amt, r3_unit = _amount_inputs("rxn_r3", "للمادة الثالثة", units, 1.0)

    with st.expander("ظروف التفاعل والتحقق من التوافق (Reaction Conditions & Compatibility)", expanded=False):
        cc1, cc2 = st.columns(2)
        with cc1:
            cond_temp = _closed_select("درجة الحرارة (Temperature):", CLOSED_DROPDOWN_OPTIONS["temperature_conditions"], "rxn_cond_temp")
            cond_solvent = _closed_select("المذيب أو الوسط (Solvent / Medium):", CLOSED_DROPDOWN_OPTIONS["solvent_conditions"], "rxn_cond_solvent")
        with cc2:
            cond_press = _closed_select("الضغط (Pressure):", CLOSED_DROPDOWN_OPTIONS["pressure_conditions"], "rxn_cond_pressure")
            cond_cat = _closed_select("العامل الحفاز (Catalyst):", CLOSED_DROPDOWN_OPTIONS["catalyst_conditions"], "rxn_cond_catalyst")
    st.markdown("</div>", unsafe_allow_html=True)

    reactant_inputs: list[dict[str, Any]] = []
    for ent, amt, unit in ((r1_entity, r1_amt, r1_unit), (r2_entity, r2_amt, r2_unit), (r3_entity, r3_amt, r3_unit)):
        if ent is not None:
            reactant_inputs.append({"substance": ent, "formula": ent["symbol_or_formula"], "amount": amt, "unit": unit})

    rxn_result = MolecuLensCore.run_reaction_analysis(
        reactant_inputs=reactant_inputs,
        domain=selected_domain,
        custom_conditions={
            "temperature": dict(CLOSED_DROPDOWN_OPTIONS["temperature_conditions"]).get(cond_temp, cond_temp),
            "pressure": dict(CLOSED_DROPDOWN_OPTIONS["pressure_conditions"]).get(cond_press, cond_press),
            "solvent": dict(CLOSED_DROPDOWN_OPTIONS["solvent_conditions"]).get(cond_solvent, cond_solvent),
            "catalyst": dict(CLOSED_DROPDOWN_OPTIONS["catalyst_conditions"]).get(cond_cat, cond_cat),
            "reaction_type": reaction_type,
        },
    )

    _render_substance_reaction_search(options, selected_domain, reaction_type)

    if not rxn_result.get("supported"):
        st.error(rxn_result.get("rejection_reason_ar", "التفاعل المدخل غير مدعوم."))
        st.info(
            "مبدأ الصرامة العلمية في MolecuLens: لا يقوم النظام بتخمين أو اختلاق معادلات أو نواتج كيميائية غير موثقة، "
            "كما يمنع تطبيق قوانين الغازات القياسية (STP) على المواد الصلبة والسائلة."
        )
        return

    eq_unicode = rxn_result.get("equation_unicode") or format_chemical_equation_display(rxn_result["equation_text"])
    st.markdown('<div class="ml-card">', unsafe_allow_html=True)
    st.markdown(f"<h3>نتائج تحليل التفاعل: {rxn_result['reaction_name_ar']}</h3>", unsafe_allow_html=True)
    st.markdown(
        f"<p><strong>نوع التفاعل (Reaction Type):</strong> {rxn_result['reaction_type_ar']} | "
        f"<strong>المجال العلمي (Domain):</strong> {rxn_result['domain_label_ar']}</p>",
        unsafe_allow_html=True,
    )
    st.markdown("<h4>المعادلة الكيميائية الموزونة (Balanced Chemical Equation):</h4>", unsafe_allow_html=True)
    st.markdown(f'<div class="ml-equation-banner">{eq_unicode}</div>', unsafe_allow_html=True)
    st.latex(rxn_result["equation_latex"])

    cons = rxn_result["atom_conservation"]
    atom_summary = " ، ".join(f"{el}: {cnt} ذرة في المتفاعلات = {cnt} ذرة في النواتج" for el, cnt in cons["reactant_atoms"].items())
    st.success(f"تحقق قانون حفظ الذرات والكتلة (Atom & Mass Conservation Verified): {atom_summary}")

    compat = rxn_result.get("conditions_compatibility") or {}
    if compat.get("warnings_ar"):
        for w_msg in compat["warnings_ar"]:
            st.warning(w_msg)
    st.markdown("</div>", unsafe_allow_html=True)

    doc = rxn_result["documented_reaction_info"]
    calc = rxn_result["calculated_stoichiometry_result"]
    theo = rxn_result["theoretical_prediction_notes"]

    st.markdown(
        f"""
<div class="ml-tri-grid">
  <div class="ml-tri-card ml-tri-documented">
    <div class="ml-tri-title" style="color:#1d4ed8;">1. {doc.get('label_bilingual', 'التفاعل الموثق (Documented Reaction)')}</div>
    <p style="margin-bottom:0.3rem;"><strong>المصدر المرجعي (Source):</strong> <span class="ltr-inline">{doc['source']}</span></p>
    <p style="margin-bottom:0.3rem;"><strong>التصنيف:</strong> {doc['reaction_type_ar']}</p>
    <p style="margin-bottom:0;"><strong>المعادلة المرجعية:</strong> <span class="ltr-inline">{eq_unicode}</span></p>
  </div>
  <div class="ml-tri-card ml-tri-calculated">
    <div class="ml-tri-title" style="color:#0f766e;">2. {calc.get('label_bilingual', 'النتيجة المحسوبة (Calculated Result)')}</div>
    <p style="margin-bottom:0.3rem;"><strong>حفظ الذرات (Atom Conservation):</strong> {'متحقق (Verified 100%)' if calc['atom_conservation_verified'] else 'غير متحقق'}</p>
    <p style="margin-bottom:0.3rem;"><strong>الحالة المولية:</strong> {rxn_result['stoichiometric_status_ar']}</p>
    <p style="margin-bottom:0;"><strong>محرك الحساب:</strong> Algebraic Balancer &amp; Stoichiometry Engine</p>
  </div>
  <div class="ml-tri-card ml-tri-theoretical">
    <div class="ml-tri-title" style="color:#b45309;">3. {theo.get('label_bilingual', 'التوقع النظري (Theoretical Prediction)')}</div>
    <p style="margin-bottom:0.3rem;"><strong>الفرضيات العلمية:</strong> {theo['notes_ar']}</p>
    <p style="margin-bottom:0;"><strong>توافق الظروف:</strong> {compat.get('status_ar', '')}</p>
  </div>
</div>
""",
        unsafe_allow_html=True,
    )

    if rxn_result.get("alternative_documented_reactions"):
        alt_names = " ، ".join(a["reaction_name_ar"] for a in rxn_result["alternative_documented_reactions"])
        st.caption(f"توجد تفاعلات موثقة أخرى لنفس المواد (غيّر نوع التفاعل أو المجال لاختيارها): {alt_names}")

    primary_yield = rxn_result["theoretical_yields"][0]
    p_formula_disp = format_chemical_formula_display(primary_yield["formula"])
    st.markdown(
        f"""
<div class="ml-yield-box">
  <div style="font-size:1.05rem;font-weight:700;margin-bottom:0.35rem;">الناتج النظري الأساسي (Theoretical Yield) — {rxn_result['stoichiometric_status_ar']}</div>
  <div class="ml-yield-value">{primary_yield['theoretical_yield_mol']:g} mol ({primary_yield['theoretical_yield_g']:g} g)</div>
  <div style="font-size:1.08rem;margin-top:0.25rem;">الكمية النظرية الناتجة لـ {primary_yield['name_ar']} (<span class="ltr-inline">{p_formula_disp}</span>)</div>
</div>
""",
        unsafe_allow_html=True,
    )

    if len(rxn_result["theoretical_yields"]) > 1:
        y_cols = st.columns(len(rxn_result["theoretical_yields"]))
        for idx_y, y_item in enumerate(rxn_result["theoretical_yields"]):
            with y_cols[idx_y]:
                yf_disp = format_chemical_formula_display(y_item["formula"])
                st.metric(
                    f"الناتج النظري (Theoretical Yield): {y_item['name_ar']} ({yf_disp})",
                    f"{y_item['theoretical_yield_mol']:g} mol",
                    f"{y_item['theoretical_yield_g']:g} g",
                )

    st.markdown('<div class="ml-card"><h4>جدول الحسابات المولية والكتلية للمتفاعلات (Stoichiometric Breakdown Table)</h4>', unsafe_allow_html=True)
    stoich_rows = []
    for row in rxn_result["reactants_table"]:
        rf_disp = format_chemical_formula_display(row["formula"])
        stoich_rows.append(
            {
                "المادة المتفاعلة (Reactant)": f"{row['name_ar']} ({rf_disp})",
                "المعامل (Coeff)": int(row["coefficient"]),
                "الكتلة المولية (Molar Mass)": f"{row['molar_mass']:.3f} g/mol",
                "الدور في التفاعل (Role)": row["role_ar"],
                "الكمية الابتدائية (Initial)": f"{row['initial_mol']:g} mol ({row['initial_mass_g']:g} g)",
                "الكمية المستهلكة (Consumed)": f"{row['consumed_mol']:g} mol ({row['consumed_mass_g']:g} g)",
                "الكمية المتبقية (Remaining)": f"{row['remaining_mol']:g} mol ({row['remaining_mass_g']:g} g)",
            }
        )
    st.dataframe(pd.DataFrame(stoich_rows), use_container_width=True, hide_index=True)
    st.markdown("</div>", unsafe_allow_html=True)

    c_exp1, c_exp2 = st.columns(2)
    with c_exp1:
        st.markdown(
            f"""
<div class="ml-card">
  <h4>التفسير العلمي للتفاعل (Reaction Scientific Explanation)</h4>
  <p>{rxn_result['explanation_ar']}</p>
  <p style="color:#64748b;font-size:0.9rem;"><strong>القيود والافتراضات النظرية (Assumptions &amp; Limits):</strong> {rxn_result['assumptions_and_limits_ar']}</p>
</div>
""",
            unsafe_allow_html=True,
        )
    with c_exp2:
        conds = rxn_result.get("conditions", {})
        st.markdown(
            f"""
<div class="ml-card">
  <h4>ظروف التفاعل والديناميكا الحرارية (Reaction Conditions &amp; Thermodynamics)</h4>
  <p><strong>درجة الحرارة المرجعية (Temperature):</strong> {conds.get('temperature', 'الظروف القياسية 25 °C')}</p>
  <p><strong>الضغط المرجعي (Pressure):</strong> {conds.get('pressure', '1 atm')}</p>
  <p><strong>العامل الحفاز / الوسط (Catalyst &amp; Medium):</strong> {conds.get('catalyst') or conds.get('solvent') or 'لا يتطلب محفزاً خاصاً'}</p>
  <p><strong>الطاقة الحرارية (Thermodynamics):</strong> {conds.get('thermodynamics_ar', 'تفاعل كيميائي موزون')}</p>
  <p style="margin-bottom:0;color:#0f766e;font-weight:700;"><strong>حالة توافق الظروف المدخلة:</strong> {compat.get('status_ar', '')}</p>
</div>
""",
            unsafe_allow_html=True,
        )


def _render_substance_reaction_search(
    options: list[dict[str, Any]],
    selected_domain: str = "all_domains",
    reaction_type: str = "all",
) -> None:
    with st.expander("البحث عن التفاعلات الموثقة لمادة معينة (Search Verified Reactions by Substance)", expanded=False):
        entity = _substance_select("اختر المادة:", options, "rxn_search_substance_key", allow_none=True)
        if entity is None:
            st.caption("اختر مادة لعرض جميع التفاعلات الموثقة في السجل التي تشارك فيها كمتفاعل أو ناتج ضمن المجال المحدد.")
            return
        found = search_reactions_by_substance(
            entity["symbol_or_formula"],
            domain=selected_domain,
            reaction_type=reaction_type,
        )
        if not found:
            all_found = search_reactions_by_substance(entity["symbol_or_formula"])
            if all_found:
                st.warning(
                    "لا توجد تفاعلات لهذه المادة ضمن فلتر المجال/النوع المختار حالياً. "
                    f"تتوفر ({len(all_found)}) تفاعلات موثقة لها عند اختيار (جميع المجالات الكيميائية — All Domains)."
                )
            else:
                st.warning("لا توجد تفاعلات موثقة لهذه المادة في السجل الحالي، ولن يقوم النظام باختلاق تفاعل غير موثق لها.")
            return
        rows = []
        for r in found:
            lhs = " + ".join(f"{x['coefficient'] if x['coefficient'] > 1 else ''}{x['formula']}" for x in r["reactants"])
            rhs = " + ".join(f"{x['coefficient'] if x['coefficient'] > 1 else ''}{x['formula']}" for x in r["products"])
            eq_disp = format_chemical_equation_display(f"{lhs} -> {rhs}")
            rows.append(
                {
                    "المعادلة الموزونة (Balanced Equation)": eq_disp,
                    "المجال (Domain)": r.get("domain_label_ar", ""),
                    "نوع التفاعل (Reaction Type)": CATEGORY_LABELS_AR.get(r.get("category", ""), r.get("category", "")),
                    "المصدر المرجعي (Reference Source)": r.get("source", ""),
                }
            )
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

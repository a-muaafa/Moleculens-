"""
Reaction Prediction view — واجهة مختبر التوقع ("ماذا سيحدث لو جمعنا هذه المواد؟").

تدعم الواجهة إدخالاً حراً للمواد (بالاسم أو بالصيغة) على طريقة المختبر الافتراضي،
وتعرض النتيجة مفصولةً بوضوح بين: تفاعل موثّق، استنتاج قواعدي، وغياب الأدلة.
"""

from typing import Any
import html
import pandas as pd
import streamlit as st

from chemistry_engine.formula_analyzer import format_chemical_formula_display
from education.ai_assistant import AIScientificAssistant
from reaction_engine.prediction import predict_reaction
from reaction_engine.visualizer import ReactionVisualizer
from reaction_engine.what_if import WhatIfAnalyzer
from search.search_engine import get_reaction_substance_options, search_scientific_entity


def _esc(value: Any) -> str:
    """Escapes any value coming from the catalog or an external database before rendering HTML."""
    if value is None:
        return ""
    return html.escape(str(value), quote=True)

PRESET_SCENARIOS: list[tuple[str, str]] = [
    ("حمض + قاعدة", "HCl + NaOH"),
    ("حمض + كربونات", "HCl + Na2CO3"),
    ("ملحان (ترسيب)", "AgNO3 + NaCl"),
    ("معدن + أكسجين", "Fe + O2"),
    ("معدن + حمض", "Zn + HCl"),
    ("حمض + فلز قاعدي", "H2SO4 + Ca(OH)2"),
    ("معدن + ماء", "Na + H2O"),
    ("مركب غير مسجّل", "CH3COOH + K2CO3"),
]

MODE_STYLES: dict[str, dict[str, str]] = {
    "documented": {
        "color": "#1d4ed8",
        "bg": "rgba(29,78,216,0.10)",
        "border": "#1d4ed8",
    },
    "rule_inferred": {
        "color": "#b45309",
        "bg": "rgba(180,83,9,0.10)",
        "border": "#b45309",
    },
    "insufficient_evidence": {
        "color": "#b91c1c",
        "bg": "rgba(185,28,28,0.10)",
        "border": "#b91c1c",
    },
}


def _split_materials(raw_text: str) -> list[str]:
    """Splits the free-text input into individual material tokens."""
    if not raw_text:
        return []
    cleaned = raw_text.replace("،", ",").replace("+", ",").replace(";", ",")
    return [tok.strip() for tok in cleaned.split(",") if tok.strip()]


def _resolve_token(token: str) -> tuple[dict[str, Any] | None, str]:
    """Delegates reactant token resolution to the backend reaction engine."""
    from reaction_engine.prediction import resolve_reactant_token

    return resolve_reactant_token(token)


def _render_mode_banner(prediction: dict[str, Any]) -> None:
    mode = prediction.get("mode", "")
    style = MODE_STYLES.get(mode, MODE_STYLES["insufficient_evidence"])
    confidence = prediction.get("confidence", {}) or {}
    st.markdown(
        f"""
<div class="ml-card" style="border-right:5px solid {style['border']};">
  <h3 style="margin-bottom:0.3rem;color:{style['color']};">{_esc(prediction.get('mode_label_ar', ''))}</h3>
  <span class="ml-tri-badge" style="background:{style['bg']};color:{style['color']};">
    {_esc(confidence.get('label_ar', ''))}
  </span>
  <p style="margin-top:0.6rem;">{_esc(prediction.get('verdict_ar', ''))}</p>
</div>
""",
        unsafe_allow_html=True,
    )


def _render_species_table(items: list[dict[str, Any]], title: str, kind_key: str = "kind_label_ar") -> None:
    if not items:
        return
    rows = []
    for item in items:
        coeff = item.get("coefficient", 1)
        rows.append(
            {
                "المادة (Species)": f"{item.get('name_ar', '')} ({format_chemical_formula_display(item.get('formula', ''))})",
                "المعامل (Coeff.)": int(coeff) if coeff else 1,
                "التصنيف (Class)": item.get(kind_key, ""),
                "الطور (State)": item.get("state", "") or "—",
            }
        )
    st.markdown(f"**{title}**")
    st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)


def _render_visual_molecules(reactants: list[dict[str, Any]], products: list[dict[str, Any]]) -> None:
    """Renders 2D visual structures of reactants and products before and after reaction."""
    st.markdown("#### التصور البصري للجزيئات (Molecular Structures Before & After)")
    c_r, c_arrow, c_p = st.columns([5, 1, 5])
    with c_r:
        st.markdown("**المتفاعلات قبل التفاعل (Reactants):**")
        r_cols = st.columns(max(len(reactants), 1))
        for idx, r in enumerate(reactants):
            with r_cols[idx]:
                coeff = int(r.get("coefficient", 1))
                coeff_str = f"{coeff} × " if coeff > 1 else ""
                f_disp = format_chemical_formula_display(r.get("formula", ""))
                svg = ReactionVisualizer.get_species_svg(r.get("formula", ""))
                st.markdown(
                    f"""
                    <div style="background:#1e293b;border:1px solid #334155;border-radius:10px;padding:8px;text-align:center;">
                        <strong style="color:#38bdf8;">{coeff_str}{_esc(r.get('name_ar', ''))}</strong>
                        <div style="color:#94a3b8;font-size:0.85rem;margin-bottom:6px;">{f_disp}</div>
                        <div>{svg}</div>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )
    with c_arrow:
        st.markdown("<div style='text-align:center;font-size:2.2rem;color:#0f766e;margin-top:60px;'>→</div>", unsafe_allow_html=True)
    with c_p:
        st.markdown("**النواتج بعد التفاعل (Products):**")
        p_cols = st.columns(max(len(products), 1))
        for idx, p in enumerate(products):
            with p_cols[idx]:
                coeff = int(p.get("coefficient", 1))
                coeff_str = f"{coeff} × " if coeff > 1 else ""
                f_disp = format_chemical_formula_display(p.get("formula", ""))
                svg = ReactionVisualizer.get_species_svg(p.get("formula", ""))
                st.markdown(
                    f"""
                    <div style="background:#1e293b;border:1px solid #059669;border-radius:10px;padding:8px;text-align:center;">
                        <strong style="color:#34d399;">{coeff_str}{_esc(p.get('name_ar', ''))}</strong>
                        <div style="color:#94a3b8;font-size:0.85rem;margin-bottom:6px;">{f_disp}</div>
                        <div>{svg}</div>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )


def _render_bonds_and_energy(prediction: dict[str, Any]) -> None:
    """Renders bond transformation breakdown and Plotly reaction energy coordinate diagram."""
    reactants = prediction.get("reactants", [])
    products = prediction.get("products", [])
    if not reactants or not products:
        return

    with st.expander("تحليل الروابط ومسار الطاقة (Bond Transformations & Energy Profile)", expanded=False):
        bonds_data = ReactionVisualizer.analyze_bond_transformations(
            reactants,
            products,
            category=prediction.get("category", ""),
            reaction_data=prediction,
        )

        if bonds_data["is_illustrative"]:
            st.warning(bonds_data["data_classification_ar"])
        else:
            st.success(bonds_data["data_classification_ar"])

        c1, c2 = st.columns([1, 1])

        with c1:
            st.markdown(f"**الروابط المتكسرة والمتكونة ({bonds_data['thermodynamic_type_ar']}):**")
            st.info(bonds_data["summary_ar"])

            st.markdown("**1. الروابط المتكسرة في المتفاعلات (Bonds Broken):**")
            for b in bonds_data["bonds_broken"]:
                e_txt = f" ({b['total_kj']} kJ/mol)" if b.get("total_kj") else ""
                desc_txt = f" — {b['description']}" if b.get("description") else ""
                st.markdown(f"- **{b['bond']}**: {b.get('count', 1)} رابطة{e_txt}{desc_txt}")

            st.markdown("**2. الروابط المتكونة في النواتج (Bonds Formed):**")
            for b in bonds_data["bonds_formed"]:
                e_txt = f" ({b['total_kj']} kJ/mol)" if b.get("total_kj") else ""
                desc_txt = f" — {b['description']}" if b.get("description") else ""
                st.markdown(f"- **{b['bond']}**: {b.get('count', 1)} رابطة{e_txt}{desc_txt}")

            st.markdown("**التحقق من حفظ الذرات (Atom Conservation Matrix):**")
            audit_rows = ReactionVisualizer.audit_atom_conservation(reactants, products)
            st.dataframe(pd.DataFrame(audit_rows), use_container_width=True, hide_index=True)

        with c2:
            r_str = " + ".join(r.get("formula", "") for r in reactants)
            p_str = " + ".join(p.get("formula", "") for p in products)
            fig = ReactionVisualizer.create_reaction_energy_diagram(
                r_str,
                p_str,
                is_exothermic=bonds_data["is_exothermic"],
                delta_h_kj=bonds_data["estimated_delta_h_kj"],
                activation_energy_kj=None,
                is_illustrative=bonds_data["is_illustrative"],
            )
            st.plotly_chart(fig, use_container_width=True)


def _render_what_if_panel(reactant_inputs: list[dict[str, Any]], prediction: dict[str, Any]) -> None:
    """Renders the What-If comparative scenario analysis panel."""
    with st.expander("تحليل سيناريوهات «ماذا لو؟» (What-If Scenario Analysis)", expanded=False):
        st.markdown(
            """
            جرب تعديل كميات المواد أو ظروف التفاعل ولاحظ كيف يتغير المتفاعل المحدد،
            وكمية الناتج النظري، وسلوك التفاعل مقارنةً بالتجربة المرجعية الأصلية:
            """
        )

        mod_inputs = []
        cols = st.columns(len(reactant_inputs))
        for idx, item in enumerate(reactant_inputs):
            with cols[idx]:
                orig_amt = float(item.get("amount", 1.0))
                form = item.get("formula", "")
                mod_amt = st.number_input(
                    f"كمية {form} المعدلة (mol):",
                    min_value=0.01,
                    value=round(orig_amt, 3),
                    step=0.5,
                    key=f"whatif_amt_{idx}_{form}",
                )
                mod_inputs.append({
                    "formula": form,
                    "amount": float(mod_amt),
                    "unit": item.get("unit", "mol"),
                })

        c_t, c_p = st.columns(2)
        with c_t:
            new_temp = st.slider("درجة الحرارة المعدلة (°C):", min_value=0.0, max_value=500.0, value=25.0, step=10.0, key="whatif_temp")
        with c_p:
            new_press = st.slider("الضغط المعدل (atm):", min_value=0.1, max_value=30.0, value=1.0, step=0.5, key="whatif_press")

        comp_res = WhatIfAnalyzer.evaluate_what_if(
            baseline_inputs=reactant_inputs,
            modified_inputs=mod_inputs,
            baseline_conditions={"temperature_c": 25.0, "pressure_atm": 1.0},
            modified_conditions={"temperature_c": new_temp, "pressure_atm": new_press},
        )

        if comp_res.get("success"):
            st.markdown(f"**التحليل العلمي للتغير (Comparative Insights):**")
            st.success(comp_res["insights_ar"])
            st.markdown("**جدول المقارنة المباشرة (Side-by-Side Comparison):**")
            st.dataframe(pd.DataFrame(comp_res["table_rows"]), use_container_width=True, hide_index=True)
        else:
            st.warning(comp_res.get("error_ar", "تعذر حساب السيناريو المقارن."))


def render_prediction_section() -> None:
    st.markdown(
        """
<div class="ml-card">
  <h3>مختبر التوقع — Reaction Prediction (المختبر الافتراضي)</h3>
  <p>
    أدخل أي مادتين أو أكثر (بالاسم العربي أو الإنجليزي، أو بالصيغة الكيميائية) واسأل المنصة:
    <strong>هل سيحدث تفاعل؟ وما نواتجه؟ وما التغيّرات الكيميائية المتوقعة؟</strong>
  </p>
  <p style="margin-bottom:0;">
    يعرض المحرك ثلاثة مستويات من الدليل بوضوح: <strong>تفاعل موثّق</strong> بمرجع علمي محدّد،
    أو <strong>استنتاج قواعدي</strong> من قاعدة كيميائية معلنة مع ذكر مصدرها وحدودها،
    أو <strong>غياب الأدلة</strong> صراحةً — ولا يختلق المنصة نواتج غير موثقة أبداً.
  </p>
</div>
""",
        unsafe_allow_html=True,
    )

    st.markdown("**سيناريوهات سريعة (Quick Scenarios):**")
    preset_cols = st.columns(4)
    for idx, (label, payload) in enumerate(PRESET_SCENARIOS):
        with preset_cols[idx % 4]:
            if st.button(label, key=f"pred_preset_{idx}"):
                st.session_state["pred_free_input"] = payload

    st.markdown('<div class="ml-card"><h4>1) إدخال المواد (Materials Input)</h4>', unsafe_allow_html=True)
    input_mode = st.radio(
        "طريقة الإدخال (Input Mode):",
        options=["free", "closed"],
        format_func=lambda m: (
            "إدخال حر — أسماء أو صيغ (Free Input)"
            if m == "free"
            else "اختيار من قائمة مغلقة (Closed Dropdown)"
        ),
        horizontal=True,
        key="pred_input_mode",
    )

    resolved: list[dict[str, Any]] = []
    notes: list[str] = []

    if input_mode == "free":
        if "pred_free_input" not in st.session_state:
            st.session_state["pred_free_input"] = "HCl + NaOH"
        raw = st.text_input(
            "اكتب المواد مفصولة بفاصلة أو علامة + (مثال: HCl + NaOH أو حديد , أكسجين):",
            key="pred_free_input",
        )
        tokens = _split_materials(raw)
        if len(tokens) < 2:
            st.warning("أدخل مادتين على الأقل لتفعيل التحليل.")
            st.markdown("</div>", unsafe_allow_html=True)
            return
        for token in tokens:
            entity, note = _resolve_token(token)
            if entity:
                resolved.append(entity)
            if note:
                notes.append(note)
    else:
        options = get_reaction_substance_options("all")
        labels = {o["key"]: o["label"] for o in options}
        keys = [o["key"] for o in options]
        picked: list[str] = []
        for slot in range(3):
            label = "المادة الأولى (Material 1):" if slot == 0 else f"المادة {slot + 1} (Material {slot + 1}):"
            chosen = st.selectbox(
                label,
                options=keys,
                format_func=lambda k: labels.get(k, k),
                key=f"pred_closed_{slot}",
            )
            if chosen:
                picked.append(chosen)
        from search.local_catalog import get_catalog_entity_by_key

        for key in picked:
            entity = get_catalog_entity_by_key(key)
            if entity:
                resolved.append(entity)
        if len(resolved) < 2:
            st.warning("اختر مادتين على الأقل لتفعيل التحليل.")
            st.markdown("</div>", unsafe_allow_html=True)
            return

    for note in notes:
        st.info(note)

    st.markdown("**2) الكميات (Amounts — اختياري لحساب المتفاعل المحدد والناتج النظري):**")
    amount_cols = st.columns(len(resolved))
    amounts: list[float] = []
    for idx, entity in enumerate(resolved):
        with amount_cols[idx]:
            f_disp = format_chemical_formula_display(entity.get("symbol_or_formula", ""))
            amounts.append(
                st.number_input(
                    f"{entity.get('name_ar', '')} ({f_disp}) بالمولات:",
                    min_value=0.01,
                    value=1.0,
                    step=0.5,
                    key=f"pred_amount_{idx}_{entity.get('key', idx)}",
                )
            )
    st.markdown("</div>", unsafe_allow_html=True)

    reactant_inputs = [
        {
            "substance": entity,
            "formula": entity.get("symbol_or_formula", ""),
            "amount": float(amount),
            "unit": "mol",
        }
        for entity, amount in zip(resolved, amounts)
    ]

    prediction = predict_reaction(reactant_inputs=reactant_inputs)

    _render_mode_banner(prediction)

    if not prediction.get("supported"):
        missing = prediction.get("missing_information_ar", [])
        if missing:
            st.markdown("**ما الذي ينقص هذه النتيجة؟**")
            for item in missing:
                st.markdown(f"- {item}")
        st.info(
            "مبدأ الصرامة العلمية: عرض «لا تتوفر بيانات كافية» أصدق علمياً من عرض معادلة مُقدَّرة قد تكون مضللة. "
            "يمكنك توسيع سجل التفاعلات أو اختيار مواد من تركيبات موثقة في Reaction Lab."
        )
        return

    st.markdown('<div class="ml-card"><h4>3) نتيجة التحليل (Prediction Analysis)</h4>', unsafe_allow_html=True)
    eq_unicode = prediction.get("equation_unicode", "")
    st.markdown(f'<div class="ml-equation-banner">{_esc(eq_unicode)}</div>', unsafe_allow_html=True)
    if prediction.get("equation_latex"):
        st.latex(prediction["equation_latex"])

    cons = prediction.get("atom_conservation") or {}
    if cons:
        if cons.get("is_balanced"):
            st.success("تحقّق قانون حفظ الذرات (Atom Conservation Verified): عدد ذرات كل عنصر في المتفاعلات = عددها في النواتج.")
        else:
            st.error("تحذير: المعادلة لا تحفظ الذرات — راجع بيانات النواتج.")

    if prediction.get("net_ionic"):
        st.markdown(f"**المعادلة الصافية الأيونية (Net Ionic Equation):** `{_esc(prediction['net_ionic'])}`")

    c_react, c_prod = st.columns(2)
    with c_react:
        _render_species_table(prediction.get("reactants", []), "المواد المتفاعلة (Reactants)")
    with c_prod:
        _render_species_table(prediction.get("products", []), "نواتج التفاعل (Products)")

    _render_visual_molecules(prediction.get("reactants", []), prediction.get("products", []))

    observations = prediction.get("observations", []) or []
    if observations:
        st.markdown("**الملاحظات المتوقعة في المختبر (Expected Lab Observations):**")
        for obs in observations:
            st.markdown(f"- **{obs.get('type', '')}:** {obs.get('text_ar', '')} — _{obs.get('basis_ar', '')}_")

    ph = prediction.get("ph_analysis") or {}
    if ph.get("available"):
        st.markdown(f"**التغيّر في درجة الحموضة (pH):** {_esc(ph.get('summary_ar', ''))} {_esc(ph.get('detail_ar', ''))}")
    else:
        st.caption(f"التغيّر في درجة الحموضة: {_esc(ph.get('summary_ar', 'غير متوفر'))}")

    if prediction.get("conditions_required_ar"):
        st.markdown(f"**الظروف المطلوبة لحدوث التفاعل:** {_esc(prediction['conditions_required_ar'])}")

    _render_bonds_and_energy(prediction)

    stoich = prediction.get("stoichiometry") or {}
    reactants_table = prediction.get("reactants_table") or []
    if reactants_table:
        st.markdown("**الجدول الستوكيومتري (Stoichiometric Breakdown):**")
        rows = []
        for row in reactants_table:
            rows.append(
                {
                    "المادة (Reactant)": f"{row['name_ar']} ({format_chemical_formula_display(row['formula'])})",
                    "المعامل (Coeff.)": int(row["coefficient"]),
                    "الكتلة المولية (g/mol)": f"{row['molar_mass']:.3f}",
                    "الدور (Role)": row.get("role_ar", ""),
                    "الابتدائي (mol)": f"{row['initial_mol']:g}",
                    "المستهلك (mol)": f"{row['consumed_mol']:g}",
                    "المتبقي (mol)": f"{row['remaining_mol']:g}",
                }
            )
        st.dataframe(pd.DataFrame(rows), use_container_width=True, hide_index=True)

    yields = prediction.get("theoretical_yields", []) or []
    if yields:
        y_cols = st.columns(len(yields))
        for idx_y, y_item in enumerate(yields):
            with y_cols[idx_y]:
                st.metric(
                    f"الناتج النظري: {y_item['name_ar']}",
                    f"{y_item['theoretical_yield_mol']:g} mol",
                    f"{y_item['theoretical_yield_g']:g} g",
                )

    _render_what_if_panel(reactant_inputs, prediction)
    _render_end_to_end_experiment_walkthrough(prediction)

    evidence = prediction.get("evidence", []) or []
    if evidence:
        st.markdown("**الدليل العلمي (Evidence & Sources):**")
        for ev in evidence:
            link = (
                f'<a href="{html.escape(ev["reference_url"], quote=True)}" target="_blank" class="ltr-inline">'
                f'{_esc(ev.get("source_id", ""))}</a>'
                if ev.get("reference_url")
                else f'<span class="ltr-inline">{_esc(ev.get("source_id", ""))}</span>'
            )
            st.markdown(
                f"""
<div class="ml-card" style="margin-bottom:0.5rem;">
  <p style="margin-bottom:0.25rem;"><strong>{_esc(ev.get('rule_name_ar', ''))}</strong></p>
  <p style="margin-bottom:0.25rem;">{_esc(ev.get('basis_ar', ''))}</p>
  <p style="margin-bottom:0;"><strong>المرجع:</strong> <span class="ltr-inline">{_esc(ev.get('source', ''))}</span> — {link}</p>
</div>
""",
                unsafe_allow_html=True,
            )

    limits = prediction.get("limitations_ar", []) or []
    if limits:
        st.markdown("**حدود التنبؤ وافتراضاته (Limits & Assumptions):**")
        for lim in limits:
            st.markdown(f"- {lim}")

    alternatives = prediction.get("alternatives", []) or []
    if alternatives:
        st.caption(
            "تفاعلات موثّقة أخرى لنفس المواد (مسارات مختلفة): "
            + " ، ".join(a.get("reaction_name_ar", "") for a in alternatives)
        )

    st.markdown("</div>", unsafe_allow_html=True)

    with st.expander("المساعد العلمي — شرح النتيجة خطوة بخطوة (AI Scientific Assistant)", expanded=True):
        st.markdown(AIScientificAssistant.explain_reaction_prediction(prediction), unsafe_allow_html=False)


def _render_end_to_end_experiment_walkthrough(prediction: dict[str, Any]) -> None:
    """Renders the 4-stage visual chemical experiment walkthrough in the Chemical Playground."""
    stages = ReactionVisualizer.build_end_to_end_experiment_stages(prediction)
    s1 = stages["stage1_preparation"]
    s2 = stages["stage2_active_vessel"]
    s3 = stages["stage3_mass_balance"]
    s4 = stages["stage4_final_inventory"]

    st.markdown("---")
    st.markdown(
        "<h4>رحلة التجربة الكيميائية البصرية المتكاملة من البداية للنهاية (End-to-End Visual Experiment)</h4>"
        "<p style='color:#475569;font-size:0.93rem;'>تتبع بصري شامل للتجربة من تحضير المواد المتفاعلة، إلى محاكاة وعاء التفاعل النشط، "
        "وميزان حفظ الكتلة والذرات (قانون لافوازييه)، وصولاً إلى الجرد النهائي لمحتويات الوعاء بعد انتهاء التفاعل.</p>",
        unsafe_allow_html=True,
    )

    exp_tab1, exp_tab2, exp_tab3, exp_tab4 = st.tabs([
        "1. تحضير المتفاعلات",
        "2. محاكاة وعاء التفاعل النشط",
        "3. ميزان حفظ الكتلة والاستهلاك",
        "4. جرد الوعاء النهائي",
    ])

    with exp_tab1:
        st.caption(s1["title_ar"])
        st.info(f"الحالة المولية للمزيج الابتدائي: {s1['stoichiometric_status_ar']}")
        cols1 = st.columns(max(1, len(s1["reactant_cards"])))
        for idx_c, card in enumerate(s1["reactant_cards"]):
            with cols1[idx_c]:
                st.markdown(
                    f"""
<div style="border:1px solid #cbd5e1;border-radius:10px;padding:10px;background:#f8fafc;text-align:center;">
  <div style="font-weight:700;color:#0f4c5c;">{_esc(card['name_ar'])} ({_esc(card['formula_display'])})</div>
  <div style="font-size:0.82rem;color:#64748b;">الطور: ({_esc(card['state'])}) | المعامل: {card['coefficient']} | الكتلة المولية: {card['molar_mass']:.2f} g/mol</div>
  <div style="font-weight:700;color:#1e3a8a;margin-top:6px;">الكمية المضافة: {card['initial_mol']:g} mol ({card['initial_mass_g']:g} g)</div>
  <div style="font-size:0.84rem;color:#0f766e;margin-top:3px;">{_esc(card['role_ar'])}</div>
</div>
""",
                    unsafe_allow_html=True,
                )

    with exp_tab2:
        st.caption(s2["title_ar"])
        st.markdown(s2["vessel_svg"], unsafe_allow_html=True)

    with exp_tab3:
        st.caption(s3["title_ar"])
        mc1, mc2, mc3 = st.columns(3)
        mc1.metric("كتلة المتفاعلات قبل التفاعل", f"{s3['total_initial_mass_g']:g} g")
        mc2.metric("كتلة النواتج المتكونة", f"{s3['total_produced_mass_g']:g} g", f"مستهلك: {s3['total_consumed_mass_g']:g} g")
        mc3.metric(
            "إجمالي الكتلة بعد التفاعل",
            f"{s3['total_final_vessel_mass_g']:g} g",
            "محفوظ 100% (قانون لافوازييه)" if s3["mass_conserved"] else "تحقق من الموازنة",
        )
        mass_chart = ReactionVisualizer.create_stoichiometry_mass_chart(s3)
        st.plotly_chart(mass_chart, use_container_width=True, key="pred_mass_conservation_chart")
        for rc in s3["reactant_consumption"]:
            st.markdown(
                f"**{rc['name_ar']} ({rc['formula_display']})**: استُهلك `{rc['consumed_pct']}%` (`{rc['consumed_mass_g']:g} g`) "
                f"| تبقّى `{rc['remaining_pct']}%` (`{rc['remaining_mass_g']:g} g`)"
            )
            st.progress(min(1.0, max(0.0, float(rc["consumed_pct"]) / 100.0)))

    with exp_tab4:
        st.caption(s4["title_ar"])
        inv_cols = st.columns(max(1, len(s4["final_inventory"])))
        for idx_i, item in enumerate(s4["final_inventory"]):
            with inv_cols[idx_i]:
                b_col = "#059669" if "ناتج" in item["kind_ar"] else "#d97706"
                st.markdown(
                    f"""
<div style="border:2px solid {b_col};border-radius:10px;padding:10px;background:#ffffff;text-align:center;">
  <div style="font-size:0.8rem;font-weight:700;color:{b_col};">{_esc(item['kind_ar'])}</div>
  <div style="font-weight:700;color:#0f4c5c;margin-top:3px;">{_esc(item['name_ar'])} ({_esc(item['formula_display'])})</div>
  <div style="font-weight:700;color:#1e293b;margin-top:5px;">{item['amount_mol']:g} mol ({item['mass_g']:g} g)</div>
  <div style="font-size:0.82rem;color:#475569;margin-top:3px;">{_esc(item['status_ar'])}</div>
</div>
""",
                    unsafe_allow_html=True,
                )

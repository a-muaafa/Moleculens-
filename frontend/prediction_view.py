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
    """Resolves a single token (name, formula or SMILES) into a catalog entity."""
    res = search_scientific_entity(token, search_type="auto", allow_network=True)
    if res.get("found") and res.get("selected"):
        entity = res["selected"]
        return entity, ""
    if res.get("is_ambiguous") and res.get("candidates"):
        entity = res["candidates"][0]
        return entity, f"المادة «{token}» لها عدة احتمالات؛ استُخدم الاحتمال الأول: {entity.get('name_ar')}"
    return None, f"تعذّر التعرف على المادة «{token}» (استخدم الاسم أو الصيغة الكيميائية أو SMILES)."


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
        raw = st.text_input(
            "اكتب المواد مفصولة بفاصلة أو علامة + (مثال: HCl + NaOH أو حديد , أكسجين):",
            value=st.session_state.get("pred_free_input", "HCl + NaOH"),
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

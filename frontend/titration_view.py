"""
Titration Lab frontend view for MolecuLens.
Interactive acid-base titration simulator with 41 acids, 30 bases, 23 indicators,
automatic optimal indicator selection, dynamic Erlenmeyer flask color SVG,
polyprotic stepwise equivalence points, and endpoint vs. equivalence point detection.
"""

from __future__ import annotations

import html
import pandas as pd
import streamlit as st
from titration.titration_engine import TitrationEngine


def render_titration_view(lang: str = "ar") -> None:
    """Renders the Titration Virtual Laboratory interface."""
    is_ar = lang == "ar"

    header_title = "مختبر المعايرة الافتراضي الموسّع (Virtual Titration Lab)" if is_ar else "Virtual Titration Laboratory"
    header_desc = (
        "محاكاة ديناميكية لتفاعلات المعايرة للأحماض والقواعد القوية والضعيفة ومتعددة البروتونات، "
        "مع اختيار الكاشف اللوني تلقائياً ومحاكاة تغير اللون داخل الدورق المخروطي ورسم نقطة التكافؤ ونقطة النهاية."
        if is_ar
        else "Dynamic simulation of strong, weak, and polyprotic acid-base titrations with automatic indicator selection, flask color simulation, and equivalence/endpoint tracking."
    )

    st.markdown(
        f"""
        <div style="background: rgba(30, 41, 59, 0.7); border: 1px solid rgba(255, 255, 255, 0.1); border-radius: 12px; padding: 20px; margin-bottom: 24px;">
            <h2 style="color: #60a5fa; margin: 0 0 8px 0; font-size: 1.5rem;">{html.escape(header_title)}</h2>
            <p style="color: #94a3b8; margin: 0; font-size: 0.95rem;">{html.escape(header_desc)}</p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    acids = TitrationEngine.ACID_DATABASE
    bases = TitrationEngine.BASE_DATABASE
    indicators = TitrationEngine.INDICATORS

    col1, col2 = st.columns([1, 1])

    with col1:
        st.subheader("1. " + ("تحديد الحمض والقاعدة وبياناتهما الكيميائية" if is_ar else "Select Acid & Base"))
        acid_options = list(acids.keys())
        acid_labels = {
            k: f"{acids[k]['name_ar' if is_ar else 'name_en']} ({k}) — {'قوي' if acids[k]['type'] == 'strong' else 'ضعيف'} (التكافؤ={acids[k]['protons']})"
            for k in acid_options
        }
        selected_acid_key = st.selectbox(
            "الحمض (Acid)" if is_ar else "Acid",
            options=acid_options,
            format_func=lambda k: acid_labels[k],
            index=0,
        )
        a_spec = acids[selected_acid_key]
        st.caption(
            f"**{a_spec['name_ar']} ({selected_acid_key}):** {a_spec.get('strength_ar', '')} | "
            f"التصنيف: {a_spec.get('category_ar', '')} | "
            f"الكتلة المولية: `{a_spec.get('molar_mass', 0)} g/mol` | "
            f"ثوابت التأين pKa: `{', '.join(str(x) for x in a_spec.get('pKa_list', [a_spec['pKa']]))}` | "
            f"الذوبانية: {a_spec.get('solubility_ar', '')} | "
            f"المصدر: *{a_spec.get('source', 'CRC Handbook / Harris QCA')}*"
        )

        base_options = list(bases.keys())
        base_labels = {
            k: f"{bases[k]['name_ar' if is_ar else 'name_en']} ({k}) — {'قوية' if bases[k]['type'] == 'strong' else 'ضعيفة'} (التكافؤ={bases[k]['hydroxides']})"
            for k in base_options
        }
        selected_base_key = st.selectbox(
            "القاعدة (Base)" if is_ar else "Base",
            options=base_options,
            format_func=lambda k: base_labels[k],
            index=0,
        )
        b_spec = bases[selected_base_key]
        st.caption(
            f"**{b_spec['name_ar']} ({selected_base_key}):** {b_spec.get('strength_ar', '')} | "
            f"التصنيف: {b_spec.get('category_ar', '')} | "
            f"الكتلة المولية: `{b_spec.get('molar_mass', 0)} g/mol` | "
            f"ثوابت التأين pKb: `{', '.join(str(x) for x in b_spec.get('pKb_list', [b_spec['pKb']]))}` | "
            f"الذوبانية: {b_spec.get('solubility_ar', '')} | "
            f"المصدر: *{b_spec.get('source', 'CRC Handbook / Harris QCA')}*"
        )

        mode_options = [
            ("acid_with_base", "معايرة حمض في الدورق بواسطة قاعدة في السحاحة (Analyte: Acid, Titrant: Base)" if is_ar else "Titrate Acid with Base"),
            ("base_with_acid", "معايرة قاعدة في الدورق بواسطة حمض في السحاحة (Analyte: Base, Titrant: Acid)" if is_ar else "Titrate Base with Acid"),
        ]
        titration_mode = st.selectbox(
            "مسار المعايرة" if is_ar else "Titration Mode",
            options=[m[0] for m in mode_options],
            format_func=lambda val: next(m[1] for m in mode_options if m[0] == val),
            index=0,
        )

    with col2:
        st.subheader("2. " + ("الحجوم والتراكيز واختيار الكاشف اللوني" if is_ar else "Volumes, Concentrations & Indicator"))
        analyte_label = "الحمض في الدورق" if titration_mode == "acid_with_base" else "القاعدة في الدورق"
        analyte_v_label = f"حجم {analyte_label} (مل mL)" if is_ar else "Analyte Volume (mL)"
        analyte_volume_ml = st.number_input(
            analyte_v_label, min_value=5.0, max_value=200.0, value=25.0, step=5.0
        )

        analyte_c_label = f"تركيز {analyte_label} (مولار M)" if is_ar else "Analyte Concentration (M)"
        analyte_conc_m = st.number_input(
            analyte_c_label, min_value=0.001, max_value=5.0, value=0.100, step=0.01
        )

        titrant_label = "القاعدة في السحاحة" if titration_mode == "acid_with_base" else "الحمض في السحاحة"
        titrant_c_label = f"تركيز {titrant_label} (مولار M)" if is_ar else "Titrant Concentration (M)"
        titrant_conc_m = st.number_input(
            titrant_c_label, min_value=0.001, max_value=5.0, value=0.100, step=0.01
        )

        auto_select_ind = st.checkbox(
            "اختيار الكاشف اللوني الأنسب تلقائياً حسب نقطة التكافؤ (Auto-Select Optimal Indicator)"
            if is_ar
            else "Auto-Select Optimal Indicator by Equivalence Point",
            value=True,
        )

    # Perform Titration Simulation
    if titration_mode == "acid_with_base":
        sim_data = TitrationEngine.simulate_titration(
            acid_key=selected_acid_key,
            base_key=selected_base_key,
            acid_conc_m=analyte_conc_m,
            acid_volume_ml=analyte_volume_ml,
            base_conc_m=titrant_conc_m,
            titration_mode="acid_with_base",
        )
    else:
        sim_data = TitrationEngine.simulate_titration(
            acid_key=selected_acid_key,
            base_key=selected_base_key,
            acid_conc_m=titrant_conc_m,
            acid_volume_ml=analyte_volume_ml,
            base_conc_m=analyte_conc_m,
            titration_mode="base_with_acid",
            base_volume_ml=analyte_volume_ml,
        )

    rec_ind_key = sim_data.get("recommended_indicator_key", "bromothymol_blue")
    ind_keys = list(indicators.keys())
    with col2:
        if auto_select_ind:
            selected_indicator = rec_ind_key
            rec_ind = indicators[selected_indicator]
            st.success(
                f"الكاشف المختار تلقائياً: **{rec_ind['name_ar']}** "
                f"(مدى التحول pH {rec_ind['transition_range'][0]}–{rec_ind['transition_range'][1]} | "
                f"حمضي: {rec_ind['color_acid']} -> قاعدي: {rec_ind['color_base']})"
            )
        else:
            selected_indicator = st.selectbox(
                "الكاشف اللوني (Chemical Indicator)" if is_ar else "Chemical Indicator",
                options=ind_keys,
                index=ind_keys.index(rec_ind_key) if rec_ind_key in ind_keys else 0,
                format_func=lambda k: (
                    f"{indicators[k]['name_ar' if is_ar else 'name_en']} "
                    f"(pH {indicators[k]['transition_range'][0]}-{indicators[k]['transition_range'][1]}: "
                    f"{indicators[k]['color_acid']} -> {indicators[k]['color_base']})"
                ),
            )

    v_eq = sim_data["v_eq_ml"]
    ph_eq = sim_data["ph_eq"]
    half_eq = sim_data.get("half_v_eq_ml")
    half_ph = sim_data.get("half_ph")

    # Compute endpoint volume and pH for the active indicator
    active_ind = indicators[selected_indicator]
    ind_mid = 0.5 * (active_ind["transition_range"][0] + active_ind["transition_range"][1])
    curve_pts = sim_data["sim_result"]["curve_points"]
    endpt_row = min(curve_pts, key=lambda pt: abs(pt["ph"] - ind_mid))

    st.markdown("---")

    st.subheader("3. " + ("محاكاة إضافة المحلول وتغير لون الكاشف داخل الدورق" if is_ar else "Interactive Addition & Flask Color Simulation"))
    max_slider_vol = float(sim_data["volumes_ml"][-1])
    current_vol = st.slider(
        "حجم المحلول المضاف من السحاحة (مل mL)" if is_ar else "Volume Added from Burette (mL)",
        min_value=0.0,
        max_value=max_slider_vol,
        value=float(min(v_eq, max_slider_vol)),
        step=0.2,
    )

    current_state = TitrationEngine.calculate_point_ph(
        v_add=current_vol,
        v_eq=v_eq,
        v_acid=analyte_volume_ml,
        ca=analyte_conc_m if titration_mode == "acid_with_base" else titrant_conc_m,
        cb=titrant_conc_m if titration_mode == "acid_with_base" else analyte_conc_m,
        acid_info=acids[selected_acid_key],
        base_info=bases[selected_base_key],
        titration_mode=titration_mode,
    )
    current_ph = current_state["ph"]

    m1, m2, m3, m4, m5 = st.columns(5)
    m1.metric("حجم التكافؤ الكلي V_eq", f"{v_eq:.2f} mL")
    m2.metric("pH نقطة التكافؤ", f"{ph_eq:.2f}")
    m3.metric("نقطة نهاية الكاشف (Endpoint)", f"{endpt_row['volume_added_ml']:.2f} mL", f"pH {endpt_row['ph']:.2f}")
    m4.metric("الحجم المضاف حالياً", f"{current_vol:.2f} mL")
    m5.metric("pH الدورق الحالي", f"{current_ph:.2f}", delta=current_state["region"])

    if sim_data.get("is_polyprotic") and len(sim_data.get("step_equivalence_volumes_ml", [])) > 1:
        steps_str = " ، ".join(
            f"التكافؤ {i + 1}: {v:.2f} mL"
            for i, v in enumerate(sim_data["step_equivalence_volumes_ml"])
        )
        st.info(f"معايرة متعددة البروتونات/التكافؤ — نقاط التكافؤ المرحلية: {steps_str}")

    viz_col1, viz_col2 = st.columns([2, 3])
    with viz_col1:
        flask_svg = TitrationEngine.create_flask_svg(
            current_ph=current_ph,
            indicator_key=selected_indicator,
            current_vol_ml=current_vol,
            v_eq_ml=v_eq,
        )
        st.markdown(flask_svg, unsafe_allow_html=True)

    with viz_col2:
        fig = TitrationEngine.create_plotly_curve(
            sim_data=sim_data,
            current_v=current_vol,
            indicator_key=selected_indicator,
        )
        st.plotly_chart(fig, use_container_width=True)

    with st.expander("دليل الكواشف اللونية المعتمدة (الألوان في الوسط الحمضي والقاعدي ومجالات التحول)", expanded=False):
        ind_rows = [
            {
                "الكاشف اللوني (Indicator)": v["name_ar"],
                "مجال تغير اللون (pH Range)": f"{v['pH_range'][0]} – {v['pH_range'][1]}",
                "اللون في الوسط الحمضي (Acidic Color)": v["color_acid"],
                "اللون في الوسط القاعدي (Basic Color)": v["color_base"],
                "الملائمة التحليلية": v["suitability_ar"],
                "المصدر العلمي المرجعي": v.get("source", "Harris QCA / CRC Handbook"),
            }
            for v in indicators.values()
        ]
        st.dataframe(pd.DataFrame(ind_rows), use_container_width=True, hide_index=True)

    # Indicator Suitability Assessment
    ind_range = indicators[selected_indicator]["transition_range"]
    is_suitable = ind_range[0] <= ph_eq <= ind_range[1] or abs(ph_eq - ((ind_range[0] + ind_range[1]) / 2.0)) <= 1.5

    suit_color = "#10b981" if is_suitable else "#f59e0b"
    suit_text = (
        ("الكاشف ملائم لهذه المعايرة لأن مدى تغير لونه يتقاطع مع القفزة الهيدروجينية لنقطة التكافؤ." if is_ar
         else "Selected indicator is suitable: its transition interval encompasses the sharp equivalence pH inflection.")
        if is_suitable
        else ("الكاشف غير مثالي؛ يفضل تفعيل الاختيار التلقائي أو اختيار كاشف يتطابق مدى تحوله اللوني بدقة مع pH نقطة التكافؤ." if is_ar
              else "Suboptimal indicator choice; select an indicator whose transition range tightly matches the equivalence pH.")
    )

    st.markdown(
        f"""
        <div style="background: rgba(15, 23, 42, 0.6); border-left: 4px solid {suit_color}; border-radius: 8px; padding: 14px 18px; margin-top: 16px;">
            <strong style="color: {suit_color}; font-size: 1rem;">
                {("تقييم ملائمة الكاشف: " if is_ar else "Indicator Assessment: ")}
                {indicators[selected_indicator]['name_ar' if is_ar else 'name_en']}
            </strong>
            <p style="color: #cbd5e1; margin: 6px 0 0 0; font-size: 0.9rem;">{html.escape(suit_text)}</p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    st.markdown("---")
    st.subheader("4. " + ("التفسير العلمي ومراحل المنحنى" if is_ar else "Scientific Breakdown"))

    desc_col1, desc_col2 = st.columns(2)
    with desc_col1:
        st.markdown(
            f"""
            **{('التفاعل الأساسي:' if is_ar else 'Chemical Reaction:')}**
            `{acids[selected_acid_key]['name_en']} ({selected_acid_key}) + {bases[selected_base_key]['name_en']} ({selected_base_key}) -> Salt + H2O`

            - **{('طبيعة الحمض:' if is_ar else 'Acid Type:')}** {acids[selected_acid_key]['strength_ar']}
            - **{('طبيعة القاعدة:' if is_ar else 'Base Type:')}** {bases[selected_base_key]['strength_ar']}
            - **{('الرقم الهيدروجيني الأولي:' if is_ar else 'Initial pH:')}** {sim_data['initial_ph']:.2f}
            """
        )

    with desc_col2:
        if half_eq is not None and half_ph is not None:
            st.markdown(
                f"""
                **{('نقطة نصف التكافؤ (Half-Equivalence):' if is_ar else 'Half-Equivalence Point:')}**
                - **V_half:** {half_eq:.2f} mL
                - **pH_half:** {half_ph:.2f} *(pH = pKa1)*
                - **{('الأهمية الكيميائية:' if is_ar else 'Significance:')}** {('عند هذا الحجم يتساوى تركيز الحمض الضعيف مع قاعدته المرافقة، مما يعطي أقصى سعة للمحلول المنظم (Buffer Capacity).' if is_ar else 'Equal concentrations of weak acid and conjugate base maximize buffer capacity.')}
                """
            )
        else:
            st.markdown(
                f"""
                **{('نقطة التكافؤ (Equivalence Point):' if is_ar else 'Equivalence Point:')}**
                - **V_eq:** {v_eq:.2f} mL
                - **pH_eq:** {ph_eq:.2f}
                - **{('التعليل:' if is_ar else 'Explanation:')}** {('ينتج ملح متعادل تماماً ولا يحدث تميؤ ملحي، لذلك pH = 7.0 عند 25 درجة مئوية.' if is_ar else 'Neutral salt formed with no salt hydrolysis; pH = 7.0 at 25 °C.')}
                """
            )

    buf_info = sim_data.get("buffer_analysis") or {}
    err_info = sim_data.get("titration_error_analysis") or {}
    if buf_info.get("summary_ar"):
        st.info(f"**تحليل المحلول المنظم وسعة التنظيم (Buffer Capacity β):** {buf_info['summary_ar']}")
    if err_info.get("summary_ar"):
        st.caption(f"**التفريق بين نقطة التكافؤ ونقطة النهاية اللونية:** {err_info['distinction_explanation_ar']} — {err_info['summary_ar']}")

    if sim_data.get("step_equivalence_points") and len(sim_data["step_equivalence_points"]) > 1:
        with st.expander("جدول نقاط التكافؤ المرحلية للأحماض والقواعد متعددة البروتونات (Stepwise Equivalence Points)", expanded=True):
            step_rows = [
                {
                    "المرحلة (Step)": s["step"],
                    "حجم التكافؤ V_eq (mL)": s["volume_ml"],
                    "pH نقطة التكافؤ المرحلية": s["ph"],
                    "حجم نصف التكافؤ V_half (mL)": s["half_equivalence_volume_ml"],
                    "pH نصف التكافؤ (pKa)": s["half_equivalence_ph"] if s["half_equivalence_ph"] is not None else "—",
                    "الكاشف الأنسب للمرحلة": s["recommended_indicator_name_ar"],
                }
                for s in sim_data["step_equivalence_points"]
            ]
            st.dataframe(pd.DataFrame(step_rows), use_container_width=True, hide_index=True)

    with st.expander("الحسابات الرياضية خطوة بخطوة (Step-by-Step Titration Calculations)", expanded=False):
        for step_line in sim_data.get("step_by_step_calculation", []):
            st.markdown(f"- {step_line}")

    with st.expander("حدود نموذج المعايرة المبسط والحالات التي لا ينطبق عليها (Model Limitations)", expanded=False):
        for lim_line in sim_data.get("model_limitations_ar", []):
            st.markdown(f"- {lim_line}")

    from reaction_engine.report_generator import ExperimentReportGenerator
    rep_tit = ExperimentReportGenerator.build_titration_report_data(sim_data)
    pdf_tit = ExperimentReportGenerator.generate_experiment_pdf_bytes(rep_tit)
    html_tit = ExperimentReportGenerator.generate_printable_html_report(rep_tit)
    t_dl1, t_dl2 = st.columns(2)
    with t_dl1:
        st.download_button(
            "تحميل تقرير المعايرة بصيغة PDF (Download Titration PDF Report)",
            data=pdf_tit,
            file_name=f"moleculens_titration_{selected_acid_key}_{selected_base_key}.pdf",
            mime="application/pdf",
            key="tit_download_pdf_btn",
        )
    with t_dl2:
        st.download_button(
            "تحميل تقرير المعايرة للطباعة (Printable HTML Report)",
            data=html_tit.encode("utf-8"),
            file_name=f"moleculens_titration_{selected_acid_key}_{selected_base_key}.html",
            mime="text/html",
            key="tit_download_html_btn",
        )

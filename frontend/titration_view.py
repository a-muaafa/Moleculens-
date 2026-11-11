"""
Titration Lab frontend view for MolecuLens.
Interactive acid-base titration simulator with Plotly curves, buffer modeling,
equivalence point detection, and indicator transition bands.
"""

import streamlit as st
import html
from titration.titration_engine import TitrationEngine


def render_titration_view(lang: str = "ar"):
    """Renders the Titration Virtual Laboratory interface."""
    is_ar = lang == "ar"

    header_title = "مختبر المعايرة الافتراضي (Titration Lab)" if is_ar else "Virtual Titration Laboratory"
    header_desc = (
        "محاكاة حركية وديناميكية لتفاعلات المعايرة حمض-قاعدة مع الحسابات الترموديناميكية لنقاط التكافؤ ومنحنيات pH."
        if is_ar
        else "Thermodynamic simulation of acid-base titration reactions with pH curves, buffer capacity, and equivalence calculations."
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
        st.subheader("1. " + ("تحديد الحمض والقاعدة" if is_ar else "Select Acid & Base"))
        acid_options = list(acids.keys())
        acid_labels = {
            k: f"{acids[k]['name_ar' if is_ar else 'name_en']} ({k}) - {'قوي' if acids[k]['type'] == 'strong' else 'ضعيف'}"
            for k in acid_options
        }
        selected_acid_key = st.selectbox(
            "الحمض (Acid)" if is_ar else "Acid",
            options=acid_options,
            format_func=lambda k: acid_labels[k],
            index=0,
        )

        base_options = list(bases.keys())
        base_labels = {
            k: f"{bases[k]['name_ar' if is_ar else 'name_en']} ({k}) - {'قوي' if bases[k]['type'] == 'strong' else 'ضعيف'}"
            for k in base_options
        }
        selected_base_key = st.selectbox(
            "القاعدة (Base)" if is_ar else "Base",
            options=base_options,
            format_func=lambda k: base_labels[k],
            index=0,
        )

        mode_options = [
            ("acid_with_base", "معايرة حمض بواسطة قاعدة (Analyte: Acid, Titrant: Base)" if is_ar else "Titrate Acid with Base"),
            ("base_with_acid", "معايرة قاعدة بواسطة حمض (Analyte: Base, Titrant: Acid)" if is_ar else "Titrate Base with Acid"),
        ]
        titration_mode = st.selectbox(
            "مسار المعايرة" if is_ar else "Titration Mode",
            options=[m[0] for m in mode_options],
            format_func=lambda val: next(m[1] for m in mode_options if m[0] == val),
            index=0,
        )

    with col2:
        st.subheader("2. " + ("الحجوم والتراكيز" if is_ar else "Volumes & Concentrations"))
        analyte_label = "الحمض المجهول/المعاير" if titration_mode == "acid_with_base" else "القاعدة المجهولة/المعايرة"
        analyte_v_label = f"حجم {analyte_label} (مل)" if is_ar else "Analyte Volume (mL)"
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

        ind_keys = list(indicators.keys())
        selected_indicator = st.selectbox(
            "الكاشف اللوني (Chemical Indicator)" if is_ar else "Chemical Indicator",
            options=ind_keys,
            format_func=lambda k: f"{indicators[k]['name_ar' if is_ar else 'name_en']} (pH {indicators[k]['transition_range'][0]}-{indicators[k]['transition_range'][1]})",
            index=0,
        )

    # Perform Titration Simulation
    if titration_mode == "acid_with_base":
        sim_data = TitrationEngine.simulate_titration(
            acid_key=selected_acid_key,
            base_key=selected_base_key,
            acid_conc_m=analyte_conc_m,
            acid_volume_ml=analyte_volume_ml,
            base_conc_m=titrant_conc_m,
        )
    else:
        sim_data = TitrationEngine.simulate_titration(
            acid_key=selected_acid_key,
            base_key=selected_base_key,
            acid_conc_m=titrant_conc_m,
            acid_volume_ml=25.0,
            base_conc_m=analyte_conc_m,
        )

    v_eq = sim_data["v_eq_ml"]
    ph_eq = sim_data["ph_eq"]
    half_eq = sim_data.get("half_v_eq_ml")
    half_ph = sim_data.get("half_ph")

    st.markdown("---")

    # Step-by-step Interactive Titration addition
    st.subheader("3. " + ("محاكاة إضافة المحلول خطوة بخطوة" if is_ar else "Interactive Progressive Addition"))
    max_slider_vol = float(sim_data["volumes_ml"][-1])
    current_vol = st.slider(
        "حجم المحلول المضاف من السحاحة (مل)" if is_ar else "Volume Added from Burette (mL)",
        min_value=0.0,
        max_value=max_slider_vol,
        value=float(min(v_eq, max_slider_vol)),
        step=0.2,
    )

    # Calculate current pH at slider volume
    current_state = TitrationEngine.calculate_point_ph(
        v_add=current_vol,
        v_eq=v_eq,
        v_acid=analyte_volume_ml,
        ca=analyte_conc_m,
        cb=titrant_conc_m,
        acid_info=acids[selected_acid_key],
        base_info=bases[selected_base_key],
    )
    current_ph = current_state["ph"]

    # Metrics Display
    m1, m2, m3, m4 = st.columns(4)
    with m1:
        st.metric(
            label="نقطة التكافؤ V_eq" if is_ar else "Equivalence Vol",
            value=f"{v_eq:.2f} mL",
        )
    with m2:
        st.metric(
            label="درجة حموضة التكافؤ pH_eq" if is_ar else "Equivalence pH",
            value=f"{ph_eq:.2f}",
        )
    with m3:
        st.metric(
            label="الحجم المضاف حالياً" if is_ar else "Current Added",
            value=f"{current_vol:.2f} mL",
        )
    with m4:
        st.metric(
            label="الرقم الهيدروجيني pH الحالي" if is_ar else "Current pH",
            value=f"{current_ph:.2f}",
            delta=f"{current_state['region']}",
        )

    # Build and render Plotly Chart
    fig = TitrationEngine.create_plotly_curve(
        sim_data=sim_data,
        current_v=current_vol,
        indicator_key=selected_indicator,
    )
    st.plotly_chart(fig, use_container_width=True)

    # Indicator Suitability Assessment
    ind_range = indicators[selected_indicator]["transition_range"]
    is_suitable = ind_range[0] <= ph_eq <= ind_range[1] or abs(ph_eq - ((ind_range[0] + ind_range[1]) / 2.0)) <= 1.5

    suit_color = "#10b981" if is_suitable else "#f59e0b"
    suit_text = (
        ("الكاشف ملائم لهذه المعايرة لأن مدى تغير لونه يتقاطع مع القفزة الهيدروجينية لنقطة التكافؤ." if is_ar
         else "Selected indicator is suitable: its transition interval encompasses the sharp equivalence pH inflection.")
        if is_suitable
        else ("الكاشف غير مثالي؛ يفضل اختيار كاشف يتطابق مدى تحوله اللوني بدقة مع pH نقطة التكافؤ." if is_ar
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

    # Scientific Breakdown of the Titration
    st.markdown("---")
    st.subheader("4. " + ("التفسير العلمي ومراحل المنحنى" if is_ar else "Scientific Breakdown"))

    desc_col1, desc_col2 = st.columns(2)
    with desc_col1:
        st.markdown(
            f"""
            **{('المعادلة الكيميائية:' if is_ar else 'Chemical Reaction:')}**
            `{acids[selected_acid_key]['name_en']} + {bases[selected_base_key]['name_en']} -> Salt + H2O`

            - **{('نوع الحمض:' if is_ar else 'Acid Type:')}** {acids[selected_acid_key]['type'].capitalize()} ({'حمض قوي' if acids[selected_acid_key]['type'] == 'strong' else 'حمض ضعيف'})
            - **{('نوع القاعدة:' if is_ar else 'Base Type:')}** {bases[selected_base_key]['type'].capitalize()} ({'قاعدة قوية' if bases[selected_base_key]['type'] == 'strong' else 'قاعدة ضعيفة'})
            - **{('الرقم الهيدروجيني الأولي:' if is_ar else 'Initial pH:')}** {sim_data['initial_ph']:.2f}
            """
        )

    with desc_col2:
        if half_eq is not None and half_ph is not None:
            st.markdown(
                f"""
                **{('نقطة نصف التكافؤ (Half-Equivalence):' if is_ar else 'Half-Equivalence Point:')}**
                - **V_half:** {half_eq:.2f} mL
                - **pH_half:** {half_ph:.2f} *(pH = pKa)*
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


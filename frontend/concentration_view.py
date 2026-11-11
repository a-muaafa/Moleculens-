"""
Concentration and solution preparation calculator frontend view for MolecuLens.
Supports Molarity, Normality, ppm/ppb, Percentage (% w/w, % w/v, % v/v), and Dilution.
"""

import streamlit as st
import html
from calculators.concentration_calculator import ConcentrationCalculator


def render_concentration_view(lang: str = "ar"):
    """Renders the comprehensive solution concentration calculator view."""
    is_ar = lang == "ar"

    header_title = "حاسبة التراكيز وإعداد المحاليل (Concentration Lab)" if is_ar else "Solution Concentration & Dilution Calculator"
    header_desc = (
        "أداة معملية دقيقة لحساب المولارية، العيارية، النسب المئوية، أجزاء المليون، وتخفيف المحاليل مع بروتوكولات التحضير المعملية."
        if is_ar
        else "Laboratory-grade calculator for Molarity, Normality, %, ppm/ppb, and Dilution with step-by-step preparation protocols."
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

    tabs = st.tabs([
        "المولارية (M)" if is_ar else "Molarity (M)",
        "العيارية (N)" if is_ar else "Normality (N)",
        "أجزاء المليون (ppm/ppb)" if is_ar else "ppm / ppb",
        "النسبة المئوية (%)" if is_ar else "Percentage (%)",
        "قانون التخفيف (C1V1 = C2V2)" if is_ar else "Dilution (C1V1 = C2V2)",
    ])

    # Tab 1: Molarity
    with tabs[0]:
        st.subheader("حساب المولارية وتحضير المحاليل" if is_ar else "Molarity & Solution Preparation")
        col_m1, col_m2 = st.columns([1, 1])

        with col_m1:
            solve_target = st.selectbox(
                "المتغير المراد حسابه" if is_ar else "Target to calculate",
                options=["mass", "molarity", "volume"],
                format_func=lambda x: {
                    "mass": "كتلة المادة المذابة (Mass in g)" if is_ar else "Solute Mass (g)",
                    "molarity": "التركيز المولاري (Molarity in M)" if is_ar else "Molarity (M)",
                    "volume": "حجم المحلول (Volume in L)" if is_ar else "Solution Volume (L)",
                }[x],
            )
            mw = st.number_input("الوزن الجزيئي للمذاب (g/mol)" if is_ar else "Molar Mass (g/mol)", min_value=1.0, value=58.44, step=0.01)

            if solve_target == "mass":
                vol_in = st.number_input("حجم المحلول المراد تحضيره (L)" if is_ar else "Target Volume (L)", min_value=0.001, value=0.500, step=0.05)
                molar_in = st.number_input("التركيز المطلوب (M)" if is_ar else "Target Molarity (M)", min_value=0.0001, value=0.100, step=0.01)
                res = ConcentrationCalculator.calculate_molarity(molar_mass_g_mol=mw, volume_l=vol_in, molarity_m=molar_in)
            elif solve_target == "molarity":
                mass_in = st.number_input("كتلة المذاب الموزونة (g)" if is_ar else "Weighed Solute Mass (g)", min_value=0.0001, value=5.844, step=0.1)
                vol_in = st.number_input("حجم المحلول الإجمالي (L)" if is_ar else "Solution Volume (L)", min_value=0.001, value=1.000, step=0.1)
                res = ConcentrationCalculator.calculate_molarity(mass_g=mass_in, molar_mass_g_mol=mw, volume_l=vol_in)
            else:
                mass_in = st.number_input("كتلة المذاب (g)" if is_ar else "Solute Mass (g)", min_value=0.0001, value=5.844, step=0.1)
                molar_in = st.number_input("التركيز المولاري المطلوب (M)" if is_ar else "Desired Molarity (M)", min_value=0.0001, value=0.100, step=0.01)
                res = ConcentrationCalculator.calculate_molarity(mass_g=mass_in, molar_mass_g_mol=mw, molarity_m=molar_in)

        with col_m2:
            st.markdown(f"**{('النتيجة العلمية:' if is_ar else 'Result:')}**")
            if res.get("success"):
                if solve_target == "mass":
                    st.metric(label="الكتلة المطلوبة" if is_ar else "Required Mass", value=f"{res['mass_g']} g")
                elif solve_target == "molarity":
                    st.metric(label="التركيز المولاري" if is_ar else "Calculated Molarity", value=f"{res['molarity_m']} M")
                else:
                    st.metric(label="الحجم النهائي" if is_ar else "Calculated Volume", value=f"{res['volume_l']} L ({res['volume_ml']} mL)")

                st.info(res["recipe_ar"] if is_ar else res["recipe_en"])
                st.caption(f"القانون المستخدم: `{res['formula']}`")
            else:
                st.error(res.get("error", "Error"))

    # Tab 2: Normality
    with tabs[1]:
        st.subheader("حساب العيارية (Normality = M * n)" if is_ar else "Normality Calculator")
        col_n1, col_n2 = st.columns(2)
        with col_n1:
            n_molarity = st.number_input("التركيز المولاري (M)" if is_ar else "Molarity (M)", min_value=0.0001, value=0.100, step=0.01, key="norm_m")
            n_equiv = st.number_input(
                "عامل التكافؤ (Valence factor / n_eq)" if is_ar else "Equivalence Factor (n)",
                min_value=1.0,
                max_value=10.0,
                value=2.0,
                step=1.0,
                help="عدد أيونات الهيدروجين القابلة للاستبدال للحمض (مثل H2SO4 = 2) أو الهيدروكسيد للقاعدة (مثل Ca(OH)2 = 2) أو الشحنة للتأكسد/اختزال.",
            )
            norm_res = ConcentrationCalculator.calculate_normality(molarity_m=n_molarity, equivalents_per_mole=n_equiv)

        with col_n2:
            if norm_res.get("success"):
                st.metric(label="العيارية المحسوبة" if is_ar else "Calculated Normality", value=f"{norm_res['normality_n']} N")
                st.success(norm_res["explanation_ar"] if is_ar else norm_res["explanation_en"])
                st.caption(f"القانون: `{norm_res['formula']}`")

    # Tab 3: ppm / ppb
    with tabs[2]:
        st.subheader("أجزاء في المليون وأجزاء في المليار" if is_ar else "ppm & ppb Calculator")
        col_p1, col_p2 = st.columns(2)
        with col_p1:
            solute_mg = st.number_input("كتلة المذاب بالمليجرام (mg)" if is_ar else "Solute Mass (mg)", min_value=0.001, value=15.0, step=1.0)
            sol_vol_l = st.number_input("حجم المحلول المائي باللتر (L)" if is_ar else "Solution Volume (L)", min_value=0.01, value=1.0, step=0.1, key="ppm_vol")
            ppm_res = ConcentrationCalculator.calculate_ppm_ppb(solute_mg=solute_mg, solution_volume_l=sol_vol_l)
        with col_p2:
            if ppm_res.get("success"):
                st.metric(label="جزء في المليون (ppm)" if is_ar else "ppm", value=f"{ppm_res['ppm']} ppm")
                st.metric(label="جزء في المليار (ppb)" if is_ar else "ppb", value=f"{ppm_res['ppb']} ppb")
                st.info(ppm_res["explanation_ar"] if is_ar else ppm_res["explanation_en"])

    # Tab 4: Percentage
    with tabs[3]:
        st.subheader("التركيز بالنسبة المئوية" if is_ar else "Percentage Concentration")
        col_pct1, col_pct2 = st.columns(2)
        with col_pct1:
            pct_type = st.selectbox(
                "نوع النسبة المئوية" if is_ar else "Percentage Type",
                options=["w/v", "w/w", "v/v"],
                format_func=lambda x: {
                    "w/v": "% w/v (كتلة مذاب جم / حجم محلول مل)" if is_ar else "% w/v (g solute / mL solution)",
                    "w/w": "% w/w (كتلة مذاب جم / كتلة محلول جم)" if is_ar else "% w/w (g solute / g solution)",
                    "v/v": "% v/v (حجم مذاب مل / حجم محلول مل)" if is_ar else "% v/v (mL solute / mL solution)",
                }[x],
            )
            solute_q = st.number_input("كمية المذاب" if is_ar else "Solute Quantity", min_value=0.001, value=5.0, step=0.5)
            solution_q = st.number_input("كمية المحلول الكلي" if is_ar else "Total Solution Quantity", min_value=0.01, value=100.0, step=10.0)
            pct_res = ConcentrationCalculator.calculate_percentage(percent_type=pct_type, solute_amount=solute_q, solution_amount=solution_q)
        with col_pct2:
            if pct_res.get("success"):
                st.metric(label="النسبة المئوية" if is_ar else "Percentage", value=f"{pct_res['percentage']:.2f} %")
                st.info(pct_res["explanation_ar"] if is_ar else pct_res["explanation_en"])
                st.caption(f"القانون: `{pct_res['formula']}`")
            else:
                st.error(pct_res.get("error", "Error"))

    # Tab 5: Dilution
    with tabs[4]:
        st.subheader("قانون التخفيف (C1 * V1 = C2 * V2)" if is_ar else "Dilution Formula (C1 * V1 = C2 * V2)")
        col_d1, col_d2 = st.columns(2)
        with col_d1:
            dil_target = st.selectbox(
                "المتغير المراد حسابه" if is_ar else "Variable to Solve",
                options=["v1", "c2", "v2"],
                format_func=lambda x: {
                    "v1": "الحجم المسحوب من المحلول المركز (V1)" if is_ar else "Stock Volume to Take (V1)",
                    "c2": "تركيز المحلول المخفف النهائي (C2)" if is_ar else "Final Concentration (C2)",
                    "v2": "الحجم النهائي المطلوب الوصول إليه (V2)" if is_ar else "Final Total Volume (V2)",
                }[x],
            )
            c1_val = st.number_input("التركيز الأولي C1 (Stock Conc)" if is_ar else "Initial Stock Conc (C1)", min_value=0.0001, value=5.0, step=0.5)
            if dil_target == "v1":
                c2_val = st.number_input("التركيز المخفف المطلوب C2" if is_ar else "Target Conc (C2)", min_value=0.0001, value=0.5, step=0.1)
                v2_val = st.number_input("الحجم النهائي المراد تحضيره V2 (مل)" if is_ar else "Final Volume V2 (mL)", min_value=0.1, value=250.0, step=10.0)
                dil_res = ConcentrationCalculator.calculate_dilution(c1=c1_val, c2=c2_val, v2=v2_val)
            elif dil_target == "c2":
                v1_val = st.number_input("الحجم المسحوب V1 (مل)" if is_ar else "Stock Volume Taken V1 (mL)", min_value=0.1, value=25.0, step=5.0)
                v2_val = st.number_input("الحجم النهائي بعد إضافة المذيب V2 (مل)" if is_ar else "Final Volume V2 (mL)", min_value=0.1, value=250.0, step=10.0)
                dil_res = ConcentrationCalculator.calculate_dilution(c1=c1_val, v1=v1_val, v2=v2_val)
            else:
                v1_val = st.number_input("الحجم المسحوب V1 (مل)" if is_ar else "Stock Volume Taken V1 (mL)", min_value=0.1, value=50.0, step=5.0)
                c2_val = st.number_input("التركيز المخفف المطلوب C2" if is_ar else "Target Conc (C2)", min_value=0.0001, value=1.0, step=0.2)
                dil_res = ConcentrationCalculator.calculate_dilution(c1=c1_val, v1=v1_val, c2=c2_val)

        with col_d2:
            if dil_res.get("success"):
                if dil_target == "v1":
                    st.metric(label="الحجم المطلوب سحبه (V1)" if is_ar else "Required Stock Volume (V1)", value=f"{dil_res['v1']} mL")
                    st.metric(label="حجم المذيب المضاف" if is_ar else "Solvent to Add", value=f"{dil_res['solvent_to_add']} mL")
                elif dil_target == "c2":
                    st.metric(label="التركيز النهائي (C2)" if is_ar else "Final Concentration (C2)", value=f"{dil_res['c2']}")
                else:
                    st.metric(label="الحجم النهائي (V2)" if is_ar else "Final Volume (V2)", value=f"{dil_res['v2']} mL")
                    st.metric(label="حجم المذيب المضاف" if is_ar else "Solvent to Add", value=f"{dil_res['solvent_to_add']} mL")

                st.success(dil_res["recipe_ar"] if is_ar else dil_res["recipe_en"])
                st.caption(f"القانون: `{dil_res['formula']}`")
            else:
                st.error(dil_res.get("error", "Error"))


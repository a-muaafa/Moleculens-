"""
Drug-Drug Interaction (DDI) frontend view for MolecuLens.
Provides pairwise drug screening, polypharmacy checking, severity categorization,
pharmacological mechanisms, evidence citations, and mandatory clinical disclaimers.
"""

import streamlit as st
import html
from pharmacy.ddi_engine import DDIEngine


def render_ddi_view(lang: str = "ar"):
    """Renders the Drug-Drug Interaction screening module."""
    is_ar = lang == "ar"

    header_title = "فاحص التداخلات الدوائية (Drug-Drug Interactions)" if is_ar else "Drug-Drug Interaction (DDI) Checker"
    header_desc = (
        "تحليل دوائي متقدم للتداخلات بين الأدوية، آلياتها الحركية والديناميكية، ومستوى خطورتها مع توصيات التدبير السريري."
        if is_ar
        else "Pharmacological screening of drug interactions, pharmacokinetic/dynamic mechanisms, severity tiers, and management guidelines."
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

    # Mandatory Academic / Clinical Disclaimer Box
    st.markdown(
        f"""
        <div style="background: rgba(239, 68, 68, 0.1); border: 1px solid rgba(239, 68, 68, 0.3); border-radius: 8px; padding: 14px 18px; margin-bottom: 20px;">
            <strong style="color: #f87171; font-size: 0.95rem;">
                {("تنبيه أكاديمي وسريري إلزامي:" if is_ar else "Mandatory Clinical Disclaimer:")}
            </strong>
            <p style="color: #fca5a5; margin: 6px 0 0 0; font-size: 0.85rem; line-height: 1.5;">
                {html.escape(DDIEngine.DISCLAIMER_AR if is_ar else DDIEngine.DISCLAIMER_EN)}
            </p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    supported_drugs = DDIEngine.get_supported_drugs()
    drug_display_map = {
        "aspirin": "Aspirin (حمض أسيتيل ساليسيليك)",
        "ibuprofen": "Ibuprofen (إيبوبروفين)",
        "paracetamol": "Paracetamol / Acetaminophen (باراسيتامول)",
        "sitagliptin": "Sitagliptin (سيتاغليبتين - Januvia)",
        "metformin": "Metformin (ميتفورمين - Glucophage)",
        "caffeine": "Caffeine (كافيين)",
        "amoxicillin": "Amoxicillin (أموكسيسيلين)",
        "omeprazole": "Omeprazole (أوميبرازول)",
        "warfarin": "Warfarin (وارفارين)",
        "digoxin": "Digoxin (ديجوكسين)",
        "amiodarone": "Amiodarone (أميودارون)",
    }

    tab_pair, tab_multi = st.tabs([
        "فحص ثنائي (دواء مع دواء)" if is_ar else "Pairwise Drug Check",
        "فحص الوصفات المتعددة (Polypharmacy)" if is_ar else "Polypharmacy Multi-Drug Check",
    ])

    with tab_pair:
        col1, col2 = st.columns(2)
        with col1:
            drug_a = st.selectbox(
                "الدواء الأول (Drug A)" if is_ar else "Drug A",
                options=supported_drugs,
                format_func=lambda d: drug_display_map.get(d, d.capitalize()),
                index=0,
                key="ddi_drug_a",
            )
        with col2:
            default_b_idx = 1 if len(supported_drugs) > 1 else 0
            drug_b = st.selectbox(
                "الدواء الثاني (Drug B)" if is_ar else "Drug B",
                options=supported_drugs,
                format_func=lambda d: drug_display_map.get(d, d.capitalize()),
                index=default_b_idx,
                key="ddi_drug_b",
            )

        if st.button("تحليل التداخل الدوائي" if is_ar else "Analyze Interaction", type="primary"):
            res = DDIEngine.check_interaction(drug_a, drug_b)

            if res.get("has_interaction"):
                severity = res["severity"]
                sev_color = {
                    "Major": "#ef4444",
                    "Moderate": "#f59e0b",
                    "Minor": "#3b82f6",
                }.get(severity, "#94a3b8")

                st.markdown(
                    f"""
                    <div style="background: rgba(30, 41, 59, 0.8); border: 2px solid {sev_color}; border-radius: 10px; padding: 20px; margin-top: 16px;">
                        <div style="display: flex; justify-content: space-between; align-items: center; margin-bottom: 14px;">
                            <h3 style="color: #f1f5f9; margin: 0;">
                                {html.escape(drug_a.capitalize())} + {html.escape(drug_b.capitalize())}
                            </h3>
                            <span style="background: {sev_color}; color: white; padding: 6px 14px; border-radius: 20px; font-weight: bold; font-size: 0.9rem;">
                                {html.escape(res['severity_ar'] if is_ar else res['severity'])}
                            </span>
                        </div>
                        <div style="margin-bottom: 12px;">
                            <strong style="color: #60a5fa;">{("الآلية الدوائية (Pharmacological Mechanism):" if is_ar else "Mechanism:")}</strong>
                            <p style="color: #cbd5e1; margin: 4px 0 0 0;">{html.escape(res['mechanism_ar'] if is_ar else res['mechanism_en'])}</p>
                        </div>
                        <div style="margin-bottom: 12px;">
                            <strong style="color: #60a5fa;">{("الأثر السريري (Clinical Effect):" if is_ar else "Clinical Effect:")}</strong>
                            <p style="color: #cbd5e1; margin: 4px 0 0 0;">{html.escape(res['clinical_effect_ar'] if is_ar else res['clinical_effect_en'])}</p>
                        </div>
                        <div style="margin-bottom: 12px;">
                            <strong style="color: #10b981;">{("التدبير السريري الموصى به (Clinical Management):" if is_ar else "Management:")}</strong>
                            <p style="color: #cbd5e1; margin: 4px 0 0 0;">{html.escape(res['management_ar'] if is_ar else res['management_en'])}</p>
                        </div>
                        <div>
                            <strong style="color: #94a3b8; font-size: 0.85rem;">{("المراجع العلمية المعتمدة (Citations):" if is_ar else "Citations:")}</strong>
                            <ul style="color: #94a3b8; font-size: 0.85rem; margin: 4px 0 0 16px; padding: 0;">
                                {''.join(f'<li>{html.escape(c)}</li>' for c in res.get('citations', []))}
                            </ul>
                        </div>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )

            elif res.get("is_same_drug"):
                st.warning(
                    f"{'تحذير: لقد اخترت نفس الدواء في الخانتين (' if is_ar else 'Notice: Same drug selected ('}{drug_a}). "
                    f"{res['mechanism_ar'] if is_ar else res['mechanism_en']}"
                )
            else:
                st.markdown(
                    f"""
                    <div style="background: rgba(30, 41, 59, 0.5); border: 1px solid #10b981; border-radius: 10px; padding: 18px; margin-top: 16px;">
                        <h4 style="color: #34d399; margin: 0 0 8px 0;">
                            {("لا يوجد تداخل مسجل بقاعدة البيانات" if is_ar else "No Interaction Documented")}
                        </h4>
                        <p style="color: #cbd5e1; margin: 0; font-size: 0.9rem;">
                            {html.escape(res['message_ar'] if is_ar else res['message_en'])}
                        </p>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )

    with tab_multi:
        st.write(
            "اختر الأدوية التي يتناولها المريض في نفس الوقت لفحص جميع الأزواج المحتملة:"
            if is_ar
            else "Select multiple medications to screen all concurrent combinations:"
        )
        selected_multi = st.multiselect(
            "قائمة الأدوية الموصوفة" if is_ar else "Prescribed Drugs",
            options=supported_drugs,
            default=["aspirin", "ibuprofen", "paracetamol"],
            format_func=lambda d: drug_display_map.get(d, d.capitalize()),
        )

        if st.button("فحص الوصفة الشاملة" if is_ar else "Screen Polypharmacy", type="primary"):
            multi_res = DDIEngine.screen_polypharmacy(selected_multi)

            c_total, c_inter, c_maj = st.columns(3)
            with c_total:
                st.metric("عدد الأدوية المفحوصة" if is_ar else "Total Drugs", multi_res["total_drugs"])
            with c_inter:
                st.metric("التداخلات المكتشفة" if is_ar else "Interactions Found", multi_res["total_interactions_found"])
            with c_maj:
                st.metric("تداخلات عالية الخطورة" if is_ar else "Major Alerts", multi_res["severity_summary"]["Major"])

            if multi_res["interactions"]:
                for item in multi_res["interactions"]:
                    sev = item["severity"]
                    col_bar = "#ef4444" if sev == "Major" else ("#f59e0b" if sev == "Moderate" else "#3b82f6")
                    st.markdown(
                        f"""
                        <div style="background: rgba(15, 23, 42, 0.7); border-left: 5px solid {col_bar}; border-radius: 6px; padding: 12px 16px; margin: 10px 0;">
                            <strong style="color: #f1f5f9;">{html.escape(item['drug_a'].capitalize())} + {html.escape(item['drug_b'].capitalize())}</strong>
                            <span style="color: {col_bar}; font-size: 0.85rem; margin-right: 12px;">[{html.escape(item['severity_ar'] if is_ar else item['severity'])}]</span>
                            <p style="color: #cbd5e1; margin: 4px 0 0 0; font-size: 0.88rem;">{html.escape(item['clinical_effect_ar'] if is_ar else item['clinical_effect_en'])}</p>
                        </div>
                        """,
                        unsafe_allow_html=True,
                    )
            else:
                st.success("لم يتم العثور على تداخلات مسجلة بين هذه الأدوية في قاعدة البيانات." if is_ar else "No documented interactions found between the selected drugs.")


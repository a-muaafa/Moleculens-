"""
AI Scientific Assistant frontend view for MolecuLens.
Provides pedagogical explanations, step-by-step breakdowns of chemical properties,
Lipinski rule interpretation, and interactive student Q&A grounded strictly in calculated data.
"""

import streamlit as st
import html
from education.ai_assistant import AIScientificAssistant


def render_ai_assistant_component(compound_data: dict = None, lang: str = "ar"):
    """
    Renders the AI Scientific Assistant component either embedded inside Explore View
    or as a standalone section.
    """
    is_ar = lang == "ar"

    box_title = "المساعد العلمي الذكي (AI Scientific Assistant)" if is_ar else "AI Scientific Assistant"
    box_desc = (
        "تفسير علمي مدعوم بالذكاء الاصطناعي يشرح الخصائص الدوائية والمفاهيم الكيميائية خطوة بخطوة باللغة العربية دون تخمين أو اختراع بيانات."
        if is_ar
        else "AI-guided scientific explanations breaking down physicochemical metrics and drug-likeness without data hallucination."
    )

    st.markdown(
        f"""
        <div class="ml-module-header" style="margin-bottom: 18px;">
            <div style="display: flex; align-items: center; gap: 10px; margin-bottom: 6px;">
                <span class="ml-ai-tag">AI</span>
                <h3 style="margin: 0; font-size: 1.28rem;">{html.escape(box_title)}</h3>
            </div>
            <p>{html.escape(box_desc)}</p>
        </div>
        """,
        unsafe_allow_html=True,
    )

    assistant = AIScientificAssistant()

    tab_explain, tab_breakdown, tab_qa = st.tabs([
        "تفسير الخصائص الدوائية" if is_ar else "Property Explanations",
        "التحليل التفصيلي للمركب" if is_ar else "Step-by-Step Breakdown",
        "اسأل المساعد العلمي" if is_ar else "Ask the Assistant",
    ])

    with tab_explain:
        prop_options = [
            ("mw", "الوزن الجزيئي (Molecular Weight)" if is_ar else "Molecular Weight (MW)"),
            ("logp", "معامل التقسيم (LogP - Lipophilicity)" if is_ar else "Partition Coefficient (LogP)"),
            ("tpsa", "مساحة السطح القطبية (TPSA)" if is_ar else "Topological Polar Surface Area (TPSA)"),
            ("hbd", "مانحات روابط الهيدروجين (H-Bond Donors)" if is_ar else "H-Bond Donors (HBD)"),
            ("hba", "مستقبلات روابط الهيدروجين (H-Bond Acceptors)" if is_ar else "H-Bond Acceptors (HBA)"),
            ("rotatable_bonds", "الروابط القابلة للدوران (Rotatable Bonds)" if is_ar else "Rotatable Bonds"),
        ]
        selected_prop = st.selectbox(
            "اختر الخاصية لفهم معناها العلمي وأثرها الدوائي:" if is_ar else "Select property to inspect pharmacological impact:",
            options=[p[0] for p in prop_options],
            format_func=lambda k: next(p[1] for p in prop_options if p[0] == k),
            index=0,
            key="ai_prop_select",
        )

        curr_val = None
        if compound_data and "properties" in compound_data:
            props = compound_data["properties"]
            prop_key_map = {
                "mw": "molecular_weight",
                "logp": "logp",
                "tpsa": "tpsa",
                "hbd": "hbd",
                "hba": "hba",
                "rotatable_bonds": "rotatable_bonds",
            }
            curr_val = props.get(prop_key_map.get(selected_prop))

        explanation = assistant.explain_property(selected_prop, value=curr_val, lang=lang)

        st.markdown(
            f"""
            <div class="ml-ai-card">
                <h4>{html.escape(explanation['property_name'])}</h4>
                <p style="color: #334155; font-size: 0.95rem; line-height: 1.7; margin-bottom: 12px;">
                    {html.escape(explanation['definition'])}
                </p>
                <div class="ml-ai-impact-box">
                    <strong style="color: #0E5A6D; font-size: 0.86rem;">{("الأهمية الدوائية والحركية:" if is_ar else "Pharmacological Impact:")}</strong>
                    <p style="color: #1E293B; margin: 4px 0 0 0; font-size: 0.92rem; line-height: 1.65;">{html.escape(explanation['pharmacological_impact'])}</p>
                </div>
                <div style="color: #52657A; font-size: 0.88rem; margin-top: 8px;">
                    <strong>{("الحد المثالي لقاعدة ليبينسكي:" if is_ar else "Lipinski Optimal Threshold:")}</strong>
                    <span style="color: #B45309; font-weight: 700;">{html.escape(explanation['optimal_range'])}</span>
                </div>
                {f'<div style="margin-top: 8px; font-weight: 700; color: #047857;">{html.escape(explanation.get("context_evaluation", ""))}</div>' if explanation.get("context_evaluation") else ''}
            </div>
            """,
            unsafe_allow_html=True,
        )

    with tab_breakdown:
        if compound_data:
            steps = assistant.generate_step_by_step_breakdown(compound_data, lang=lang)
            for idx, step in enumerate(steps, 1):
                st.markdown(
                    f"""
                    <div class="ml-ai-step">
                        <span style="background: #0E5A6D; color: white; padding: 2px 9px; border-radius: 6px; font-size: 0.76rem; font-weight: 700;">
                            {("خطوة" if is_ar else "Step")} {idx}
                        </span>
                        <strong style="color: #0D2B42; margin-right: 8px; font-size: 0.96rem;">{html.escape(step['title'])}</strong>
                        <p style="color: #475569; margin: 6px 0 0 0; font-size: 0.92rem; line-height: 1.65;">{html.escape(step['detail'])}</p>
                    </div>
                    """,
                    unsafe_allow_html=True,
                )
        else:
            st.info(
                "يرجى تحديد أو البحث عن مركب في المختبر لعرض التحليل التدريجي الكامل لخواصه وتركيبه."
                if is_ar
                else "Please search for or select a compound to view its complete step-by-step breakdown."
            )

    with tab_qa:
        preset_questions = [
            "لماذا يُعتبر معامل التقسيم LogP مهماً في تصميم الدواء؟" if is_ar else "Why is LogP critical in drug design?",
            "كيف تؤثر الروابط القابلة للدوران على التوافر الحيوي؟" if is_ar else "How do rotatable bonds affect oral bioavailability?",
            "ما أهمية مساحة السطح القطبية TPSA في عبور الحاجز الدموي الدماغي؟" if is_ar else "Why is TPSA vital for Blood-Brain Barrier penetration?",
            "ما هي قاعدة ليبينسكي الخماسية (Rule of Five) ولماذا سُميت بذلك؟" if is_ar else "What is Lipinski's Rule of Five and why was it named so?",
            "ما الفرق بين مانحات روابط الهيدروجين (HBD) ومستقبلاتها (HBA)؟" if is_ar else "What is the difference between H-bond donors (HBD) and acceptors (HBA)?",
            "كيف يؤثر الوزن الجزيئي على عبور المركب للأغشية الحيوية؟" if is_ar else "How does molecular weight influence biological membrane permeability?",
        ]

        chosen_preset = st.selectbox(
            "اختر سؤالاً علمياً لاستكشاف الإجابة:" if is_ar else "Select a scientific question to inspect explanation:",
            options=["-- اختر سؤالاً --" if is_ar else "-- Choose a question --"] + preset_questions,
            key="ai_preset_q",
        )

        active_question = chosen_preset if chosen_preset != ("-- اختر سؤالاً --" if is_ar else "-- Choose a question --") else ""

        if active_question:
            answer = assistant.answer_student_question(active_question, compound_data=compound_data, lang=lang)
            st.markdown(
                f"""
                <div class="ml-ai-answer">
                    <div style="color: #0E5A6D; font-size: 0.88rem; margin-bottom: 6px;">
                        <strong>{("السؤال:" if is_ar else "Question:")}</strong> {html.escape(active_question)}
                    </div>
                    <div style="color: #1E293B; font-size: 0.95rem; line-height: 1.75;">
                        {html.escape(answer)}
                    </div>
                </div>
                """,
                unsafe_allow_html=True,
            )

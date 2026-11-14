"""
Molecular Intelligence frontend view for MolecuLens.
Features Morgan Fingerprint similarity, structure modification with valence validation,
and transparent LogP prediction evaluation with benchmark metrics.
"""

from __future__ import annotations

from typing import Any
import html
import pandas as pd
import streamlit as st

from chemistry_engine.formula_analyzer import format_chemical_formula_display
from chemistry_engine.molecular_intelligence import MolecularIntelligence
from search.local_catalog import LOCAL_SCIENTIFIC_CATALOG


def _esc(val: Any) -> str:
    if val is None:
        return ""
    return html.escape(str(val), quote=True)


def render_molecular_intelligence_section() -> None:
    """Renders the comprehensive Molecular Intelligence interface."""
    st.markdown(
        """
<div class="ml-card">
  <h3>الذكاء الجزيئي — Molecular Intelligence &amp; Structure Studio</h3>
  <p>
    وحدة حاسوبية متقدمة مبنية على خوارزميات الكيمياء المعلوماتية (Cheminformatics):
    <strong>مقارنة التشابه الجزيئي (Molecular Similarity)</strong> عبر بصمات مورغان وبُعد تانيموتو،
    <strong>تعديل البنية الكيميائية (Structure Modification)</strong> مع التحقق الصارم من التكافؤ،
    و<strong>تقييم شفافية تنبؤ LogP</strong> مع فواصل الثقة ومقاييس أداء النموذج.
  </p>
</div>
""",
        unsafe_allow_html=True,
    )

    # Filter catalog entities that have canonical_smiles
    compounds = [
        c for c in LOCAL_SCIENTIFIC_CATALOG
        if c.get("canonical_smiles") and len(c.get("canonical_smiles", "")) > 1
        and c.get("entity_type") != "element"
    ]

    if not compounds:
        st.warning("لا توجد مركبات مسجلة بصيغ بنيوية في الكتالوج.")
        return

    comp_options = {c["key"]: f"{c['name_ar']} | {c['name_en']} ({format_chemical_formula_display(c['symbol_or_formula'])})" for c in compounds}
    comp_keys = list(comp_options.keys())

    tab_sim, tab_mod, tab_logp = st.tabs([
        "1. التشابه الجزيئي (Molecular Similarity)",
        "2. تعديل البنية الجزيئية (Structure Modification)",
        "3. تقييم نموذج تنبؤ LogP (LogP Evaluation)",
    ])

    # Tab 1: Similarity
    with tab_sim:
        st.markdown(
            """
<div class="ml-card">
  <h4>مقارنة التشابه الجزيئي عبر بصمات مورغان (Morgan Fingerprints &amp; Tanimoto)</h4>
  <p>
    تعتمد الخوارزمية على توليد البصمات الجزيئية الدائرية (Circular Fingerprints بقطر نصفه 2 وبطول 2048 بت)،
    ثم حساب معامل تشابه تانيموتو (Tanimoto Coefficient) مقارنة بجميع مركبات الكتالوج العلمي.
  </p>
</div>
""",
            unsafe_allow_html=True,
        )

        chosen_sim_key = st.selectbox(
            "اختر المركب المراد البحث عن أشباهه البنيوية في الكتالوج:",
            options=comp_keys,
            format_func=lambda k: comp_options.get(k, k),
            key="sim_select_key",
        )

        chosen_comp = next((c for c in compounds if c["key"] == chosen_sim_key), None)
        if chosen_comp:
            smiles = chosen_comp["canonical_smiles"]
            st.caption(f"التمثيل البنيوي SMILES: `{smiles}`")

            sims = MolecularIntelligence.find_most_similar(smiles, top_n=6)
            if sims:
                st.markdown(f"**أقرب المركبات بنيوياً إلى «{_esc(chosen_comp['name_ar'])}» في الكتالوج المعتمد:**")
                cols = st.columns(3)
                for idx, item in enumerate(sims):
                    with cols[idx % 3]:
                        pct = item["similarity_pct"]
                        bar_color = "#059669" if pct >= 50 else ("#0f4c5c" if pct >= 25 else "#475569")
                        st.markdown(
                            f"""
<div class="ml-card" style="border-top:4px solid {bar_color};">
  <div style="display:flex;justify-content:space-between;align-items:center;">
    <strong style="color:#0f4c5c;font-size:1.05rem;">{_esc(item['name_ar'])}</strong>
    <span class="ml-tri-badge" style="background:rgba(15,76,92,0.12);color:{bar_color};font-weight:bold;">
      {pct}% تشابه
    </span>
  </div>
  <p style="margin:0.25rem 0;color:#64748b;font-size:0.88rem;">{_esc(item['name_en'])} ({format_chemical_formula_display(item['formula'])})</p>
  <span class="ml-badge" style="background:#334155;margin-bottom:0.4rem;">{_esc(item['subclass_label_ar'])}</span>
  <p style="margin:0;font-size:0.8rem;color:#94a3b8;font-family:monospace;word-break:break-all;">{_esc(item['smiles'])}</p>
</div>
""",
                            unsafe_allow_html=True,
                        )
            else:
                st.info("لم يتم العثور على مركبات تشترك في بصمات قريبة كافية مع هذا المركب في الكتالوج.")

    # Tab 2: Structure Modification
    with tab_mod:
        st.markdown(
            """
<div class="ml-card">
  <h4>محرر تعديل البنية الجزيئية واستبدال المجموعات الوظيفية</h4>
  <p>
    اختر مركباً أساسياً ونوع التحويل الكيميائي المستهدف (مثل استبدال الهيدروكسيل بهالوجين أو أسترة الكربوكسيل)،
    ليقوم المحرك بتطبيق التعديل والتحقق من صحة التكافؤ واستقرار البنية، مع مقارنة الخصائص الفيزيوكيميائية قبل وبعد التعديل.
  </p>
</div>
""",
            unsafe_allow_html=True,
        )

        c_m1, c_m2 = st.columns(2)
        with c_m1:
            mod_comp_key = st.selectbox(
                "المركب الأساسي للتعديل:",
                options=comp_keys,
                format_func=lambda k: comp_options.get(k, k),
                key="mod_select_comp",
            )
        with c_m2:
            mod_keys = list(MolecularIntelligence.MODIFICATIONS.keys())
            mod_labels = {k: MolecularIntelligence.MODIFICATIONS[k]["name_ar"] for k in mod_keys}
            chosen_mod_key = st.selectbox(
                "نوع التعديل البنيوي المستهدف:",
                options=mod_keys,
                format_func=lambda k: mod_labels.get(k, k),
                key="mod_select_type",
            )

        base_entity = next((c for c in compounds if c["key"] == mod_comp_key), None)
        if base_entity:
            mod_result = MolecularIntelligence.modify_structure(base_entity["canonical_smiles"], chosen_mod_key)
            if mod_result["success"]:
                st.success(mod_result["chemical_validity"]["status_ar"])

                b = mod_result["before"]
                a = mod_result["after"]
                d = mod_result["deltas"]

                col_b, col_a = st.columns(2)
                with col_b:
                    st.markdown(
                        f"""
<div class="ml-card" style="border-right:4px solid #64748b;">
  <h4>المركب قبل التعديل (Original Molecule)</h4>
  <p><strong>الاسم:</strong> {_esc(base_entity['name_ar'])} ({_esc(base_entity['name_en'])})</p>
  <p><strong>الصيغة الجزيئية:</strong> {format_chemical_formula_display(b['formula'])}</p>
  <p><strong>الكتلة الجزيئية (MW):</strong> {b['molecular_weight']} g/mol</p>
  <p><strong>معامل LogP:</strong> {b['logp']}</p>
  <p><strong>مساحة السطح القطبي (TPSA):</strong> {b['tpsa']} Å²</p>
  <p><strong>مانحات / مستقبلات H:</strong> {b['hbd']} HBD / {b['hba']} HBA</p>
  <p><strong>روابط الدوران:</strong> {b['rotatable_bonds']}</p>
  <p style="font-family:monospace;font-size:0.8rem;color:#94a3b8;">{_esc(b['smiles'])}</p>
</div>
""",
                        unsafe_allow_html=True,
                    )

                with col_a:
                    st.markdown(
                        f"""
<div class="ml-card" style="border-right:4px solid #0f766e;">
  <h4>المركب بعد التعديل (Modified Molecule)</h4>
  <p><strong>التحويل:</strong> {_esc(mod_result['modification_name_ar'])}</p>
  <p><strong>الصيغة الجزيئية:</strong> {format_chemical_formula_display(a['formula'])}</p>
  <p><strong>الكتلة الجزيئية (MW):</strong> {a['molecular_weight']} g/mol <span style="color:#0f766e;">({'+' if d['mw'] > 0 else ''}{d['mw']})</span></p>
  <p><strong>معامل LogP:</strong> {a['logp']} <span style="color:#0f766e;">({'+' if d['logp'] > 0 else ''}{d['logp']})</span></p>
  <p><strong>مساحة السطح القطبي (TPSA):</strong> {a['tpsa']} Å² <span style="color:#0f766e;">({'+' if d['tpsa'] > 0 else ''}{d['tpsa']})</span></p>
  <p><strong>مانحات / مستقبلات H:</strong> {a['hbd']} HBD / {a['hba']} HBA</p>
  <p><strong>روابط الدوران:</strong> {a['rotatable_bonds']}</p>
  <p style="font-family:monospace;font-size:0.8rem;color:#0f766e;">{_esc(a['smiles'])}</p>
</div>
""",
                        unsafe_allow_html=True,
                    )

                st.markdown(
                    f"""
<div class="ml-card" style="border-left:4px solid #3b82f6;">
  <strong style="color:#1d4ed8;">التحليل البيداغوجي للأثر الدوائي والفيزيوكيميائي:</strong>
  <p style="margin-top:0.35rem;margin-bottom:0;">{_esc(mod_result['impact_summary_ar'])}</p>
</div>
""",
                    unsafe_allow_html=True,
                )
            else:
                st.warning(mod_result.get("error_ar", "تعذر تنفيذ التعديل."))

    # Tab 3: LogP Prediction Evaluation
    with tab_logp:
        st.markdown(
            """
<div class="ml-card">
  <h4>تقييم نموذج التنبؤ بمعامل التوزع (LogP Prediction Evaluation &amp; Transparency)</h4>
  <p>
    تعرض هذه اللوحة تقييماً علمياً لنموذج التجزئة الذرية المعتمد (Wildman-Crippen SlogP)،
    مع إيضاح أن القيمة المحسوبة هي قيمة تقديرية حاسوبية مصحوبة بفاصل ثقة ومقاييس أداء إحصائية معلنة.
  </p>
</div>
""",
            unsafe_allow_html=True,
        )

        logp_comp_key = st.selectbox(
            "اختر المركب المراد تقييم تنبؤ LogP له:",
            options=comp_keys,
            format_func=lambda k: comp_options.get(k, k),
            key="logp_select_comp",
        )

        target_comp = next((c for c in compounds if c["key"] == logp_comp_key), None)
        if target_comp:
            eval_res = MolecularIntelligence.evaluate_logp_prediction(target_comp["canonical_smiles"])
            if eval_res["available"]:
                m1, m2, m3, m4 = st.columns(4)
                m1.metric("القيمة التقديرية (LogP)", str(eval_res["predicted_logp"]))
                m2.metric("فاصل الثقة 95%", eval_res["confidence_interval_95"]["range_text"])
                m3.metric("معامل التحديد R²", str(eval_res["r2_score"]))
                m4.metric("جذر متوسط مربع الخطأ RMSE", f"±{eval_res['rmse']}")

                st.info(eval_res["pedagogical_note_ar"])

                st.markdown(
                    f"""
<div class="ml-card">
  <h5>معلومات النموذج والحدود العلمية (Model Methodology &amp; Limitations)</h5>
  <p><strong>اسم النموذج والمنهجية:</strong> {_esc(eval_res['methodology_ar'])}</p>
  <p><strong>المرجع الأكاديمي القياسي:</strong> <code>{_esc(eval_res['citation'])}</code></p>
  <p><strong>حدود الاستخدام العلمية:</strong> {_esc(eval_res['limitations_ar'])}</p>
</div>
""",
                    unsafe_allow_html=True,
                )

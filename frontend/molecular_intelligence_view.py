"""
Molecular Intelligence frontend view for MolecuLens.
Features Morgan Fingerprint similarity, virtual SMARTS structure modification with
chemical verification and explicit caveats, and strict separation between deterministic
Wildman-Crippen cLogP calculation and an actually trained/tested QSAR LogP model.
"""

from __future__ import annotations

import html
from typing import Any
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
    وحدة حاسوبية مبنية على الكيمياء المعلوماتية (Cheminformatics):
    <strong>مقارنة التشابه الجزيئي (Molecular Similarity)</strong> عبر بصمات مورغان وبُعد تانيموتو،
    <strong>محرر التعديل البنيوي الافتراضي (In-Silico Virtual Editing)</strong> مع الفحص الكيميائي للنواتج،
    و<strong>حساب cLogP وتقييم نموذج QSAR المختبر فعلياً</strong> بفصل صارم بين الحساب الحتمي والتنبؤ الإحصائي.
  </p>
</div>
""",
        unsafe_allow_html=True,
    )

    compounds = [
        c
        for c in LOCAL_SCIENTIFIC_CATALOG
        if c.get("canonical_smiles")
        and len(c.get("canonical_smiles", "")) > 1
        and c.get("entity_type") != "element"
    ]

    if not compounds:
        st.warning("لا توجد مركبات مسجلة بصيغ بنيوية في الكتالوج.")
        return

    comp_options = {
        c["key"]: f"{c['name_ar']} | {c['name_en']} ({format_chemical_formula_display(c['symbol_or_formula'])})"
        for c in compounds
    }
    comp_keys = list(comp_options.keys())

    tab_sim, tab_mod, tab_logp = st.tabs([
        "1. التشابه الجزيئي (Molecular Similarity)",
        "2. التعديل البنيوي الافتراضي (Virtual Editing)",
        "3. حساب cLogP وتقييم نموذج التنبؤ (LogP Evaluation)",
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

    # Tab 2: Virtual Structure Modification
    with tab_mod:
        st.markdown(
            """
<div class="ml-card">
  <h4>محرر التعديل البنيوي الافتراضي (In-Silico Virtual Structure Modification)</h4>
  <p>
    يحاكي هذا المحرر استبدال المجموعات الوظيفية على الرسم الجزيئي عبر قواعد SMARTS، ثم يفحص التكافؤ والروابط وتنبيهات الثبات البنيوي ومطابقة السجل المرجعي.
    التعديلات هنا <strong>افتراضية حاسوبياً</strong> ولا تمثل بالضرورة تفاعلاً تخليقياً واقعياً بخطوة واحدة أو ضماناً لأي تحسن دوائي.
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
                "نوع التعديل الافتراضي المستهدف:",
                options=mod_keys,
                format_func=lambda k: mod_labels.get(k, k),
                key="mod_select_type",
            )

        base_entity = next((c for c in compounds if c["key"] == mod_comp_key), None)
        if base_entity:
            mod_result = MolecularIntelligence.modify_structure(base_entity["canonical_smiles"], chosen_mod_key)
            if mod_result["success"]:
                validity = mod_result["chemical_validity"]
                cat_match = validity.get("catalog_match", {})
                alerts = validity.get("stability_alerts", [])

                st.info(f"{mod_result['modification_nature_ar']}: {validity['status_ar']}")

                if alerts:
                    for al in alerts:
                        st.warning(f"تنبيه ثبات كيميائي على الناتج الافتراضي: {al}")

                if cat_match.get("is_documented_compound"):
                    st.success(
                        f"مطابقة مرجعية: البنية الناتجة تطابق مركباً موثقاً في السجل المحلي: "
                        f"{cat_match.get('name_ar')} ({cat_match.get('name_en')})."
                    )
                else:
                    st.caption("التحقق المرجعي: البنية الناتجة هي مركب افتراضي غير مدرج في السجل المحلي المعتمد (241 مادة).")

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
  <p><strong>معامل cLogP (حسابي):</strong> {b['logp']}</p>
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
  <h4>البنية الافتراضية بعد التعديل (Virtual Modified Molecule)</h4>
  <p><strong>التعديل الافتراضي:</strong> {_esc(mod_result['modification_name_ar'])}</p>
  <p><strong>الصيغة الجزيئية:</strong> {format_chemical_formula_display(a['formula'])}</p>
  <p><strong>الكتلة الجزيئية (MW):</strong> {a['molecular_weight']} g/mol <span style="color:#0f766e;">({'+' if d['mw'] > 0 else ''}{d['mw']})</span></p>
  <p><strong>معامل cLogP (حسابي):</strong> {a['logp']} <span style="color:#0f766e;">({'+' if d['logp'] > 0 else ''}{d['logp']})</span></p>
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
  <strong style="color:#1d4ed8;">تحليل الفروق الحسابية للواصفات الجزيئية (Descriptor Deltas Analysis):</strong>
  <p style="margin-top:0.35rem;margin-bottom:0;">{_esc(mod_result['impact_summary_ar'])}</p>
</div>
""",
                    unsafe_allow_html=True,
                )

                site_info = mod_result.get("site_info")
                if site_info:
                    st.info(site_info["note_ar"])
                    if site_info["alternative_products_smiles"]:
                        with st.expander("نواتج المواقع الأخرى المحتملة (SMILES)", expanded=False):
                            for alt in site_info["alternative_products_smiles"]:
                                st.code(alt)

                caveats = mod_result.get("scientific_caveats_ar", [])
                if caveats:
                    with st.expander("القيود العلمية وعدم الجزم بالنتائج الدوائية (Scientific Caveats)", expanded=False):
                        for cv in caveats:
                            st.markdown(f"- {cv}")
            else:
                st.warning(mod_result.get("error_ar", "تعذر تنفيذ التعديل الافتراضي."))

    # Tab 3: Deterministic cLogP vs Trained QSAR Model Evaluation
    with tab_logp:
        st.markdown(
            """
<div class="ml-card">
  <h4>حساب cLogP الحاسوبي مقابل تقييم نموذج التنبؤ المختبر فعلياً</h4>
  <p>
    تفصل هذه اللوحة بوضوح تام بين مسارين مختلفين علمياً:
    <strong>(1) الحساب الوصفي الحتمي (Wildman-Crippen cLogP)</strong> الذي يجمع مساهمات الذرات وفق معادلة ثابتة دون فاصل ثقة إحصائي،
    و<strong>(2) نموذج الانحدار الإحصائي (Ridge QSAR)</strong> الذي جرى تدريبه واختباره فعلياً على مجموعة بيانات تجريبية موثقة لقياس الأداء الحقيقي (R², RMSE, MAE).
  </p>
</div>
""",
            unsafe_allow_html=True,
        )

        logp_comp_key = st.selectbox(
            "اختر المركب المراد تحليل LogP له:",
            options=comp_keys,
            format_func=lambda k: comp_options.get(k, k),
            key="logp_select_comp",
        )

        target_comp = next((c for c in compounds if c["key"] == logp_comp_key), None)
        if target_comp:
            eval_res = MolecularIntelligence.evaluate_logp_prediction(target_comp["canonical_smiles"])
            if eval_res["available"]:
                det = eval_res["deterministic_clogp"]
                qsar = eval_res["qsar_model_evaluation"]

                st.markdown("### القسم الأول: الحساب الحاسوبي الحتمي (Deterministic Wildman-Crippen cLogP)")
                d1, d2, d3, d4 = st.columns(4)
                d1.metric("cLogP المحسوب حتمياً", str(det["calculated_clogp"]))
                d2.metric("مجموع المساهمات الدهنية (+)", f"+{det['positive_contrib_sum']}")
                d3.metric("مجموع المساهمات القطبية (-)", str(det["negative_contrib_sum"]))
                d4.metric(
                    "LogP التجريبي الموثق" if det["experimental_logp"] is not None else "المرجع التجريبي",
                    str(det["experimental_logp"]) if det["experimental_logp"] is not None else "غير مسجل محلياً",
                )
                st.caption(det["confidence_interval_note_ar"])

                st.divider()
                st.markdown("### القسم الثاني: نموذج التنبؤ المدرَّب والمختبر فعلياً (Trained & Holdout-Tested Ridge QSAR)")
                if qsar.get("model_available"):
                    q1, q2, q3, q4 = st.columns(4)
                    q1.metric("تنبؤ نموذج QSAR", str(qsar["qsar_predicted_logp"]))
                    q2.metric("R² (اختبار محجوز / حجب الفئة الكيميائية)", f"{qsar['test_r2']} / {qsar['group_cv_r2']}")
                    q3.metric("RMSE (اختبار محجوز / حجب الفئة الكيميائية)", f"±{qsar['test_rmse']} / ±{qsar['group_cv_rmse']}")
                    q4.metric(
                        "حجم التدريب / الاختبار (الإجمالي)",
                        f"{qsar['n_train']} / {qsar['n_test']} ({qsar['n_total']})",
                    )

                    if qsar["in_applicability_domain"]:
                        st.success(qsar["applicability_note_ar"])
                    else:
                        st.warning(qsar["applicability_note_ar"])

                    st.info(eval_res["pedagogical_note_ar"])

                    st.markdown(
                        f"**نتائج مجموعة الاختبار المستقلة الفعلية ({qsar['n_test']} مركباً محتجزاً — {qsar['split_method']}):**"
                    )
                    st.caption(
                        f"اختبار محجوز ({qsar['n_test']} مركباً) — Ridge QSAR: R² = {qsar['test_r2']} | RMSE = {qsar['test_rmse']} | MAE = {qsar['test_mae']} "
                        f"|| Wildman-Crippen على نفس العينة: R² = {qsar['crippen_test_r2']} | RMSE = {qsar['crippen_test_rmse']} | MAE = {qsar['crippen_test_mae']}"
                    )
                    st.caption(
                        f"تحقق متقاطع خماسي عشوائي: R² = {qsar['cv_r2']} | RMSE = {qsar['cv_rmse']} "
                        f"|| تحقق بحجب فئة كيميائية كاملة: R² = {qsar['group_cv_r2']} | RMSE = {qsar['group_cv_rmse']} | MAE = {qsar['group_cv_mae']} "
                        f"|| Crippen على كل المركبات: R² = {qsar['crippen_all_r2']} | RMSE = {qsar['crippen_all_rmse']}"
                    )
                    st.caption(
                        "التقسيم العشوائي متفائل لأن المجموعة تضم سلاسل متجانسة قريبة بنيوياً؛ "
                        "التحقق بحجب الفئة أصدق لأنه يختبر النموذج على فئة لم يرها."
                    )
                    with st.expander("الأداء حسب الفئة الكيميائية (Leave-One-Class-Out)", expanded=False):
                        st.dataframe(pd.DataFrame(qsar["group_cv_by_class"]), use_container_width=True, hide_index=True)
                    st.dataframe(pd.DataFrame(qsar["test_set_table"]), use_container_width=True, hide_index=True)
                else:
                    st.warning("نموذج QSAR الإحصائي غير متاح في البيئة الحالية.")

                st.markdown(
                    f"""
<div class="ml-card">
  <h5>المنهجية والمراجع والحدود العلمية (Methodology, Citations &amp; Limitations)</h5>
  <p><strong>المنهجية:</strong> {_esc(eval_res['methodology_ar'])}</p>
  <p><strong>مرجع Wildman-Crippen:</strong> <code>{_esc(det['citation'])}</code></p>
  <p><strong>مصادر البيانات التجريبية للاختبار:</strong> <code>{_esc(MolecularIntelligence.DATASET_CITATION)}</code></p>
  <p><strong>حدود الاستخدام العلمية:</strong> {_esc(eval_res['limitations_ar'])}</p>
</div>
""",
                    unsafe_allow_html=True,
                )

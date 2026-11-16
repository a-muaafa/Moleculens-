"""
Expanded Drug-Drug Interaction (DDI) frontend view for MolecuLens.
Supports 10 therapeutic groups, 45+ medications, brand/generic/Arabic name resolution,
pairwise and multi-drug polypharmacy analysis with separate pair-by-pair results,
update timestamps, source license documentation, and live RxNorm/openFDA synchronization.
"""

from __future__ import annotations

import html
import pandas as pd
import streamlit as st
from pharmacy.ddi_engine import DDIEngine


def _format_drug_option(key: str, is_ar: bool = True) -> str:
    meta = DDIEngine.DRUG_CATALOG.get(key)
    if not meta:
        return key.capitalize()
    brands = ", ".join(meta.get("brand_names", [])[:2])
    brand_str = f" [{brands}]" if brands else ""
    if is_ar:
        return f"{meta['name_ar']} — {meta['name_en']}{brand_str}"
    return f"{meta['name_en']} ({meta['name_ar']}){brand_str}"


def render_ddi_view(lang: str = "ar") -> None:
    """Renders the expanded Drug-Drug Interaction screening module."""
    is_ar = lang == "ar"

    header_title = "فاحص التداخلات الدوائية الموسّع (Drug-Drug Interactions — DDI)" if is_ar else "Expanded Drug-Drug Interaction (DDI) Checker"
    header_desc = (
        "تحليل دوائي شامل للتداخلات بين الأدوية عبر 10 مجموعات علاجية (أكثر من 45 دواءً و180 اسماً تجارياً وعلمياً)، "
        "مع آلية التحديث من مصادر RxNorm و openFDA و PubChem، وتاريخ تحديث كل سجل."
        if is_ar
        else "Comprehensive pharmacological screening across 10 therapeutic classes with brand/generic resolution, RxNorm/openFDA synchronization, and update provenance."
    )

    st.markdown(
        f"""
        <div style="background: rgba(30, 41, 59, 0.7); border: 1px solid rgba(255, 255, 255, 0.1); border-radius: 12px; padding: 20px; margin-bottom: 20px;">
            <h2 style="color: #60a5fa; margin: 0 0 8px 0; font-size: 1.5rem;">{html.escape(header_title)}</h2>
            <p style="color: #94a3b8; margin: 0; font-size: 0.95rem;">{html.escape(header_desc)}</p>
        </div>
        """,
        unsafe_allow_html=True,
    )

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

    all_drugs = DDIEngine.get_supported_drugs()

    tab_pair, tab_multi, tab_catalog = st.tabs([
        "1. فحص ثنائي (بالاسم العلمي أو التجاري)" if is_ar else "1. Pairwise Drug Check",
        "2. فحص الوصفات المتعددة (دواءان أو أكثر — تحليل كل زوج منفصلاً)" if is_ar else "2. Multi-Drug Polypharmacy Analysis",
        "3. قاعدة بيانات الأدوية والمجموعات العلاجية والتحديث الدوري" if is_ar else "3. Drug Catalog & External Sync",
    ])

    with tab_pair:
        f_col1, f_col2 = st.columns(2)
        class_options = [("all", "جميع المجموعات العلاجية العشر (All 10 Therapeutic Classes)")] + [
            (k, v["name_ar"] if is_ar else v["name_en"])
            for k, v in DDIEngine.THERAPEUTIC_CLASSES.items()
        ]
        with f_col1:
            cls_a = st.selectbox(
                "تصفية المجموعة العلاجية للدواء الأول:",
                options=[c[0] for c in class_options],
                format_func=lambda k: dict(class_options).get(k, k),
                key="ddi_cls_a",
            )
            drugs_a = DDIEngine.get_supported_drugs(cls_a) or all_drugs
            drug_a = st.selectbox(
                "الدواء الأول (الاسم العلمي / التجاري — Drug A):" if is_ar else "Drug A",
                options=drugs_a,
                format_func=lambda d: _format_drug_option(d, is_ar),
                index=0,
                key="ddi_drug_a",
            )
            meta_a = DDIEngine.get_drug_metadata(drug_a) or {}
            if meta_a:
                st.caption(
                    f"المجموعة: {meta_a.get('therapeutic_class_ar', '')} | "
                    f"الأسماء التجارية: {', '.join(meta_a.get('brand_names', []))} | "
                    f"RxCUI: `{meta_a.get('rxcui', '')}` | آخر تحديث: `{meta_a.get('last_updated', '')}`"
                )

        with f_col2:
            cls_b = st.selectbox(
                "تصفية المجموعة العلاجية للدواء الثاني:",
                options=[c[0] for c in class_options],
                format_func=lambda k: dict(class_options).get(k, k),
                key="ddi_cls_b",
            )
            drugs_b = DDIEngine.get_supported_drugs(cls_b) or all_drugs
            default_b_idx = 1 if len(drugs_b) > 1 else 0
            drug_b = st.selectbox(
                "الدواء الثاني (الاسم العلمي / التجاري — Drug B):" if is_ar else "Drug B",
                options=drugs_b,
                format_func=lambda d: _format_drug_option(d, is_ar),
                index=default_b_idx,
                key="ddi_drug_b",
            )
            meta_b = DDIEngine.get_drug_metadata(drug_b) or {}
            if meta_b:
                st.caption(
                    f"المجموعة: {meta_b.get('therapeutic_class_ar', '')} | "
                    f"الأسماء التجارية: {', '.join(meta_b.get('brand_names', []))} | "
                    f"RxCUI: `{meta_b.get('rxcui', '')}` | آخر تحديث: `{meta_b.get('last_updated', '')}`"
                )

        if st.button("تحليل التداخل الدوائي" if is_ar else "Analyze Interaction", type="primary"):
            res = DDIEngine.check_interaction(drug_a, drug_b)
            _render_single_pair_result(res, is_ar)

    with tab_multi:
        st.write(
            "أضف دواءين أو أكثر (بالاسم العلمي أو التجاري) لتحليل جميع أزواج الأدوية المدخلة وعرض نتيجة كل زوج بشكل منفصل:"
            if is_ar
            else "Select two or more medications to analyze all pairwise combinations separately:"
        )
        selected_multi = st.multiselect(
            "قائمة الأدوية الموصوفة (اختر دواءين أو أكثر)" if is_ar else "Prescribed Drugs (Select 2 or more)",
            options=all_drugs,
            default=["warfarin", "aspirin", "ibuprofen", "omeprazole"],
            format_func=lambda d: _format_drug_option(d, is_ar),
        )

        if st.button("فحص الوصفة الشاملة وتحليل جميع الأزواج" if is_ar else "Screen Polypharmacy Pairs", type="primary"):
            if len(selected_multi) < 2:
                st.warning("يرجى اختيار دواءين على الأقل لإجراء تحليل التداخلات الزوجية.")
            else:
                multi_res = DDIEngine.screen_polypharmacy(selected_multi)

                c1, c2, c3, c4, c5 = st.columns(5)
                c1.metric("الأدوية المفحوصة", multi_res["total_drugs"])
                c2.metric("إجمالي الأزواج المفحوصة", multi_res["total_pairs_evaluated"])
                c3.metric("خطورة عالية (Major)", multi_res["severity_summary"]["Major"])
                c4.metric("خطورة متوسطة (Moderate)", multi_res["severity_summary"]["Moderate"])
                c5.metric("خطورة طفيفة (Minor)", multi_res["severity_summary"]["Minor"])

                st.markdown("#### تفصيل نتائج جميع أزواج الأدوية المدخلة (Pair-by-Pair Analysis)")
                for pair_res in multi_res.get("all_pairs_results", []):
                    _render_single_pair_result(pair_res, is_ar)

    with tab_catalog:
        st.subheader("موسوعة الأدوية المعتمدة (10 مجموعات علاجية) والتحديث الدوري من المصادر الموثوقة")
        cat_rows = [
            {
                "المادة الفعالة (Generic)": v["name_en"],
                "الاسم العربي": v["name_ar"],
                "الأسماء التجارية الشائعة (Brands)": "، ".join(v.get("brand_names", [])),
                "المجموعة العلاجية": v["therapeutic_class_ar"],
                "معرف RxNorm (RxCUI)": v.get("rxcui", ""),
                "معرف DrugBank": v.get("drugbank_id", ""),
                "معرف PubChem CID": str(v.get("pubchem_cid", "")),
                "تاريخ التحديث": v.get("last_updated", ""),
            }
            for v in DDIEngine.DRUG_CATALOG.values()
        ]
        st.dataframe(pd.DataFrame(cat_rows), use_container_width=True, hide_index=True)

        st.markdown("#### التحديث المباشر من المصادر الدوائية المعتمدة (RxNorm + openFDA + PubChem)")
        sync_col1, sync_col2 = st.columns([2, 1])
        with sync_col1:
            sync_drug_key = st.selectbox(
                "اختر الدواء لتحديث بياناته ونشرته الرسمية من RxNorm و openFDA:",
                options=all_drugs,
                format_func=lambda d: _format_drug_option(d, is_ar),
                key="ddi_sync_select",
            )
        with sync_col2:
            st.write("")
            if st.button("تحديث بيانات الدواء من المصادر المعتمدة", key="ddi_sync_btn"):
                sync_out = DDIEngine.sync_drug_from_external_sources(sync_drug_key)
                if sync_out.get("ok"):
                    rec = sync_out["record"]
                    st.success(
                        f"تم التحقق وتوحيد السجل للدواء ({rec['name_ar']} — RxCUI: {rec.get('rxcui', 'N/A')}) "
                        f"دون تكرار | تاريخ التحديث: {rec.get('last_updated')}"
                    )
                    if rec.get("openfda_interactions_excerpt"):
                        st.info(f"**مقتطف التداخلات من نشرة openFDA الرسمية:** {rec['openfda_interactions_excerpt']}")

        with st.expander("توثيق المصادر الدوائية وتراخيص استخدامها (Data Sources & Licensing Provenance)", expanded=False):
            lic_rows = [
                {
                    "المصدر العلمي": k,
                    "الجهة الرسمية": v["name"],
                    "نطاق الاستخدام في المنصة": v["scope_ar"],
                    "الترخيص القانوني (License)": v["license"],
                }
                for k, v in DDIEngine.DATA_SOURCE_LICENSES.items()
            ]
            st.dataframe(pd.DataFrame(lic_rows), use_container_width=True, hide_index=True)


def _render_single_pair_result(res: dict, is_ar: bool = True) -> None:
    meta_a = res.get("drug_a_metadata") or {}
    meta_b = res.get("drug_b_metadata") or {}
    label_a = f"{meta_a.get('name_ar', res.get('drug_a', ''))} ({meta_a.get('name_en', res.get('drug_a', ''))})" if meta_a else str(res.get("drug_a", "")).capitalize()
    label_b = f"{meta_b.get('name_ar', res.get('drug_b', ''))} ({meta_b.get('name_en', res.get('drug_b', ''))})" if meta_b else str(res.get("drug_b", "")).capitalize()

    if res.get("has_interaction"):
        severity = res["severity"]
        sev_color = {
            "Major": "#ef4444",
            "Moderate": "#f59e0b",
            "Minor": "#3b82f6",
        }.get(severity, "#94a3b8")

        st.markdown(
            f"""
            <div style="background: rgba(30, 41, 59, 0.85); border: 2px solid {sev_color}; border-radius: 10px; padding: 18px; margin-top: 12px;">
                <div style="display: flex; justify-content: space-between; align-items: center; flex-wrap: wrap; gap: 8px; margin-bottom: 12px;">
                    <h4 style="color: #f1f5f9; margin: 0;">
                        {html.escape(label_a)} + {html.escape(label_b)}
                    </h4>
                    <span style="background: {sev_color}; color: white; padding: 5px 14px; border-radius: 20px; font-weight: bold; font-size: 0.88rem;">
                        {html.escape(str(res['severity']))} — {html.escape(str(res['severity_ar']))}
                    </span>
                </div>
                <div style="margin-bottom: 10px;">
                    <strong style="color: #60a5fa;">{("الآلية الدوائية (Mechanism):" if is_ar else "Mechanism:")}</strong>
                    <p style="color: #cbd5e1; margin: 3px 0 0 0;">{html.escape(res['mechanism_ar'] if is_ar else res['mechanism_en'])}</p>
                </div>
                <div style="margin-bottom: 10px;">
                    <strong style="color: #60a5fa;">{("الأثر السريري (Clinical Effect):" if is_ar else "Clinical Effect:")}</strong>
                    <p style="color: #cbd5e1; margin: 3px 0 0 0;">{html.escape(res['clinical_effect_ar'] if is_ar else res['clinical_effect_en'])}</p>
                </div>
                <div style="margin-bottom: 10px;">
                    <strong style="color: #10b981;">{("التوصيات والتدبير السريري (Management):" if is_ar else "Management:")}</strong>
                    <p style="color: #cbd5e1; margin: 3px 0 0 0;">{html.escape(res['management_ar'] if is_ar else res['management_en'])}</p>
                </div>
                <div style="display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;color:#94a3b8;font-size:0.82rem;border-top:1px solid rgba(255,255,255,0.1);padding-top:8px;">
                    <div><strong>المصدر العلمي:</strong> {' | '.join(html.escape(c) for c in res.get('citations', []))}</div>
                    <div><strong>تاريخ تحديث البيانات:</strong> <code>{html.escape(str(res.get('last_updated', '2026-11-15')))}</code></div>
                </div>
            </div>
            """,
            unsafe_allow_html=True,
        )
    elif res.get("is_same_drug"):
        st.warning(
            f"تحذير تكرار المادة الفعالة: ({res.get('drug_a')}) و ({res.get('drug_b')}) يشيران لنفس الدواء ({res.get('canonical_a')}). "
            f"{res.get('mechanism_ar', '')}"
        )
    else:
        st.markdown(
            f"""
            <div style="background: rgba(30, 41, 59, 0.55); border: 1px solid #64748b; border-radius: 10px; padding: 14px 18px; margin-top: 12px;">
                <div style="display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;">
                    <strong style="color: #e2e8f0;">
                        {html.escape(label_a)} + {html.escape(label_b)}
                    </strong>
                    <span style="background:#475569;color:#f8fafc;padding:3px 10px;border-radius:14px;font-size:0.8rem;">
                        غير مسجل في السجل المحلي (تاريخ التحديث: {html.escape(str(res.get('last_updated', '2026-11-15')))})
                    </span>
                </div>
                <p style="color: #cbd5e1; margin: 6px 0 0 0; font-size: 0.88rem;">
                    {html.escape(res.get('message_ar' if is_ar else 'message_en', ''))}
                </p>
            </div>
            """,
            unsafe_allow_html=True,
        )

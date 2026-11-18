"""
Pharmacology – علم الأدوية (Coming Soon – قيد التطوير) frontend view for MolecuLens.
Provides a structured, bilingual (Arabic/English) preview architecture of the upcoming
Pharmacology & Clinical Therapeutics module without enabling incomplete functionality.
"""

from __future__ import annotations

import html
import streamlit as st
from pharmacy.pharmacology_module import (
    PHARMACOLOGY_EDUCATIONAL_DISCLAIMER_AR,
    PHARMACOLOGY_PILLARS,
)


def render_pharmacology_view(lang: str = "ar") -> None:
    """Renders the Pharmacology (Coming Soon) module preview interface."""
    is_ar = lang == "ar"

    st.markdown(
        f"""
<div class="ml-card" style="border-right:5px solid #0284c7;background:linear-gradient(135deg,#0f4c5c 0%,#1e293b 100%);color:#f8fafc;">
  <div style="display:flex;justify-content:space-between;align-items:center;flex-wrap:wrap;gap:10px;">
    <div>
      <h2 style="color:#38bdf8;margin:0 0 6px 0;">
        {("Pharmacology – علم الأدوية" if is_ar else "Pharmacology & Clinical Therapeutics")}
      </h2>
      <p style="color:#e2e8f0;margin:0;font-size:0.96rem;">
        {("وحدة علم الأدوية الجزيئي والسريري لطلاب الصيدلة والعلوم الصحية — واجهة مبدئية مهيأة للربط مع النواة العلمية."
          if is_ar else
          "Molecular and Clinical Pharmacology module for pharmacy and health-science students — preview architecture.")}
      </p>
    </div>
    <span style="background:#f59e0b;color:#0f172a;padding:8px 18px;border-radius:20px;font-weight:800;font-size:0.95rem;">
      Coming Soon – قيد التطوير
    </span>
  </div>
</div>
""",
        unsafe_allow_html=True,
    )

    st.info(
        "هذا القسم في حالة (Coming Soon – قيد التطوير): تم تصميم الهيكل المعماري للأقسام الفرعية السبعة أدناه "
        "وتجهيز نقاط الربط البرمجية دون تفعيل وظائف غير مكتملة أو عرض بيانات وهمية، التزاماً بمعيار الصرامة العلمية في MolecuLens. "
        "يمكنك حالياً استكشاف الأدوية ومحددات ليبينسكي وبطاقات السلامة في قسم (Explore & Understand) وقسم (الذكاء الجزيئي — Molecular Intelligence)."
        if is_ar
        else "This module is currently in 'Coming Soon – Under Development' status. The 7 architectural pillars below are structured and ready for backend integration without exposing incomplete features or placeholder data."
    )

    st.warning(PHARMACOLOGY_EDUCATIONAL_DISCLAIMER_AR)

    st.subheader(
        "المحاور العلمية السبعة المخطط تفعيلها في قسم علم الأدوية (Planned Pharmacology Modules)"
        if is_ar
        else "Seven Planned Pharmacology & Therapeutics Pillars"
    )

    cols = st.columns(2)
    for idx, pillar in enumerate(PHARMACOLOGY_PILLARS):
        with cols[idx % 2]:
            title = pillar["title_ar"] if is_ar else pillar["title_en"]
            desc = pillar["desc_ar"] if is_ar else pillar["desc_en"]
            sub_title = pillar["title_en"] if is_ar else pillar["title_ar"]
            st.markdown(
                f"""
<div class="ml-card" style="border-top:4px solid #0284c7;min-height:190px;display:flex;flex-direction:column;justify-content:space-between;">
  <div>
    <div style="display:flex;justify-content:space-between;align-items:flex-start;gap:8px;">
      <h4 style="color:#0f4c5c;margin:0 0 4px 0;">{html.escape(title)}</h4>
      <span style="background:rgba(2,132,199,0.12);color:#0369a1;padding:2px 10px;border-radius:12px;font-size:0.78rem;font-weight:700;white-space:nowrap;">
        قيد التطوير (Coming Soon)
      </span>
    </div>
    <div style="color:#64748b;font-size:0.82rem;margin-bottom:8px;">{html.escape(sub_title)}</div>
    <p style="color:#334155;font-size:0.9rem;margin-bottom:8px;">{html.escape(desc)}</p>
  </div>
  <div style="font-size:0.8rem;color:#0f766e;border-top:1px solid #e2e8f0;padding-top:6px;">
    <strong>{("المصادر المعتمدة للربط:" if is_ar else "Planned Data Sources:")}</strong>
    <span class="ltr-inline">{html.escape(pillar["planned_sources"])}</span>
  </div>
</div>
""",
                unsafe_allow_html=True,
            )

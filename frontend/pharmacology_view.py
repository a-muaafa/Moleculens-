"""
Pharmacology – علم الأدوية (Coming Soon – قيد التطوير) frontend view for MolecuLens.
Provides a structured, bilingual (Arabic/English) preview architecture of the upcoming
Pharmacology & Clinical Therapeutics module without enabling incomplete functionality.
"""

from __future__ import annotations

import html
import streamlit as st

PHARMACOLOGY_PILLARS: list[dict[str, str]] = [
    {
        "code": "classification",
        "title_ar": "1. تصنيف الأدوية والمجموعات العلاجية",
        "title_en": "1. Drug Classification & Therapeutic Classes",
        "desc_ar": (
            "تصنيف منهجي للأدوية وفق نظام ATC العالمي والمجموعات العلاجية والدوائية "
            "(المسكنات، المضادات الحيوية، أدوية القلب والضغط، السكري، الجهاز الهضمي، والتنفسي)."
        ),
        "desc_en": (
            "Systematic organization of medications by WHO ATC classification and therapeutic families "
            "(analgesics, antibiotics, cardiovascular, antidiabetics, GI, and respiratory agents)."
        ),
        "planned_sources": "WHO ATC Index, RxNorm, ChEMBL",
    },
    {
        "code": "moa",
        "title_ar": "2. آليات عمل الأدوية (Mechanism of Action)",
        "title_en": "2. Mechanism of Action (MoA)",
        "desc_ar": (
            "تفسير جزيئي وخلوي لكيفية إحداث الدواء لتأثيره العلاجي: تثبيط الإنزيمات، "
            "تنشيط أو حصر المستقبلات، تعديل القنوات الأيونية، أو التداخل مع المسارات الأيضية."
        ),
        "desc_en": (
            "Molecular and cellular explanation of therapeutic action: enzyme inhibition, "
            "receptor agonism/antagonism, ion channel modulation, and metabolic pathway targeting."
        ),
        "planned_sources": "ChEMBL Mechanism API, IUPHAR/BPS Guide to Pharmacology",
    },
    {
        "code": "targets",
        "title_ar": "3. المستقبلات والأهداف الحيوية (Receptors & Biological Targets)",
        "title_en": "3. Receptors & Biological Targets",
        "desc_ar": (
            "ربط المادة الدوائية بالبروتينات والمستقبلات المستهدفة (GPCRs، القنوات الأيونية، "
            "الإنزيمات، والمستقبلات النووية) مع معرفات UniProt و PDB."
        ),
        "desc_en": (
            "Linking active pharmaceutical ingredients to biological macromolecular targets "
            "(GPCRs, ion channels, enzymes, nuclear receptors) with UniProt and PDB cross-references."
        ),
        "planned_sources": "UniProtKB, RCSB PDB, Open Targets",
    },
    {
        "code": "adme",
        "title_ar": "4. الحرائك الدوائية (Pharmacokinetics — ADME)",
        "title_en": "4. Pharmacokinetics (ADME Profile)",
        "desc_ar": (
            "دراسة رحلة الدواء داخل الجسم: الامتصاص (Absorption) والتوافر الحيوي، التوزيع (Distribution) "
            "والارتباط ببروتينات البلازما، الأيض الكبدي (Metabolism / CYP450)، والإطراح (Excretion) وعمر النصف."
        ),
        "desc_en": (
            "Quantitative profiling of Absorption (bioavailability), Distribution (volume of distribution, "
            "plasma protein binding), Metabolism (CYP450 isoenzymes), and Excretion (clearance, half-life)."
        ),
        "planned_sources": "FDA Clinical Pharmacology Monographs, ChEMBL ADMET",
    },
    {
        "code": "dose_response",
        "title_ar": "5. العلاقة بين الجرعة والاستجابة (Dose–Response Relationship)",
        "title_en": "5. Dose–Response Relationship & Therapeutic Index",
        "desc_ar": (
            "نمذجة منحنيات التركيز والاستجابة (EC50, IC50, Emax)، الفعالية القصوى (Efficacy) مقابل القوة الدوائية (Potency)، "
            "والنافذة العلاجية (Therapeutic Window) ومؤشر الأمان الدوائي."
        ),
        "desc_en": (
            "Modeling graded and quantal dose-response curves (EC50, IC50, Emax), potency vs. efficacy, "
            "and therapeutic index safety margins."
        ),
        "planned_sources": "Goodman & Gilman's Pharmacological Basis of Therapeutics",
    },
    {
        "code": "adverse_effects",
        "title_ar": "6. الآثار الجانبية وموانع الاستخدام (Adverse Effects & Contraindications)",
        "title_en": "6. Adverse Effects & Contraindications",
        "desc_ar": (
            "توثيق التفاعلات الضارة الشائعة والخطيرة، التحذيرات المؤطرة (Boxed Warnings)، "
            "موانع الاستخدام، والاحتياطات في القصور الكلوي أو الكبدي والحمل والرضاعة."
        ),
        "desc_en": (
            "Structured catalog of adverse drug reactions, boxed warnings, contraindications, "
            "and special population precautions (renal/hepatic impairment, pregnancy)."
        ),
        "planned_sources": "openFDA Drug Labeling API, DailyMed",
    },
    {
        "code": "references",
        "title_ar": "7. المعلومات الدوائية والمراجع العلمية (Drug Monographs & References)",
        "title_en": "7. Drug Monographs & Scientific References",
        "desc_ar": (
            "ربط كل بطاقة دوائية بمراجعها الدستورية والسريرية المعتمدة (FDA، EMA، BNF، USP-NF) "
            "مع توثيق تاريخ آخر تحديث وتراخيص البيانات."
        ),
        "desc_en": (
            "Direct linkage of each drug monograph to authoritative regulatory and pharmacopeial "
            "references (FDA, EMA, BNF, USP-NF) with update timestamps and license provenance."
        ),
        "planned_sources": "openFDA (CC0), RxNorm (NLM), PubChem (Public Domain)",
    },
]


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
        "يمكنك حالياً استخدام قسم (التداخلات الدوائية — Drug Interactions) وقسم (الذكاء الجزيئي — Molecular Intelligence) المتاحين بالكامل."
        if is_ar
        else "This module is currently in 'Coming Soon – Under Development' status. The 7 architectural pillars below are structured and ready for backend integration without exposing incomplete features or placeholder data."
    )

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

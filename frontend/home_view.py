"""
Home overview and hero header components for MolecuLens.
Features modular cards for Explore & Understand, Reaction Lab, Titration Lab,
Concentration Lab, and Drug-Drug Interactions with bilingual support.
"""

from typing import Callable, Any
import streamlit as st
import html
from config.settings import APP_NAME, APP_SLOGAN, APP_SUBTITLE


def render_hero_header(lang: str = "ar") -> None:
    """Renders the top branding hero banner."""
    is_ar = lang == "ar"
    title = APP_NAME
    subtitle = (
        "المنصة العلمية الذكية لاستكشاف المواد والتفاعلات الكيميائية والحسابات الدوائية"
        if is_ar
        else "Intelligent Scientific Platform for Molecular Exploration, Reactions & Clinical Calculations"
    )
    desc = (
        "منصة علمية تفاعلية للتعلم والاستكشاف الكيميائي والحيوي والدوائي وتهيئة طلاب الصيدلة والعلوم قبل دخول المختبر (Pre-Lab Simulation). "
        "ابحث عن أي عنصر أو مركب أو دواء أو بروتين، أو ادخل إلى مختبرات المعايرة والتفاعلات والتراكيز والتداخلات الدوائية، "
        "وستتولى النواة العلمية للمنصة جلب البيانات الموثقة وإجراء الحسابات بدقة وموثوقية عالية."
        if is_ar
        else "Interactive scientific learning platform for pharmacy and chemistry students before entering the laboratory (Pre-Lab Simulation). "
        "Explore elements, compounds, drugs, and proteins, or simulate titrations, chemical reactions, concentration recipes, and drug-drug interactions with scientific rigor."
    )

    st.markdown(
        f"""
<div class="ml-hero">
  <div class="ml-hero-top">
    <div>
      <h1 class="ml-brand-title">{title}</h1>
      <div class="ml-brand-sub">{html.escape(subtitle)}</div>
    </div>
    <span class="ml-slogan-pill">{APP_SLOGAN}</span>
  </div>
  <p class="ml-hero-desc">
    {html.escape(desc)}
  </p>
</div>
""",
        unsafe_allow_html=True,
    )


def render_home_overview(on_navigate: Callable[[str], None], lang: str = "ar") -> None:
    """Renders the grid of 5 core modules on the home overview."""
    is_ar = lang == "ar"

    c1, c2 = st.columns(2)
    with c1:
        st.markdown(
            f"""
<div class="ml-section-card">
  <h3>1. Explore &amp; Understand</h3>
  <p><strong>{('استكشاف المادة وفهم تركيبها وخصائصها ومساعد الذكاء الاصطناعي' if is_ar else 'Molecular Exploration, Properties & AI Assistant')}</strong></p>
  <p>
    {('بحث موحد ذكي عن العناصر والمركبات والأدوية والبروتينات مع تحليل البنية الجزيئية ومحددات ليبينسكي والتفسير البيداغوجي خطوة بخطوة.' if is_ar
      else 'Unified intelligent search for elements, compounds, drugs, and proteins with 3D structure visualization, Lipinski druglikeness, and AI explanations.')}
  </p>
</div>
""",
            unsafe_allow_html=True,
        )
        if st.button("الدخول إلى Explore & Understand" if is_ar else "Open Explore & Understand", key="home_go_explore"):
            on_navigate("explore")

    with c2:
        st.markdown(
            f"""
<div class="ml-section-card">
  <h3>2. Reaction Lab</h3>
  <p><strong>{('مختبر تحليل التفاعلات والموازنة والحسابات الستوكيومترية' if is_ar else 'Reaction Analysis & Stoichiometric Yields')}</strong></p>
  <p>
    {('موازنة المعادلات الكيميائية جبرياً، وتحديد المتفاعل المحدد (Limiting Reagent)، وحساب المردود النظري مع التحقق من شروط التفاعل وقواعد بيانات ORD.' if is_ar
      else 'Algebraic equation balancing, limiting reagent identification, theoretical yields calculation, and Open Reaction Database (ORD) integration.')}
  </p>
</div>
""",
            unsafe_allow_html=True,
        )
        if st.button("الدخول إلى Reaction Lab" if is_ar else "Open Reaction Lab", key="home_go_reaction"):
            on_navigate("reaction_lab")

    st.markdown("<br/>", unsafe_allow_html=True)
    c3, c4, c5 = st.columns(3)

    with c3:
        st.markdown(
            f"""
<div class="ml-section-card">
  <h3>3. Titration Lab</h3>
  <p><strong>{('مختبر المعايرة الافتراضي' if is_ar else 'Virtual Titration Lab')}</strong></p>
  <p>
    {('محاكاة ديناميكية لمنحنيات pH لتفاعلات حمض-قاعدة قوية وضعيفة، ونقاط التكافؤ ونصف التكافؤ واختيار الكواشف اللونية.' if is_ar
      else 'Thermodynamic pH curve simulation for strong/weak acid-base pairs, buffer capacity, equivalence detection, and indicator matching.')}
  </p>
</div>
""",
            unsafe_allow_html=True,
        )
        if st.button("الدخول إلى Titration Lab" if is_ar else "Open Titration Lab", key="home_go_titration"):
            on_navigate("titration")

    with c4:
        st.markdown(
            f"""
<div class="ml-section-card">
  <h3>4. Concentration Lab</h3>
  <p><strong>{('حاسبة التراكيز وإعداد المحاليل' if is_ar else 'Concentration Calculator')}</strong></p>
  <p>
    {('حسابات المولارية، العيارية، النسب المئوية، أجزاء المليون، وقانون التخفيف C1V1=C2V2 مع بروتوكولات التحضير المعملية خطوة بخطوة.' if is_ar
      else 'Molarity, Normality, %, ppm/ppb, and Dilution calculations with step-by-step laboratory preparation protocols.')}
  </p>
</div>
""",
            unsafe_allow_html=True,
        )
        if st.button("الدخول إلى Concentration Lab" if is_ar else "Open Concentration Lab", key="home_go_concentration"):
            on_navigate("concentration")

    with c5:
        st.markdown(
            f"""
<div class="ml-section-card">
  <h3>5. Drug Interactions</h3>
  <p><strong>{('فاحص التداخلات الدوائية (DDI)' if is_ar else 'Drug-Drug Interactions')}</strong></p>
  <p>
    {('فحص التداخلات الدوائية الثنائية والوصفات المتعددة، وتحديد الآليات الحركية والديناميكية ودرجة الخطورة مع التنبيه السريري الإلزامي.' if is_ar
      else 'Pairwise and polypharmacy interaction screening, mechanism analysis, severity categorization, citations, and safety disclaimers.')}
  </p>
</div>
""",
            unsafe_allow_html=True,
        )
        if st.button("الدخول إلى Drug Interactions" if is_ar else "Open Drug Interactions", key="home_go_ddi"):
            on_navigate("ddi")

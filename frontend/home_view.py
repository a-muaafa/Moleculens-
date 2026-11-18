"""
Home overview and hero header components for MolecuLens.
Features modular cards for Explore & Understand, Reaction Lab, Titration Lab,
and Concentration Lab with bilingual support.
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
        "ابحث عن أي عنصر أو مركب أو دواء أو بروتين، أو ادخل إلى مختبر التفاعلات الكيميائية (Reaction Lab) ومختبرات المعايرة والتراكيز والذكاء الجزيئي، "
        "وستتولى النواة العلمية للمنصة جلب البيانات الموثقة وإجراء الحسابات بدقة وموثوقية عالية."
        if is_ar
        else "Interactive scientific learning platform for pharmacy and chemistry students before entering the laboratory (Pre-Lab Simulation). "
        "Explore elements, compounds, drugs, and proteins, or simulate titrations, chemical reactions and concentration recipes with scientific rigor."
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
  <h3>2. Reaction Lab – مختبر التفاعلات الكيميائية</h3>
  <p><strong>{('المختبر الموحد للتوقع التفاعلي، الموازنة، والحسابات الستوكيومترية (152 تفاعلاً موثقاً)' if is_ar else 'Unified Reaction Prediction, Stoichiometry & Visual Lab (152 Verified Reactions)')}</strong></p>
  <p>
    {('أدخل مادتين أو أكثر (عنصر + عنصر، عنصر + مركب، مركب + مركب، أو عدة مواد معاً) لتوقع النواتج بمستويات الدليل الثلاثة، '
      'وحساب المتفاعل المحدد والفائض والمردود النظري، والتمثيل البصري للجزيئات قبل التفاعل وبعده، وتجربة What-If دون اختلاق نواتج.' if is_ar
      else 'Enter two or more substances for 3-tier reaction prediction, limiting/excess stoichiometry, 2D molecular visualization before/after reaction, and What-If analysis.')}
  </p>
</div>
""",
            unsafe_allow_html=True,
        )
        if st.button("الدخول إلى Reaction Lab – مختبر التفاعلات الكيميائية" if is_ar else "Open Reaction Lab", key="home_go_reaction"):
            on_navigate("reaction_lab")

    st.markdown("<br/>", unsafe_allow_html=True)
    c3, c4 = st.columns(2)

    with c3:
        st.markdown(
            f"""
<div class="ml-section-card">
  <h3>3. Molecular Intelligence</h3>
  <p><strong>{('الذكاء الجزيئي والتعديل البنيوي ومقارنة التشابه' if is_ar else 'Molecular Intelligence & Structure Modification')}</strong></p>
  <p>
    {('حساب بصمات مورغان (Morgan Fingerprints) ومعامل تانيموتو، التعديل البنيوي للمجموعات الوظيفية مع دراسة الفروق الدوائية، وتقييم نموذج QSAR لتنبؤ LogP على 120 مركباً مرجعياً (90 تدريب / 30 اختبار + تحقق متقاطع 5-Fold CV).' if is_ar
      else 'Morgan Fingerprints and Tanimoto similarity, functional group structure editing with delta ADMET metrics, and QSAR LogP evaluation on 120 experimental reference compounds (90 train / 30 test + 5-Fold CV).')}
  </p>
</div>
""",
            unsafe_allow_html=True,
        )
        if st.button("الدخول إلى Molecular Intelligence" if is_ar else "Open Molecular Intelligence", key="home_go_molecular_intelligence"):
            on_navigate("molecular_intelligence")

    with c4:
        st.markdown(
            f"""
<div class="ml-section-card">
  <h3>4. Titration Lab</h3>
  <p><strong>{('مختبر المعايرة الافتراضي الموسّع (41 حمضاً و30 قاعدة و23 كاشفاً)' if is_ar else 'Expanded Virtual Titration Lab')}</strong></p>
  <p>
    {('محاكاة ديناميكية لمنحنيات pH للأحماض والقواعد القوية والضعيفة ومتعددة البروتونات، مع اختيار الكاشف اللوني تلقائياً ومحاكاة تغير اللون داخل الدورق ونقطة التكافؤ ونقطة النهاية.' if is_ar
      else 'Dynamic pH curve simulation for strong, weak, and polyprotic acid-base pairs, automatic indicator selection, and flask color simulation.')}
  </p>
</div>
""",
            unsafe_allow_html=True,
        )
        if st.button("الدخول إلى Titration Lab" if is_ar else "Open Titration Lab", key="home_go_titration"):
            on_navigate("titration")

    st.markdown("<br/>", unsafe_allow_html=True)
    c5, c6 = st.columns(2)

    with c5:
        st.markdown(
            f"""
<div class="ml-section-card">
  <h3>5. Concentration Lab</h3>
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

    with c6:
        st.markdown(
            f"""
<div class="ml-section-card">
  <h3>6. Pharmacology – علم الأدوية</h3>
  <p><strong>{('علم الأدوية الجزيئي والسريري (Coming Soon – قيد التطوير)' if is_ar else 'Pharmacology & Clinical Therapeutics (Coming Soon)')}</strong></p>
  <p>
    {('واجهة مبدئية مهيكلة لـ 7 محاور: تصنيف الأدوية، آليات العمل (MoA)، المستقبلات، الحرائك الدوائية ADME، الجرعة والاستجابة، الآثار الجانبية، والمراجع العلمية.' if is_ar
      else 'Structured preview of 7 pharmacology pillars: Drug Classification, MoA, Receptors, ADME, Dose-Response, Adverse Effects, and Monographs.')}
  </p>
</div>
""",
            unsafe_allow_html=True,
        )
        if st.button("معاينة Pharmacology (Coming Soon)" if is_ar else "Preview Pharmacology", key="home_go_pharmacology"):
            on_navigate("pharmacology")

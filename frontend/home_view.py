"""
Home overview and hero header components for MolecuLens.
Features a concise 3-step welcome guide and modular cards for Explore & Understand,
Reaction Lab, Molecular Intelligence, Titration Lab, and Concentration Lab.
"""

from typing import Callable
import streamlit as st
import html
from config.settings import (
    APP_NAME,
    APP_SLOGAN,
    SCIENTIFIC_DISCLAIMER_AR,
    SHOW_EXPERIMENTAL_PHARMACOLOGY,
)


def render_hero_header(lang: str = "ar") -> None:
    """Renders the top branding hero banner."""
    is_ar = lang == "ar"
    title = APP_NAME
    subtitle = (
        "المنصة العلمية الذكية لاستكشاف المواد والتفاعلات الكيميائية والتحليل الجزيئي"
        if is_ar
        else "Intelligent Scientific Platform for Molecular Exploration, Reactions & Titration"
    )
    desc = (
        "منصة علمية تفاعلية للتعلم والاستكشاف الكيميائي وتهيئة طلاب الصيدلة والعلوم قبل دخول المختبر (Pre-Lab Simulation). "
        "ابحث عن أي عنصر أو مركب أو دواء أو بروتين، أو ادخل إلى مختبر التفاعلات الكيميائية (Reaction Lab) ومختبرات المعايرة والذكاء الجزيئي والتراكيز، "
        "مع عمل كامل للوظائف الأساسية دون إنترنت ومصادر علمية موثقة."
        if is_ar
        else "Interactive scientific learning platform for pharmacy and chemistry students before entering the laboratory (Pre-Lab Simulation). "
        "Explore elements, compounds, drugs, and proteins, or simulate chemical reactions, titrations, and molecular analysis offline-first."
    )
    chips_html = (
        """
  <div class="ml-hero-stats">
    <span class="ml-hero-chip">241 مادة علمية موثقة محلياً</span>
    <span class="ml-hero-chip">152 تفاعلاً معتمداً بمرجع علمي</span>
    <span class="ml-hero-chip">41 حمضاً و 30 قاعدة و 23 كاشفاً</span>
    <span class="ml-hero-chip">توثيق علمي صريح دون اختلاق</span>
  </div>
"""
        if is_ar
        else """
  <div class="ml-hero-stats">
    <span class="ml-hero-chip">241 Offline Substances</span>
    <span class="ml-hero-chip">152 Verified Reactions</span>
    <span class="ml-hero-chip">41 Acids / 30 Bases / 23 Indicators</span>
    <span class="ml-hero-chip">Evidence-Grounded Core</span>
  </div>
"""
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
  {chips_html}
</div>
""",
        unsafe_allow_html=True,
    )


def render_home_overview(on_navigate: Callable[[str], None], lang: str = "ar") -> None:
    """Renders the 3-step welcome card and the grid of completed core modules."""
    is_ar = lang == "ar"

    if is_ar:
        st.markdown(
            """
<div class="ml-card" style="border-right:4px solid var(--ml-accent);margin-bottom:1.25rem;">
  <h3 style="margin-top:0;">ابدأ رحلتك في MolecuLens بثلاث خطوات</h3>
  <div class="ml-step-grid">
    <div class="ml-step-tile">
      <div style="display:flex;align-items:center;margin-bottom:0.35rem;">
        <span class="ml-step-num">01</span>
        <strong>اختر مركباً أو تجربة جاهزة</strong>
      </div>
      <p style="margin:0;font-size:0.9rem;color:var(--ml-text-muted);">ابحث عن أي مادة بالعربية أو الإنجليزية أو الصيغة الكيميائية، أو اختر تجربة جاهزة تعمل دون إنترنت.</p>
    </div>
    <div class="ml-step-tile">
      <div style="display:flex;align-items:center;margin-bottom:0.35rem;">
        <span class="ml-step-num">02</span>
        <strong>شغّل التحليل أو المحاكاة</strong>
      </div>
      <p style="margin:0;font-size:0.9rem;color:var(--ml-text-muted);">نفّذ محاكاة التفاعل الكيميائي أو المعايرة الحجمية أو تحليل الخصائص الجزيئية وQSAR.</p>
    </div>
    <div class="ml-step-tile">
      <div style="display:flex;align-items:center;margin-bottom:0.35rem;">
        <span class="ml-step-num">03</span>
        <strong>استعرض النتائج والتفسير العلمي</strong>
      </div>
      <p style="margin:0;font-size:0.9rem;color:var(--ml-text-muted);">اطّلع على المعادلات الموزونة والرسوم التفاعلية ومصدر البيانات وحدود النموذج العلمي.</p>
    </div>
  </div>
</div>
""",
            unsafe_allow_html=True,
        )
    else:
        st.markdown(
            """
<div class="ml-card" style="border-left:4px solid var(--ml-accent);margin-bottom:1.25rem;">
  <h3 style="margin-top:0;">Start with MolecuLens in 3 Steps</h3>
  <ol style="margin:0.5rem 0 0 1.2rem;">
    <li><strong>Choose a compound or offline preset.</strong></li>
    <li><strong>Run the analysis or simulation.</strong></li>
    <li><strong>Review the results and scientific explanation.</strong></li>
  </ol>
</div>
""",
            unsafe_allow_html=True,
        )

    c1, c2 = st.columns(2)
    with c1:
        st.markdown(
            f"""
<div class="ml-section-card">
  <span class="ml-module-pill">{('المختبر 01 • التفاعلات والستوكيومتري' if is_ar else 'Lab 01 • Reactions')}</span>
  <h3>1. Reaction Lab – مختبر التفاعلات الكيميائية</h3>
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

    with c2:
        st.markdown(
            f"""
<div class="ml-section-card">
  <span class="ml-module-pill">{('المختبر 02 • المعايرة الحجمية' if is_ar else 'Lab 02 • Titration')}</span>
  <h3>2. Titration Lab – مختبر المعايرة الافتراضي</h3>
  <p><strong>{('محاكاة منحنيات pH والكواشف اللونية (41 حمضاً و30 قاعدة و23 كاشفاً)' if is_ar else 'Expanded Virtual Titration Lab (41 Acids, 30 Bases, 23 Indicators)')}</strong></p>
  <p>
    {('محاكاة ديناميكية لمنحنيات pH للأحماض والقواعد القوية والضعيفة ومتعددة البروتونات، مع اختيار الكاشف اللوني تلقائياً، '
      'وحساب سعة المحلول المنظم، والتمييز بين نقطة التكافؤ ونقطة النهاية، وتجارب جاهزة بنقرة واحدة.' if is_ar
      else 'Dynamic pH curve simulation for strong, weak, and polyprotic acid-base pairs, automatic indicator selection, and one-click presets.')}
  </p>
</div>
""",
            unsafe_allow_html=True,
        )
        if st.button("الدخول إلى Titration Lab" if is_ar else "Open Titration Lab", key="home_go_titration"):
            on_navigate("titration")

    st.markdown("<div style='height:0.65rem;'></div>", unsafe_allow_html=True)
    c3, c4 = st.columns(2)

    with c3:
        st.markdown(
            f"""
<div class="ml-section-card">
  <span class="ml-module-pill">{('المختبر 03 • الكيمياء المعلوماتية وQSAR' if is_ar else 'Lab 03 • Cheminformatics')}</span>
  <h3>3. Molecular Intelligence – الذكاء الجزيئي</h3>
  <p><strong>{('التشابه الجزيئي والتعديل البنيوي ونماذج QSAR الشفافة' if is_ar else 'Molecular Intelligence, Virtual Editing & Transparent QSAR')}</strong></p>
  <p>
    {('حساب بصمات مورغان (Morgan Fingerprints) ومعامل تانيموتو، التعديل البنيوي الافتراضي للمجموعات الوظيفية، '
      'وتقييم نماذج Ridge وHybrid Residual QSAR مقابل Wildman-Crippen عبر Scaffold Split وLOCO مع نطاق عدم اليقين.' if is_ar
      else 'Morgan Fingerprints and Tanimoto similarity, virtual functional group editing, and QSAR LogP evaluation with Scaffold Split, LOCO, and uncertainty.')}
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
  <span class="ml-module-pill">{('المختبر 04 • الاستكشاف ثنائي اللغة' if is_ar else 'Lab 04 • Exploration')}</span>
  <h3>4. Explore &amp; Understand – استكشاف المادة</h3>
  <p><strong>{('بحث موحد بالعربية والإنجليزية والصيغ الكيميائية (241 مادة محلية)' if is_ar else 'Unified Multilingual Chemical Search (241 Offline Substances)')}</strong></p>
  <p>
    {('استكشاف العناصر والمركبات والأدوية والبروتينات مع تحليل البنية الجزيئية ثنائية وثلاثية الأبعاد ومحددات ليبينسكي '
      'والتفسير العلمي المتدرج مع عمل كامل دون إنترنت.' if is_ar
      else 'Unified intelligent search for elements, compounds, drugs, and proteins with 2D/3D visualization, Lipinski druglikeness, and offline-first support.')}
  </p>
</div>
""",
            unsafe_allow_html=True,
        )
        if st.button("الدخول إلى Explore & Understand" if is_ar else "Open Explore & Understand", key="home_go_explore"):
            on_navigate("explore")

    st.markdown("<div style='height:0.65rem;'></div>", unsafe_allow_html=True)
    c5, c6 = st.columns(2)

    with c5:
        st.markdown(
            f"""
<div class="ml-section-card">
  <span class="ml-module-pill">{('المختبر 05 • المحاليل والتخفيف' if is_ar else 'Lab 05 • Solutions')}</span>
  <h3>5. Concentration Lab – حاسبة التراكيز</h3>
  <p><strong>{('حسابات التراكيز وإعداد المحاليل المخبرية' if is_ar else 'Concentration Calculator & Lab Recipes')}</strong></p>
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

    if SHOW_EXPERIMENTAL_PHARMACOLOGY:
        with c6:
            st.markdown(
                f"""
<div class="ml-section-card">
  <span class="ml-module-pill">{('خطة التطوير' if is_ar else 'Roadmap')}</span>
  <h3>6. Pharmacology – علم الأدوية</h3>
  <p><strong>{('علم الأدوية الجزيئي والسريري (قيد التطوير المستقبلي)' if is_ar else 'Pharmacology & Clinical Therapeutics (Experimental)')}</strong></p>
  <p>
    {('واجهة مبدئية مهيكلة لـ 7 محاور في علم الأدوية مخصصة للمراحل المستقبلية.' if is_ar
      else 'Structured preview of 7 pharmacology pillars for future roadmap phases.')}
  </p>
</div>
""",
                unsafe_allow_html=True,
            )
            if st.button("معاينة Pharmacology" if is_ar else "Preview Pharmacology", key="home_go_pharmacology"):
                on_navigate("pharmacology")

    st.info(SCIENTIFIC_DISCLAIMER_AR)

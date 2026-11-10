import streamlit as st
from config.settings import APP_NAME, APP_SLOGAN, APP_SUBTITLE


def render_hero_header() -> None:
    st.markdown(
        f"""
<div class="ml-hero">
  <div class="ml-hero-top">
    <div>
      <h1 class="ml-brand-title">{APP_NAME}</h1>
      <div class="ml-brand-sub">{APP_SUBTITLE}</div>
    </div>
    <span class="ml-slogan-pill">{APP_SLOGAN}</span>
  </div>
  <p class="ml-hero-desc">
    منصة علمية تفاعلية للتعلم والاستكشاف الكيميائي والحيوي والدوائي. ابحث عن أي عنصر أو مركب أو دواء أو بروتين،
    أو ادخل إلى مختبر التفاعلات لتحليل المعادلات والحسابات الكيميائية، وستتولى النواة العلمية للمنصة جلب البيانات
    وتصنيف الكيانات وإجراء الحسابات بدقة وموثوقية.
  </p>
</div>
""",
        unsafe_allow_html=True,
    )


def render_home_overview(on_navigate) -> None:
    col1, col2 = st.columns(2)
    with col1:
        st.markdown(
            """
<div class="ml-section-card">
  <h3>1. Explore &amp; Understand</h3>
  <p><strong>استكشاف المادة وفهم تركيبها وخصائصها</strong></p>
  <p>
    ابحث بالاسم العربي أو الإنجليزي أو الاسم العلمي (IUPAC) أو الصيغة الكيميائية أو الاختصار أو الاسم التجاري.
    يصنف النظام المادة تلقائيا (عنصر كيميائي، مركب كيميائي، مادة دوائية، بروتين/جزيء حيوي) ويعرض ملفها العلمي الكامل
    مع التفسير المبسط والجامعي المتقدم.
  </p>
</div>
""",
            unsafe_allow_html=True,
        )
        if st.button("الدخول إلى قسم Explore & Understand", key="home_go_explore"):
            on_navigate("explore")

    with col2:
        st.markdown(
            """
<div class="ml-section-card">
  <h3>2. Reaction Lab</h3>
  <p><strong>تحليل التفاعلات والحسابات والنواتج</strong></p>
  <p>
    اختر المواد المتفاعلة وحدد كمياتها ووحداتها وظروف التفاعل، ليقوم المحرك العلمي بالتحقق من المواد، وموازنة المعادلة
    الكيميائية، وحساب المتفاعل المحدد والفائض، وتحديد الناتج النظري مع جدول التغير المولي الكامل والشرح العلمي.
  </p>
</div>
""",
            unsafe_allow_html=True,
        )
        if st.button("الدخول إلى قسم Reaction Lab", key="home_go_reaction"):
            on_navigate("reaction_lab")


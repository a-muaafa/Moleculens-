import streamlit as st
import streamlit.components.v1 as components

LOCK_SELECT_SCRIPT = """
<script>
(function () {
  var doc = window.parent.document;
  function lock() {
    doc.querySelectorAll('div[data-baseweb="select"] input').forEach(function (el) {
      el.setAttribute('readonly', 'true');
      el.setAttribute('inputmode', 'none');
      el.setAttribute('autocomplete', 'off');
    });
  }
  lock();
  new MutationObserver(lock).observe(doc.body, { childList: true, subtree: true });
})();
</script>
"""


def inject_platform_styles(lang: str = "ar") -> None:
    is_ar = lang == "ar"
    dir_css = "rtl" if is_ar else "ltr"
    align_css = "right" if is_ar else "left"
    border_side = "border-right" if is_ar else "border-left"

    raw_css = f"""
@import url('https://fonts.googleapis.com/css2?family=Cairo:wght@400;500;600;700;800&family=JetBrains+Mono:wght@400;600;700&display=swap');
:root {{
    --ml-bg: #F4F7FA;
    --ml-surface: #FFFFFF;
    --ml-surface-soft: #EBF2F7;
    --ml-surface-tint: #F8FBFD;
    --ml-primary: #0E5A6D;
    --ml-primary-dark: #0D2B42;
    --ml-secondary: #244C73;
    --ml-accent: #0D9488;
    --ml-accent-soft: rgba(13, 148, 136, 0.10);
    --ml-amber: #B45309;
    --ml-amber-bg: #FFFBEB;
    --ml-emerald: #047857;
    --ml-emerald-bg: #ECFDF5;
    --ml-crimson: #B91C1C;
    --ml-crimson-bg: #FEF2F2;
    --ml-text: #1E293B;
    --ml-text-muted: #52657A;
    --ml-border: rgba(14, 90, 109, 0.12);
    --ml-border-strong: rgba(14, 90, 109, 0.24);
    --ml-radius-xl: 18px;
    --ml-radius-lg: 15px;
    --ml-radius-md: 13px;
    --ml-radius-sm: 10px;
    --ml-shadow-sm: 0 1px 4px rgba(13, 43, 66, 0.04), 0 4px 12px -2px rgba(13, 43, 66, 0.04);
    --ml-shadow: 0 4px 20px -4px rgba(13, 43, 66, 0.07), 0 2px 6px -1px rgba(13, 43, 66, 0.03);
    --ml-shadow-hover: 0 12px 28px -6px rgba(14, 90, 109, 0.14), 0 4px 10px -2px rgba(13, 43, 66, 0.05);
    --ml-ease: cubic-bezier(0.22, 1, 0.36, 1);
}}
@keyframes mlFadeSlideUp {{
    from {{ opacity: 0; transform: translateY(8px); }}
    to {{ opacity: 1; transform: translateY(0); }}
}}
html, body, .stApp {{
    font-family: 'Cairo', sans-serif !important;
    background-color: var(--ml-bg) !important;
    color: var(--ml-text) !important;
    direction: {dir_css} !important;
    text-align: {align_css} !important;
}}
[data-testid="stAppViewContainer"],
[data-testid="stMain"],
[data-testid="stMarkdownContainer"],
[data-testid="stMarkdownContainer"] p,
[data-testid="stMarkdownContainer"] li,
[data-testid="stMarkdownContainer"] h1,
[data-testid="stMarkdownContainer"] h2,
[data-testid="stMarkdownContainer"] h3,
[data-testid="stMarkdownContainer"] h4,
[data-testid="stMarkdownContainer"] h5,
[data-testid="stWidgetLabel"],
[data-testid="stWidgetLabel"] p,
[data-testid="stCaptionContainer"],
[data-testid="stCaptionContainer"] p,
[data-testid="stAlertContainer"],
[data-testid="stAlertContainer"] p,
.ml-hero, .ml-hero h1, .ml-hero p, .ml-hero div,
.ml-card, .ml-card h3, .ml-card h4, .ml-card h5, .ml-card p,
.ml-module-header, .ml-module-header h2, .ml-module-header h3, .ml-module-header p,
.ml-section-card, .ml-section-card h3, .ml-section-card p,
.ml-ai-card, .ml-safety-card, .ml-no-reaction-card,
.ml-tri-card, .ml-suitability-card, .ml-footer, .ml-step-tile {{
    direction: {dir_css} !important;
    text-align: {align_css} !important;
}}
.stApp > header {{
    background: transparent !important;
}}
.block-container {{
    padding-top: 1rem !important;
    padding-bottom: 2.75rem !important;
    max-width: 1320px !important;
}}
#MainMenu, footer {{
    visibility: hidden;
}}
h1, h2, h3, h4, h5 {{
    font-family: 'Cairo', sans-serif !important;
    letter-spacing: 0 !important;
    color: var(--ml-primary-dark);
}}
a {{
    color: var(--ml-primary);
    text-decoration: none;
    font-weight: 600;
    transition: color 0.18s ease;
}}
a:hover {{
    color: var(--ml-accent);
    text-decoration: underline;
}}
.ltr-inline, code, .ml-formula-pill, .ml-equation-banner,
.ml-slogan-pill, .ml-structure-2d-stage, .ml-cpk-legend,
.katex, .katex-display {{
    direction: ltr !important;
    unicode-bidi: embed;
    text-align: left !important;
}}
[data-testid="stSlider"] div[data-baseweb="slider"] {{
    direction: ltr !important;
}}
.ml-equation-banner, .ml-slogan-pill {{
    text-align: center !important;
}}
.ml-hero {{
    position: relative;
    overflow: hidden;
    background: linear-gradient(130deg, #0D2B42 0%, #144966 55%, #0D7377 100%);
    color: #FFFFFF;
    padding: 1.45rem 1.85rem;
    border-radius: var(--ml-radius-xl);
    margin-bottom: 0.9rem;
    border: 1px solid rgba(255, 255, 255, 0.16);
    box-shadow: 0 12px 32px -10px rgba(13, 43, 66, 0.28);
    animation: mlFadeSlideUp 0.38s var(--ml-ease) both;
}}
.ml-hero::before {{
    content: "";
    position: absolute;
    inset: 0;
    background-image:
        radial-gradient(circle at 12% 20%, rgba(20, 184, 166, 0.22) 0%, transparent 44%),
        radial-gradient(circle at 86% 80%, rgba(56, 189, 248, 0.16) 0%, transparent 45%),
        radial-gradient(rgba(255, 255, 255, 0.08) 1px, transparent 1px);
    background-size: 100% 100%, 100% 100%, 22px 22px;
    pointer-events: none;
}}
.ml-hero > * {{
    position: relative;
    z-index: 1;
}}
.ml-hero-top {{
    display: flex;
    justify-content: space-between;
    align-items: center;
    flex-wrap: wrap;
    gap: 0.85rem;
    margin-bottom: 0.55rem;
}}
.ml-brand-title {{
    font-size: 1.95rem;
    font-weight: 800;
    color: #FFFFFF !important;
    margin: 0;
    line-height: 1.2;
}}
.ml-brand-sub {{
    font-size: 0.96rem;
    color: #CBE4EE !important;
    font-weight: 600;
    margin-top: 0.15rem;
}}
.ml-slogan-pill {{
    background: rgba(255, 255, 255, 0.14);
    backdrop-filter: blur(8px);
    border: 1px solid rgba(255, 255, 255, 0.28);
    color: #F8FAFC !important;
    padding: 0.35rem 1rem;
    border-radius: 999px;
    font-weight: 700;
    font-size: 0.86rem;
    white-space: nowrap;
}}
.ml-hero-desc {{
    font-size: 0.95rem;
    line-height: 1.75;
    color: #E6F2F6 !important;
    margin: 0;
}}
.ml-hero-stats {{
    display: flex;
    flex-wrap: wrap;
    gap: 0.5rem;
    margin-top: 0.85rem;
    justify-content: flex-start;
}}
.ml-hero-chip {{
    display: inline-flex;
    align-items: center;
    padding: 0.22rem 0.72rem;
    border-radius: 999px;
    font-size: 0.78rem;
    font-weight: 700;
    color: #E0F7FA !important;
    background: rgba(255, 255, 255, 0.11);
    border: 1px solid rgba(255, 255, 255, 0.18);
}}
div[data-testid="stRadio"]:has([aria-label="الأقسام الرئيسية لمنصة MolecuLens:"]) {{
    background: var(--ml-surface);
    border: 1px solid var(--ml-border);
    border-radius: var(--ml-radius-lg);
    padding: 0.42rem;
    margin-bottom: 1.1rem;
    box-shadow: var(--ml-shadow-sm);
    animation: mlFadeSlideUp 0.34s var(--ml-ease) both;
}}
div[data-testid="stRadio"]:has([aria-label="الأقسام الرئيسية لمنصة MolecuLens:"]) > label[data-testid="stWidgetLabel"] {{
    display: none !important;
}}
[data-testid="stRadioGroup"][aria-label="الأقسام الرئيسية لمنصة MolecuLens:"] {{
    display: flex !important;
    flex-direction: row !important;
    flex-wrap: wrap !important;
    gap: 6px !important;
    width: 100% !important;
    background: var(--ml-surface-soft) !important;
    padding: 5px !important;
    border-radius: 12px !important;
    border: 1px solid rgba(14, 90, 109, 0.08) !important;
}}
[data-testid="stRadioGroup"][aria-label="الأقسام الرئيسية لمنصة MolecuLens:"] > div {{
    flex: 1 1 auto !important;
    display: flex !important;
    margin: 0 !important;
}}
[data-testid="stRadioGroup"][aria-label="الأقسام الرئيسية لمنصة MolecuLens:"] label[data-testid="stRadioOption"] {{
    width: 100% !important;
    display: flex !important;
    align-items: center !important;
    justify-content: center !important;
    margin: 0 !important;
    padding: 0.56rem 0.9rem !important;
    border-radius: 9px !important;
    background: transparent !important;
    border: 1px solid transparent !important;
    cursor: pointer !important;
    transition: all 0.22s var(--ml-ease) !important;
}}
[data-testid="stRadioGroup"][aria-label="الأقسام الرئيسية لمنصة MolecuLens:"] label[data-testid="stRadioOption"] > div {{
    width: 100% !important;
    justify-content: center !important;
}}
[data-testid="stRadioGroup"][aria-label="الأقسام الرئيسية لمنصة MolecuLens:"] label[data-testid="stRadioOption"] > div > div:not([data-testid="stMarkdownContainer"]) {{
    display: none !important;
}}
[data-testid="stRadioGroup"][aria-label="الأقسام الرئيسية لمنصة MolecuLens:"] label[data-testid="stRadioOption"] [data-testid="stMarkdownContainer"] {{
    width: 100% !important;
    text-align: center !important;
}}
[data-testid="stRadioGroup"][aria-label="الأقسام الرئيسية لمنصة MolecuLens:"] label[data-testid="stRadioOption"] p {{
    font-size: 0.86rem !important;
    font-weight: 700 !important;
    color: var(--ml-text-muted) !important;
    margin: 0 !important;
    text-align: center !important;
    white-space: nowrap !important;
    transition: color 0.2s ease !important;
}}
[data-testid="stRadioGroup"][aria-label="الأقسام الرئيسية لمنصة MolecuLens:"] label[data-testid="stRadioOption"]:hover {{
    background: rgba(255, 255, 255, 0.92) !important;
    border-color: var(--ml-border) !important;
    transform: translateY(-1px);
}}
[data-testid="stRadioGroup"][aria-label="الأقسام الرئيسية لمنصة MolecuLens:"] label[data-testid="stRadioOption"]:hover p {{
    color: var(--ml-primary-dark) !important;
}}
[data-testid="stRadioGroup"][aria-label="الأقسام الرئيسية لمنصة MolecuLens:"] label[data-testid="stRadioOption"][data-selected="true"] {{
    background: linear-gradient(135deg, #0E5A6D 0%, #0D9488 100%) !important;
    border-color: transparent !important;
    box-shadow: 0 4px 14px -2px rgba(14, 90, 109, 0.30) !important;
    transform: translateY(-1px);
}}
[data-testid="stRadioGroup"][aria-label="الأقسام الرئيسية لمنصة MolecuLens:"] label[data-testid="stRadioOption"][data-selected="true"] p {{
    color: #FFFFFF !important;
}}
.ml-card {{
    background: var(--ml-surface);
    border: 1px solid var(--ml-border);
    border-radius: var(--ml-radius-md);
    padding: 1.35rem 1.6rem;
    margin-bottom: 1.05rem;
    box-shadow: var(--ml-shadow);
    transition: box-shadow 0.24s var(--ml-ease), border-color 0.24s var(--ml-ease), transform 0.24s var(--ml-ease);
    animation: mlFadeSlideUp 0.35s var(--ml-ease) both;
}}
.ml-card:hover {{
    border-color: var(--ml-border-strong);
    box-shadow: var(--ml-shadow-hover);
}}
.ml-card h3 {{
    margin-top: 0;
    margin-bottom: 0.5rem;
    font-size: 1.28rem;
    font-weight: 800;
    color: var(--ml-primary-dark) !important;
    line-height: 1.4;
}}
.ml-card h4 {{
    margin-top: 0;
    margin-bottom: 0.55rem;
    font-size: 1.08rem;
    font-weight: 800;
    color: var(--ml-primary) !important;
    line-height: 1.45;
}}
.ml-card h5 {{
    margin-top: 0;
    margin-bottom: 0.45rem;
    font-size: 0.98rem;
    font-weight: 800;
    color: var(--ml-secondary) !important;
}}
.ml-card p {{
    margin-bottom: 0.45rem;
    line-height: 1.75;
    color: var(--ml-text);
}}
.ml-module-header {{
    background: linear-gradient(135deg, #FFFFFF 0%, #EDF6F8 100%);
    border: 1px solid var(--ml-border);
    {border_side}: 5px solid var(--ml-primary);
    border-radius: var(--ml-radius-lg);
    padding: 1.25rem 1.55rem;
    margin-bottom: 1.1rem;
    box-shadow: var(--ml-shadow);
    animation: mlFadeSlideUp 0.34s var(--ml-ease) both;
}}
.ml-module-header h2, .ml-module-header h3 {{
    color: var(--ml-primary-dark) !important;
    margin: 0 0 0.35rem 0;
    font-size: 1.34rem;
    font-weight: 800;
}}
.ml-module-header p {{
    color: var(--ml-text-muted);
    margin: 0;
    font-size: 0.94rem;
    line-height: 1.72;
}}
.ml-ai-tag {{
    background: linear-gradient(135deg, var(--ml-primary) 0%, var(--ml-accent) 100%);
    color: #FFFFFF;
    padding: 0.2rem 0.68rem;
    border-radius: 7px;
    font-size: 0.76rem;
    font-weight: 800;
    letter-spacing: 0.04em;
    direction: ltr;
}}
.ml-section-card {{
    background: var(--ml-surface);
    border: 1px solid var(--ml-border);
    border-top: 4px solid var(--ml-primary);
    border-radius: var(--ml-radius-lg);
    padding: 1.35rem 1.5rem;
    height: 100%;
    margin-bottom: 0.65rem;
    box-shadow: var(--ml-shadow);
    transition: transform 0.24s var(--ml-ease), box-shadow 0.24s var(--ml-ease), border-color 0.24s var(--ml-ease);
    animation: mlFadeSlideUp 0.38s var(--ml-ease) both;
}}
.ml-section-card:hover {{
    transform: translateY(-3px);
    box-shadow: var(--ml-shadow-hover);
    border-color: var(--ml-border-strong);
    border-top-color: var(--ml-accent);
}}
.ml-section-card h3 {{
    margin-top: 0;
    margin-bottom: 0.45rem;
    font-size: 1.18rem;
    font-weight: 800;
    color: var(--ml-primary-dark) !important;
}}
.ml-section-card p {{
    color: var(--ml-text-muted);
    line-height: 1.7;
    font-size: 0.92rem;
}}
.ml-module-pill {{
    display: inline-block;
    padding: 0.16rem 0.65rem;
    border-radius: 999px;
    font-size: 0.75rem;
    font-weight: 700;
    background: var(--ml-accent-soft);
    color: var(--ml-primary);
    border: 1px solid rgba(13, 148, 136, 0.20);
    margin-bottom: 0.5rem;
}}
.ml-step-grid {{
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(220px, 1fr));
    gap: 0.85rem;
    margin-top: 0.65rem;
}}
.ml-step-tile {{
    background: var(--ml-surface-tint);
    border: 1px solid var(--ml-border);
    border-radius: 12px;
    padding: 0.95rem 1.05rem;
    transition: transform 0.22s var(--ml-ease), border-color 0.22s ease, background 0.22s ease;
}}
.ml-step-tile:hover {{
    transform: translateY(-2px);
    background: #FFFFFF;
    border-color: var(--ml-accent);
}}
.ml-step-num {{
    display: inline-flex;
    align-items: center;
    justify-content: center;
    width: 26px;
    height: 26px;
    border-radius: 8px;
    background: linear-gradient(135deg, var(--ml-primary) 0%, var(--ml-accent) 100%);
    color: #FFFFFF;
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.78rem;
    font-weight: 700;
    margin-left: 0.45rem;
}}
.ml-safety-card {{
    background: linear-gradient(135deg, #FFFDF5 0%, #FFFBEB 100%);
    border: 1px solid rgba(217, 119, 6, 0.26);
    {border_side}: 5px solid #D97706;
    border-radius: var(--ml-radius-md);
    padding: 1.2rem 1.45rem;
    margin-bottom: 1.05rem;
    box-shadow: var(--ml-shadow-sm);
    animation: mlFadeSlideUp 0.35s var(--ml-ease) both;
}}
.ml-no-reaction-card {{
    background: linear-gradient(135deg, #FFF8F6 0%, #FEF2F2 100%);
    border: 1px solid rgba(194, 65, 12, 0.26);
    {border_side}: 5px solid #C2410C;
    border-radius: var(--ml-radius-md);
    padding: 1.2rem 1.45rem;
    margin-bottom: 1.05rem;
    box-shadow: var(--ml-shadow-sm);
    animation: mlFadeSlideUp 0.35s var(--ml-ease) both;
}}
.ml-no-reaction-card h4 {{
    color: #9A3412;
    margin: 0 0 0.5rem 0;
    font-weight: 800;
}}
.ml-ai-card {{
    background: var(--ml-surface-tint);
    border: 1px solid rgba(14, 90, 109, 0.16);
    {border_side}: 4px solid var(--ml-primary);
    border-radius: 12px;
    padding: 1.1rem 1.35rem;
    margin-top: 0.75rem;
    box-shadow: var(--ml-shadow-sm);
}}
.ml-ai-card h4 {{
    color: var(--ml-primary-dark);
    margin: 0 0 0.5rem 0;
    font-weight: 800;
}}
.ml-ai-impact-box {{
    background: #FFFFFF;
    border: 1px solid var(--ml-border);
    border-radius: 10px;
    padding: 0.75rem 1rem;
    margin: 0.65rem 0;
}}
.ml-ai-step {{
    background: var(--ml-surface);
    border: 1px solid var(--ml-border);
    border-radius: 10px;
    padding: 0.9rem 1.15rem;
    margin-bottom: 0.65rem;
    box-shadow: var(--ml-shadow-sm);
    transition: border-color 0.2s ease, transform 0.2s var(--ml-ease);
}}
.ml-ai-step:hover {{
    border-color: var(--ml-accent);
    transform: translateX(-2px);
}}
.ml-ai-answer {{
    background: linear-gradient(135deg, #F0FDF9 0%, #FFFFFF 100%);
    border: 1px solid rgba(13, 148, 136, 0.24);
    {border_side}: 4px solid var(--ml-accent);
    border-radius: 12px;
    padding: 1.1rem 1.3rem;
    margin-top: 0.85rem;
    box-shadow: var(--ml-shadow-sm);
}}
.ml-species-card {{
    background: var(--ml-surface);
    border: 1px solid var(--ml-border);
    border-top: 3px solid var(--ml-primary);
    border-radius: 12px;
    padding: 0.85rem;
    text-align: center !important;
    box-shadow: var(--ml-shadow-sm);
    transition: transform 0.2s var(--ml-ease), box-shadow 0.2s var(--ml-ease);
}}
.ml-species-card:hover {{
    transform: translateY(-2px);
    box-shadow: var(--ml-shadow);
}}
.ml-species-card-product {{
    border-top-color: var(--ml-emerald);
    background: linear-gradient(180deg, #F4FBF8 0%, #FFFFFF 100%);
}}
.ml-stage-card {{
    background: var(--ml-surface-tint);
    border: 1px solid var(--ml-border);
    border-radius: 12px;
    padding: 0.95rem;
    text-align: center !important;
    box-shadow: var(--ml-shadow-sm);
}}
.ml-suitability-card {{
    background: var(--ml-surface);
    border: 1px solid var(--ml-border);
    border-radius: 12px;
    padding: 1rem 1.25rem;
    margin-top: 1rem;
    box-shadow: var(--ml-shadow-sm);
}}
.ml-formula-pill {{
    display: inline-block;
    font-family: 'JetBrains Mono', monospace;
    font-size: 1.06rem;
    font-weight: 700;
    color: var(--ml-primary);
    background: rgba(14, 90, 109, 0.08);
    border: 1px solid rgba(14, 90, 109, 0.22);
    padding: 0.14rem 0.72rem;
    border-radius: 8px;
    margin: 0 0.25rem;
}}
.ml-structure-2d-stage {{
    overflow: auto;
    background: #FFFFFF;
    border: 1.5px solid #D7E3EC;
    border-radius: var(--ml-radius-md);
    padding: 1.25rem;
    margin: 0.75rem 0 1.1rem 0;
    text-align: center !important;
    min-height: 360px;
    display: flex;
    align-items: center;
    justify-content: center;
    box-shadow: inset 0 2px 10px rgba(13, 43, 66, 0.03);
}}
.ml-structure-2d-stage svg {{
    max-width: 100%;
    height: auto;
    display: block;
    margin: 0 auto;
}}
.ml-structure-3d-stage {{
    background: #FFFFFF;
    border: 1.5px solid #D7E3EC;
    border-radius: var(--ml-radius-md);
    padding: 0.65rem;
    margin: 0.65rem 0 1rem 0;
    box-shadow: inset 0 2px 10px rgba(13, 43, 66, 0.03);
}}
.ml-cpk-legend {{
    display: flex;
    flex-wrap: wrap;
    gap: 0.5rem;
    align-items: center;
    margin: 0.5rem 0 0.75rem 0;
}}
.ml-cpk-pill {{
    display: inline-flex;
    align-items: center;
    gap: 0.35rem;
    padding: 0.2rem 0.65rem;
    border-radius: 999px;
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.8rem;
    font-weight: 700;
    background: #FFFFFF;
    border: 1px solid #D7E3EC;
    color: var(--ml-text);
}}
.ml-cpk-dot {{
    width: 10px;
    height: 10px;
    border-radius: 50%;
    display: inline-block;
}}
.ml-tri-grid {{
    display: grid;
    grid-template-columns: repeat(3, minmax(260px, 1fr));
    gap: 1rem;
    margin-bottom: 1.2rem;
}}
.ml-tri-card {{
    border-radius: var(--ml-radius-md);
    padding: 1.2rem 1.35rem;
    box-shadow: var(--ml-shadow-sm);
    height: 100%;
    display: flex;
    flex-direction: column;
    justify-content: space-between;
    transition: transform 0.22s var(--ml-ease), box-shadow 0.22s var(--ml-ease);
    animation: mlFadeSlideUp 0.38s var(--ml-ease) both;
}}
.ml-tri-card:hover {{
    transform: translateY(-2px);
    box-shadow: var(--ml-shadow);
}}
.ml-tri-documented {{
    background: linear-gradient(180deg, #F5F8FF 0%, #FFFFFF 100%);
    border: 1px solid rgba(37, 99, 235, 0.20);
    border-top: 4px solid #2563EB;
}}
.ml-tri-calculated {{
    background: linear-gradient(180deg, #F2FBF9 0%, #FFFFFF 100%);
    border: 1px solid rgba(13, 148, 136, 0.22);
    border-top: 4px solid #0D9488;
}}
.ml-tri-theoretical {{
    background: linear-gradient(180deg, #FFFBF2 0%, #FFFFFF 100%);
    border: 1px solid rgba(217, 119, 6, 0.22);
    border-top: 4px solid #D97706;
}}
.ml-tri-title {{
    font-size: 1.02rem;
    font-weight: 800;
    margin-bottom: 0.5rem;
}}
.ml-tri-badge {{
    display: inline-block;
    padding: 0.16rem 0.65rem;
    border-radius: 999px;
    font-size: 0.76rem;
    font-weight: 700;
    margin-bottom: 0.5rem;
}}
.ml-badges-row {{
    display: flex;
    flex-wrap: wrap;
    gap: 0.45rem;
    margin-bottom: 0.65rem;
}}
.ml-badge {{
    display: inline-flex;
    align-items: center;
    padding: 0.22rem 0.78rem;
    border-radius: 999px;
    font-size: 0.79rem;
    font-weight: 700;
    color: var(--ml-primary-dark);
    background: #EEF4F8;
    border: 1px solid #CFDDE7;
}}
.ml-badge-element {{ background: #EFF6FF; color: #1D4ED8; border-color: #BFDBFE; }}
.ml-badge-compound {{ background: #ECFDF5; color: #047857; border-color: #A7F3D0; }}
.ml-badge-drug {{ background: #F5F3FF; color: #6D28D9; border-color: #DDD6FE; }}
.ml-badge-protein {{ background: #FFFBEB; color: #B45309; border-color: #FDE68A; }}
.ml-badge-neutral {{ background: #F1F5F9; color: #334155; border-color: #CBD5E1; }}
.ml-badge-local {{ background: #F0FDF9; color: #0F766E; border-color: #99F6E4; }}
.ml-equation-banner {{
    font-family: 'JetBrains Mono', monospace;
    font-size: 1.32rem;
    font-weight: 700;
    background: linear-gradient(135deg, #EDF6F8 0%, #F8FBFD 100%);
    border: 1.5px solid rgba(14, 90, 109, 0.22);
    border-radius: 14px;
    padding: 0.95rem 1.25rem;
    margin: 0.65rem 0 0.9rem 0;
    color: var(--ml-primary-dark);
    overflow-x: auto;
}}
.ml-yield-box {{
    background: linear-gradient(135deg, #0D2B42 0%, #0E5A6D 55%, #0D9488 100%);
    color: #FFFFFF;
    padding: 1.35rem 1.6rem;
    border-radius: var(--ml-radius-lg);
    text-align: center !important;
    margin-bottom: 1.1rem;
    box-shadow: 0 10px 26px -6px rgba(14, 90, 109, 0.28);
    animation: mlFadeSlideUp 0.38s var(--ml-ease) both;
}}
.ml-yield-value {{
    font-size: 1.95rem;
    font-weight: 800;
    color: #FFFFFF;
    direction: ltr;
    unicode-bidi: embed;
    font-family: 'JetBrains Mono', monospace;
}}
[data-testid="stElementContainer"]:has([data-testid="stButton"]),
[data-testid="stElementContainer"]:has([data-testid="stDownloadButton"]),
[data-testid="stButton"],
[data-testid="stDownloadButton"],
.stButton,
.stDownloadButton {{
    width: 100% !important;
}}
.stButton > button,
[data-testid="stButton"] button {{
    width: 100% !important;
    min-height: 2.58rem !important;
    border-radius: 10px !important;
    font-family: 'Cairo', sans-serif !important;
    font-weight: 700 !important;
    font-size: 0.9rem !important;
    padding: 0.48rem 0.95rem !important;
    background: #FFFFFF !important;
    color: var(--ml-primary-dark) !important;
    border: 1.5px solid rgba(14, 90, 109, 0.22) !important;
    transition: all 0.22s var(--ml-ease) !important;
    box-shadow: var(--ml-shadow-sm) !important;
}}
.stButton > button:hover,
[data-testid="stButton"] button:hover {{
    background: linear-gradient(135deg, #0E5A6D 0%, #0D9488 100%) !important;
    color: #FFFFFF !important;
    border-color: transparent !important;
    transform: translateY(-2px);
    box-shadow: 0 6px 16px -3px rgba(14, 90, 109, 0.28) !important;
}}
.stButton > button:active,
.stButton > button:focus,
[data-testid="stButton"] button:active,
[data-testid="stButton"] button:focus {{
    background: #0D2B42 !important;
    color: #FFFFFF !important;
    border-color: #0E5A6D !important;
    outline: none !important;
    transform: translateY(0);
}}
.stDownloadButton > button,
[data-testid="stDownloadButton"] button {{
    width: 100% !important;
    min-height: 2.7rem !important;
    border-radius: 10px !important;
    font-family: 'Cairo', sans-serif !important;
    font-weight: 700 !important;
    font-size: 0.9rem !important;
    background: linear-gradient(135deg, #0D2B42 0%, #0E5A6D 100%) !important;
    color: #FFFFFF !important;
    border: none !important;
    transition: all 0.22s var(--ml-ease) !important;
    box-shadow: 0 4px 14px -2px rgba(13, 43, 66, 0.22) !important;
}}
.stDownloadButton > button:hover,
[data-testid="stDownloadButton"] button:hover {{
    background: linear-gradient(135deg, #0E5A6D 0%, #0D9488 100%) !important;
    color: #FFFFFF !important;
    transform: translateY(-2px);
    box-shadow: 0 8px 20px -4px rgba(14, 90, 109, 0.32) !important;
}}
[data-testid="stMetric"] {{
    background: var(--ml-surface);
    border: 1px solid var(--ml-border);
    border-radius: var(--ml-radius-md);
    padding: 0.95rem 1.1rem;
    box-shadow: var(--ml-shadow-sm);
    transition: transform 0.22s var(--ml-ease), border-color 0.22s ease, box-shadow 0.22s ease;
    animation: mlFadeSlideUp 0.34s var(--ml-ease) both;
}}
[data-testid="stMetric"]:hover {{
    transform: translateY(-2px);
    border-color: var(--ml-border-strong);
    box-shadow: var(--ml-shadow);
}}
[data-testid="stMetricLabel"] p {{
    font-weight: 700 !important;
    font-size: 0.84rem !important;
    color: var(--ml-text-muted) !important;
}}
[data-testid="stMetricValue"] {{
    color: var(--ml-primary-dark) !important;
    font-weight: 800 !important;
    font-size: 1.15rem !important;
    line-height: 1.35 !important;
    white-space: normal !important;
    overflow: visible !important;
}}
[data-testid="stMetricValue"] > div {{
    white-space: normal !important;
    overflow: visible !important;
    text-overflow: unset !important;
}}
.stTabs [data-baseweb="tab-list"] {{
    gap: 0.4rem;
    background: var(--ml-surface-soft);
    padding: 0.35rem;
    border-radius: 12px;
    border: 1px solid var(--ml-border);
}}
.stTabs [data-baseweb="tab"] {{
    border-radius: 9px;
    padding: 0.5rem 1rem;
    font-family: 'Cairo', sans-serif !important;
    font-weight: 700 !important;
    font-size: 0.9rem;
    color: var(--ml-text-muted);
    background: transparent;
    border: none !important;
    transition: all 0.2s var(--ml-ease);
}}
.stTabs [data-baseweb="tab"]:hover {{
    background: rgba(255, 255, 255, 0.78);
    color: var(--ml-primary-dark);
}}
.stTabs [aria-selected="true"] {{
    background: #FFFFFF !important;
    color: var(--ml-primary-dark) !important;
    box-shadow: 0 2px 8px rgba(13, 43, 66, 0.08);
}}
.stTabs [data-baseweb="tab-highlight"] {{
    display: none !important;
}}
[data-baseweb="tab-panel"] {{
    animation: mlFadeSlideUp 0.28s var(--ml-ease) both;
    padding-top: 0.85rem !important;
}}
[data-testid="stExpander"] {{
    background: var(--ml-surface);
    border: 1px solid var(--ml-border) !important;
    border-radius: var(--ml-radius-md) !important;
    box-shadow: var(--ml-shadow-sm);
    overflow: hidden;
    margin-bottom: 0.85rem;
}}
[data-testid="stExpander"] summary {{
    font-weight: 700 !important;
    color: var(--ml-primary-dark) !important;
    padding: 0.8rem 1.1rem !important;
    transition: background 0.18s ease;
}}
[data-testid="stExpander"] summary:hover {{
    background: var(--ml-surface-tint);
}}
[data-testid="stTextInput"] > div > div {{
    border-radius: 12px !important;
    border: 1.5px solid #CBD8E3 !important;
    background: #FFFFFF !important;
    box-shadow: var(--ml-shadow-sm) !important;
    min-height: 46px !important;
    transition: border-color 0.2s ease, box-shadow 0.2s ease !important;
}}
[data-testid="stTextInput"] > div > div:focus-within {{
    border-color: var(--ml-accent) !important;
    box-shadow: 0 0 0 3.5px rgba(13, 148, 136, 0.14) !important;
}}
[data-testid="stTextInput"] input {{
    font-family: 'Cairo', sans-serif !important;
    font-size: 1rem !important;
    font-weight: 600 !important;
    color: var(--ml-primary-dark) !important;
    padding: 0.55rem 0.9rem !important;
}}
[data-testid="stNumberInput"] > div > div {{
    border-radius: 10px !important;
    border: 1.5px solid #CBD8E3 !important;
    background: #FFFFFF !important;
    transition: border-color 0.2s ease, box-shadow 0.2s ease !important;
}}
[data-testid="stNumberInput"] > div > div:focus-within {{
    border-color: var(--ml-accent) !important;
    box-shadow: 0 0 0 3px rgba(13, 148, 136, 0.14) !important;
}}
.ml-highlight {{
    background: #FEF08A;
    color: #111827;
    padding: 0 4px;
    border-radius: 4px;
    font-weight: 800;
}}
.ml-suggestion {{
    background: var(--ml-surface);
    border: 1px solid var(--ml-border);
    {border_side}: 4px solid var(--ml-primary);
    border-radius: var(--ml-radius-sm);
    padding: 0.58rem 1rem;
    display: flex;
    gap: 0.7rem;
    align-items: center;
    flex-wrap: wrap;
    box-shadow: var(--ml-shadow-sm);
    transition: border-color 0.2s ease, transform 0.2s var(--ml-ease), box-shadow 0.2s ease;
}}
.ml-suggestion:hover {{
    border-color: var(--ml-accent);
    transform: translateX(-3px);
    box-shadow: var(--ml-shadow);
}}
.ml-sug-name {{ font-weight: 700; font-size: 0.96rem; color: var(--ml-primary-dark); }}
.ml-sug-sep {{ color: #94A3B8; }}
.ml-sug-formula {{
    background: rgba(14, 90, 109, 0.09);
    border: 1px solid rgba(14, 90, 109, 0.20);
    color: var(--ml-primary);
    padding: 0.08rem 0.6rem;
    border-radius: 6px;
    font-family: 'JetBrains Mono', monospace;
    font-weight: 700;
}}
.ml-sug-class {{
    margin-right: auto;
    font-size: 0.8rem;
    color: var(--ml-accent);
    background: var(--ml-accent-soft);
    padding: 0.12rem 0.65rem;
    border-radius: 999px;
    font-weight: 700;
}}
div[data-baseweb="select"] > div {{
    border-radius: var(--ml-radius-sm) !important;
    border: 1.5px solid #CBD8E3 !important;
    background: #FFFFFF !important;
    cursor: pointer !important;
    transition: border-color 0.2s ease;
}}
div[data-baseweb="select"] > div:hover {{
    border-color: var(--ml-accent) !important;
}}
div[data-baseweb="select"] input {{
    caret-color: transparent !important;
    pointer-events: none !important;
    user-select: none !important;
    cursor: pointer !important;
}}
[data-testid="stDataFrame"] {{
    border-radius: var(--ml-radius-sm);
    overflow: hidden;
    border: 1px solid var(--ml-border);
    box-shadow: var(--ml-shadow-sm);
}}
.ml-footer {{
    background: var(--ml-surface);
    border: 1px solid var(--ml-border);
    border-radius: var(--ml-radius-md);
    padding: 0.9rem 1.25rem;
    margin-top: 1.35rem;
    color: var(--ml-text-muted);
    font-size: 0.85rem;
    line-height: 1.7;
    box-shadow: var(--ml-shadow-sm);
}}
@media (prefers-reduced-motion: reduce) {{
    *, *::before, *::after {{
        animation-duration: 0.01ms !important;
        transition-duration: 0.01ms !important;
    }}
}}
@media (max-width: 992px) {{
    .ml-tri-grid {{ grid-template-columns: 1fr; }}
}}
@media (max-width: 768px) {{
    .block-container {{
        padding-top: 0.75rem !important;
        padding-left: 0.65rem !important;
        padding-right: 0.65rem !important;
    }}
    .ml-hero {{ padding: 1.15rem; }}
    .ml-brand-title {{ font-size: 1.55rem; }}
    .ml-card, .ml-module-header, .ml-section-card {{ padding: 1rem; }}
    .ml-structure-2d-stage {{ padding: 0.65rem; min-height: 260px; }}
    .ml-sug-class {{ margin-right: 0; }}
}}
"""
    st.markdown(f"<style>{raw_css}</style>", unsafe_allow_html=True)
    components.html(LOCK_SELECT_SCRIPT, height=0)

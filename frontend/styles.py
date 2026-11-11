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


def inject_platform_styles() -> None:
    st.markdown(
        """
<style>
@import url('https://fonts.googleapis.com/css2?family=Cairo:wght@400;500;600;700;800&family=JetBrains+Mono:wght@400;600;700&display=swap');

:root {
    --ml-primary: #0f4c5c;
    --ml-primary-dark: #0b192c;
    --ml-secondary: #1e3e62;
    --ml-accent: #0f766e;
    --ml-radius-lg: 16px;
    --ml-radius-md: 14px;
    --ml-radius-sm: 10px;
    --ml-border: rgba(148, 163, 184, 0.30);
    --ml-shadow: 0 4px 16px rgba(15, 23, 42, 0.06);
}

html, body, [class*="css"], .stApp {
    font-family: 'Cairo', sans-serif;
    direction: rtl;
    text-align: right;
}

.block-container {
    padding-top: 1.25rem;
    padding-bottom: 3rem;
    max-width: 1340px;
}

#MainMenu, footer {
    visibility: hidden;
}

/* Unified Heading Typography Scale */
h1, h2, h3, h4, h5 {
    font-family: 'Cairo', sans-serif !important;
    letter-spacing: 0 !important;
}

.ml-hero {
    background: linear-gradient(135deg, #0b192c 0%, #1e3e62 58%, #0f4c5c 100%);
    color: #ffffff;
    padding: 2rem 2.25rem;
    border-radius: var(--ml-radius-lg);
    margin-bottom: 1.35rem;
    border: 1px solid rgba(255, 255, 255, 0.16);
    box-shadow: 0 10px 28px rgba(11, 25, 44, 0.18);
}

.ml-hero-top {
    display: flex;
    justify-content: space-between;
    align-items: center;
    flex-wrap: wrap;
    gap: 1rem;
    margin-bottom: 0.75rem;
}

.ml-brand-title {
    font-size: 2.15rem;
    font-weight: 800;
    color: #ffffff !important;
    margin: 0;
    line-height: 1.25;
}

.ml-brand-sub {
    font-size: 1rem;
    color: #cbd5e1;
    font-weight: 600;
}

.ml-slogan-pill {
    background: rgba(255, 255, 255, 0.14);
    border: 1px solid rgba(255, 255, 255, 0.25);
    color: #f8fafc;
    padding: 0.35rem 1.05rem;
    border-radius: 999px;
    font-weight: 700;
    font-size: 0.92rem;
    direction: ltr;
}

.ml-hero-desc {
    font-size: 1.03rem;
    line-height: 1.85;
    color: #e2e8f0 !important;
    margin: 0;
}

/* Unified Card Design */
.ml-card {
    background: var(--secondary-background-color);
    border: 1px solid var(--ml-border);
    border-radius: var(--ml-radius-md);
    padding: 1.45rem 1.65rem;
    margin-bottom: 1.15rem;
    box-shadow: var(--ml-shadow);
    transition: box-shadow 0.2s ease, border-color 0.2s ease;
}

.ml-card:hover {
    border-color: rgba(15, 76, 92, 0.38);
}

.ml-card h3 {
    margin-top: 0;
    margin-bottom: 0.6rem;
    font-size: 1.38rem;
    font-weight: 800;
    line-height: 1.4;
}

.ml-card h4 {
    margin-top: 0;
    margin-bottom: 0.65rem;
    font-size: 1.15rem;
    font-weight: 800;
    color: #0f4c5c;
    line-height: 1.45;
}

.ml-card p {
    margin-bottom: 0.45rem;
    line-height: 1.75;
}

.ml-section-card {
    background: var(--secondary-background-color);
    border: 1px solid rgba(15, 76, 92, 0.3);
    border-top: 5px solid #0f4c5c;
    border-radius: var(--ml-radius-md);
    padding: 1.5rem 1.6rem;
    height: 100%;
    margin-bottom: 1rem;
    box-shadow: var(--ml-shadow);
}

.ml-safety-card {
    background: linear-gradient(135deg, rgba(245, 158, 11, 0.08) 0%, rgba(15, 76, 92, 0.06) 100%);
    border: 1px solid rgba(217, 119, 6, 0.35);
    border-right: 5px solid #d97706;
    border-radius: var(--ml-radius-md);
    padding: 1.25rem 1.5rem;
    margin-bottom: 1.15rem;
    box-shadow: var(--ml-shadow);
}

/* Prominent Chemical Formula Display */
.ml-formula-pill {
    display: inline-block;
    direction: ltr;
    unicode-bidi: embed;
    font-family: 'JetBrains Mono', monospace;
    font-size: 1.18rem;
    font-weight: 700;
    color: #0f4c5c;
    background: rgba(15, 76, 92, 0.10);
    border: 1px solid rgba(15, 76, 92, 0.26);
    padding: 0.18rem 0.78rem;
    border-radius: 8px;
    margin: 0 0.25rem;
}

/* Prominent 2D & 3D Molecular Structure Stages */
.ml-structure-2d-stage {
    overflow: auto;
    direction: ltr;
    background: #ffffff;
    border: 2px solid #cbd5e1;
    border-radius: 14px;
    padding: 1.25rem;
    margin: 0.75rem 0 1.1rem 0;
    text-align: center;
    min-height: 360px;
    display: flex;
    align-items: center;
    justify-content: center;
    box-shadow: inset 0 2px 8px rgba(15, 23, 42, 0.04);
}

.ml-structure-2d-stage svg {
    max-width: 100%;
    height: auto;
    display: block;
    margin: 0 auto;
}

.ml-structure-3d-stage {
    background: #ffffff;
    border: 2px solid #cbd5e1;
    border-radius: 14px;
    padding: 0.65rem;
    margin: 0.65rem 0 1rem 0;
    box-shadow: inset 0 2px 8px rgba(15, 23, 42, 0.04);
}

.ml-cpk-legend {
    display: flex;
    flex-wrap: wrap;
    gap: 0.5rem;
    align-items: center;
    margin: 0.5rem 0 0.75rem 0;
    direction: ltr;
}

.ml-cpk-pill {
    display: inline-flex;
    align-items: center;
    gap: 0.35rem;
    padding: 0.18rem 0.6rem;
    border-radius: 999px;
    font-family: 'JetBrains Mono', monospace;
    font-size: 0.8rem;
    font-weight: 700;
    background: #f8fafc;
    border: 1px solid #cbd5e1;
    color: #1e293b;
}

.ml-cpk-dot {
    width: 10px;
    height: 10px;
    border-radius: 50%;
    display: inline-block;
}

/* Three-Way Visual Separation Cards (Documented vs Calculated vs Theoretical) */
.ml-tri-grid {
    display: grid;
    grid-template-columns: repeat(3, minmax(260px, 1fr));
    gap: 1.1rem;
    margin-bottom: 1.25rem;
}

.ml-tri-card {
    border-radius: var(--ml-radius-md);
    padding: 1.25rem 1.4rem;
    box-shadow: var(--ml-shadow);
    height: 100%;
    display: flex;
    flex-direction: column;
    justify-content: space-between;
}

.ml-tri-documented {
    background: rgba(29, 78, 216, 0.06);
    border: 1px solid rgba(29, 78, 216, 0.30);
    border-top: 5px solid #1d4ed8;
}

.ml-tri-calculated {
    background: rgba(15, 118, 110, 0.06);
    border: 1px solid rgba(15, 118, 110, 0.30);
    border-top: 5px solid #0f766e;
}

.ml-tri-theoretical {
    background: rgba(180, 83, 9, 0.06);
    border: 1px solid rgba(180, 83, 9, 0.30);
    border-top: 5px solid #b45309;
}

.ml-tri-title {
    font-size: 1.06rem;
    font-weight: 800;
    margin-bottom: 0.55rem;
}

.ml-tri-badge {
    display: inline-block;
    padding: 0.16rem 0.65rem;
    border-radius: 999px;
    font-size: 0.78rem;
    font-weight: 700;
    margin-bottom: 0.55rem;
}

/* Unified Badges */
.ml-badge {
    display: inline-block;
    padding: 0.25rem 0.85rem;
    border-radius: 999px;
    font-size: 0.84rem;
    font-weight: 700;
    color: #ffffff;
    background: #0f4c5c;
    margin-left: 0.45rem;
    margin-bottom: 0.25rem;
}

.ml-badge-element { background: #1d4ed8; }
.ml-badge-compound { background: #0f766e; }
.ml-badge-drug { background: #6d28d9; }
.ml-badge-protein { background: #b45309; }

/* Unified Equation & Yield Display */
.ml-equation-banner {
    direction: ltr;
    unicode-bidi: embed;
    text-align: center;
    font-family: 'JetBrains Mono', monospace;
    font-size: 1.38rem;
    font-weight: 700;
    background: rgba(15, 76, 92, 0.08);
    border: 1.5px solid rgba(15, 76, 92, 0.28);
    border-radius: 12px;
    padding: 1rem 1.25rem;
    margin: 0.65rem 0 0.9rem 0;
    color: #0f4c5c;
    overflow-x: auto;
}

.ml-yield-box {
    background: linear-gradient(135deg, #0f4c5c 0%, #1e3e62 100%);
    color: #ffffff;
    padding: 1.4rem 1.65rem;
    border-radius: var(--ml-radius-md);
    text-align: center;
    margin-bottom: 1.15rem;
    box-shadow: 0 6px 18px rgba(15, 76, 92, 0.18);
}

.ml-yield-value {
    font-size: 2rem;
    font-weight: 800;
    color: #ffffff;
    direction: ltr;
    unicode-bidi: embed;
    font-family: 'JetBrains Mono', monospace;
}

.ltr-inline, code {
    direction: ltr;
    unicode-bidi: embed;
    font-family: 'JetBrains Mono', monospace;
    display: inline-block;
    text-align: left;
}

/* Unified Buttons & Interactive States */
.stButton > button {
    width: 100%;
    min-height: 2.68rem;
    border-radius: var(--ml-radius-sm);
    font-family: 'Cairo', sans-serif;
    font-weight: 700;
    font-size: 0.94rem;
    padding: 0.55rem 0.95rem;
    background: #1e3e62;
    color: #ffffff !important;
    border: 1px solid #0f4c5c;
    transition: all 0.18s ease-in-out;
    box-shadow: 0 2px 5px rgba(15, 23, 42, 0.08);
}

.stButton > button:hover {
    background: #0f4c5c;
    color: #ffffff !important;
    border-color: #0f766e;
    transform: translateY(-1px);
    box-shadow: 0 4px 10px rgba(15, 76, 92, 0.2);
}

.stButton > button:active, .stButton > button:focus {
    background: #0b192c;
    color: #ffffff !important;
    border-color: #0f4c5c;
    outline: none;
}

/* Yellow Highlight & Live Suggestions */
.ml-highlight {
    background: #fef08a;
    color: #111827;
    padding: 0 4px;
    border-radius: 4px;
    font-weight: 800;
}

.ml-suggestion {
    background: var(--secondary-background-color);
    border: 1px solid var(--ml-border);
    border-right: 4px solid #0f4c5c;
    border-radius: var(--ml-radius-sm);
    padding: 0.58rem 1rem;
    margin-bottom: 0.42rem;
    display: flex;
    gap: 0.7rem;
    align-items: center;
    flex-wrap: wrap;
    transition: border-color 0.15s ease, background 0.15s ease;
}

.ml-suggestion:hover {
    border-color: #0f766e;
}

.ml-sug-name { font-weight: 700; font-size: 0.98rem; }
.ml-sug-sep { color: #94a3b8; }
.ml-sug-formula {
    background: rgba(15, 76, 92, 0.13);
    border: 1px solid rgba(15, 76, 92, 0.22);
    padding: 0.08rem 0.6rem;
    border-radius: 6px;
    font-family: 'JetBrains Mono', monospace;
    font-weight: 700;
    direction: ltr;
    unicode-bidi: embed;
}
.ml-sug-class {
    margin-right: auto;
    font-size: 0.82rem;
    color: #0f766e;
    background: rgba(15, 118, 110, 0.10);
    padding: 0.1rem 0.6rem;
    border-radius: 999px;
    font-weight: 700;
}

/* Unified Closed Dropdown Styling */
div[data-baseweb="select"] > div {
    border-radius: var(--ml-radius-sm) !important;
    border-color: var(--ml-border) !important;
    cursor: pointer !important;
}

div[data-baseweb="select"] input {
    caret-color: transparent !important;
    pointer-events: none !important;
    user-select: none !important;
    cursor: pointer !important;
}

/* Unified Table Styling */
[data-testid="stDataFrame"] {
    border-radius: var(--ml-radius-sm);
    overflow: hidden;
    border: 1px solid var(--ml-border);
}

/* Responsive Mobile & Tablet Layout */
@media (max-width: 992px) {
    .ml-tri-grid {
        grid-template-columns: 1fr;
    }
}

@media (max-width: 768px) {
    .block-container {
        padding-top: 0.75rem;
        padding-left: 0.65rem;
        padding-right: 0.65rem;
        padding-bottom: 2rem;
    }
    .ml-hero {
        padding: 1.25rem 1.15rem;
    }
    .ml-brand-title {
        font-size: 1.6rem;
    }
    .ml-hero-desc {
        font-size: 0.94rem;
        line-height: 1.65;
    }
    .ml-card {
        padding: 1.05rem 1.1rem;
    }
    .ml-card h3 {
        font-size: 1.18rem;
    }
    .ml-card h4 {
        font-size: 1.04rem;
    }
    .ml-equation-banner {
        font-size: 1.08rem;
        padding: 0.75rem 0.85rem;
    }
    .ml-yield-value {
        font-size: 1.48rem;
    }
    .ml-structure-2d-stage {
        padding: 0.65rem;
        min-height: 260px;
    }
    .ml-sug-class {
        margin-right: 0;
    }
}
</style>
""",
        unsafe_allow_html=True,
    )
    components.html(LOCK_SELECT_SCRIPT, height=0)

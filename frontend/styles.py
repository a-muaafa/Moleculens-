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
@import url('https://fonts.googleapis.com/css2?family=Cairo:wght@400;500;600;700;800&family=JetBrains+Mono:wght@400;600&display=swap');

html, body, [class*="css"], .stApp {
    font-family: 'Cairo', sans-serif;
    direction: rtl;
    text-align: right;
}

.block-container {
    padding-top: 1.2rem;
    padding-bottom: 3rem;
    max-width: 1320px;
}

#MainMenu, footer {
    visibility: hidden;
}

.ml-hero {
    background: linear-gradient(135deg, #0b192c 0%, #1e3e62 60%, #0f4c5c 100%);
    color: #ffffff;
    padding: 2rem 2.2rem;
    border-radius: 16px;
    margin-bottom: 1.25rem;
    border: 1px solid rgba(255, 255, 255, 0.14);
    box-shadow: 0 10px 28px rgba(0, 0, 0, 0.16);
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
    line-height: 1.2;
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
    padding: 0.35rem 1rem;
    border-radius: 999px;
    font-weight: 700;
    font-size: 0.92rem;
    direction: ltr;
}

.ml-hero-desc {
    font-size: 1.03rem;
    line-height: 1.8;
    color: #e2e8f0 !important;
    margin: 0;
}

.ml-card {
    background: var(--secondary-background-color);
    border: 1px solid rgba(128, 128, 128, 0.22);
    border-radius: 14px;
    padding: 1.35rem 1.5rem;
    margin-bottom: 1.1rem;
    box-shadow: 0 3px 10px rgba(0, 0, 0, 0.04);
}

.ml-card h3, .ml-card h4 {
    margin-top: 0;
    font-weight: 800;
}

.ml-section-card {
    background: var(--secondary-background-color);
    border: 1px solid rgba(15, 76, 92, 0.3);
    border-top: 5px solid #0f4c5c;
    border-radius: 14px;
    padding: 1.5rem;
    height: 100%;
    margin-bottom: 1rem;
}

.ml-badge {
    display: inline-block;
    padding: 0.25rem 0.85rem;
    border-radius: 999px;
    font-size: 0.84rem;
    font-weight: 700;
    color: #ffffff;
    background: #0f4c5c;
    margin-left: 0.45rem;
}

.ml-badge-element { background: #1d4ed8; }
.ml-badge-compound { background: #0f766e; }
.ml-badge-drug { background: #6d28d9; }
.ml-badge-protein { background: #b45309; }

.ml-yield-box {
    background: linear-gradient(135deg, #0f4c5c 0%, #1e3e62 100%);
    color: #ffffff;
    padding: 1.3rem 1.5rem;
    border-radius: 12px;
    text-align: center;
    margin-bottom: 1rem;
}

.ml-yield-value {
    font-size: 2rem;
    font-weight: 800;
    color: #ffffff;
    direction: ltr;
}

.ltr-inline, code {
    direction: ltr;
    unicode-bidi: embed;
    font-family: 'JetBrains Mono', monospace;
    display: inline-block;
    text-align: left;
}

.stButton > button {
    width: 100%;
    border-radius: 10px;
    font-weight: 700;
    padding: 0.55rem 0.9rem;
    background: #1e3e62;
    color: #ffffff;
    border: 1px solid #0f4c5c;
}

.stButton > button:hover {
    background: #0f4c5c;
    color: #ffffff;
    border-color: #0f4c5c;
}
.ml-highlight {
    background: #fef08a;
    color: #111827;
    padding: 0 2px;
    border-radius: 3px;
    font-weight: 800;
}

.ml-suggestion {
    background: var(--secondary-background-color);
    border: 1px solid rgba(128, 128, 128, 0.25);
    border-radius: 10px;
    padding: 0.5rem 0.9rem;
    margin-bottom: 0.35rem;
    display: flex;
    gap: 0.6rem;
    align-items: center;
    flex-wrap: wrap;
}

.ml-sug-name { font-weight: 700; }
.ml-sug-sep { color: #94a3b8; }
.ml-sug-formula {
    background: rgba(15, 76, 92, 0.12);
    padding: 0 0.5rem;
    border-radius: 6px;
    font-family: 'JetBrains Mono', monospace;
}
.ml-sug-class {
    margin-right: auto;
    font-size: 0.82rem;
    color: #0f766e;
    font-weight: 700;
}

div[data-baseweb="select"] input {
    caret-color: transparent !important;
    pointer-events: none !important;
    user-select: none !important;
    cursor: pointer !important;
}

div[data-baseweb="select"] > div { cursor: pointer !important; }
</style>
""",
        unsafe_allow_html=True,
    )
    components.html(LOCK_SELECT_SCRIPT, height=0)


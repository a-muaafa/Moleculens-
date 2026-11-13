"""
Main entry point for the MolecuLens scientific platform.
Integrates Explore & Understand, Reaction Lab, Titration Lab,
Concentration Lab, and Drug-Drug Interactions with bilingual support.
"""

from typing import Any
import streamlit as st
from config.settings import APP_NAME, APP_SLOGAN
from frontend.explore_view import render_explore_section
from frontend.home_view import render_hero_header, render_home_overview
from frontend.prediction_view import render_prediction_section
from frontend.reaction_lab_view import render_reaction_lab_section
from frontend.titration_view import render_titration_view
from frontend.concentration_view import render_concentration_view
from frontend.ddi_view import render_ddi_view
from frontend.styles import inject_platform_styles

st.set_page_config(page_title=f"{APP_NAME} | {APP_SLOGAN}", layout="wide")

if "app_lang" not in st.session_state:
    st.session_state["app_lang"] = "ar"

if "active_section" not in st.session_state:
    st.session_state["active_section"] = "explore"

# Language Selector in Sidebar
with st.sidebar:
    st.markdown("### Language / اللغة")
    lang_choice = st.radio(
        "اختر لغة الواجهة / Select Language:",
        options=["ar", "en"],
        format_func=lambda x: "العربية (Arabic)" if x == "ar" else "English (الإنجليزية)",
        index=0 if st.session_state["app_lang"] == "ar" else 1,
        key="sidebar_lang_radio",
    )
    if lang_choice != st.session_state["app_lang"]:
        st.session_state["app_lang"] = lang_choice
        st.rerun()

app_lang = st.session_state["app_lang"]
is_ar = app_lang == "ar"

inject_platform_styles(lang=app_lang)


def _navigate_to(section_key: str) -> None:
    st.session_state["active_section"] = section_key
    st.rerun()


def _open_reaction_lab_with_entity(entity_info: dict[str, Any]) -> None:
    st.session_state["passed_reaction_entity"] = entity_info
    st.session_state["active_section"] = "reaction_lab"
    st.rerun()


def _open_prediction_with_entity(entity_info: dict[str, Any]) -> None:
    formula = entity_info.get("symbol_or_formula", "")
    if formula:
        st.session_state["pred_free_input"] = formula
    st.session_state["active_section"] = "prediction"
    st.rerun()


render_hero_header(lang=app_lang)

nav_options_ar = {
    "home": "الصفحة الرئيسية (Home)",
    "explore": "1. Explore & Understand (استكشاف المادة)",
    "prediction": "2. Reaction Prediction (توقع التفاعل)",
    "reaction_lab": "3. Reaction Lab (مختبر التفاعلات)",
    "titration": "4. Titration Lab (مختبر المعايرة)",
    "concentration": "5. Concentration Lab (حاسبة التراكيز)",
    "ddi": "6. Drug Interactions (التداخلات الدوائية)",
}

nav_options_en = {
    "home": "Home Overview",
    "explore": "1. Explore & Understand",
    "prediction": "2. Reaction Prediction",
    "reaction_lab": "3. Reaction Lab",
    "titration": "4. Titration Lab",
    "concentration": "5. Concentration Lab",
    "ddi": "6. Drug Interactions (DDI)",
}

nav_options = nav_options_ar if is_ar else nav_options_en
nav_keys = list(nav_options.keys())
current_idx = nav_keys.index(st.session_state["active_section"]) if st.session_state["active_section"] in nav_keys else 1

chosen_nav = st.radio(
    "الأقسام الرئيسية لمنصة MolecuLens:" if is_ar else "Main Sections of MolecuLens:",
    options=nav_keys,
    index=current_idx,
    format_func=lambda k: nav_options[k],
    horizontal=True,
)

if chosen_nav != st.session_state["active_section"]:
    st.session_state["active_section"] = chosen_nav
    st.rerun()

active_tab = st.session_state["active_section"]

if active_tab == "home":
    render_home_overview(_navigate_to, lang=app_lang)
elif active_tab == "explore":
    render_explore_section(_open_reaction_lab_with_entity, _open_prediction_with_entity, lang=app_lang)
elif active_tab == "prediction":
    render_prediction_section()
elif active_tab == "reaction_lab":
    render_reaction_lab_section()
elif active_tab == "titration":
    render_titration_view(lang=app_lang)
elif active_tab == "concentration":
    render_concentration_view(lang=app_lang)
elif active_tab == "ddi":
    render_ddi_view(lang=app_lang)

st.divider()
st.caption(
    "تنبيه علمي وأكاديمي: المعلومات والحسابات هنا لأغراض التعلم والتهيئة المخبرية فقط، "
    "وليست بديلاً عن الإرشاد الطبي أو الصيدلاني أو أوراق بيانات السلامة (SDS) الرسمية. "
    "Educational use only - not medical advice; always consult the official SDS and your laboratory supervisor."
)

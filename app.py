"""
Main entry point for the MolecuLens scientific platform.
Integrates Explore & Understand, the unified Reaction Lab (Reaction Prediction & Stoichiometry),
Molecular Intelligence, Titration Lab, Concentration Lab, and Pharmacology (Coming Soon).
"""

from typing import Any
import streamlit as st
from config.settings import APP_NAME, APP_SLOGAN
from frontend.explore_view import render_explore_section
from frontend.home_view import render_hero_header, render_home_overview
from frontend.reaction_lab_view import render_reaction_lab_section
from frontend.molecular_intelligence_view import render_molecular_intelligence_section
from frontend.titration_view import render_titration_view
from frontend.concentration_view import render_concentration_view
from frontend.pharmacology_view import render_pharmacology_view
from frontend.styles import inject_platform_styles

st.set_page_config(page_title=f"{APP_NAME} | {APP_SLOGAN}", layout="wide")

app_lang = "ar"
is_ar = True

if "active_section" not in st.session_state:
    st.session_state["active_section"] = "explore"

# Redirect legacy "prediction" navigation key into the unified "reaction_lab"
if st.session_state.get("active_section") == "prediction":
    st.session_state["active_section"] = "reaction_lab"

inject_platform_styles(lang=app_lang)


def _navigate_to(section_key: str) -> None:
    st.session_state["active_section"] = "reaction_lab" if section_key == "prediction" else section_key
    st.rerun()


def _open_reaction_lab_with_entity(entity_info: dict[str, Any]) -> None:
    st.session_state["passed_reaction_entity"] = entity_info
    formula = entity_info.get("symbol_or_formula", "")
    if formula:
        st.session_state["pred_free_input"] = formula
    st.session_state["active_section"] = "reaction_lab"
    st.rerun()


def _open_prediction_with_entity(entity_info: dict[str, Any]) -> None:
    _open_reaction_lab_with_entity(entity_info)


render_hero_header(lang=app_lang)

nav_options = {
    "home": "الصفحة الرئيسية (Home)",
    "explore": "1. Explore & Understand (استكشاف المادة)",
    "reaction_lab": "2. Reaction Lab – مختبر التفاعلات الكيميائية",
    "molecular_intelligence": "3. Molecular Intelligence (الذكاء الجزيئي والتعديل البنيوي)",
    "titration": "4. Titration Lab (مختبر المعايرة)",
    "concentration": "5. Concentration Lab (حاسبة التراكيز)",
    "pharmacology": "6. Pharmacology – علم الأدوية (Coming Soon)",
}

nav_keys = list(nav_options.keys())
current_idx = nav_keys.index(st.session_state["active_section"]) if st.session_state["active_section"] in nav_keys else 1

chosen_nav = st.radio(
    "الأقسام الرئيسية لمنصة MolecuLens:",
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
elif active_tab in {"reaction_lab", "prediction"}:
    render_reaction_lab_section()
elif active_tab == "molecular_intelligence":
    render_molecular_intelligence_section()
elif active_tab == "titration":
    render_titration_view(lang=app_lang)
elif active_tab == "concentration":
    render_concentration_view(lang=app_lang)
elif active_tab == "pharmacology":
    render_pharmacology_view(lang=app_lang)

st.divider()
st.caption(
    "تنبيه علمي وأكاديمي: المعلومات والحسابات هنا لأغراض التعلم والتهيئة المخبرية فقط، "
    "وليست بديلاً عن الإرشاد الطبي أو الصيدلاني أو أوراق بيانات السلامة (SDS) الرسمية."
)

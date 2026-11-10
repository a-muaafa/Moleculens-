from typing import Any
import streamlit as st
from config.settings import APP_NAME, APP_SLOGAN
from frontend.explore_view import render_explore_section
from frontend.home_view import render_hero_header, render_home_overview
from frontend.reaction_lab_view import render_reaction_lab_section
from frontend.styles import inject_platform_styles

st.set_page_config(page_title=f"{APP_NAME} | {APP_SLOGAN}", layout="wide")
inject_platform_styles()

if "active_section" not in st.session_state:
    st.session_state["active_section"] = "explore"


def _navigate_to(section_key: str) -> None:
    st.session_state["active_section"] = section_key
    st.rerun()


def _open_reaction_lab_with_entity(entity_info: dict[str, Any]) -> None:
    st.session_state["passed_reaction_entity"] = entity_info
    st.session_state["active_section"] = "reaction_lab"
    st.rerun()


render_hero_header()

nav_options = {
    "home": "الصفحة الرئيسية (Home Overview)",
    "explore": "1. Explore & Understand (استكشاف المادة وفهم خصائصها)",
    "reaction_lab": "2. Reaction Lab (مختبر التفاعلات والحسابات الكيميائية)",
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

if st.session_state["active_section"] == "home":
    render_home_overview(_navigate_to)
elif st.session_state["active_section"] == "explore":
    render_explore_section(_open_reaction_lab_with_entity)
else:
    render_reaction_lab_section()


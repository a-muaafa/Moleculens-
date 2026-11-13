"""UI smoke tests for the Reaction Prediction lab (Streamlit AppTest, no browser needed)."""

import os

from streamlit.testing.v1 import AppTest

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
APP_FILE = os.path.join(ROOT, "app.py")


def _open_prediction_tab() -> AppTest:
    at = AppTest.from_file(APP_FILE, default_timeout=180)
    at.session_state["active_section"] = "prediction"
    at.run()
    return at


def test_prediction_section_renders_without_errors():
    at = _open_prediction_tab()
    assert not at.exception
    assert any("مختبر التوقع" in md.value for md in at.markdown)


def test_prediction_accepts_free_text_materials():
    at = _open_prediction_tab()
    assert not at.exception
    at.text_input(key="pred_free_input").set_value("Mn + HCl")
    at.run()
    assert not at.exception
    body = " ".join(md.value for md in at.markdown)
    assert "Mn" in body and "HCl" in body


def test_prediction_reports_insufficient_evidence_honestly():
    at = _open_prediction_tab()
    at.text_input(key="pred_free_input").set_value("HCl + HNO3")
    at.run()
    assert not at.exception
    body = " ".join(md.value for md in at.markdown)
    assert "لا تتوفر بيانات كافية" in body


def test_prediction_explains_documented_reaction_with_source():
    at = _open_prediction_tab()
    at.text_input(key="pred_free_input").set_value("HCl + NaOH")
    at.run()
    assert not at.exception
    body = " ".join(md.value for md in at.markdown)
    assert "موثّق" in body
    assert "USP" in body or "المرجع" in body

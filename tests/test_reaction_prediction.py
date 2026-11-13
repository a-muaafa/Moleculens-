"""Tests for the Reaction Prediction engine (prediction tiers, rules, and honesty guarantees)."""

import pytest

from education.ai_assistant import AIScientificAssistant
from reaction_engine.prediction import (
    classify_species,
    is_salt_insoluble,
    predict_reaction,
)


def _inputs(formulas, amounts=None):
    return [
        {
            "formula": f,
            "amount": float(amounts[i]) if amounts else 1.0,
            "unit": "mol",
        }
        for i, f in enumerate(formulas)
    ]


# --------------------------------------------------------------------------------------
# Species classification
# --------------------------------------------------------------------------------------


def test_classify_strong_acid_and_base():
    acid = classify_species("HCl")
    assert acid["kind"] == "strong_acid"
    assert acid["cation"] == "H"
    assert acid["anion"] == "Cl"

    base = classify_species("NaOH")
    assert base["kind"] == "strong_base"
    assert base["cation"] == "Na"
    assert base["cation_charge"] == 1


def test_classify_bicarbonate_is_not_hydroxide():
    """NaHCO3 must be read as sodium bicarbonate, not as a hydroxide/carbonate mix."""
    spec = classify_species("NaHCO3")
    assert spec["kind"] == "salt"
    assert spec["anion"] == "HCO3"
    assert spec["cation"] == "Na"


def test_classify_sulfuric_acid_is_diprotic_strong_acid():
    spec = classify_species("H2SO4")
    assert spec["kind"] == "strong_acid"
    assert spec["proticity"] == 2
    assert spec["anion"] == "SO4"


def test_classify_hydrated_salt_strips_water():
    spec = classify_species("CuSO4.5H2O")
    assert spec["formula"] == "CuSO4"
    assert spec["cation"] == "Cu"
    assert spec["cation_charge"] == 2


def test_classify_hydrocarbon_and_ammonia():
    assert classify_species("CH4")["kind"] == "hydrocarbon"
    assert classify_species("NH3")["kind"] == "ammonia"


def test_classify_invalid_formula():
    assert classify_species("XYZ!!")["is_valid_formula"] is False


# --------------------------------------------------------------------------------------
# Solubility rules
# --------------------------------------------------------------------------------------


def test_solubility_rules():
    assert is_salt_insoluble("Ag", "Cl") is True
    assert is_salt_insoluble("Ba", "SO4") is True
    assert is_salt_insoluble("Ca", "CO3") is True
    assert is_salt_insoluble("Na", "Cl") is False
    assert is_salt_insoluble("Ag", "NO3") is False
    assert is_salt_insoluble("K", "SO4") is False


# --------------------------------------------------------------------------------------
# Tier 1 — documented reactions
# --------------------------------------------------------------------------------------


def test_documented_tier_for_hcl_naoh():
    res = predict_reaction(_inputs(["HCl", "NaOH"]))
    assert res["supported"] is True
    assert res["mode"] == "documented"
    assert res["equation_unicode"].replace(" ", "") == "HCl+NaOH→NaCl+H₂O"
    assert res["evidence"][0]["source_id"]
    assert res["evidence"][0]["reference_url"]


def test_documented_tier_reports_net_ionic_and_ph():
    res = predict_reaction(_inputs(["HCl", "NaOH"]))
    assert "H⁺" in res["net_ionic"]
    assert res["ph_analysis"]["available"] is True


def test_documented_tier_detects_precipitate_observation():
    res = predict_reaction(_inputs(["AgNO3", "NaCl"]))
    types = [o["type"] for o in res["observations"]]
    assert "precipitate" in types


# --------------------------------------------------------------------------------------
# Tier 2 — rule based inference
# --------------------------------------------------------------------------------------


def test_rule_inference_acid_base():
    res = predict_reaction(_inputs(["H2SO4", "Ca(OH)2"]))
    assert res["supported"] is True
    assert res["mode"] == "rule_inferred"
    assert "CaSO₄" in res["equation_unicode"]
    assert res["evidence"][0]["rule_id"] == "acid_base_neutralization"
    assert res["evidence"][0]["source"]


def test_documented_tier_takes_precedence_over_rule_inference():
    """Zn + HCl exists in the verified registry, so the documented tier must win."""
    res = predict_reaction(_inputs(["Zn", "HCl"]))
    assert res["mode"] == "documented"
    assert "H₂" in res["equation_unicode"]
    assert res["evidence"][0]["source_id"]


def test_rule_inference_metal_plus_acid_releases_hydrogen():
    res = predict_reaction(_inputs(["Mn", "HCl"]))
    assert res["mode"] == "rule_inferred"
    assert "H₂" in res["equation_unicode"]
    assert "MnCl₂" in res["equation_unicode"]
    assert any(o["type"] == "gas" for o in res["observations"])


def test_rule_inference_salt_hydroxide_precipitation():
    res = predict_reaction(_inputs(["MgSO4", "KOH"]))
    assert res["supported"] is True
    assert res["mode"] == "rule_inferred"
    assert "Mg(OH)₂" in res["equation_unicode"]
    assert any(o["type"] == "precipitate" for o in res["observations"])


def test_rule_inference_salt_salt_precipitation():
    res = predict_reaction(_inputs(["BaCl2", "K2SO4"]))
    assert res["supported"] is True
    assert "BaSO₄" in res["equation_unicode"]
    assert any(o["type"] == "precipitate" for o in res["observations"])


def test_rule_inference_combustion():
    res = predict_reaction(_inputs(["C6H6", "O2"]))
    assert res["supported"] is True
    assert "CO₂" in res["equation_unicode"]
    assert any(o["type"] == "heat" for o in res["observations"])


def test_rule_inference_alkali_metal_with_water():
    res = predict_reaction(_inputs(["Na", "H2O"]))
    assert res["supported"] is True
    assert "NaOH" in res["equation_unicode"]
    assert "H₂" in res["equation_unicode"]


def test_rule_inference_stoichiometry_and_limiting_reagent():
    res = predict_reaction(_inputs(["H2SO4", "Ca(OH)2"], amounts=[1.0, 0.5]))
    stoich = res["stoichiometry"]
    assert stoich["is_exact_stoichiometric_ratio"] is False
    assert stoich["limiting_reagent"]["formula"] == "Ca(OH)2"
    yields = {y["formula"]: y["theoretical_yield_mol"] for y in res["theoretical_yields"]}
    assert yields["CaSO4"] == pytest.approx(0.5)
    assert yields["H2O"] == pytest.approx(1.0)


def test_rule_inference_preserves_atom_conservation():
    res = predict_reaction(_inputs(["H2SO4", "Ca(OH)2"]))
    assert res["atom_conservation"]["is_balanced"] is True


# --------------------------------------------------------------------------------------
# Tier 3 — insufficient evidence (never invent products)
# --------------------------------------------------------------------------------------


def test_insufficient_evidence_when_no_rule_applies():
    res = predict_reaction(_inputs(["HCl", "HNO3"]))
    assert res["supported"] is False
    assert res["mode"] == "insufficient_evidence"
    assert res["products"] == []
    assert res["missing_information_ar"]


def test_insufficient_evidence_requires_at_least_two_materials():
    res = predict_reaction(_inputs(["HCl"]))
    assert res["supported"] is False
    assert res["mode"] == "insufficient_evidence"


def test_invalid_formula_is_rejected_without_invention():
    res = predict_reaction(_inputs(["HCl", "????"]))
    assert res["supported"] is False
    assert res["mode"] == "insufficient_evidence"
    assert res["products"] == []


# --------------------------------------------------------------------------------------
# AI assistant grounding
# --------------------------------------------------------------------------------------


def test_ai_assistant_explains_prediction_in_arabic():
    res = predict_reaction(_inputs(["H2SO4", "Ca(OH)2"]))
    text = AIScientificAssistant.explain_reaction_prediction(res)
    assert "مستوى الدليل" in text
    assert "المعادلة الموزونة" in text
    assert "حدود التنبؤ" in text
    assert "المرجع" in text


def test_ai_assistant_handles_insufficient_evidence():
    res = predict_reaction(_inputs(["HCl", "HNO3"]))
    text = AIScientificAssistant.explain_reaction_prediction(res)
    assert "البيانات الناقصة" in text

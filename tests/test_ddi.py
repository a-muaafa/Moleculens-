"""
Unit tests for the DDIEngine module in MolecuLens.
Verifies pairwise interaction evaluations, severity tiers, clinical disclaimers,
polypharmacy screening, and synonym normalization.
"""

from pharmacy.ddi_engine import DDIEngine


def test_aspirin_ibuprofen_interaction():
    res = DDIEngine.check_interaction("aspirin", "ibuprofen")
    assert res["success"] is True
    assert res["has_interaction"] is True
    assert res["severity"] == "Major"
    assert "COX-1" in res["mechanism_en"]
    assert len(res["citations"]) > 0
    assert len(res["disclaimer_ar"]) > 20
    assert len(res["disclaimer_en"]) > 20


def test_sitagliptin_metformin_interaction():
    res = DDIEngine.check_interaction("sitagliptin", "metformin")
    assert res["success"] is True
    assert res["has_interaction"] is True
    assert res["severity"] == "Minor"


def test_synonym_resolution():
    res1 = DDIEngine.check_interaction("tylenol", "aspirin")
    res2 = DDIEngine.check_interaction("paracetamol", "acetylsalicylic acid")
    assert res1["severity"] == res2["severity"]


def test_no_documented_interaction_has_disclaimer():
    res = DDIEngine.check_interaction("sitagliptin", "amoxicillin")
    assert res["success"] is True
    assert res["has_interaction"] is False
    assert res["severity"] == "No Interaction Documented"
    assert "disclaimer_ar" in res and "disclaimer_en" in res
    assert "عدم وجود تفاعل" in res["disclaimer_ar"]
    assert "does NOT prove safety" in res["disclaimer_en"]


def test_duplicate_therapy_check():
    res = DDIEngine.check_interaction("panadol", "paracetamol")
    assert res["success"] is True
    assert res.get("is_same_drug") is True
    assert res["severity"] == "Duplicate Therapy"


def test_polypharmacy_screening():
    drugs = ["aspirin", "ibuprofen", "warfarin"]
    res = DDIEngine.screen_polypharmacy(drugs)
    assert res["success"] is True
    assert res["total_drugs"] == 3
    assert res["total_pairs_evaluated"] == 3
    assert res["total_interactions_found"] >= 2
    assert res["severity_summary"]["Major"] >= 2


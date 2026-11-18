"""
Comprehensive tests for:
1. Diverse reaction scenarios (non-reacting pairs, unknown substances, invalid formulas,
   3-product reactions, special catalytic/thermal conditions, limiting/excess shifts,
   and zero fabricated products).
2. Security, input validation, XSS sanitization, resilient HTTP retry/timeout, and
   clean environment configuration (.env.example, Dockerfile, no hardcoded user paths).
3. Standardized Backend-Frontend API Contract (MolecuLensAPI status codes, envelope,
   canonical units, and UI states).
"""

from __future__ import annotations

import os
from pathlib import Path
from unittest.mock import MagicMock

import pytest
import requests

from core.api_contract import CANONICAL_UNITS, UI_STATE_MESSAGES, MolecuLensAPI
from core.security import (
    contains_unsafe_markup,
    is_safe_https_url,
    sanitize_external_record,
    sanitize_external_text,
    validate_safe_chemical_input,
)
from data_integrations.http_client import request_with_retry
from reaction_engine.prediction import predict_reaction
from reaction_engine.reaction_service import analyze_chemical_reaction

ROOT = Path(__file__).resolve().parents[1]


# ==============================================================================
# 1. Diverse Reaction Engine Tests (Item 4 & Item 3)
# ==============================================================================


def test_diverse_reactions_non_reacting_substances_produce_no_fabricated_products():
    # Two soluble salts that form no precipitate (NaCl + KNO3), noble metal + acid (Au + HCl), noble gas + water
    for pair in [("NaCl", "KNO3"), ("Au", "HCl"), ("He", "H2O")]:
        res = predict_reaction([
            {"substance": {"symbol_or_formula": pair[0]}},
            {"substance": {"symbol_or_formula": pair[1]}},
        ])
        assert res["supported"] is False
        assert res["mode"] == "insufficient_evidence"
        assert res["products"] == []


def test_diverse_reactions_unknown_substance_and_invalid_formula():
    # Unknown substance via API
    api_unknown = MolecuLensAPI.predict_reaction("UnobtainiumXYZ + HCl")
    assert api_unknown["ok"] is False
    assert api_unknown["status_code"] == 400
    assert api_unknown["error"]["code"] == "INSUFFICIENT_REACTANTS"

    # Invalid chemical formula passed directly to predict_reaction
    bad_formula_res = predict_reaction([
        {"substance": {"symbol_or_formula": "H2O???###"}},
        {"substance": {"symbol_or_formula": "NaOH"}},
    ])
    assert bad_formula_res["supported"] is False
    assert bad_formula_res["products"] == []


def test_diverse_reactions_three_products_and_charge_conservation():
    # Reaction producing 3 distinct products: CH3COOH + NaHCO3 -> CH3COONa + H2O + CO2
    res = predict_reaction([
        {"substance": {"symbol_or_formula": "CH3COOH"}, "amount": 1.0, "unit": "mol"},
        {"substance": {"symbol_or_formula": "NaHCO3"}, "amount": 1.0, "unit": "mol"},
    ])
    assert res["supported"] is True
    assert len(res["products"]) == 3
    prod_formulas = {p["formula"] for p in res["products"]}
    assert prod_formulas == {"CH3COONa", "H2O", "CO2"}
    assert res["atom_conservation"]["is_balanced"] is True


def test_diverse_reactions_special_conditions_and_limiting_reagent_shift():
    # Haber-Bosch N2 + 3H2 -> 2NH3 requires high T (380-550 C), high P (150-300 atm), and fe_haber catalyst
    rxn_bad_cond = analyze_chemical_reaction(
        reactant_inputs=[
            {"substance": {"symbol_or_formula": "N2"}, "amount": 1.0, "unit": "mol"},
            {"substance": {"symbol_or_formula": "H2"}, "amount": 2.0, "unit": "mol"},
        ],
        domain="general_chemistry",
        custom_conditions={
            "temperature": "25 °C",
            "pressure": "1 atm",
            "solvent": "aqueous",
            "catalyst": "none",
        },
    )
    assert rxn_bad_cond["supported"] is True
    compat = rxn_bad_cond["conditions_compatibility"]
    assert compat["is_compatible"] is False
    assert len(compat["warnings_ar"]) >= 2

    # Limiting reagent: 1 mol N2 needs 3 mol H2; with 2 mol H2, H2 is limiting (remaining H2 = 0, remaining N2 = 1/3 mol)
    assert rxn_bad_cond["limiting_reagent"]["formula"] == "H2"
    n2_row = next(r for r in rxn_bad_cond["reactants_table"] if r["formula"] == "N2")
    assert n2_row["remaining_mol"] > 0.3


# ==============================================================================
# 2. Security, Resilient External HTTP, and Deployment Hygiene (Items 6, 7, 8)
# ==============================================================================


def test_security_xss_and_markup_sanitization():
    malicious = "<script>alert('xss')</script><img src=x onerror=alert(1)>SafeTitle"
    assert contains_unsafe_markup(malicious) is True

    cleaned = sanitize_external_text(malicious)
    assert "<script" not in cleaned.lower()
    assert "onerror=" not in cleaned.lower()
    assert "SafeTitle" in cleaned

    nested = {"title": "<iframe src='javascript:alert(1)'></iframe>Aspirin", "synonyms": ["<b>ASA</b>"]}
    san_nested = sanitize_external_record(nested)
    assert "<iframe" not in san_nested["title"].lower()
    assert "&lt;b&gt;ASA&lt;/b&gt;" in san_nested["synonyms"][0]

    ok, _, err = validate_safe_chemical_input("<script>evil()</script>")
    assert ok is False
    assert "رُفض" in err

    assert is_safe_https_url("https://pubchem.ncbi.nlm.nih.gov/rest/pug") is True
    assert is_safe_https_url("http://pubchem.ncbi.nlm.nih.gov") is False
    assert is_safe_https_url("https://localhost:8501/secret") is False


def test_http_client_retry_on_transient_error(monkeypatch):
    calls = {"count": 0}

    def fake_get(url, params=None, headers=None, timeout=5):
        calls["count"] += 1
        if calls["count"] == 1:
            raise requests.exceptions.Timeout("Transient timeout")
        resp = MagicMock()
        resp.status_code = 200
        resp.json.return_value = {"ok": True}
        return resp

    monkeypatch.setattr(requests, "get", fake_get)
    resp = request_with_retry("https://pubchem.ncbi.nlm.nih.gov/rest/pug", max_retries=2, backoff_seconds=0.01)
    assert resp.status_code == 200
    assert calls["count"] == 2


def test_deployment_files_and_no_hardcoded_developer_paths():
    assert (ROOT / ".env.example").is_file()
    assert (ROOT / "Dockerfile").is_file()
    assert (ROOT / "requirements.txt").is_file()
    assert (ROOT / "docs" / "API_CONTRACT.md").is_file()
    assert (ROOT / "results" / "logp_eval.json").is_file()

    env_text = (ROOT / ".env.example").read_text(encoding="utf-8")
    assert "MOLECULENS_ENV" in env_text
    assert "EXTERNAL_API_TIMEOUT_SECONDS" in env_text

    # Ensure no Python source file contains hardcoded developer machine paths
    for py_file in ROOT.rglob("*.py"):
        if "venv" in py_file.parts or ".git" in py_file.parts or "tests" in py_file.parts:
            continue
        content = py_file.read_text(encoding="utf-8", errors="ignore")
        assert "C:\\Users\\" not in content and "c:\\Users\\" not in content, f"Hardcoded path in {py_file}"


# ==============================================================================
# 3. Unified Backend-Frontend API Contract (Item 9)
# ==============================================================================


def test_unified_api_contract_endpoints_and_status_codes():
    # 1. Explore endpoint: 200 OK, 400 Bad Request, 404 Not Found
    ok_exp = MolecuLensAPI.explore_entity("H2O", allow_network=False)
    assert ok_exp["ok"] is True
    assert ok_exp["status_code"] == 200
    assert ok_exp["ui_state"] == "ready"
    assert ok_exp["units"] == CANONICAL_UNITS
    assert ok_exp["meta"]["computed_in_backend"] is True
    assert ok_exp["data"]["entity"]["symbol_or_formula"] == "H2O"

    bad_exp = MolecuLensAPI.explore_entity("")
    assert bad_exp["ok"] is False
    assert bad_exp["status_code"] == 400
    assert bad_exp["ui_state"] == "error"
    assert bad_exp["error"]["code"] == "INVALID_SEARCH_QUERY"

    not_found_exp = MolecuLensAPI.explore_entity("NonExistentCompound999999", allow_network=False)
    assert not_found_exp["ok"] is False
    assert not_found_exp["status_code"] == 404
    assert not_found_exp["ui_state"] == "empty"

    # 2. Predict reaction endpoint: 200 OK (with visual stages), 422 Insufficient Evidence
    ok_pred = MolecuLensAPI.predict_reaction("CH3COOH + K2CO3")
    assert ok_pred["ok"] is True
    assert ok_pred["status_code"] == 200
    assert "visual_experiment_stages" in ok_pred["data"]

    insuff_pred = MolecuLensAPI.predict_reaction("Au + H2O")
    assert insuff_pred["ok"] is False
    assert insuff_pred["status_code"] == 422
    assert insuff_pred["error"]["code"] == "INSUFFICIENT_REACTION_EVIDENCE"
    assert insuff_pred["data"]["products"] == []

    # 3. Analyze reaction endpoint: 422 on L_STP with solid NaCl
    stp_err = MolecuLensAPI.analyze_reaction([
        {"substance": {"symbol_or_formula": "NaCl", "name_ar": "كلوريد الصوديوم"}, "amount": 5.0, "unit": "L_STP"}
    ])
    assert stp_err["ok"] is False
    assert stp_err["status_code"] == 422
    assert stp_err["error"]["code"] == "INVALID_STP_PHASE"

    # 4. LogP & Structure Modification endpoints
    ok_logp = MolecuLensAPI.evaluate_logp("CCO")
    assert ok_logp["ok"] is True
    assert ok_logp["status_code"] == 200
    assert ok_logp["data"]["qsar_model_evaluation"]["loco_zero_leakage_verified"] is True

    ok_mod = MolecuLensAPI.modify_structure("c1ccccc1O", "hydroxyl_to_chloro")
    assert ok_mod["ok"] is True
    assert ok_mod["status_code"] == 200
    assert ok_mod["data"]["pharmacological_efficacy_claimed"] is False

    # 5. Scientific Audit endpoint
    audit_res = MolecuLensAPI.scientific_audit()
    assert audit_res["ok"] is True
    assert audit_res["status_code"] == 200
    assert audit_res["data"]["manually_verified_references"] >= 30
    assert audit_res["data"]["total_reactions"] == 152

    # 6. Pharmacology Coming Soon endpoint & Titration API endpoints
    pharm_res = MolecuLensAPI.pharmacology_overview()
    assert pharm_res["ok"] is True
    assert pharm_res["data"]["status"] == "coming_soon"
    assert len(pharm_res["data"]["pillars"]) == 7

    tit_res = MolecuLensAPI.simulate_titration("H3PO4", "NaOH", 0.1, 25.0, 0.1)
    assert tit_res["ok"] is True
    assert tit_res["data"]["is_polyprotic"] is True


def test_experiment_report_pdf_and_ai_architecture_api_endpoints():
    from reaction_engine.report_generator import ExperimentReportGenerator

    rxn_pred = MolecuLensAPI.predict_reaction("HCl + NaOH")
    assert rxn_pred["ok"] is True

    rep_res = MolecuLensAPI.generate_experiment_report(rxn_pred["data"], report_type="reaction_experiment")
    assert rep_res["ok"] is True
    assert rep_res["status_code"] == 200
    assert rep_res["data"]["pdf_available"] is True
    assert rep_res["data"]["pdf_size_bytes"] > 200
    assert "<html" in rep_res["data"]["printable_html"]

    pdf_bytes = ExperimentReportGenerator.generate_experiment_pdf_bytes(rep_res["data"]["report"])
    assert pdf_bytes.startswith(b"%PDF-1.")

    tit_sim = MolecuLensAPI.simulate_titration("CH3COOH", "NaOH", 0.1, 25.0, 0.1)
    tit_rep = MolecuLensAPI.generate_experiment_report(tit_sim["data"], report_type="titration_experiment")
    assert tit_rep["ok"] is True
    assert tit_rep["data"]["pdf_available"] is True

    ai_mat = MolecuLensAPI.ai_architecture_matrix()
    assert ai_mat["ok"] is True
    assert ai_mat["status_code"] == 200
    assert len(ai_mat["data"]["architecture_matrix"]) >= 5


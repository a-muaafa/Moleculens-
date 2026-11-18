"""
Quality & security gates for MolecuLens.

Covers: secret scanning, XSS / input validation, SSRF-safe URL filter, API contract for every
endpoint (valid and invalid inputs), missing data, external-service failure, reaction mass and
charge conservation (including negative tests), titration, concentration, and report-generation
calculations.

Tests that need optional heavy libraries (plotly, streamlit, RDKit) skip cleanly when the library
is missing; in CI (requirements.txt installed) they all run.
"""

from __future__ import annotations

import copy
import re
from pathlib import Path

import pytest
import requests

from core.security import (
    contains_unsafe_markup,
    is_safe_https_url,
    sanitize_external_record,
    sanitize_external_text,
    validate_safe_chemical_input,
)
from core.validators import validate_search_query

ROOT = Path(__file__).resolve().parents[1]

XSS_PAYLOADS = [
    "<script>alert(1)</script>",
    "<img src=x onerror=alert(1)>",
    "<svg/onload=alert(1)>",
    "javascript:alert(1)",
    "<iframe src=javascript:alert(1)></iframe>",
    "\"><script>alert(document.cookie)</script>",
]

ENVELOPE_KEYS = {"ok", "status_code", "ui_state", "endpoint", "message_ar", "message_en", "units", "data", "error", "meta"}


# ==============================================================================
# 1. Secrets must never be committed
# ==============================================================================

_SECRET_PATTERNS = {
    "Anthropic/OpenAI style key": re.compile(r"sk-[A-Za-z0-9_\-]{20,}"),
    "Google API key": re.compile(r"AIza[0-9A-Za-z_\-]{30,}"),
    "GitHub token": re.compile(r"gh[pousr]_[A-Za-z0-9]{30,}"),
    "AWS access key id": re.compile(r"AKIA[0-9A-Z]{16}"),
    "Slack token": re.compile(r"xox[baprs]-[A-Za-z0-9\-]{10,}"),
    "Private key block": re.compile(r"-----BEGIN (?:RSA |EC |OPENSSH |DSA )?PRIVATE KEY-----"),
    "Hardcoded credential assignment": re.compile(
        r"(?i)\b(api[_-]?key|secret|passwd|password|access[_-]?token)\b\s*[:=]\s*['\"][A-Za-z0-9_\-\/+=]{12,}['\"]"
    ),
}
_SCAN_SUFFIXES = {".py", ".md", ".yml", ".yaml", ".toml", ".json", ".txt", ".cfg", ".ini", ".example", ".env", ".sh"}
_SKIP_DIRS = {".git", "__pycache__", "venv", ".venv", "node_modules", ".pytest_cache", "htmlcov"}


def _iter_repo_text_files():
    for path in ROOT.rglob("*"):
        if not path.is_file() or any(part in _SKIP_DIRS for part in path.parts):
            continue
        if path.name == "test_quality_gates.py":
            continue  # this file contains the detection patterns themselves
        if path.suffix.lower() in _SCAN_SUFFIXES or path.name.startswith(".env"):
            yield path


def test_repository_contains_no_api_keys_or_secrets():
    findings = []
    for path in _iter_repo_text_files():
        try:
            text = path.read_text(encoding="utf-8", errors="ignore")
        except OSError:
            continue
        for label, pattern in _SECRET_PATTERNS.items():
            if pattern.search(text):
                findings.append(f"{path.relative_to(ROOT)}: {label}")
    assert not findings, "Possible secrets committed:\n" + "\n".join(findings)


def test_env_files_are_gitignored_and_example_has_no_values_that_look_secret():
    gitignore = (ROOT / ".gitignore").read_text(encoding="utf-8")
    assert re.search(r"(?m)^\.env\s*$", gitignore), ".env must be listed in .gitignore"
    example = (ROOT / ".env.example").read_text(encoding="utf-8")
    for line in example.splitlines():
        if "=" in line and not line.lstrip().startswith("#"):
            key, _, value = line.partition("=")
            assert not re.search(r"(?i)key|secret|token|password", key) or value.strip() == "", line


# ==============================================================================
# 2. XSS and user-input validation
# ==============================================================================


@pytest.mark.parametrize("payload", XSS_PAYLOADS)
def test_xss_payloads_are_detected_rejected_and_neutralised(payload):
    assert contains_unsafe_markup(payload) is True

    ok, cleaned, err = validate_safe_chemical_input(payload)
    assert ok is False and cleaned == "" and err

    ok_q, cleaned_q, _ = validate_search_query(payload)
    assert ok_q is False and cleaned_q == ""

    safe = sanitize_external_text(payload)
    assert "<" not in safe and ">" not in safe
    assert "onerror=" not in safe.lower() and "onload=" not in safe.lower()

    record = sanitize_external_record({"title": payload, "items": [payload, {"deep": payload}]})
    assert "<" not in record["title"]
    assert "<" not in record["items"][0]
    assert "<" not in record["items"][1]["deep"]


@pytest.mark.parametrize("payload", XSS_PAYLOADS)
def test_unknown_reactant_text_is_not_echoed_as_raw_html(payload):
    from reaction_engine.prediction import parse_and_resolve_reactant_text

    _items, notes = parse_and_resolve_reactant_text(f"{payload} + HCl")
    assert all("<" not in note and ">" not in note.replace("»", "").replace("«", "") for note in notes), notes


@pytest.mark.parametrize(
    "bad",
    [None, 123, [], {}, "", "   ", "x" * 301, "<script>x</script>"],
)
def test_validate_safe_chemical_input_rejects_missing_or_malformed_values(bad):
    ok, cleaned, err = validate_safe_chemical_input(bad)
    assert ok is False
    assert cleaned == ""
    assert err


def test_validate_safe_chemical_input_strips_control_characters_but_keeps_valid_text():
    ok, cleaned, _ = validate_safe_chemical_input("C\x00C\x07O")
    assert ok is True
    assert cleaned == "CCO"


@pytest.mark.parametrize("bad", [None, 5, "", "   ", "y" * 251])
def test_validate_search_query_rejects_missing_or_oversized_input(bad):
    ok, cleaned, err = validate_search_query(bad)
    assert ok is False and cleaned == "" and err


@pytest.mark.parametrize(
    "url,expected",
    [
        ("https://pubchem.ncbi.nlm.nih.gov/rest/pug", True),
        ("https://www.ebi.ac.uk/chembl/api/data", True),
        ("http://pubchem.ncbi.nlm.nih.gov", False),  # not https
        ("ftp://example.com/file", False),
        ("https://localhost:8501/secret", False),
        ("https://127.0.0.1/admin", False),
        ("https://0.0.0.0/", False),
        ("https://[::1]/", False),
        ("https://[::ffff:127.0.0.1]/", False),  # IPv4-mapped loopback
        ("https://10.0.0.5/internal", False),  # RFC1918
        ("https://192.168.1.1/router", False),
        ("https://172.16.0.9/", False),
        ("https://169.254.169.254/latest/meta-data", False),  # cloud metadata endpoint
        ("https://service.internal/x", False),
        ("https://printer.local/", False),
        ("https://trusted.example.com@evil.example.net/", False),  # credential-in-URL trick
        ("https://", False),
        ("", False),
        (None, False),
    ],
)
def test_https_url_filter_blocks_non_https_and_internal_targets(url, expected):
    assert is_safe_https_url(url) is expected


# ==============================================================================
# 3. API contract: every endpoint, valid and invalid inputs
# ==============================================================================


def _assert_envelope(res, *, allowed_status):
    assert ENVELOPE_KEYS <= set(res), set(res)
    assert res["status_code"] in allowed_status, (res["endpoint"], res["status_code"], res.get("message_ar"))
    assert res["ok"] is (res["status_code"] == 200)
    assert res["meta"]["computed_in_backend"] is True
    if res["ok"]:
        assert res["error"] is None
    else:
        assert res["error"]["code"]


@pytest.fixture(scope="module")
def api():
    pytest.importorskip("plotly")
    from core.api_contract import MolecuLensAPI

    return MolecuLensAPI


@pytest.mark.parametrize("query", ["", "   ", "<script>alert(1)</script>", "x" * 400])
def test_api_explore_rejects_invalid_queries_with_400(api, query):
    res = api.explore_entity(query, allow_network=False)
    _assert_envelope(res, allowed_status={400})
    assert "<script" not in res["message_ar"].lower()


def test_api_explore_unknown_entity_returns_404_without_network(api):
    res = api.explore_entity("zzzqqq-not-a-real-substance-123", allow_network=False)
    _assert_envelope(res, allowed_status={404})
    assert res["ui_state"] == "empty"


def test_api_predict_reaction_valid_and_invalid_inputs(api):
    good = api.predict_reaction("HCl + NaOH")
    _assert_envelope(good, allowed_status={200})
    assert good["data"]["atom_conservation"]["is_balanced"] is True

    for bad_input in ["", "<script>x</script> + HCl", "HCl", 12345, None]:
        res = api.predict_reaction(bad_input)
        _assert_envelope(res, allowed_status={400})

    unsupported = api.predict_reaction([
        {"substance": {"symbol_or_formula": "He"}},
        {"substance": {"symbol_or_formula": "H2O"}},
    ])
    _assert_envelope(unsupported, allowed_status={422})
    assert unsupported["data"]["products"] == []


def test_api_analyze_reaction_valid_and_invalid_inputs(api):
    good = api.analyze_reaction([
        {"substance": {"symbol_or_formula": "HCl"}, "amount": 1.0, "unit": "mol"},
        {"substance": {"symbol_or_formula": "NaOH"}, "amount": 1.0, "unit": "mol"},
    ])
    _assert_envelope(good, allowed_status={200})

    empty = api.analyze_reaction([])
    _assert_envelope(empty, allowed_status={400, 422})


def test_api_titration_valid_and_invalid_parameters(api):
    good = api.simulate_titration("HCl", "NaOH", acid_conc_m=0.1, acid_volume_ml=25.0, base_conc_m=0.1)
    _assert_envelope(good, allowed_status={200})
    assert good["data"]["v_eq_ml"] == pytest.approx(25.0, abs=0.01)

    for kwargs in (
        {"acid_key": "NotAnAcid", "base_key": "NaOH"},
        {"acid_key": "HCl", "base_key": "NotABase"},
        {"acid_key": "HCl", "base_key": "NaOH", "acid_conc_m": -0.1},
        {"acid_key": "HCl", "base_key": "NaOH", "acid_volume_ml": 0},
        {"acid_key": "HCl", "base_key": "NaOH", "base_conc_m": 0},
    ):
        res = api.simulate_titration(**kwargs)
        _assert_envelope(res, allowed_status={400})


def test_api_scientific_audit_endpoint(api):
    res = api.scientific_audit()
    _assert_envelope(res, allowed_status={200})
    assert res["data"]["overall_passed"] is True
    assert res["data"]["total_reactions"] == res["data"]["passed_reactions"]


def test_api_logp_and_modify_endpoints_with_valid_and_invalid_smiles(api):
    pytest.importorskip("rdkit")
    _assert_envelope(api.evaluate_logp("CCO"), allowed_status={200})
    for bad in ["", "<script>", "not_a_smiles!!", None]:
        _assert_envelope(api.evaluate_logp(bad), allowed_status={400})

    _assert_envelope(api.modify_structure("c1ccccc1O", "hydroxyl_to_chloro"), allowed_status={200})
    _assert_envelope(api.modify_structure("CCO", "no_such_modification"), allowed_status={400, 422})
    _assert_envelope(api.modify_structure("", "hydroxyl_to_chloro"), allowed_status={400})


def test_api_pharmacology_overview_is_marked_coming_soon(api):
    pytest.importorskip("streamlit")
    res = api.pharmacology_overview()
    _assert_envelope(res, allowed_status={200})
    assert res["data"]["status"] == "coming_soon"


# ==============================================================================
# 4. Missing data and wrong inputs in the scientific engines
# ==============================================================================


def test_predict_reaction_with_missing_or_malformed_reactants_never_fabricates_products():
    from reaction_engine.prediction import predict_reaction

    for bad in ([], [{}], [{"substance": {}}], [{"substance": {"symbol_or_formula": ""}}, {"substance": {"symbol_or_formula": "HCl"}}]):
        res = predict_reaction(bad)
        assert res["supported"] is False
        assert res["products"] == []


def test_concentration_calculator_rejects_missing_and_invalid_values():
    from calculators.concentration_calculator import ConcentrationCalculator as C

    assert C.calculate_molarity()["success"] is False
    assert C.calculate_molarity(mass_g=-1, molar_mass_g_mol=58.44, volume_l=1)["success"] is False
    assert C.calculate_molarity(mass_g=1, molar_mass_g_mol=0, volume_l=1)["success"] is False
    assert C.calculate_dilution(c1=1, v1=1, c2=1, v2=1)["success"] is False  # over-specified
    assert C.calculate_dilution(c1=1, v1=1)["success"] is False  # under-specified


def test_titration_engine_rejects_invalid_parameters():
    pytest.importorskip("plotly")
    from titration.titration_engine import TitrationEngine

    for kwargs in (
        {"acid_key": "Nope", "base_key": "NaOH", "acid_conc_m": 0.1, "acid_volume_ml": 25.0, "base_conc_m": 0.1},
        {"acid_key": "HCl", "base_key": "Nope", "acid_conc_m": 0.1, "acid_volume_ml": 25.0, "base_conc_m": 0.1},
        {"acid_key": "HCl", "base_key": "NaOH", "acid_conc_m": 0.0, "acid_volume_ml": 25.0, "base_conc_m": 0.1},
        {"acid_key": "HCl", "base_key": "NaOH", "acid_conc_m": 0.1, "acid_volume_ml": -5.0, "base_conc_m": 0.1},
    ):
        with pytest.raises(ValueError):
            TitrationEngine.simulate_titration(**kwargs)


# ==============================================================================
# 5. External service failure must degrade gracefully
# ==============================================================================


def test_pubchem_connection_failure_returns_error_dict_instead_of_raising(monkeypatch):
    from data_integrations.pubchem_client import fetch_pubchem_by_name

    def boom(*_a, **_k):
        raise requests.exceptions.ConnectionError("network down")

    monkeypatch.setattr(requests, "get", boom)
    monkeypatch.setattr("time.sleep", lambda *_a, **_k: None)
    res = fetch_pubchem_by_name("zz-unique-offline-query-001")
    assert res["ok"] is False
    assert res["records"] == []


def test_pubchem_timeout_returns_error_dict_instead_of_raising(monkeypatch):
    from data_integrations.pubchem_client import fetch_pubchem_by_name

    def slow(*_a, **_k):
        raise requests.exceptions.Timeout("too slow")

    monkeypatch.setattr(requests, "get", slow)
    monkeypatch.setattr("time.sleep", lambda *_a, **_k: None)
    res = fetch_pubchem_by_name("zz-unique-timeout-query-002")
    assert res["ok"] is False


def test_search_falls_back_cleanly_when_network_is_unavailable(monkeypatch):
    from search.search_engine import search_scientific_entity

    def boom(*_a, **_k):
        raise requests.exceptions.ConnectionError("network down")

    monkeypatch.setattr(requests, "get", boom)
    monkeypatch.setattr("time.sleep", lambda *_a, **_k: None)
    res = search_scientific_entity("zz-unique-offline-search-003", search_type="auto", allow_network=True)
    assert res["found"] is False
    assert res.get("message")


def test_http_client_retries_on_503_then_returns_and_rejects_unsafe_urls(monkeypatch):
    from data_integrations.http_client import request_with_retry

    calls = {"n": 0}

    class _Resp:
        def __init__(self, code):
            self.status_code = code

    def flaky(url, params=None, headers=None, timeout=5):
        calls["n"] += 1
        return _Resp(503 if calls["n"] < 3 else 200)

    monkeypatch.setattr(requests, "get", flaky)
    monkeypatch.setattr("time.sleep", lambda *_a, **_k: None)
    resp = request_with_retry("https://pubchem.ncbi.nlm.nih.gov/rest/pug", max_retries=3, backoff_seconds=0)
    assert resp.status_code == 200 and calls["n"] == 3

    with pytest.raises(ValueError):
        request_with_retry("https://169.254.169.254/latest/meta-data")


def test_http_client_raises_after_all_retries_exhausted(monkeypatch):
    from data_integrations.http_client import request_with_retry

    def boom(*_a, **_k):
        raise requests.exceptions.ConnectionError("down")

    monkeypatch.setattr(requests, "get", boom)
    monkeypatch.setattr("time.sleep", lambda *_a, **_k: None)
    with pytest.raises(requests.exceptions.ConnectionError):
        request_with_retry("https://pubchem.ncbi.nlm.nih.gov/rest/pug", max_retries=2, backoff_seconds=0)


# ==============================================================================
# 6. Reactions: mass and charge conservation (positive and negative tests)
# ==============================================================================


def test_every_registry_reaction_passes_mass_charge_and_atom_audit():
    from reaction_engine.reaction_registry import VERIFIED_REACTIONS_REGISTRY
    from reaction_engine.scientific_audit import audit_single_reaction

    failures = []
    for rxn in VERIFIED_REACTIONS_REGISTRY:
        report = audit_single_reaction(rxn)
        for check in ("atom_and_mass_conservation", "charge_conservation"):
            if not report["checks"][check]:
                failures.append((rxn["reaction_id"], check))
    assert not failures, failures


def test_every_registry_reaction_conserves_mass_independently_of_the_audit_module():
    from chemistry_engine.formula_analyzer import calculate_molar_mass
    from reaction_engine.reaction_registry import VERIFIED_REACTIONS_REGISTRY

    for rxn in VERIFIED_REACTIONS_REGISTRY:
        def side(species):
            return sum(int(s["coefficient"]) * calculate_molar_mass(s["formula"]) for s in species)

        assert side(rxn["reactants"]) == pytest.approx(side(rxn["products"]), abs=0.05), rxn["reaction_id"]


def test_audit_detects_broken_atom_balance_and_broken_charge_balance():
    from reaction_engine.reaction_registry import VERIFIED_REACTIONS_REGISTRY
    from reaction_engine.scientific_audit import audit_single_reaction

    base = VERIFIED_REACTIONS_REGISTRY[0]
    assert audit_single_reaction(copy.deepcopy(base))["passed"] is True

    wrong_coeff = copy.deepcopy(base)
    wrong_coeff["reactants"][0]["coefficient"] += 1
    broken = audit_single_reaction(wrong_coeff)
    assert broken["passed"] is False and broken["checks"]["atom_and_mass_conservation"] is False

    wrong_charge = copy.deepcopy(base)
    wrong_charge["products"][0]["charge"] = 1
    broken_q = audit_single_reaction(wrong_charge)
    assert broken_q["passed"] is False and broken_q["checks"]["charge_conservation"] is False


def test_manual_verification_count_is_consistent_with_registry():
    from reaction_engine.reaction_registry import VERIFIED_REACTIONS_REGISTRY
    from reaction_engine.scientific_audit import MANUALLY_VERIFIED_REACTION_SOURCES, audit_registry_and_rules

    registry_ids = {r["reaction_id"] for r in VERIFIED_REACTIONS_REGISTRY}
    assert set(MANUALLY_VERIFIED_REACTION_SOURCES) <= registry_ids, "manual-verification entry for an unknown reaction id"
    report = audit_registry_and_rules()
    assert report["manually_verified_references"] == len(MANUALLY_VERIFIED_REACTION_SOURCES)
    assert str(report["total_reactions"]) in report["reference_check_scope_ar"]


# ==============================================================================
# 7. Titration, concentration and pharmacology calculations
# ==============================================================================


@pytest.mark.parametrize(
    "c_acid,v_acid,c_base",
    [(0.1, 25.0, 0.1), (0.05, 50.0, 0.1), (0.2, 10.0, 0.1), (0.1, 20.0, 0.25)],
)
def test_strong_acid_strong_base_equivalence_volume_follows_c1v1_equals_c2v2(c_acid, v_acid, c_base):
    pytest.importorskip("plotly")
    from titration.titration_engine import TitrationEngine

    res = TitrationEngine.simulate_titration(
        acid_key="HCl", base_key="NaOH", acid_conc_m=c_acid, acid_volume_ml=v_acid, base_conc_m=c_base
    )
    assert res["success"] is True
    assert res["v_eq_ml"] == pytest.approx(c_acid * v_acid / c_base, abs=0.05)
    assert res["ph_eq"] == pytest.approx(7.0, abs=0.1)


def test_weak_acid_strong_base_equivalence_point_is_basic():
    pytest.importorskip("plotly")
    from titration.titration_engine import TitrationEngine

    res = TitrationEngine.simulate_titration(
        acid_key="CH3COOH", base_key="NaOH", acid_conc_m=0.1, acid_volume_ml=25.0, base_conc_m=0.1
    )
    assert res["v_eq_ml"] == pytest.approx(25.0, abs=0.05)
    assert res["ph_eq"] > 8.0
    assert 1.0 < res["initial_ph"] < 4.0


def test_concentration_calculator_known_values():
    from calculators.concentration_calculator import ConcentrationCalculator as C

    mol = C.calculate_molarity(mass_g=5.844, molar_mass_g_mol=58.44, volume_l=1.0)
    assert mol["success"] is True
    assert mol["molarity_m"] == pytest.approx(0.1, rel=1e-3)

    dil = C.calculate_dilution(c1=1.0, v1=10.0, c2=0.1)
    assert dil["success"] is True
    assert dil["v2"] == pytest.approx(100.0)
    assert dil["solvent_to_add"] == pytest.approx(90.0)


# ==============================================================================
# 8. Regression tests for consistency fixes
# ==============================================================================


def _equation_key_from_registry(rxn):
    return (
        sorted((s["formula"], int(s["coefficient"])) for s in rxn["reactants"]),
        sorted((s["formula"], int(s["coefficient"])) for s in rxn["products"]),
    )


def _equation_key_from_text(equation: str):
    def side(text):
        items = []
        for token in text.split(" + "):
            m = re.match(r"\s*(\d*)\s*(\S+)", token)
            items.append((re.sub(r"\((?:s|l|g|aq)\)$", "", m.group(2)), int(m.group(1) or 1)))
        return sorted(items)

    left, right = equation.split("->")
    return side(left), side(right)


def test_every_manual_verification_entry_matches_its_registry_reaction_exactly():
    """A 'manually verified' badge must never be attached to a different reaction (e.g. after renumbering)."""
    from reaction_engine.reaction_registry import VERIFIED_REACTIONS_REGISTRY
    from reaction_engine.scientific_audit import MANUALLY_VERIFIED_REACTION_SOURCES

    by_id = {r["reaction_id"]: r for r in VERIFIED_REACTIONS_REGISTRY}
    mismatches = []
    for rid, meta in MANUALLY_VERIFIED_REACTION_SOURCES.items():
        assert rid in by_id, f"unknown reaction id {rid}"
        if _equation_key_from_registry(by_id[rid]) != _equation_key_from_text(meta["verified_equation"]):
            mismatches.append((rid, meta["verified_equation"]))
    assert not mismatches, mismatches


def test_no_runtime_database_or_env_file_is_tracked_by_gitignore_rules():
    gitignore = (ROOT / ".gitignore").read_text(encoding="utf-8")
    assert "data/*.db" in gitignore
    assert re.search(r"(?m)^\.env\s*$", gitignore)


def test_readme_does_not_hardcode_a_passing_test_count_or_static_success_badge():
    readme = (ROOT / "README.md").read_text(encoding="utf-8")
    assert "img.shields.io/badge/tests-" not in readme
    assert not re.search(r"\b\d+\s*(?:passed|اختباراً ناجحاً)", readme)


def test_no_ddi_references_remain_in_repository():
    """Verifies that Drug-Drug Interactions (DDI) has been completely removed from code and docs."""
    ddi_pattern = re.compile(r"(\bddi\b|drug[-_ ]drug|التداخلات الدوائية)", re.IGNORECASE)
    offenders = []
    for path in _iter_repo_text_files():
        text = path.read_text(encoding="utf-8", errors="ignore")
        if ddi_pattern.search(text):
            offenders.append(str(path.relative_to(ROOT)))
    assert not offenders, f"Found remaining DDI mentions in: {offenders}"


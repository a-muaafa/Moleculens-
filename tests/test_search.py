from unittest.mock import MagicMock, patch

from search.autocomplete import (
    build_instant_search_catalog,
    format_suggestion_label,
    get_autocomplete_suggestions,
)
from search.compound_search import (
    SEARCH_ONLINE_STATE_KEY,
    SEARCH_QUERY_STATE_KEY,
    SEARCH_TYPE_STATE_KEY,
    _compute_fuzzy_name_score,
    clear_pubchem_cache,
    query_pubchem_by_formula,
    search_compound,
    update_search_session_state,
    verify_compound_identity,
)
from search.local_dictionary import (
    LOCAL_COMPOUND_DICTIONARY,
    formulas_match,
    is_chemical_formula_query,
    normalize_chemical_formula,
)


def test_formula_co2_search():
    for q in ["CO2", "co2", "CO₂"]:
        res = search_compound(q, allow_network=False)
        assert res["found"] is True
        assert res["is_ambiguous"] is False
        assert res["selected"] is not None
        assert res["selected"]["name_ar"] == "ثاني أكسيد الكربون"
        assert res["selected"]["name_en"] == "Carbon Dioxide"
        assert res["selected"]["molecular_formula"] == "CO2"
        assert res["selected"]["matched_by"] == "exact_formula"
        assert res["selected"]["match_category"] == "Exact Match"
        assert "CO2" in res["selected"]["match_reason_ar"]

        ver = verify_compound_identity(res["selected"])
        assert ver["is_verified"] is True
        assert ver["compound"]["rdkit_formula"] == "CO2"


def test_formula_h2o_o2_h2o2_exact_separation():
    cases = [
        ("H2O", "الماء", "Water", "H2O"),
        ("h2o", "الماء", "Water", "H2O"),
        ("H₂O", "الماء", "Water", "H2O"),
        ("O2", "الأكسجين الجزيئي", "Molecular Oxygen", "O2"),
        ("O₂", "الأكسجين الجزيئي", "Molecular Oxygen", "O2"),
        ("H2O2", "فوق أكسيد الهيدروجين", "Hydrogen Peroxide", "H2O2"),
        ("H₂O₂", "فوق أكسيد الهيدروجين", "Hydrogen Peroxide", "H2O2"),
    ]
    for query, exp_ar, exp_en, exp_formula in cases:
        res = search_compound(query, allow_network=False)
        assert res["found"] is True
        assert res["is_ambiguous"] is False
        assert len(res["candidates"]) == 1
        assert res["selected"]["name_ar"] == exp_ar
        assert res["selected"]["name_en"] == exp_en
        assert res["selected"]["molecular_formula"] == exp_formula
        assert res["selected"]["matched_by"] == "exact_formula"


def test_formula_c6h12o6_and_ambiguous_queries_never_pick_random():
    res_hex = search_compound("C6H12O6", allow_network=False)
    assert res_hex["found"] is True
    assert res_hex["is_ambiguous"] is True
    assert res_hex["selected"] is None
    names_en = {c["name_en"] for c in res_hex["candidates"]}
    assert {"D-Glucose", "D-Fructose", "D-Galactose"}.issubset(names_en)
    for c in res_hex["candidates"]:
        assert c["match_reason_ar"]

    res_c2h6o = search_compound("C2H6O", allow_network=False)
    assert res_c2h6o["found"] is True
    assert res_c2h6o["is_ambiguous"] is True
    assert res_c2h6o["selected"] is None
    isomer_names = {c["name_en"] for c in res_c2h6o["candidates"]}
    assert "Ethanol" in isomer_names and "Dimethyl Ether" in isomer_names

    res_o = search_compound("O", allow_network=False)
    assert res_o["found"] is True
    assert res_o["is_ambiguous"] is True
    assert res_o["selected"] is None
    assert len(res_o["candidates"]) >= 2


def test_no_fuzzy_matching_on_chemical_formulas():
    unrelated_formulas = ["H2S", "SO2", "NO2", "C2H4O2", "C3H8O", "H3PO4", "C7H6O2"]
    for formula_q in unrelated_formulas:
        assert is_chemical_formula_query(formula_q) is True
        res = search_compound(formula_q, allow_network=False)
        assert res["found"] is False, f"Formula {formula_q} should not fuzzy-match any unrelated compound!"
        assert res["candidates"] == []
        assert res["selected"] is None

        for entry in LOCAL_COMPOUND_DICTIONARY:
            score, _, _ = _compute_fuzzy_name_score(formula_q.lower(), entry, search_type="auto")
            assert score == 0.0
            score_f, _, _ = _compute_fuzzy_name_score(formula_q.lower(), entry, search_type="formula")
            assert score_f == 0.0

        suggs = get_autocomplete_suggestions(formula_q, search_type="auto")
        assert all(s["matched_by"] != "approximate_name" for s in suggs)


def test_search_arabic_and_english_names_and_iupac():
    res_ar = search_compound("ثاني أكسيد الكربون", allow_network=False)
    assert res_ar["found"] is True
    assert res_ar["selected"]["molecular_formula"] == "CO2"
    assert res_ar["selected"]["matched_by"] == "exact_name_ar"

    res_water = search_compound("Water", allow_network=False)
    assert res_water["found"] is True
    assert res_water["selected"]["name_ar"] == "الماء"
    assert res_water["selected"]["matched_by"] == "exact_name_en"

    res_h2o2 = search_compound("Hydrogen Peroxide", allow_network=False)
    assert res_h2o2["found"] is True
    assert res_h2o2["selected"]["name_ar"] == "فوق أكسيد الهيدروجين"

    res_iupac = search_compound("oxidane", allow_network=False)
    assert res_iupac["found"] is True
    assert res_iupac["selected"]["name_en"] == "Water"
    assert res_iupac["selected"]["matched_by"] == "exact_iupac"


def test_search_abbreviations_synonyms_and_brands():
    res_atp = search_compound("ATP", allow_network=False)
    assert res_atp["found"] is True
    assert res_atp["selected"]["name_en"] == "Adenosine Triphosphate"
    assert res_atp["selected"]["matched_by"] == "exact_abbreviation"

    res_nadh = search_compound("NADH", allow_network=False)
    assert res_nadh["found"] is True
    assert res_nadh["selected"]["name_en"] == "Nicotinamide Adenine Dinucleotide"

    res_syn = search_compound("ماء الأكسجين", allow_network=False)
    assert res_syn["found"] is True
    assert res_syn["selected"]["molecular_formula"] == "H2O2"
    assert res_syn["selected"]["matched_by"] == "exact_synonym_ar"

    res_brand = search_compound("Januvia", allow_network=False)
    assert res_brand["found"] is True
    assert res_brand["selected"]["name_en"] == "Sitagliptin"
    assert res_brand["selected"]["matched_by"] == "exact_brand"


def test_search_fuzzy_and_typo_correction_for_names():
    res = search_compound("Sitaglptin", allow_network=False)
    assert res["found"] is True
    assert res["match_type"] == "approximate"
    assert len(res["candidates"]) >= 1
    assert any(c["name_en"] == "Sitagliptin" for c in res["candidates"])
    assert "تشابه" in res["candidates"][0]["match_reason_ar"]


def test_pubchem_formula_search_and_caching():
    clear_pubchem_cache()
    mock_cids_resp = MagicMock()
    mock_cids_resp.status_code = 200
    mock_cids_resp.json.return_value = {"IdentifierList": {"CID": [1118]}}

    mock_props_resp = MagicMock()
    mock_props_resp.status_code = 200
    mock_props_resp.json.return_value = {
        "PropertyTable": {
            "Properties": [
                {
                    "CID": 1118,
                    "Title": "Sulfuric Acid",
                    "MolecularFormula": "H2O4S",
                    "IUPACName": "sulfuric acid",
                    "CanonicalSMILES": "OS(=O)(=O)O",
                }
            ]
        }
    }

    with patch("search.compound_search.requests.get", side_effect=[mock_cids_resp, mock_props_resp]) as mock_get:
        res = search_compound("H2SO4", allow_network=True)
        assert res["found"] is True
        assert res["is_ambiguous"] is False
        assert res["selected"]["name_en"] == "Sulfuric Acid"
        assert res["selected"]["pubchem_cid"] == 1118
        assert res["data_source"] == "pubchem_online"
        assert mock_get.call_count == 2

        res_cached = query_pubchem_by_formula("H2SO4")
        assert res_cached["ok"] is True
        assert res_cached["cache_hit"] is True
        assert mock_get.call_count == 2

    clear_pubchem_cache()


def test_pubchem_failure_and_offline_fallback():
    clear_pubchem_cache()
    with patch("search.compound_search.requests.get", side_effect=ConnectionError("PubChem unreachable")):
        res_exact = search_compound("فيلداجليبتين", allow_network=True)
        assert res_exact["found"] is True
        assert res_exact["selected"]["name_en"] == "Vildagliptin"

        res_fuzzy = search_compound("Vildaglptin", allow_network=True)
        assert res_fuzzy["found"] is True
        assert "المسار البديل" in res_fuzzy["connection_status"]


def test_autocomplete_interactive_fields_and_typo_rules():
    suggs_ar = get_autocomplete_suggestions("سيتا")
    assert len(suggs_ar) >= 1
    assert suggs_ar[0]["name_en"] == "Sitagliptin"
    assert "display_label" in suggs_ar[0]
    assert "match_reason_ar" in suggs_ar[0]

    suggs_co2 = get_autocomplete_suggestions("CO2")
    assert len(suggs_co2) >= 1
    assert suggs_co2[0]["molecular_formula"] == "CO2"
    assert suggs_co2[0]["name_ar"] == "ثاني أكسيد الكربون"

    suggs_typo = get_autocomplete_suggestions("Sitaglptin")
    assert len(suggs_typo) >= 1
    assert suggs_typo[0]["name_en"] == "Sitagliptin"
    assert suggs_typo[0]["match_category"] == "Did You Mean"

    catalog = build_instant_search_catalog()
    assert len(catalog) == len(LOCAL_COMPOUND_DICTIONARY)
    assert any("CO2" in item["instant_label"] for item in catalog)


def test_compound_identity_verification_strict_checks():
    res_co2 = search_compound("CO2", allow_network=False)
    ver_ok = verify_compound_identity(res_co2["selected"])
    assert ver_ok["is_verified"] is True
    assert all(ver_ok["verification_steps"].values())

    mismatched_formula_entry = dict(res_co2["selected"])
    mismatched_formula_entry["molecular_formula"] = "H2O"
    ver_bad_formula = verify_compound_identity(mismatched_formula_entry)
    assert ver_bad_formula["is_verified"] is False
    assert ver_bad_formula["verification_steps"]["formula_verified"] is False

    missing_ref_entry = dict(res_co2["selected"])
    missing_ref_entry["pubchem_cid"] = None
    missing_ref_entry["chembl_id"] = None
    ver_no_ref = verify_compound_identity(missing_ref_entry)
    assert ver_no_ref["is_verified"] is False
    assert ver_no_ref["verification_steps"]["has_reference_id"] is False

    invalid_smiles_entry = dict(res_co2["selected"])
    invalid_smiles_entry["canonical_smiles"] = "INVALID_SMILES_###"
    ver_bad_smi = verify_compound_identity(invalid_smiles_entry)
    assert ver_bad_smi["is_verified"] is False
    assert ver_bad_smi["verification_steps"]["smiles_valid"] is False


def test_search_ui_session_state_preserves_query_when_changing_type():
    state = {}
    update_search_session_state(state, default_query="سيتاجليبتين")
    assert state[SEARCH_QUERY_STATE_KEY] == "سيتاجليبتين"
    assert state[SEARCH_TYPE_STATE_KEY] == "auto"
    assert state[SEARCH_ONLINE_STATE_KEY] is True

    update_search_session_state(state, new_query="CO2")
    assert state[SEARCH_QUERY_STATE_KEY] == "CO2"

    for new_type in ["formula", "name_ar", "name_en", "iupac", "abbreviation", "brand", "auto"]:
        update_search_session_state(state, new_search_type=new_type)
        assert state[SEARCH_QUERY_STATE_KEY] == "CO2"
        assert state[SEARCH_TYPE_STATE_KEY] == new_type

    update_search_session_state(state, new_allow_online=False)
    assert state[SEARCH_QUERY_STATE_KEY] == "CO2"
    assert state[SEARCH_ONLINE_STATE_KEY] is False


def test_formula_normalization_helpers():
    assert normalize_chemical_formula("CO2") == "CO2"
    assert normalize_chemical_formula("co2") == "CO2"
    assert normalize_chemical_formula("H₂O") == "H2O"
    assert normalize_chemical_formula("C₆H₁₂O₆") == "C6H12O6"
    assert formulas_match("NH3", "H3N") is True
    assert formulas_match("CO2", "H2O") is False

from unittest.mock import patch
from search.autocomplete import get_autocomplete_suggestions
from search.compound_search import search_compound, verify_compound_identity


def test_search_arabic_name():
    res = search_compound("الماء", allow_network=False)
    assert res["found"] is True
    assert res["is_ambiguous"] is False
    assert res["selected"]["name_en"] == "Water"
    assert res["selected"]["molecular_formula"] == "H2O"


def test_search_english_name():
    res = search_compound("Water", allow_network=False)
    assert res["found"] is True
    assert res["selected"]["name_ar"] == "الماء"

    res_h2o2 = search_compound("Hydrogen Peroxide", allow_network=False)
    assert res_h2o2["found"] is True
    assert res_h2o2["selected"]["name_ar"] == "فوق أكسيد الهيدروجين"


def test_search_iupac_name():
    res = search_compound("oxidane", allow_network=False)
    assert res["found"] is True
    assert res["selected"]["name_en"] == "Water"


def test_search_chemical_formula():
    for formula, expected_ar in [("H2O", "الماء"), ("O2", "الأكسجين الجزيئي"), ("H2O2", "فوق أكسيد الهيدروجين")]:
        res = search_compound(formula, allow_network=False)
        assert res["found"] is True
        assert res["selected"]["name_ar"] == expected_ar


def test_search_abbreviations():
    res_atp = search_compound("ATP", allow_network=False)
    assert res_atp["found"] is True
    assert res_atp["selected"]["name_en"] == "Adenosine Triphosphate"

    res_nadh = search_compound("NADH", allow_network=False)
    assert res_nadh["found"] is True
    assert res_nadh["selected"]["name_en"] == "Nicotinamide Adenine Dinucleotide"


def test_search_synonyms_and_brands():
    res_syn = search_compound("ماء الأكسجين", allow_network=False)
    assert res_syn["found"] is True
    assert res_syn["selected"]["molecular_formula"] == "H2O2"

    res_brand = search_compound("Januvia", allow_network=False)
    assert res_brand["found"] is True
    assert res_brand["selected"]["name_en"] == "Sitagliptin"


def test_search_fuzzy_and_typo_correction():
    res = search_compound("Sitaglptin", allow_network=False)
    assert res["found"] is True
    assert res["match_type"] == "approximate"
    assert len(res["candidates"]) >= 1
    assert any(c["name_en"] == "Sitagliptin" for c in res["candidates"])


def test_search_ambiguous_queries():
    res_o = search_compound("O", allow_network=False)
    assert res_o["found"] is True
    assert res_o["is_ambiguous"] is True
    assert res_o["selected"] is None
    assert len(res_o["candidates"]) >= 2

    res_hex = search_compound("C6H12O6", allow_network=False)
    assert res_hex["found"] is True
    assert res_hex["is_ambiguous"] is True
    assert res_hex["selected"] is None
    names = {c["name_en"] for c in res_hex["candidates"]}
    assert "D-Glucose" in names and "D-Fructose" in names


def test_search_nonexistent_compound():
    res = search_compound("ZZZZ_UNKNOWN_CHEMICAL_9999", allow_network=False)
    assert res["found"] is False
    assert res["match_type"] == "none"
    assert res["selected"] is None


def test_search_fallback_when_pubchem_fails():
    with patch("search.compound_search.requests.get", side_effect=ConnectionError("PubChem unreachable")):
        res_exact = search_compound("فيلداجليبتين", allow_network=True)
        assert res_exact["found"] is True
        assert res_exact["selected"]["name_en"] == "Vildagliptin"

        res_fuzzy = search_compound("Vildaglptin", allow_network=True)
        assert res_fuzzy["found"] is True
        assert "المسار البديل" in res_fuzzy["connection_status"]


def test_autocomplete_and_identity_verification():
    suggs = get_autocomplete_suggestions("سيتا")
    assert len(suggs) >= 1
    assert suggs[0]["name_en"] == "Sitagliptin"

    ver = verify_compound_identity(suggs[0])
    assert ver["is_verified"] is True
    assert all(ver["verification_steps"].values())


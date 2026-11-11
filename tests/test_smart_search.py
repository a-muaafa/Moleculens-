from config.settings import (
    ENTITY_BIOMOLECULE_PROTEIN,
    ENTITY_CHEMICAL_COMPOUND,
    ENTITY_CHEMICAL_ELEMENT,
    ENTITY_DRUG_SUBSTANCE,
)
from search.search_engine import get_live_suggestions, search_scientific_entity


def test_water_multilingual_and_formula_resolution():
    for q in ["الماء", "ماء", "Water", "H2O", "h2o", "H₂O", "oxidane"]:
        res = search_scientific_entity(q, allow_network=False)
        assert res["found"] is True
        assert res["is_ambiguous"] is False
        assert res["selected"]["key"] == "water_h2o"
        assert res["selected"]["symbol_or_formula"] == "H2O"
        assert res["selected"]["entity_type"] == ENTITY_CHEMICAL_COMPOUND


def test_element_search_and_periodic_table_fallback():
    for q in ["Fe", "Iron", "الحديد", "حديد"]:
        res = search_scientific_entity(q, allow_network=False)
        assert res["found"] is True
        assert res["selected"]["symbol_or_formula"] == "Fe"
        assert res["selected"]["entity_type"] == ENTITY_CHEMICAL_ELEMENT

    res_au = search_scientific_entity("Au", allow_network=False)
    assert res_au["found"] is True
    assert res_au["selected"]["entity_type"] == ENTITY_CHEMICAL_ELEMENT
    assert res_au["selected"]["identifiers"]["atomic_number"] == 79


def test_chemical_formula_search_and_no_fuzzy_pollution():
    for formula_q, expected_key in [("CO2", "carbon_dioxide_co2"), ("CO₂", "carbon_dioxide_co2"), ("O2", "oxygen_gas_o2"), ("H2", "hydrogen_gas_h2"), ("H2O2", "hydrogen_peroxide_h2o2")]:
        res = search_scientific_entity(formula_q, allow_network=False)
        assert res["found"] is True
        assert res["selected"]["key"] == expected_key

    for unknown_formula in ["C99H2", "Zz9"]:
        res_unk = search_scientific_entity(unknown_formula, allow_network=False)
        assert res_unk["found"] is False
        assert res_unk["candidates"] == []


def test_ambiguous_formula_c6h12o6():
    res = search_scientific_entity("C6H12O6", allow_network=False)
    assert res["found"] is True
    assert res["is_ambiguous"] is True
    assert res["selected"] is None
    keys = {c["key"] for c in res["candidates"]}
    assert "glucose_c6h12o6" in keys and "fructose_c6h12o6" in keys


def test_drug_and_protein_search():
    for q in ["Sitagliptin", "سيتاجليبتين", "Januvia", "MK-0431"]:
        res = search_scientific_entity(q, allow_network=False)
        assert res["found"] is True
        assert res["selected"]["key"] == "sitagliptin"
        assert res["selected"]["entity_type"] == ENTITY_DRUG_SUBSTANCE

    for q_p in ["Insulin", "الإنسولين", "INS"]:
        res_p = search_scientific_entity(q_p, allow_network=False)
        assert res_p["found"] is True
        assert res_p["selected"]["key"] == "human_insulin"
        assert res_p["selected"]["entity_type"] == ENTITY_BIOMOLECULE_PROTEIN


def test_approximate_name_typo_and_suggestions():
    res = search_scientific_entity("Sitaglptin", allow_network=False)
    assert res["found"] is True
    assert any(c["key"] == "sitagliptin" for c in res["candidates"])

    suggs = get_live_suggestions("الما")
    assert len(suggs) >= 1
    assert suggs[0]["key"] == "water_h2o"



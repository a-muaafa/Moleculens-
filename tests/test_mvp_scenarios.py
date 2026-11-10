from config.settings import (
    ENTITY_CHEMICAL_COMPOUND,
    ENTITY_CHEMICAL_ELEMENT,
    ENTITY_DRUG_SUBSTANCE,
)
from core.orchestrator import MolecuLensCore


def test_mvp_scenario_1_search_fe_element():
    res = MolecuLensCore.explore_entity("Fe", allow_network=False)
    assert res["found"] is True
    assert res["entity"]["entity_type"] == ENTITY_CHEMICAL_ELEMENT
    assert res["profile"]["basic_info"]["symbol"] == "Fe"
    assert res["profile"]["basic_info"]["atomic_number"] == 26


def test_mvp_scenario_2_search_h2o_structure_and_composition():
    res = MolecuLensCore.explore_entity("H2O", allow_network=False)
    assert res["found"] is True
    assert res["entity"]["entity_type"] == ENTITY_CHEMICAL_COMPOUND
    assert res["profile"]["basic_info"]["formula"] == "H2O"
    assert "ذرتان هيدروجين" in res["profile"]["composition"]["summary_ar"]
    assert "104.45" in res["profile"]["structure"]["molecular_geometry"]


def test_mvp_scenario_3_search_co2_by_chemical_formula():
    res = MolecuLensCore.explore_entity("CO2", allow_network=False)
    assert res["found"] is True
    assert res["entity"]["entity_type"] == ENTITY_CHEMICAL_COMPOUND
    assert res["profile"]["basic_info"]["formula"] == "CO2"
    assert res["profile"]["basic_info"]["name_ar"] == "ثاني أكسيد الكربون"


def test_mvp_scenario_4_search_sitagliptin_pharmacology():
    res = MolecuLensCore.explore_entity("Sitagliptin", allow_network=False)
    assert res["found"] is True
    assert res["entity"]["entity_type"] == ENTITY_DRUG_SUBSTANCE
    assert "pharmacology" in res["profile"]
    assert "CHEMBL1422" in str(res["profile"]["identifiers"].get("chembl_id"))


def test_mvp_scenario_5_reaction_h2_and_o2():
    rxn = MolecuLensCore.run_reaction_analysis(
        [
            {"formula": "H2", "amount": 2.0, "unit": "mol"},
            {"formula": "O2", "amount": 1.0, "unit": "mol"},
        ]
    )
    assert rxn["supported"] is True
    assert rxn["equation_text"] == "2H2 + O2 -> 2H2O"
    assert rxn["atom_conservation"]["is_balanced"] is True
    assert rxn["theoretical_yields"][0]["theoretical_yield_mol"] == 2.0


def test_mvp_scenario_6_unsupported_reaction_rejection():
    rxn = MolecuLensCore.run_reaction_analysis(
        [
            {"formula": "H2O", "amount": 1.0, "unit": "mol"},
            {"formula": "O2", "amount": 1.0, "unit": "mol"},
        ]
    )
    assert rxn["supported"] is False
    assert rxn["rejection_reason_ar"] != ""


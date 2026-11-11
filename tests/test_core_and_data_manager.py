from config.settings import (
    ENTITY_BIOMOLECULE_PROTEIN,
    ENTITY_CHEMICAL_COMPOUND,
    ENTITY_CHEMICAL_ELEMENT,
    ENTITY_DRUG_SUBSTANCE,
    OP_SCIENTIFIC_CALCULATION,
    OP_SCIENTIFIC_RETRIEVAL,
)
from core.schemas import build_provenance_tag
from core.validators import validate_reaction_inputs, validate_search_query
from data_integrations.cache_store import clear_all_caches, get_cached_item, set_cached_item
from data_integrations.data_manager import DataManager
from data_integrations.ord_client import search_ord_reactions
from search.entity_classifier import classify_scientific_entity
from search.local_catalog import get_catalog_entity_by_key


def test_entity_classification_all_four_types():
    fe = get_catalog_entity_by_key("iron_fe")
    c_fe = classify_scientific_entity(fe)
    assert c_fe["entity_type"] == ENTITY_CHEMICAL_ELEMENT
    assert c_fe["label_ar"] == "عنصر كيميائي"

    h2o = get_catalog_entity_by_key("water_h2o")
    c_h2o = classify_scientific_entity(h2o)
    assert c_h2o["entity_type"] == ENTITY_CHEMICAL_COMPOUND
    assert c_h2o["label_ar"] == "مركب كيميائي"

    sita = get_catalog_entity_by_key("sitagliptin")
    c_sita = classify_scientific_entity(sita)
    assert c_sita["entity_type"] == ENTITY_DRUG_SUBSTANCE
    assert c_sita["label_ar"] == "مادة دوائية"

    ins = get_catalog_entity_by_key("human_insulin")
    c_ins = classify_scientific_entity(ins)
    assert c_ins["entity_type"] == ENTITY_BIOMOLECULE_PROTEIN
    assert c_ins["label_ar"] == "بروتين أو جزيء حيوي"


def test_data_manager_unification_and_sources():
    sita = get_catalog_entity_by_key("sitagliptin")
    unified = DataManager.unify_entity_record(sita, allow_network_enrichment=False)
    assert unified["key"] == "sitagliptin"
    assert unified["entity_type"] == ENTITY_DRUG_SUBSTANCE
    assert "ChEMBL" in unified["sources"]
    assert "PubChem" in unified["sources"]
    assert unified["identifiers"]["chembl_id"] == "CHEMBL1422"
    assert unified["identity_verification"]["is_verified"] is True


def test_structure_and_identifier_record_matching():
    primary = {
        "name_en": "Acetaminophen",
        "symbol_or_formula": "C8H9NO2",
        "canonical_smiles": "CC(=O)Nc1ccc(O)cc1",
        "inchikey": "RZVAJINKPMORJF-UHFFFAOYSA-N",
    }
    chembl_cands = [
        {
            "chembl_id": "CHEMBL112",
            "name_en": "PARACETAMOL",
            "symbol_or_formula": "C8H9NO2",
            "canonical_smiles": "CC(=O)Nc1ccc(O)cc1",
            "inchikey": "RZVAJINKPMORJF-UHFFFAOYSA-N",
            "max_phase": 4,
        }
    ]
    matched = DataManager.match_records_by_structure_and_identifiers(primary, chembl_cands)
    assert matched is not None
    assert matched["chembl_id"] == "CHEMBL112"


def test_ord_reaction_search_integration():
    ord_res = search_ord_reactions("CH3COOH", allow_network=False)
    assert ord_res["ok"] is True
    assert len(ord_res["records"]) >= 2
    assert all(r["atom_conservation_verified"] is True for r in ord_res["records"])

    ext_rxns = DataManager.search_documented_reactions_external("CH3COOH", domain="organic_chemistry", allow_network=False)
    assert ext_rxns["ok"] is True
    assert len(ext_rxns["records"]) >= 2


def test_cache_store_and_validators():
    clear_all_caches()
    set_cached_item("test_key", {"val": 42})
    hit, val = get_cached_item("test_key")
    assert hit is True
    assert val == {"val": 42}
    clear_all_caches()

    ok_q, clean_q, _ = validate_search_query("  H2O  ")
    assert ok_q is True
    assert clean_q == "H2O"

    bad_q, _, err = validate_search_query("   ")
    assert bad_q is False
    assert err != ""

    ok_r, _ = validate_reaction_inputs([{"substance": {"symbol_or_formula": "H2"}, "amount": 2.0, "unit": "mol"}])
    assert ok_r is True

    ok_stp_gas, _ = validate_reaction_inputs([{"substance": {"symbol_or_formula": "H2"}, "amount": 22.414, "unit": "L_STP"}])
    assert ok_stp_gas is True

    bad_stp_liquid, err_stp = validate_reaction_inputs([{"substance": {"symbol_or_formula": "H2O"}, "amount": 22.414, "unit": "L_STP"}])
    assert bad_stp_liquid is False
    assert "L_STP" in err_stp

    bad_r, _ = validate_reaction_inputs([{"substance": {"symbol_or_formula": "H2"}, "amount": -1.0, "unit": "mol"}])
    assert bad_r is False


def test_provenance_tags_separation():
    ret_tag = build_provenance_tag(OP_SCIENTIFIC_RETRIEVAL, ["PubChem"])
    calc_tag = build_provenance_tag(OP_SCIENTIFIC_CALCULATION, ["RDKit"])
    assert ret_tag["operation_type"] == OP_SCIENTIFIC_RETRIEVAL
    assert calc_tag["operation_type"] == OP_SCIENTIFIC_CALCULATION
    assert ret_tag["operation_label_ar"] != calc_tag["operation_label_ar"]

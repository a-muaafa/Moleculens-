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

    bad_r, _ = validate_reaction_inputs([{"substance": {"symbol_or_formula": "H2"}, "amount": -1.0, "unit": "mol"}])
    assert bad_r is False


def test_provenance_tags_separation():
    ret_tag = build_provenance_tag(OP_SCIENTIFIC_RETRIEVAL, ["PubChem"])
    calc_tag = build_provenance_tag(OP_SCIENTIFIC_CALCULATION, ["RDKit"])
    assert ret_tag["operation_type"] == OP_SCIENTIFIC_RETRIEVAL
    assert calc_tag["operation_type"] == OP_SCIENTIFIC_CALCULATION
    assert ret_tag["operation_label_ar"] != calc_tag["operation_label_ar"]


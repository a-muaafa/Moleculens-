import requests
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
from data_integrations.chebi_client import search_chebi_entities
from data_integrations.chembl_client import fetch_chembl_drug_info, search_chembl_molecules
from data_integrations.data_manager import DataManager
from data_integrations.ord_client import (
    is_authentic_ord_id,
    search_ord_reactions,
    search_rhea_reactions,
    verify_external_equation_balance,
)
from data_integrations.pdb_client import fetch_pdb_structure_metadata
from data_integrations.pubchem_client import fetch_pubchem_by_formula, fetch_pubchem_by_name
from data_integrations.rhea_client import fetch_rhea_reaction_info
from data_integrations.uniprot_client import fetch_uniprot_protein_info
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


def test_ord_reaction_search_integration_and_separation(monkeypatch):
    clear_all_caches()
    # 1. Local registry search must NOT fabricate ord-rxn_... IDs or claim to be from ORD
    ord_res = search_ord_reactions("CH3COOH", allow_network=False)
    assert ord_res["ok"] is True
    assert len(ord_res["local_organic_records"]) >= 2
    assert ord_res["ord_records"] == []
    for r in ord_res["local_organic_records"]:
        assert r["record_origin"] == "local_curated_registry"
        assert r["is_verified_ord_record"] is False
        assert r["ord_id"] is None
        assert r["source_id"] != ""
        assert r["reference_url"].startswith("https://")
        assert r["atom_conservation_verified"] is True

    # 2. Verify external equation balance checker on balanced and unbalanced equations & reaction SMILES
    assert verify_external_equation_balance("CH3COOH + C2H5OH -> CH3COOCH2CH3 + H2O")["is_balanced"] is True
    assert verify_external_equation_balance("CH4 + O2 -> CO2 + H2O")["is_balanced"] is False
    assert verify_external_equation_balance("CC(=O)O.CCO>>CC(=O)OCC.O")["is_balanced"] is True
    assert verify_external_equation_balance("CC(=O)O.CCO>>CC(=O)OCC")["is_balanced"] is False

    # 3. Verify authentic ORD ID discriminator
    assert is_authentic_ord_id("ord-1234567890abcdef1234567890abcdef") is True
    assert is_authentic_ord_id("ord-rxn_065_organic_reaction") is False

    # 4. Simulate live ORD API returning one balanced reaction, one unbalanced reaction, and one non-ORD ID
    class _MockOrdResponse:
        status_code = 200
        headers = {"content-type": "application/json"}

        @staticmethod
        def json():
            return {
                "reactions": [
                    {
                        "reaction_id": "ord-1234567890abcdef1234567890abcdef",
                        "name": "Fischer Esterification",
                        "reaction_smiles": "CC(=O)O.CCO>>CC(=O)OCC.O",
                    },
                    {
                        "reaction_id": "ord-abcdef1234567890abcdef1234567890",
                        "name": "Unbalanced Incomplete Entry",
                        "reaction_smiles": "CC(=O)O.CCO>>CC(=O)OCC",
                    },
                    {
                        "reaction_id": "custom_local_99",
                        "name": "Unverified Source Entry",
                        "equation": "C2H4 + H2 -> C2H6",
                    },
                ]
            }

    clear_all_caches()
    monkeypatch.setattr(requests, "get", lambda *a, **kw: _MockOrdResponse())
    live_res = search_ord_reactions("ester", allow_network=True)
    assert len(live_res["ord_records"]) == 3

    bal_rec = live_res["ord_records"][0]
    assert bal_rec["ord_id"] == "ord-1234567890abcdef1234567890abcdef"
    assert bal_rec["is_verified_ord_record"] is True
    assert bal_rec["record_origin"] == "ord_external"
    assert bal_rec["atom_conservation_verified"] is True
    assert bal_rec["reference_url"].endswith("/ord-1234567890abcdef1234567890abcdef")

    unbal_rec = live_res["ord_records"][1]
    assert unbal_rec["ord_id"] == "ord-abcdef1234567890abcdef1234567890"
    assert unbal_rec["is_verified_ord_record"] is True
    assert unbal_rec["atom_conservation_verified"] is False

    unver_rec = live_res["ord_records"][2]
    assert unver_rec["ord_id"] is None
    assert unver_rec["is_verified_ord_record"] is False
    assert unver_rec["record_origin"] == "external_unverified"
    assert unver_rec["atom_conservation_verified"] is True


def test_external_sources_network_failures_and_incomplete_data(monkeypatch):
    clear_all_caches()

    def _raise_timeout(*args, **kwargs):
        raise requests.exceptions.ConnectTimeout("Simulated network timeout")

    monkeypatch.setattr(requests, "get", _raise_timeout)

    assert fetch_pubchem_by_name("Aspirin")["ok"] is False
    assert fetch_pubchem_by_formula("C9H8O4")["ok"] is False
    assert fetch_chembl_drug_info("CHEMBL25")["ok"] is False
    assert search_chembl_molecules("Aspirin")["ok"] is False
    assert search_chebi_entities("Aspirin")["ok"] is False
    assert fetch_uniprot_protein_info("P01308")["ok"] is False
    assert fetch_pdb_structure_metadata("4INS")["ok"] is False
    assert fetch_rhea_reaction_info("RHEA:15889")["ok"] is False
    assert search_rhea_reactions("catalase")["ok"] is False

    # Test incomplete / mismatched PubChem payload (e.g., formula C6H12O6 with SMILES for water 'O')
    class _MockBadPubChem:
        status_code = 200

        @staticmethod
        def json():
            return {
                "PropertyTable": {
                    "Properties": [
                        {
                            "CID": 999999,
                            "Title": "Corrupted Record",
                            "MolecularFormula": "C6H12O6",
                            "CanonicalSMILES": "O",
                        },
                        {
                            "CID": 888888,
                            "Title": "Empty Record",
                        },
                    ]
                }
            }

    clear_all_caches()
    monkeypatch.setattr(requests, "get", lambda *a, **kw: _MockBadPubChem())
    bad_pc = fetch_pubchem_by_name("Corrupted")
    assert bad_pc["ok"] is False
    assert bad_pc["records"] == []


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

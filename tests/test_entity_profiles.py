from config.settings import (
    ENTITY_BIOMOLECULE_PROTEIN,
    ENTITY_CHEMICAL_COMPOUND,
    ENTITY_CHEMICAL_ELEMENT,
    ENTITY_DRUG_SUBSTANCE,
)
from core.orchestrator import MolecuLensCore


def test_iron_element_profile_differs_from_molecule_profile():
    res = MolecuLensCore.explore_entity("Fe", allow_network=False)
    assert res["found"] is True
    profile = res["profile"]
    assert profile["profile_kind"] == ENTITY_CHEMICAL_ELEMENT
    assert profile["has_molecular_structure"] is False
    assert profile["basic_info"]["atomic_number"] == 26
    assert profile["atomic_structure"]["electron_configuration"] == "[Ar] 3d6 4s2"
    assert "composition" not in profile


def test_water_and_co2_compound_profiles():
    res_w = MolecuLensCore.explore_entity("H2O", allow_network=False)
    prof_w = res_w["profile"]
    assert prof_w["profile_kind"] == ENTITY_CHEMICAL_COMPOUND
    assert prof_w["has_molecular_structure"] is True
    assert prof_w["basic_info"]["formula"] == "H2O"
    assert prof_w["basic_info"]["molar_mass_value"] == 18.015
    assert "ذرتان هيدروجين" in prof_w["composition"]["summary_ar"]
    assert "ذرة أكسجين" in prof_w["composition"]["summary_ar"]
    assert prof_w["structure"]["image_2d_png"] is not None
    assert prof_w["explanation"]["student_level_ar"] != ""
    assert prof_w["explanation"]["advanced_level_ar"] != ""

    res_co2 = MolecuLensCore.explore_entity("CO2", allow_network=False)
    prof_co2 = res_co2["profile"]
    assert prof_co2["basic_info"]["formula"] == "CO2"
    assert prof_co2["basic_info"]["molar_mass_value"] == 44.009
    assert "ذرة كربون" in prof_co2["composition"]["summary_ar"]
    assert "ذرتان أكسجين" in prof_co2["composition"]["summary_ar"]


def test_sitagliptin_drug_profile():
    res = MolecuLensCore.explore_entity("Sitagliptin", allow_network=False)
    prof = res["profile"]
    assert prof["profile_kind"] == ENTITY_DRUG_SUBSTANCE
    assert "pharmacology" in prof
    assert "DPP4" in prof["pharmacology"]["target_protein_ar"]
    assert prof["pharmacology"]["calculated_drug_descriptors"]["lipinski_compliant"] is True


def test_insulin_protein_profile():
    res = MolecuLensCore.explore_entity("Insulin", allow_network=False)
    prof = res["profile"]
    assert prof["profile_kind"] == ENTITY_BIOMOLECULE_PROTEIN
    assert prof["biology"]["uniprot_id"] == "P01308"
    assert len(prof["biology"]["pdb_links"]) >= 1


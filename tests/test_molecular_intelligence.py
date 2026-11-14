"""
Unit tests for Molecular Intelligence module (Tanimoto similarity, Morgan fingerprints,
chemical structure modification, and LogP benchmark evaluation).
"""

import pytest
from chemistry_engine.molecular_intelligence import MolecularIntelligence


def test_compute_fingerprint_valid_and_invalid():
    # Valid SMILES (Ethanol)
    fp = MolecularIntelligence.compute_fingerprint("CCO")
    assert fp is not None

    # Invalid SMILES
    fp_invalid = MolecularIntelligence.compute_fingerprint("INVALID_SMILES_XYZ")
    assert fp_invalid is None

    # Empty string
    fp_empty = MolecularIntelligence.compute_fingerprint("")
    assert fp_empty is None


def test_calculate_similarity_identical_and_different():
    # Identical molecules should have similarity of 1.0
    sim_identical = MolecularIntelligence.calculate_similarity("CC(=O)O", "CC(=O)O")
    assert pytest.approx(sim_identical, 0.001) == 1.0

    # Structurally related molecules (Acetic acid and Ethanol)
    sim_related = MolecularIntelligence.calculate_similarity("CC(=O)O", "CCO")
    assert 0.0 < sim_related < 1.0

    # Dissimilar molecules
    sim_distant = MolecularIntelligence.calculate_similarity("C1CCCCC1", "c1ccccc1")
    assert 0.0 <= sim_distant <= 1.0

    # Invalid SMILES returns None
    sim_bad = MolecularIntelligence.calculate_similarity("CCO", "BAD_SMILES")
    assert sim_bad is None


def test_find_most_similar_compounds():
    # Query aspirin SMILES
    aspirin_smiles = "CC(=O)Oc1ccccc1C(=O)O"
    results = MolecularIntelligence.find_most_similar(aspirin_smiles, top_n=5)
    assert len(results) > 0
    assert len(results) <= 5

    # Should be sorted in descending order of similarity
    scores = [r["similarity_score"] for r in results]
    assert scores == sorted(scores, reverse=True)

    # Top hit or high hit should be related (e.g. salicylic acid or aspirin derivatives)
    top_keys = [r["key"] for r in results]
    assert any("salicylic" in k or "aspirin" in k for k in top_keys)


def test_structure_modification_hydroxyl_to_chloro():
    # Phenol (c1ccccc1O) with -OH -> -Cl
    phenol_smiles = "c1ccccc1O"
    result = MolecularIntelligence.modify_structure(phenol_smiles, "hydroxyl_to_chloro")
    assert result["success"] is True
    assert "Cl" in result["after"]["smiles"]
    assert result["deltas"]["mw"] > 0
    assert result["deltas"]["logp"] > 0
    assert result["impact_summary_ar"] != ""
    assert result["chemical_validity"]["is_sanitized"] is True


def test_structure_modification_unmatched_substructure():
    # Methane (C) has no hydroxyl group to replace with chloro
    methane_smiles = "C"
    result = MolecularIntelligence.modify_structure(methane_smiles, "hydroxyl_to_chloro")
    assert result["success"] is False
    assert "لم يحتوي" in result["error_ar"] or "لا يحتوي" in result["error_ar"]


def test_structure_modification_invalid_smiles():
    result = MolecularIntelligence.modify_structure("NOT_A_SMILES", "hydroxyl_to_chloro")
    assert result["success"] is False
    assert "error_ar" in result


def test_evaluate_logp_prediction_known_molecule():
    # Ibuprofen: CC(C)Cc1ccc(cc1)C(C)C(=O)O
    ibuprofen_smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    eval_result = MolecularIntelligence.evaluate_logp_prediction(ibuprofen_smiles)
    assert eval_result["available"] is True
    assert eval_result["predicted_logp"] > 2.0
    assert eval_result["confidence_interval_95"]["lower"] < eval_result["predicted_logp"] < eval_result["confidence_interval_95"]["upper"]
    assert pytest.approx(eval_result["r2_score"], 0.001) == 0.884
    assert pytest.approx(eval_result["rmse"], 0.001) == 0.652
    assert eval_result["methodology_ar"] != ""
    assert "pedagogical_note_ar" in eval_result


def test_evaluate_logp_prediction_invalid_molecule():
    eval_result = MolecularIntelligence.evaluate_logp_prediction("INVALID_XYZ")
    assert eval_result["available"] is False
    assert "message_ar" in eval_result

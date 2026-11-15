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
    assert result["is_virtual_modification"] is True
    assert result["modification_nature_ar"] != ""
    assert isinstance(result["scientific_caveats_ar"], list) and len(result["scientific_caveats_ar"]) > 0
    assert "Cl" in result["after"]["smiles"]
    assert result["deltas"]["mw"] > 0
    assert result["deltas"]["logp"] > 0
    assert result["impact_summary_ar"] != ""
    assert result["chemical_validity"]["is_sanitized"] is True
    assert result["chemical_validity"]["is_single_fragment"] is True
    assert result["chemical_validity"]["has_zero_radicals"] is True
    assert "catalog_match" in result["chemical_validity"]


def test_structure_modification_catalog_match_and_alerts():
    # Toluene -> benzoic acid (matches local catalog)
    tol_res = MolecularIntelligence.modify_structure("Cc1ccccc1", "methyl_to_carboxyl")
    assert tol_res["success"] is True
    assert tol_res["chemical_validity"]["catalog_match"]["is_documented_compound"] is True
    assert tol_res["chemical_validity"]["catalog_match"]["key"] == "benzoic_acid"

    # Acetic acid -OH -> -Cl produces acetyl chloride (acyl halide stability alert)
    acid_res = MolecularIntelligence.modify_structure("CC(=O)O", "hydroxyl_to_chloro")
    assert acid_res["success"] is True
    assert len(acid_res["chemical_validity"]["stability_alerts"]) > 0

    # Benzene H -> F succeeds without Kekulize errors
    bz_res = MolecularIntelligence.modify_structure("c1ccccc1", "hydrogen_to_fluoro")
    assert bz_res["success"] is True
    assert "F" in bz_res["after"]["smiles"]


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
    # Ibuprofen (held-out molecule with experimental LogP=3.97, MW above small training range max 194.19)
    ibuprofen_smiles = "CC(C)Cc1ccc(cc1)C(C)C(=O)O"
    eval_result = MolecularIntelligence.evaluate_logp_prediction(ibuprofen_smiles)
    assert eval_result["available"] is True
    assert "confidence_interval_95" not in eval_result
    assert pytest.approx(eval_result["experimental_logp"], 0.01) == 3.97

    det = eval_result["deterministic_clogp"]
    assert det["available"] is True
    assert det["calculated_clogp"] > 2.0
    assert det["has_confidence_interval"] is False

    qsar = eval_result["qsar_model_evaluation"]
    assert qsar["model_available"] is True
    assert qsar["n_total"] == 32
    assert qsar["n_train"] == 24
    assert qsar["n_test"] == 8
    assert qsar["test_r2"] > 0.7
    assert qsar["test_rmse"] > 0.0
    assert len(qsar["test_set_table"]) == 8
    assert qsar["in_applicability_domain"] is False
    assert eval_result["methodology_ar"] != ""
    assert "pedagogical_note_ar" in eval_result

    # Aspirin falls inside training MW range -> applicability domain True
    aspirin_eval = MolecularIntelligence.evaluate_logp_prediction("CC(=O)Oc1ccccc1C(=O)O")
    assert aspirin_eval["available"] is True
    assert aspirin_eval["qsar_model_evaluation"]["in_applicability_domain"] is True


def test_evaluate_logp_prediction_invalid_molecule():
    eval_result = MolecularIntelligence.evaluate_logp_prediction("INVALID_XYZ")
    assert eval_result["available"] is False
    assert "message_ar" in eval_result


def test_evaluate_logp_prediction_numpy_fallback_when_sklearn_missing(monkeypatch):
    import chemistry_engine.molecular_intelligence as mi
    monkeypatch.setattr(mi.MolecularIntelligence, "_CACHED_QSAR_BUNDLE", None)
    monkeypatch.setattr(mi, "SKLEARN_AVAILABLE", False)

    eval_result = mi.MolecularIntelligence.evaluate_logp_prediction("CC(=O)Oc1ccccc1C(=O)O")
    assert eval_result["available"] is True
    qsar = eval_result["qsar_model_evaluation"]
    assert qsar["model_available"] is True
    assert qsar["test_r2"] > 0.7
    assert qsar["test_rmse"] > 0.0
    assert len(qsar["test_set_table"]) == 8

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
    # Ibuprofen (experimental LogP=3.97, MW=206.3 inside expanded 90-molecule training range [31.06, 360.45])
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
    assert qsar["n_total"] == 120
    assert qsar["n_train"] == 90
    assert qsar["n_test"] == 30
    assert qsar["test_r2"] > 0.8
    assert qsar["test_rmse"] > 0.0
    assert qsar["cv_folds"] == 5
    assert qsar["cv_r2"] > 0.8
    assert qsar["cv_rmse"] > 0.0
    assert len(qsar["test_set_table"]) == 30
    assert qsar["in_applicability_domain"] is True
    assert eval_result["methodology_ar"] != ""
    assert "pedagogical_note_ar" in eval_result

    # Amoxicillin (MW = 365.4 > 360.45 training max) is flagged as outside applicability domain
    amox_eval = MolecularIntelligence.evaluate_logp_prediction("CC1(C)SC2C(NC(=O)C(N)c3ccc(O)cc3)C(=O)N2C1C(=O)O")
    assert amox_eval["available"] is True
    assert amox_eval["qsar_model_evaluation"]["in_applicability_domain"] is False


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
    assert qsar["n_total"] == 120
    assert qsar["n_train"] == 90
    assert qsar["n_test"] == 30
    assert qsar["test_r2"] > 0.8
    assert qsar["cv_r2"] > 0.8
    assert qsar["test_rmse"] > 0.0
    assert len(qsar["test_set_table"]) == 30


def test_methyl_to_carboxyl_rejects_methyl_not_bonded_to_carbon():
    # Methyl formate: the only CH3 is bonded to oxygen, so no carbon-bound methyl exists
    result = MolecularIntelligence.modify_structure("COC=O", "methyl_to_carboxyl")
    assert result["success"] is False


def test_structure_modification_site_info_and_alternatives():
    # Symmetric molecule: two matching sites but a single distinct product
    sym = MolecularIntelligence.modify_structure("Cc1ccc(C)cc1", "methyl_to_carboxyl")
    assert sym["success"] is True
    assert sym["site_info"]["match_count"] == 2
    assert sym["site_info"]["distinct_products_count"] == 1

    # Asymmetric molecule: two different hydroxyls give two distinct products
    asym = MolecularIntelligence.modify_structure("OCc1ccc(O)cc1", "hydroxyl_to_chloro")
    assert asym["success"] is True
    assert asym["site_info"]["match_count"] == 2
    assert asym["site_info"]["distinct_products_count"] == 2
    assert len(asym["site_info"]["alternative_products_smiles"]) == 1
    assert asym["after"]["smiles"] not in asym["site_info"]["alternative_products_smiles"]


def test_limitations_text_matches_actual_split_sizes():
    res = MolecularIntelligence.evaluate_logp_prediction("CCO")
    q = res["qsar_model_evaluation"]
    assert q["model_available"] is True
    assert str(q["n_train"]) in res["limitations_ar"]
    assert str(q["n_test"]) in res["limitations_ar"]
    assert str(q["n_total"]) in res["methodology_ar"]


def test_qsar_split_is_consistent_and_group_validation_reported():
    q = MolecularIntelligence.evaluate_logp_prediction("CCO")["qsar_model_evaluation"]
    assert q["n_train"] + q["n_test"] == q["n_total"]
    assert q["n_test"] == int(round(q["n_total"] * 0.25))
    for key in ("group_cv_r2", "group_cv_rmse", "group_cv_mae", "crippen_all_r2", "crippen_all_rmse"):
        assert key in q
    assert sum(row["n"] for row in q["group_cv_by_class"]) == q["n_total"]
    assert q["group_cv_rmse"] > 0 and q["cv_rmse"] > 0


def test_chemical_class_assignment():
    from rdkit import Chem

    cases = {
        "CCO": "alcohol_or_phenol",
        "CC(=O)O": "carboxylic_acid",
        "COC=O": "ester",
        "c1ccccc1": "aromatic_hydrocarbon",
        "ClC(Cl)Cl": "halogenated",
    }
    for smi, expected in cases.items():
        assert MolecularIntelligence._chemical_class(Chem.MolFromSmiles(smi)) == expected

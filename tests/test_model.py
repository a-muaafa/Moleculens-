import json
import os
import joblib
import numpy as np
import pandas as pd
from config.settings import METRICS_PATH, MODEL_PATH, TEST_PRED_PATH, TRAIN_DATA_PATH, TRAINING_REF_PATH
from model.applicability_domain import classify_applicability_domain
from model.evaluate import evaluate_predictions, verify_zero_scaffold_leakage
from model.predict import predict_dpp4_activity


def test_model_and_reference_artifacts_exist():
    assert os.path.exists(MODEL_PATH)
    assert os.path.exists(TRAINING_REF_PATH)
    assert os.path.exists(METRICS_PATH)
    assert os.path.exists(TRAIN_DATA_PATH)
    assert os.path.exists(TEST_PRED_PATH)


def test_zero_scaffold_leakage_between_train_and_test():
    train_df = pd.read_csv(TRAIN_DATA_PATH)
    test_df = pd.read_csv(TEST_PRED_PATH)
    check = verify_zero_scaffold_leakage(train_df, test_df)
    assert check["smiles_leakage_count"] == 0
    assert check["scaffold_leakage_count"] == 0
    assert check["is_clean"] is True


def test_random_vs_scaffold_split_metrics():
    with open(METRICS_PATH, "r", encoding="utf-8") as f:
        metrics = json.load(f)

    assert metrics["scaffold_leakage_count"] == 0
    assert metrics["mae_pic50"] > 0
    assert metrics["rmse_pic50"] >= metrics["mae_pic50"]
    assert 0.0 < metrics["roc_auc"] <= 1.0

    rand_comp = metrics["random_split_comparison"]
    assert rand_comp["shared_scaffolds_leaked"] > 0
    assert "mae_pic50" in rand_comp and "rmse_pic50" in rand_comp and "r2_pic50" in rand_comp


def test_evaluate_predictions_computation():
    y_true = np.array([5.0, 5.5, 6.5, 7.5, 8.0])
    y_pred = np.array([5.2, 5.4, 6.3, 7.2, 7.9])
    res = evaluate_predictions(y_true, y_pred, active_threshold=6.0)
    assert res["mae_pic50"] < 0.3
    assert res["r2_pic50"] > 0.9
    assert res["roc_auc"] == 1.0


def test_applicability_domain_classification():
    ad_in = classify_applicability_domain(
        "N#C[C@@H]1CCCN1C(=O)CNC12CC3CC(CC(O)(C3)C1)C2",
        max_tanimoto=0.85,
        mean_top5_tanimoto=0.65,
    )
    assert ad_in["status_code"] == "In-Domain"
    assert ad_in["can_show_reliable_prediction"] is True

    ad_border = classify_applicability_domain(
        "N[C@@H](CC(=O)N1CCn2c(nnc2C(F)(F)F)C1)Cc1cc(F)c(F)cc1F",
        max_tanimoto=0.42,
        mean_top5_tanimoto=0.35,
    )
    assert ad_border["status_code"] == "Borderline"
    assert ad_border["can_show_reliable_prediction"] is True

    ad_ood = classify_applicability_domain(
        "CC(=O)Oc1ccccc1C(=O)O",
        max_tanimoto=0.14,
        mean_top5_tanimoto=0.11,
    )
    assert ad_ood["status_code"] == "Out-of-Domain"
    assert ad_ood["can_show_reliable_prediction"] is False

    ad_water = classify_applicability_domain("O", max_tanimoto=0.90)
    assert ad_water["status_code"] == "Out-of-Domain"
    assert ad_water["can_show_reliable_prediction"] is False


def test_predict_dpp4_activity_pipeline():
    model = joblib.load(MODEL_PATH)
    ref_bundle = joblib.load(TRAINING_REF_PATH)

    vilda_smi = "N#C[C@@H]1CCCN1C(=O)CNC12CC3CC(CC(O)(C3)C1)C2"
    pred = predict_dpp4_activity(vilda_smi, model, ref_bundle)
    assert pred["pred_pic50"] > 6.5
    assert pred["max_tanimoto"] >= 0.60
    assert pred["applicability_domain"]["status_code"] == "In-Domain"
    assert len(pred["nearest_train_compounds"]) == 5

    asp_pred = predict_dpp4_activity("CC(=O)Oc1ccccc1C(=O)O", model, ref_bundle)
    assert asp_pred["applicability_domain"]["status_code"] == "Out-of-Domain"


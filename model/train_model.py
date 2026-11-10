import json
import os
import sys
import joblib
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, rdFingerprintGenerator
from sklearn.ensemble import RandomForestRegressor
from sklearn.model_selection import train_test_split
from config.settings import (
    ACTIVE_PIC50_THRESHOLD,
    CLEAN_DATA_PATH,
    METRICS_PATH,
    MODEL_PATH,
    MORGAN_NBITS,
    MORGAN_RADIUS,
    RANDOM_SEED,
    RF_ESTIMATORS,
    TEST_PRED_PATH,
    TRAIN_DATA_PATH,
    TRAINING_REF_PATH,
)
from model.evaluate import evaluate_predictions, scaffold_split_indices, verify_zero_scaffold_leakage


def compute_fingerprints(smiles_list: list[str]):
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=MORGAN_RADIUS, fpSize=MORGAN_NBITS)
    fps_bitvect = []
    fps_numpy = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(str(smi))
        bv = gen.GetFingerprint(mol)
        arr = gen.GetFingerprintAsNumPy(mol).astype(np.float32)
        fps_bitvect.append(bv)
        fps_numpy.append(arr)
    return fps_bitvect, np.array(fps_numpy, dtype=np.float32)


def compute_training_property_bounds(train_df: pd.DataFrame) -> dict:
    mws, logps, tpsas = [], [], []
    for smi in train_df["canonical_smiles"]:
        mol = Chem.MolFromSmiles(str(smi))
        if mol is not None:
            mws.append(float(Descriptors.MolWt(mol)))
            logps.append(float(Descriptors.MolLogP(mol)))
            tpsas.append(float(Descriptors.TPSA(mol)))
    return {
        "mw_min": round(float(np.percentile(mws, 1)), 2),
        "mw_max": round(float(np.percentile(mws, 99)), 2),
        "logp_min": round(float(np.percentile(logps, 1)), 2),
        "logp_max": round(float(np.percentile(logps, 99)), 2),
        "tpsa_min": round(float(np.percentile(tpsas, 1)), 2),
        "tpsa_max": round(float(np.percentile(tpsas, 99)), 2),
    }


def train_dpp4_model(
    input_path: str = CLEAN_DATA_PATH,
    model_path: str = MODEL_PATH,
    ref_path: str = TRAINING_REF_PATH,
    metrics_path: str = METRICS_PATH,
) -> dict:
    if not os.path.exists(input_path):
        print("خطأ حرج: ملف البيانات المعالجة غير موجود. يرجى تشغيل clean_data.py أولا.", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(input_path).reset_index(drop=True)
    if len(df) < 50 or "pIC50" not in df.columns:
        print("خطأ حرج: البيانات المعالجة غير مكتملة أو غير صالحة للتدريب.", file=sys.stderr)
        sys.exit(1)

    fps_bv, X = compute_fingerprints(df["canonical_smiles"].tolist())
    y = df["pIC50"].values.astype(np.float64)

    train_idx, test_idx = scaffold_split_indices(df, train_frac=0.8)
    train_df = df.iloc[train_idx].copy().reset_index(drop=True)
    test_df = df.iloc[test_idx].copy().reset_index(drop=True)

    leakage_check = verify_zero_scaffold_leakage(train_df, test_df)
    if not leakage_check["is_clean"]:
        print("خطأ حرج: تم رصد تسرب بيانات بين مجموعتي التدريب والاختبار.", file=sys.stderr)
        sys.exit(1)

    X_train, y_train = X[train_idx], y[train_idx]
    X_test, y_test = X[test_idx], y[test_idx]

    model = RandomForestRegressor(n_estimators=RF_ESTIMATORS, random_state=RANDOM_SEED, n_jobs=-1)
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)
    scaf_eval = evaluate_predictions(y_test, y_pred, active_threshold=ACTIVE_PIC50_THRESHOLD)

    rand_train_idx, rand_test_idx = train_test_split(
        np.arange(len(df)), test_size=0.2, random_state=RANDOM_SEED
    )
    rand_model = RandomForestRegressor(n_estimators=RF_ESTIMATORS, random_state=RANDOM_SEED, n_jobs=-1)
    rand_model.fit(X[rand_train_idx], y[rand_train_idx])
    rand_pred = rand_model.predict(X[rand_test_idx])
    rand_eval = evaluate_predictions(y[rand_test_idx], rand_pred, active_threshold=ACTIVE_PIC50_THRESHOLD)
    rand_shared_scaffolds = len(
        set(df.iloc[rand_train_idx]["murcko_scaffold"]) & set(df.iloc[rand_test_idx]["murcko_scaffold"])
    )

    train_fps_bv = [fps_bv[i] for i in train_idx]
    max_sims, mean_top5_sims, nearest_ids = [], [], []
    for i in test_idx:
        sims = np.array(DataStructs.BulkTanimotoSimilarity(fps_bv[i], train_fps_bv), dtype=np.float64)
        top_order = np.argsort(sims)[::-1]
        max_sims.append(float(sims[top_order[0]]))
        mean_top5_sims.append(float(np.mean(sims[top_order[:5]])))
        nearest_ids.append(str(train_df.loc[int(top_order[0]), "molecule_chembl_id"]))

    test_df["true_pIC50"] = y_test
    test_df["pred_pIC50"] = y_pred
    test_df["abs_error_pIC50"] = np.abs(y_test - y_pred)
    test_df["residual_pIC50"] = y_pred - y_test
    test_df["true_ic50_nm"] = 10.0 ** (9.0 - y_test)
    test_df["pred_ic50_nm"] = 10.0 ** (9.0 - y_pred)
    test_df["max_train_tanimoto"] = max_sims
    test_df["mean_top5_train_tanimoto"] = mean_top5_sims
    test_df["nearest_train_chembl_id"] = nearest_ids

    os.makedirs(os.path.dirname(model_path), exist_ok=True)
    os.makedirs(os.path.dirname(TRAIN_DATA_PATH), exist_ok=True)

    train_df.to_csv(TRAIN_DATA_PATH, index=False)
    test_df.to_csv(TEST_PRED_PATH, index=False)
    joblib.dump(model, model_path)

    train_bounds = compute_training_property_bounds(train_df)
    reference_bundle = {
        "train_df": train_df,
        "train_fps_numpy": X_train,
        "train_bounds": train_bounds,
    }
    joblib.dump(reference_bundle, ref_path)

    metrics = {
        "split_method": "Bemis-Murcko Scaffold Split",
        "random_seed": RANDOM_SEED,
        "target_id": "CHEMBL284",
        "target_organism": "Homo sapiens",
        "predicted_measurement_type": "IC50 (via logarithmic pIC50)",
        "predicted_units": "nM (and pIC50)",
        "target_variable": "pIC50 (9 - log10(IC50_nM))",
        "activity_threshold_pic50": ACTIVE_PIC50_THRESHOLD,
        "activity_threshold_nm": 1000.0,
        "train_count": int(len(train_df)),
        "test_count": int(len(test_df)),
        "train_scaffolds": int(train_df["murcko_scaffold"].nunique()),
        "test_scaffolds": int(test_df["murcko_scaffold"].nunique()),
        "smiles_leakage_count": leakage_check["smiles_leakage_count"],
        "scaffold_leakage_count": leakage_check["scaffold_leakage_count"],
        "mae_pic50": scaf_eval["mae_pic50"],
        "rmse_pic50": scaf_eval["rmse_pic50"],
        "r2_pic50": scaf_eval["r2_pic50"],
        "roc_auc": scaf_eval["roc_auc"],
        "roc_auc_note_ar": scaf_eval["roc_auc_note_ar"],
        "train_property_bounds": train_bounds,
        "random_split_comparison": {
            "train_count": int(len(rand_train_idx)),
            "test_count": int(len(rand_test_idx)),
            "mae_pic50": rand_eval["mae_pic50"],
            "rmse_pic50": rand_eval["rmse_pic50"],
            "r2_pic50": rand_eval["r2_pic50"],
            "roc_auc": rand_eval["roc_auc"],
            "shared_scaffolds_leaked": int(rand_shared_scaffolds),
        },
    }

    with open(metrics_path, "w", encoding="utf-8") as f:
        json.dump(metrics, f, ensure_ascii=False, indent=2)

    print(
        f"اكتمل التدريب والتقييم بنجاح (Scaffold Split): "
        f"MAE={metrics['mae_pic50']:.3f}, RMSE={metrics['rmse_pic50']:.3f}, "
        f"R2={metrics['r2_pic50']:.3f}, ROC-AUC={metrics['roc_auc']:.3f}"
    )
    return metrics


if __name__ == "__main__":
    train_dpp4_model()


import json
import os
import sys
import joblib
import numpy as np
import pandas as pd
from rdkit import Chem, DataStructs
from rdkit.Chem import rdFingerprintGenerator
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.ensemble import RandomForestRegressor
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, roc_auc_score
from sklearn.model_selection import train_test_split

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


def compute_fingerprints(smiles_list):
    gen = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)
    fps_bitvect = []
    fps_numpy = []
    for smi in smiles_list:
        mol = Chem.MolFromSmiles(str(smi))
        bv = gen.GetFingerprint(mol)
        arr = gen.GetFingerprintAsNumPy(mol).astype(np.float32)
        fps_bitvect.append(bv)
        fps_numpy.append(arr)
    return fps_bitvect, np.array(fps_numpy, dtype=np.float32)


def scaffold_split_indices(df, train_frac=0.8):
    scaffolds = {}
    for idx, row in df.iterrows():
        scaf = str(row.get("murcko_scaffold", ""))
        if not scaf or scaf == "nan":
            mol = Chem.MolFromSmiles(str(row["canonical_smiles"]))
            scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False) if mol else str(row["canonical_smiles"])
        scaffolds.setdefault(scaf, []).append(idx)

    scaf_groups = sorted(
        scaffolds.values(),
        key=lambda g: (len(g), str(df.loc[g[0], "canonical_smiles"])),
        reverse=True,
    )

    train_idx, test_idx = [], []
    for idx_g, group in enumerate(scaf_groups):
        if idx_g % 5 == 1:
            test_idx.extend(group)
        else:
            train_idx.extend(group)
    return train_idx, test_idx


def train_dpp4_model():
    input_path = os.path.join(BASE_DIR, "data", "processed", "dpp4_cleaned.csv")
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

    smiles_overlap = len(set(train_df["canonical_smiles"]) & set(test_df["canonical_smiles"]))
    scaffold_overlap = len(set(train_df["murcko_scaffold"]) & set(test_df["murcko_scaffold"]))
    if smiles_overlap > 0 or scaffold_overlap > 0:
        print("خطأ حرج: تم رصد تسرب بيانات بين مجموعتي التدريب والاختبار.", file=sys.stderr)
        sys.exit(1)

    X_train, y_train = X[train_idx], y[train_idx]
    X_test, y_test = X[test_idx], y[test_idx]

    model = RandomForestRegressor(n_estimators=300, random_state=42, n_jobs=-1)
    model.fit(X_train, y_train)

    y_pred = model.predict(X_test)
    mae = float(mean_absolute_error(y_test, y_pred))
    rmse = float(np.sqrt(mean_squared_error(y_test, y_pred)))
    r2 = float(r2_score(y_test, y_pred))
    y_test_active = (y_test >= 6.0).astype(int)
    roc_auc = float(roc_auc_score(y_test_active, y_pred))

    rand_train_idx, rand_test_idx = train_test_split(
        np.arange(len(df)), test_size=0.2, random_state=42
    )
    rand_model = RandomForestRegressor(n_estimators=300, random_state=42, n_jobs=-1)
    rand_model.fit(X[rand_train_idx], y[rand_train_idx])
    rand_pred = rand_model.predict(X[rand_test_idx])
    rand_mae = float(mean_absolute_error(y[rand_test_idx], rand_pred))
    rand_rmse = float(np.sqrt(mean_squared_error(y[rand_test_idx], rand_pred)))
    rand_r2 = float(r2_score(y[rand_test_idx], rand_pred))
    rand_auc = float(roc_auc_score((y[rand_test_idx] >= 6.0).astype(int), rand_pred))
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

    models_dir = os.path.join(BASE_DIR, "models")
    proc_dir = os.path.join(BASE_DIR, "data", "processed")
    os.makedirs(models_dir, exist_ok=True)
    os.makedirs(proc_dir, exist_ok=True)

    train_df.to_csv(os.path.join(proc_dir, "train_compounds.csv"), index=False)
    test_df.to_csv(os.path.join(proc_dir, "test_predictions.csv"), index=False)

    joblib.dump(model, os.path.join(models_dir, "dpp4_rf_model.pkl"))

    reference_bundle = {
        "train_df": train_df,
        "train_fps_numpy": X_train,
    }
    joblib.dump(reference_bundle, os.path.join(models_dir, "training_reference.pkl"))

    metrics = {
        "split_method": "Bemis-Murcko Scaffold Split",
        "target_variable": "pIC50 (9 - log10(IC50_nM))",
        "activity_threshold_pic50": 6.0,
        "activity_threshold_nm": 1000.0,
        "train_count": int(len(train_df)),
        "test_count": int(len(test_df)),
        "train_scaffolds": int(train_df["murcko_scaffold"].nunique()),
        "test_scaffolds": int(test_df["murcko_scaffold"].nunique()),
        "smiles_leakage_count": int(smiles_overlap),
        "scaffold_leakage_count": int(scaffold_overlap),
        "mae_pic50": round(mae, 4),
        "rmse_pic50": round(rmse, 4),
        "r2_pic50": round(r2, 4),
        "roc_auc": round(roc_auc, 4),
        "random_split_comparison": {
            "mae_pic50": round(rand_mae, 4),
            "rmse_pic50": round(rand_rmse, 4),
            "r2_pic50": round(rand_r2, 4),
            "roc_auc": round(rand_auc, 4),
            "shared_scaffolds_leaked": int(rand_shared_scaffolds),
        },
    }

    with open(os.path.join(models_dir, "model_metrics.json"), "w", encoding="utf-8") as f:
        json.dump(metrics, f, ensure_ascii=False, indent=2)

    print(
        f"اكتمل التدريب والتقييم بنجاح (Scaffold Split): "
        f"MAE={mae:.3f}, RMSE={rmse:.3f}, R2={r2:.3f}, ROC-AUC={roc_auc:.3f}"
    )


if __name__ == "__main__":
    train_dpp4_model()

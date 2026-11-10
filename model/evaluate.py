import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem.Scaffolds import MurckoScaffold
from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score, roc_auc_score
from config.settings import ACTIVE_PIC50_THRESHOLD


def scaffold_split_indices(df: pd.DataFrame, train_frac: float = 0.8) -> tuple[list[int], list[int]]:
    scaffolds = {}
    for idx, row in df.iterrows():
        scaf = str(row.get("murcko_scaffold", ""))
        if not scaf or scaf == "nan":
            mol = Chem.MolFromSmiles(str(row["canonical_smiles"]))
            scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False) if mol else str(row["canonical_smiles"])
        scaffolds.setdefault(scaf, []).append(int(idx))

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


def verify_zero_scaffold_leakage(train_df: pd.DataFrame, test_df: pd.DataFrame) -> dict:
    train_smiles = set(train_df["canonical_smiles"].astype(str))
    test_smiles = set(test_df["canonical_smiles"].astype(str))
    smiles_overlap = len(train_smiles & test_smiles)

    train_scaf = set(train_df["murcko_scaffold"].astype(str))
    test_scaf = set(test_df["murcko_scaffold"].astype(str))
    scaffold_overlap = len(train_scaf & test_scaf)

    return {
        "smiles_leakage_count": int(smiles_overlap),
        "scaffold_leakage_count": int(scaffold_overlap),
        "is_clean": (smiles_overlap == 0) and (scaffold_overlap == 0),
    }


def evaluate_predictions(y_true: np.ndarray, y_pred: np.ndarray, active_threshold: float = ACTIVE_PIC50_THRESHOLD) -> dict:
    mae = float(mean_absolute_error(y_true, y_pred))
    rmse = float(np.sqrt(mean_squared_error(y_true, y_pred)))
    r2 = float(r2_score(y_true, y_pred))

    y_bin = (y_true >= active_threshold).astype(int)
    if len(np.unique(y_bin)) > 1:
        roc_auc = float(roc_auc_score(y_bin, y_pred))
    else:
        roc_auc = float("nan")

    return {
        "mae_pic50": round(mae, 4),
        "rmse_pic50": round(rmse, 4),
        "r2_pic50": round(r2, 4),
        "roc_auc": round(roc_auc, 4) if not np.isnan(roc_auc) else None,
        "roc_auc_note_ar": (
            "مقياس ROC-AUC يقيس قدرة النموذج الترتيبية على التمييز التصنيفي بين المركبات النشطة وغير النشطة "
            f"(عند عتبة pIC50 >= {active_threshold:.1f} أي IC50 <= 1000 nM)، وليس مقياسا لدقة التنبؤ العددي بقيمة IC50."
        ),
    }


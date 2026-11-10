from functools import lru_cache
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import rdFingerprintGenerator
from chemistry.molecular_structure import validate_smiles
from config.settings import MORGAN_NBITS, MORGAN_RADIUS
from model.applicability_domain import classify_applicability_domain


@lru_cache(maxsize=1)
def get_morgan_generator():
    return rdFingerprintGenerator.GetMorganGenerator(radius=MORGAN_RADIUS, fpSize=MORGAN_NBITS)


def compute_tanimoto_against_training(mol, ref_bundle: dict) -> tuple[np.ndarray, np.ndarray, float, float, pd.DataFrame]:
    gen = get_morgan_generator()
    q_arr = gen.GetFingerprintAsNumPy(mol).astype(np.float32)

    train_df = ref_bundle["train_df"]
    train_X = ref_bundle["train_fps_numpy"]

    dot_prod = np.dot(train_X, q_arr)
    q_bits = float(np.sum(q_arr))
    t_bits = np.sum(train_X, axis=1)
    denom = q_bits + t_bits - dot_prod
    sims = np.where(denom > 0, dot_prod / denom, 0.0)

    order = np.argsort(sims)[::-1]
    top_idx = order[:5]
    max_sim = float(sims[order[0]])
    mean_top5_sim = float(np.mean(sims[top_idx]))

    nearest_df = train_df.iloc[top_idx].copy()
    nearest_df["tanimoto"] = np.round(sims[top_idx], 3)
    return q_arr, sims, max_sim, mean_top5_sim, nearest_df


def interpret_predicted_activity(pic50_val: float) -> dict:
    ic50_nm = float(10.0 ** (9.0 - pic50_val))
    if pic50_val >= 7.0:
        return {
            "category_ar": "نشاط تثبيطي مرتفع (ضمن العتبات التعليمية)",
            "pill_class": "pill-high",
            "ic50_nm": round(ic50_nm, 2),
            "pic50": round(float(pic50_val), 2),
            "explanation_ar": (
                "القيمة المتوقعة تقع تحت 100 nM (pIC50 >= 7.0). يمثل هذا الرقم تقديرا حاسوبيا أوليا "
                "لنشاط تثبيطي مرتفع تجاه إنزيم DPP4 البشري، ولا يعد إثباتا مخبريا أو سريريا."
            ),
        }
    if pic50_val >= 6.0:
        return {
            "category_ar": "نشاط تثبيطي متوسط (ضمن العتبات التعليمية)",
            "pill_class": "pill-mid",
            "ic50_nm": round(ic50_nm, 2),
            "pic50": round(float(pic50_val), 2),
            "explanation_ar": (
                "القيمة المتوقعة تقع بين 100 و 1000 nM (pIC50 بين 6.0 و 7.0). يشير التقدير الحاسوبي "
                "إلى نشاط معتدل قد يستدعي تحسينات بنائية."
            ),
        }
    return {
        "category_ar": "نشاط تثبيطي ضعيف أو غير نشط (ضمن العتبات التعليمية)",
        "pill_class": "pill-low",
        "ic50_nm": round(ic50_nm, 2),
        "pic50": round(float(pic50_val), 2),
        "explanation_ar": (
            "القيمة المتوقعة أعلى من 1000 nM (pIC50 < 6.0). يشير التقدير الحاسوبي الأولي إلى ضعف "
            "القدرة التثبيطية للمركب تجاه إنزيم DPP4 البشري."
        ),
    }


def predict_dpp4_activity(smiles_or_mol, model, ref_bundle: dict) -> dict:
    if isinstance(smiles_or_mol, str):
        mol, canon_smi, err = validate_smiles(smiles_or_mol)
        if mol is None:
            raise ValueError(err)
    else:
        mol = smiles_or_mol
        canon_smi = Chem.MolToSmiles(mol, canonical=True)

    q_arr, all_sims, max_sim, mean_top5_sim, nearest_df = compute_tanimoto_against_training(mol, ref_bundle)
    pred_pic50 = float(model.predict([q_arr])[0])
    activity_info = interpret_predicted_activity(pred_pic50)
    ad_info = classify_applicability_domain(
        mol,
        max_tanimoto=max_sim,
        mean_top5_tanimoto=mean_top5_sim,
        train_bounds=ref_bundle.get("train_bounds"),
    )

    return {
        "canonical_smiles": canon_smi,
        "measurement_type": "IC50 (مشتق من التنبؤ اللوغاريتمي pIC50)",
        "units": "nM (وحدة تركيز التثبيط النصفي) / pIC50 (مقياس لوغاريتمي)",
        "pred_pic50": activity_info["pic50"],
        "pred_ic50_nm": activity_info["ic50_nm"],
        "category_ar": activity_info["category_ar"],
        "pill_class": activity_info["pill_class"],
        "explanation_ar": activity_info["explanation_ar"],
        "max_tanimoto": round(max_sim, 3),
        "mean_top5_tanimoto": round(mean_top5_sim, 3),
        "all_train_similarities": all_sims,
        "nearest_train_compounds": nearest_df,
        "applicability_domain": ad_info,
    }


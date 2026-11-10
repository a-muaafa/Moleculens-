import json
import os
import sys
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold

BASE_DIR = os.path.dirname(os.path.abspath(__file__))


def convert_to_nm(val, unit, mol_wt):
    if pd.isna(val) or val <= 0 or pd.isna(unit):
        return None
    u = str(unit).strip()
    if u == "nM":
        return float(val)
    if u in ("uM", "umol.L-1"):
        return float(val) * 1e3
    if u in ("mM", "mmol.L-1"):
        return float(val) * 1e6
    if u in ("pM", "pmol.L-1"):
        return float(val) * 1e-3
    if u in ("M", "mol.L-1"):
        return float(val) * 1e9
    if u in ("ug.mL-1", "ug/mL") and mol_wt and mol_wt > 0:
        return (float(val) * 1e6) / float(mol_wt)
    return None


def clean_and_featurize():
    input_path = os.path.join(BASE_DIR, "data", "raw", "dpp4_raw.csv")
    if not os.path.exists(input_path):
        print("خطأ حرج: ملف البيانات الخام غير موجود.", file=sys.stderr)
        sys.exit(1)

    df = pd.read_csv(input_path)
    required_cols = {
        "molecule_chembl_id",
        "canonical_smiles",
        "standard_type",
        "standard_relation",
        "standard_value",
        "standard_units",
        "assay_type",
    }
    if not required_cols.issubset(set(df.columns)) or df["canonical_smiles"].nunique() < 20:
        print("خطأ حرج: البيانات الخام غير صالحة أو تحتوي على سجلات اصطناعية.", file=sys.stderr)
        sys.exit(1)

    report = {"initial_records": int(len(df))}

    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")
    df = df.dropna(subset=["canonical_smiles", "standard_value"]).copy()
    df = df[df["standard_value"] > 0].copy()
    report["after_valid_numeric"] = int(len(df))

    if "target_chembl_id" in df.columns:
        df = df[df["target_chembl_id"] == "CHEMBL284"].copy()
    df = df[df["standard_type"] == "IC50"].copy()
    report["after_measurement_type_ic50"] = int(len(df))

    df = df[df["standard_relation"] == "="].copy()
    report["after_exact_relation"] = int(len(df))

    if "data_validity_comment" in df.columns:
        df = df[df["data_validity_comment"].isna() | (df["data_validity_comment"].astype(str).str.strip() == "")].copy()
    report["after_data_validity"] = int(len(df))

    df = df[df["assay_type"] == "B"].copy()
    report["after_binding_assay_filter"] = int(len(df))

    curated_rows = []
    converted_ug_count = 0
    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(str(row["canonical_smiles"]))
        if mol is None:
            continue
        canon_smi = Chem.MolToSmiles(mol, canonical=True)
        mw = Descriptors.MolWt(mol)
        orig_unit = str(row["standard_units"]).strip()
        val_nm = convert_to_nm(row["standard_value"], orig_unit, mw)
        if val_nm is None or val_nm <= 0:
            continue
        pic50 = 9.0 - np.log10(val_nm)
        if pic50 < 2.0 or pic50 > 12.0:
            continue
        if orig_unit != "nM":
            converted_ug_count += 1
        scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
        row_dict = row.to_dict()
        row_dict["canonical_smiles"] = canon_smi
        row_dict["original_units"] = orig_unit
        row_dict["standard_units"] = "nM"
        row_dict["standard_value"] = float(val_nm)
        row_dict["pIC50"] = float(pic50)
        row_dict["murcko_scaffold"] = scaf if scaf else canon_smi
        curated_rows.append(row_dict)

    df_curated = pd.DataFrame(curated_rows)
    report["after_unit_and_structure_standardization"] = int(len(df_curated))
    report["converted_non_nm_units"] = int(converted_ug_count)

    grp_span = df_curated.groupby("canonical_smiles")["pIC50"].transform(lambda s: s.max() - s.min())
    discordant_count = int((grp_span > 1.0).sum())
    df_curated = df_curated[grp_span <= 1.0].copy()
    report["discordant_replicates_removed"] = discordant_count

    agg_spec = {
        "molecule_chembl_id": "first",
        "standard_value": "median",
        "pIC50": "median",
        "standard_type": "first",
        "standard_relation": "first",
        "standard_units": "first",
        "original_units": "first",
        "assay_chembl_id": "first",
        "assay_type": "first",
        "assay_description": "first",
        "document_chembl_id": "first",
        "document_journal": "first",
        "document_year": "first",
        "target_chembl_id": "first",
        "target_pref_name": "first",
        "murcko_scaffold": "first",
    }
    if "molecule_pref_name" in df_curated.columns:
        agg_spec["molecule_pref_name"] = "first"
    if "bao_label" in df_curated.columns:
        agg_spec["bao_label"] = "first"

    counts = df_curated.groupby("canonical_smiles").size().rename("replicate_count")
    df_clean = df_curated.groupby("canonical_smiles", as_index=False).agg(agg_spec)
    df_clean = df_clean.merge(counts, on="canonical_smiles", how="left")
    df_clean["standard_value"] = 10.0 ** (9.0 - df_clean["pIC50"])

    report["final_unique_compounds"] = int(len(df_clean))
    report["unique_murcko_scaffolds"] = int(df_clean["murcko_scaffold"].nunique())

    out_dir = os.path.join(BASE_DIR, "data", "processed")
    os.makedirs(out_dir, exist_ok=True)
    output_path = os.path.join(out_dir, "dpp4_cleaned.csv")
    report_path = os.path.join(out_dir, "cleaning_report.json")

    df_clean.to_csv(output_path, index=False)
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2)

    print(f"اكتمل التنظيف العلمي للبيانات. عدد المركبات الفريدة الصالحة للتدريب: {len(df_clean)}")


if __name__ == "__main__":
    clean_and_featurize()


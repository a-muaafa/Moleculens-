import json
import os
import sys
import numpy as np
import pandas as pd
from rdkit import Chem
from rdkit.Chem import Descriptors
from rdkit.Chem.Scaffolds import MurckoScaffold
from config.settings import CLEAN_DATA_PATH, CLEAN_REPORT_PATH, RAW_DATA_PATH, TARGET_INFO


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


def clean_chembl_dataframe(df_raw: pd.DataFrame, target_measurement_type: str = "IC50") -> tuple[pd.DataFrame, dict]:
    required_cols = {
        "molecule_chembl_id",
        "canonical_smiles",
        "standard_type",
        "standard_relation",
        "standard_value",
        "standard_units",
        "assay_type",
    }
    if not required_cols.issubset(set(df_raw.columns)):
        raise ValueError("ملف البيانات الخام يفتقر إلى الأعمدة الإلزامية لبيانات ChEMBL.")

    initial_count = int(len(df_raw))
    raw_type_dist = df_raw["standard_type"].fillna("Unknown").value_counts().to_dict()
    raw_type_dist = {str(k): int(v) for k, v in raw_type_dist.items()}

    df = df_raw.copy()
    df["standard_value"] = pd.to_numeric(df["standard_value"], errors="coerce")

    valid_num_mask = df["canonical_smiles"].notna() & df["standard_value"].notna() & (df["standard_value"] > 0)
    excluded_missing_or_nonpositive = int((~valid_num_mask).sum())
    df = df[valid_num_mask].copy()
    after_valid_numeric = int(len(df))

    target_id = TARGET_INFO["target_id"]
    target_org = TARGET_INFO["organism"]
    if "target_chembl_id" in df.columns:
        target_mask = df["target_chembl_id"] == target_id
    elif "target_id" in df.columns:
        target_mask = df["target_id"] == target_id
    else:
        target_mask = pd.Series(True, index=df.index)

    if "target_organism" in df.columns:
        org_mask = df["target_organism"].astype(str).str.strip() == target_org
    else:
        org_mask = pd.Series(True, index=df.index)

    combined_target_mask = target_mask & org_mask
    excluded_wrong_target_or_organism = int((~combined_target_mask).sum())
    df = df[combined_target_mask].copy()

    type_mask = df["standard_type"].astype(str).str.strip() == target_measurement_type
    excluded_other_measurement_type = int((~type_mask).sum())
    df = df[type_mask].copy()
    after_measurement_type = int(len(df))

    rel_mask = df["standard_relation"].astype(str).str.strip() == "="
    excluded_censored_relations = int((~rel_mask).sum())
    df = df[rel_mask].copy()
    after_exact_relation = int(len(df))

    if "data_validity_comment" in df.columns:
        val_mask = df["data_validity_comment"].isna() | (df["data_validity_comment"].astype(str).str.strip() == "")
    else:
        val_mask = pd.Series(True, index=df.index)
    excluded_validity_flags = int((~val_mask).sum())
    df = df[val_mask].copy()
    after_data_validity = int(len(df))

    assay_mask = df["assay_type"].astype(str).str.strip() == "B"
    excluded_non_binding_assays = int((~assay_mask).sum())
    df = df[assay_mask].copy()
    after_binding_assay = int(len(df))

    curated_rows = []
    converted_units_count = 0
    excluded_invalid_structure_or_unit = 0

    for _, row in df.iterrows():
        mol = Chem.MolFromSmiles(str(row["canonical_smiles"]))
        if mol is None:
            excluded_invalid_structure_or_unit += 1
            continue
        canon_smi = Chem.MolToSmiles(mol, canonical=True)
        mw = Descriptors.MolWt(mol)
        orig_unit = str(row["standard_units"]).strip()
        orig_val = float(row["standard_value"])
        val_nm = convert_to_nm(orig_val, orig_unit, mw)
        if val_nm is None or val_nm <= 0:
            excluded_invalid_structure_or_unit += 1
            continue
        pic50 = 9.0 - np.log10(val_nm)
        if pic50 < 2.0 or pic50 > 12.0:
            excluded_invalid_structure_or_unit += 1
            continue
        if orig_unit != "nM":
            converted_units_count += 1

        scaf = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
        r_dict = row.to_dict()
        r_dict["compound_id"] = str(row.get("molecule_chembl_id") or row.get("compound_id"))
        r_dict["molecule_chembl_id"] = r_dict["compound_id"]
        r_dict["target_id"] = str(row.get("target_chembl_id") or row.get("target_id") or target_id)
        r_dict["target_chembl_id"] = r_dict["target_id"]
        r_dict["target_organism"] = str(row.get("target_organism") or target_org)
        r_dict["assay_id"] = str(row.get("assay_chembl_id") or row.get("assay_id") or "Unknown")
        r_dict["assay_chembl_id"] = r_dict["assay_id"]
        r_dict["document_id"] = str(row.get("document_chembl_id") or row.get("document_id") or "Unknown")
        r_dict["document_chembl_id"] = r_dict["document_id"]
        r_dict["canonical_smiles"] = canon_smi
        r_dict["original_value"] = orig_val
        r_dict["original_units"] = orig_unit
        r_dict["standard_value"] = float(val_nm)
        r_dict["standard_units"] = "nM"
        r_dict["pIC50"] = float(pic50)
        r_dict["murcko_scaffold"] = scaf if scaf else canon_smi
        curated_rows.append(r_dict)

    df_curated = pd.DataFrame(curated_rows)
    after_unit_and_structure = int(len(df_curated))

    if df_curated.empty:
        raise ValueError("لم يتبق أي سجل صالح بعد تطبيق شروط التنظيف العلمية.")

    grp_span = df_curated.groupby("canonical_smiles")["pIC50"].transform(lambda s: s.max() - s.min())
    discordant_mask = grp_span > 1.0
    discordant_removed = int(discordant_mask.sum())
    df_curated = df_curated[~discordant_mask].copy()
    records_after_cleaning = int(len(df_curated))

    def _join_unique(series):
        vals = [str(v) for v in series.dropna().unique() if str(v).strip() and str(v) != "nan"]
        return "; ".join(vals) if vals else "غير محدد"

    def _build_source_ids(group_df):
        pairs = []
        for _, r in group_df.iterrows():
            a_id = str(r.get("assay_id", "Unknown"))
            d_id = str(r.get("document_id", "Unknown"))
            pair_str = f"{a_id} ({d_id})"
            if pair_str not in pairs:
                pairs.append(pair_str)
        return "; ".join(pairs)

    source_map = {
        smi: _build_source_ids(grp)
        for smi, grp in df_curated.groupby("canonical_smiles")
    }
    all_assays_map = {
        smi: _join_unique(grp["assay_id"])
        for smi, grp in df_curated.groupby("canonical_smiles")
    }
    all_docs_map = {
        smi: _join_unique(grp["document_id"])
        for smi, grp in df_curated.groupby("canonical_smiles")
    }

    agg_spec = {
        "compound_id": "first",
        "molecule_chembl_id": "first",
        "target_id": "first",
        "target_chembl_id": "first",
        "target_organism": "first",
        "target_pref_name": "first" if "target_pref_name" in df_curated.columns else "first",
        "assay_id": "first",
        "assay_chembl_id": "first",
        "assay_type": "first",
        "assay_description": "first" if "assay_description" in df_curated.columns else "first",
        "document_id": "first",
        "document_chembl_id": "first",
        "document_journal": "first" if "document_journal" in df_curated.columns else "first",
        "document_year": "first" if "document_year" in df_curated.columns else "first",
        "standard_type": "first",
        "standard_relation": "first",
        "standard_units": "first",
        "original_value": "first",
        "original_units": "first",
        "standard_value": "median",
        "pIC50": "median",
        "murcko_scaffold": "first",
    }
    agg_spec = {k: v for k, v in agg_spec.items() if k in df_curated.columns}

    counts = df_curated.groupby("canonical_smiles").size().rename("measurement_count")
    df_clean = df_curated.groupby("canonical_smiles", as_index=False).agg(agg_spec)
    df_clean = df_clean.merge(counts, on="canonical_smiles", how="left")
    df_clean["replicate_count"] = df_clean["measurement_count"]
    df_clean["median_activity"] = 10.0 ** (9.0 - df_clean["pIC50"])
    df_clean["standard_value"] = df_clean["median_activity"]
    df_clean["source_ids"] = df_clean["canonical_smiles"].map(source_map)
    df_clean["all_assay_ids"] = df_clean["canonical_smiles"].map(all_assays_map)
    df_clean["all_document_ids"] = df_clean["canonical_smiles"].map(all_docs_map)

    meas_dist = df_clean["measurement_count"].value_counts().sort_index().to_dict()
    meas_dist_str = {f"{int(k)}_measurements": int(v) for k, v in meas_dist.items()}

    report = {
        "target_id": target_id,
        "target_organism": target_org,
        "aggregation_method": "Median of compatible measurements (after excluding discordant replicates with >1.0 pIC50 span)",
        "initial_records": initial_count,
        "after_valid_numeric": after_valid_numeric,
        "after_measurement_type_ic50": after_measurement_type,
        "after_exact_relation": after_exact_relation,
        "after_data_validity": after_data_validity,
        "after_binding_assay_filter": after_binding_assay,
        "after_unit_and_structure_standardization": after_unit_and_structure,
        "records_after_cleaning": records_after_cleaning,
        "converted_non_nm_units": converted_units_count,
        "discordant_replicates_removed": discordant_removed,
        "final_unique_compounds": int(len(df_clean)),
        "unique_compounds_count": int(len(df_clean)),
        "unique_murcko_scaffolds": int(df_clean["murcko_scaffold"].nunique()),
        "unique_assays_count": int(df_curated["assay_id"].nunique()),
        "unique_documents_count": int(df_curated["document_id"].nunique()),
        "measurement_type_distribution": raw_type_dist,
        "measurements_per_compound_distribution": meas_dist_str,
        "exclusion_reasons": {
            "missing_or_nonpositive_value": excluded_missing_or_nonpositive,
            "wrong_target_or_non_human_organism": excluded_wrong_target_or_organism,
            "non_ic50_measurement_type": excluded_other_measurement_type,
            "censored_or_inequality_relation": excluded_censored_relations,
            "data_validity_warning_flagged": excluded_validity_flags,
            "non_binding_assay_type": excluded_non_binding_assays,
            "invalid_structure_or_out_of_range": excluded_invalid_structure_or_unit,
            "discordant_replicate_measurements": discordant_removed,
        },
    }
    return df_clean, report


def clean_and_featurize(
    input_path: str = RAW_DATA_PATH,
    output_path: str = CLEAN_DATA_PATH,
    report_path: str = CLEAN_REPORT_PATH,
) -> tuple[pd.DataFrame, dict]:
    if not os.path.exists(input_path):
        print("خطأ حرج: ملف البيانات الخام غير موجود.", file=sys.stderr)
        sys.exit(1)

    df_raw = pd.read_csv(input_path)
    if df_raw["canonical_smiles"].nunique() < 20:
        print("خطأ حرج: البيانات الخام غير صالحة أو تحتوي على سجلات اصطناعية.", file=sys.stderr)
        sys.exit(1)

    df_clean, report = clean_chembl_dataframe(df_raw, target_measurement_type="IC50")

    os.makedirs(os.path.dirname(output_path), exist_ok=True)
    df_clean.to_csv(output_path, index=False)
    with open(report_path, "w", encoding="utf-8") as f:
        json.dump(report, f, ensure_ascii=False, indent=2)

    print(f"اكتمل التنظيف العلمي للبيانات. عدد المركبات الفريدة الصالحة للتدريب: {len(df_clean)}")
    return df_clean, report


if __name__ == "__main__":
    clean_and_featurize()


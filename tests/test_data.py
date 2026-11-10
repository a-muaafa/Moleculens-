import json
import os
import pandas as pd
from config.settings import CLEAN_DATA_PATH, CLEAN_REPORT_PATH
from data.clean_data import clean_chembl_dataframe, convert_to_nm


def test_unit_conversion_to_nm():
    assert convert_to_nm(25.0, "nM", 300.0) == 25.0
    assert convert_to_nm(1.5, "uM", 300.0) == 1500.0
    ug_conv = convert_to_nm(0.3, "ug.mL-1", 300.0)
    assert abs(ug_conv - 1000.0) < 1e-6
    assert convert_to_nm(-5.0, "nM", 300.0) is None


def test_clean_chembl_dataframe_policy():
    sample_raw = pd.DataFrame(
        [
            {
                "molecule_chembl_id": "CHEMBL_A",
                "canonical_smiles": "N#C[C@@H]1CCCN1C(=O)CNC12CC3CC(CC(O)(C3)C1)C2",
                "target_chembl_id": "CHEMBL284",
                "target_organism": "Homo sapiens",
                "assay_chembl_id": "ASSAY_1",
                "document_chembl_id": "DOC_1",
                "assay_type": "B",
                "standard_type": "IC50",
                "standard_relation": "=",
                "standard_value": 10.0,
                "standard_units": "nM",
                "data_validity_comment": None,
            },
            {
                "molecule_chembl_id": "CHEMBL_A",
                "canonical_smiles": "N#C[C@@H]1CCCN1C(=O)CNC12CC3CC(CC(O)(C3)C1)C2",
                "target_chembl_id": "CHEMBL284",
                "target_organism": "Homo sapiens",
                "assay_chembl_id": "ASSAY_2",
                "document_chembl_id": "DOC_2",
                "assay_type": "B",
                "standard_type": "IC50",
                "standard_relation": "=",
                "standard_value": 20.0,
                "standard_units": "nM",
                "data_validity_comment": None,
            },
            {
                "molecule_chembl_id": "CHEMBL_KI",
                "canonical_smiles": "CC(=O)Oc1ccccc1C(=O)O",
                "target_chembl_id": "CHEMBL284",
                "target_organism": "Homo sapiens",
                "assay_chembl_id": "ASSAY_3",
                "document_chembl_id": "DOC_3",
                "assay_type": "B",
                "standard_type": "Ki",
                "standard_relation": "=",
                "standard_value": 50.0,
                "standard_units": "nM",
                "data_validity_comment": None,
            },
            {
                "molecule_chembl_id": "CHEMBL_RAT",
                "canonical_smiles": "CC(=O)Nc1ccc(O)cc1",
                "target_chembl_id": "CHEMBL284",
                "target_organism": "Rattus norvegicus",
                "assay_chembl_id": "ASSAY_4",
                "document_chembl_id": "DOC_4",
                "assay_type": "B",
                "standard_type": "IC50",
                "standard_relation": "=",
                "standard_value": 50.0,
                "standard_units": "nM",
                "data_validity_comment": None,
            },
            {
                "molecule_chembl_id": "CHEMBL_CENSORED",
                "canonical_smiles": "Cn1cnc2c1c(=O)n(C)c(=O)n2C",
                "target_chembl_id": "CHEMBL284",
                "target_organism": "Homo sapiens",
                "assay_chembl_id": "ASSAY_5",
                "document_chembl_id": "DOC_5",
                "assay_type": "B",
                "standard_type": "IC50",
                "standard_relation": ">",
                "standard_value": 10000.0,
                "standard_units": "nM",
                "data_validity_comment": "Outside typical range",
            },
        ]
    )

    df_clean, rep = clean_chembl_dataframe(sample_raw, target_measurement_type="IC50")
    assert len(df_clean) == 1
    row = df_clean.iloc[0]
    assert row["compound_id"] == "CHEMBL_A"
    assert row["measurement_count"] == 2
    assert 13.0 <= row["median_activity"] <= 15.0
    assert "ASSAY_1 (DOC_1)" in row["source_ids"]
    assert "ASSAY_2 (DOC_2)" in row["source_ids"]
    assert rep["exclusion_reasons"]["wrong_target_or_non_human_organism"] == 1
    assert rep["exclusion_reasons"]["non_ic50_measurement_type"] == 1
    assert rep["exclusion_reasons"]["censored_or_inequality_relation"] == 1


def test_processed_dataset_schema_and_report():
    assert os.path.exists(CLEAN_DATA_PATH)
    assert os.path.exists(CLEAN_REPORT_PATH)

    df = pd.read_csv(CLEAN_DATA_PATH)
    required_schema = {
        "compound_id",
        "canonical_smiles",
        "target_id",
        "target_organism",
        "assay_id",
        "document_id",
        "standard_type",
        "standard_value",
        "standard_units",
        "median_activity",
        "measurement_count",
        "source_ids",
    }
    assert required_schema.issubset(set(df.columns))
    assert (df["target_id"] == "CHEMBL284").all()
    assert (df["target_organism"] == "Homo sapiens").all()
    assert (df["standard_type"] == "IC50").all()
    assert (df["standard_units"] == "nM").all()

    with open(CLEAN_REPORT_PATH, "r", encoding="utf-8") as f:
        rep = json.load(f)
    assert rep["initial_records"] == 1000
    assert rep["final_unique_compounds"] == len(df)


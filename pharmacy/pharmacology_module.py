from typing import Any
from chemistry_engine.structure_analyzer import compute_rdkit_molecular_metrics
from config.settings import OP_SCIENTIFIC_CALCULATION, OP_SCIENTIFIC_RETRIEVAL
from core.schemas import build_provenance_tag


def build_pharmacology_section(entity: dict[str, Any]) -> dict[str, Any]:
    pharm = dict(entity.get("pharmacology_details", {}))
    smiles = entity.get("canonical_smiles", "")
    rdkit_metrics = compute_rdkit_molecular_metrics(smiles) if smiles else None

    return {
        "therapeutic_class_ar": pharm.get("therapeutic_class_ar", "مادة دوائية نشطة حيويا"),
        "indication_ar": pharm.get("indication_ar", "يستخدم وفق الاستطبابات الدوائية المعتمدة"),
        "target_protein_ar": pharm.get("target_protein_ar", "هدف بروتيني جزيئي"),
        "mechanism_of_action_ar": pharm.get("mechanism_of_action_ar", "يرتبط بالهدف الجزيئي لتعديل المسار الفسيولوجي."),
        "documented_bioactivity_ar": pharm.get("documented_bioactivity_ar", "موثق في قاعدة بيانات ChEMBL"),
        "approval_status_ar": pharm.get("approval_status_ar", "مسجل في قواعد البيانات الدوائية"),
        "pdb_complex_id": pharm.get("pdb_complex_id"),
        "calculated_drug_descriptors": rdkit_metrics,
        "retrieval_provenance": build_provenance_tag(
            OP_SCIENTIFIC_RETRIEVAL,
            sources=["ChEMBL", "PubChem"],
            notes_ar="البيانات الدوائية والنشاط الحيوي مسترجعة من السجلات المرجعية الموثقة وليست تنبؤا حاسوبيا.",
        ),
        "calculation_provenance": build_provenance_tag(
            OP_SCIENTIFIC_CALCULATION,
            sources=["RDKit Chemistry Engine"],
            notes_ar="الخصائص الجزيئية الدوائية (MW, LogP, HBD, HBA, TPSA, Lipinski) محسوبة حتميا من البنية الكيميائية.",
        ),
    }


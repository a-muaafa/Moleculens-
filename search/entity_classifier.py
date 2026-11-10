from typing import Any
from chemistry_engine.formula_analyzer import parse_chemical_formula
from chemistry_engine.periodic_table import VALID_ELEMENTS
from config.settings import (
    ENTITY_BIOMOLECULE_PROTEIN,
    ENTITY_CHEMICAL_COMPOUND,
    ENTITY_CHEMICAL_ELEMENT,
    ENTITY_DRUG_SUBSTANCE,
    ENTITY_TYPE_LABELS_AR,
    ENTITY_TYPE_LABELS_EN,
)


def classify_scientific_entity(record: dict[str, Any]) -> dict[str, str]:
    explicit_type = record.get("entity_type")
    if explicit_type in ENTITY_TYPE_LABELS_AR:
        return {
            "entity_type": explicit_type,
            "label_ar": ENTITY_TYPE_LABELS_AR[explicit_type],
            "label_en": ENTITY_TYPE_LABELS_EN[explicit_type],
            "classification_reason_ar": record.get(
                "classification_reason_ar",
                f"تم تصنيف الكيان كـ «{ENTITY_TYPE_LABELS_AR[explicit_type]}» بناء على خصائصه ومعرفاته المرجعية.",
            ),
        }

    ids = record.get("identifiers", {})
    if record.get("uniprot_id") or ids.get("uniprot_id") or record.get("pdb_ids") or ids.get("pdb_ids"):
        return {
            "entity_type": ENTITY_BIOMOLECULE_PROTEIN,
            "label_ar": ENTITY_TYPE_LABELS_AR[ENTITY_BIOMOLECULE_PROTEIN],
            "label_en": ENTITY_TYPE_LABELS_EN[ENTITY_BIOMOLECULE_PROTEIN],
            "classification_reason_ar": "يحتوي السجل على معرفات بروتينية/تركيبية مرجعية (UniProt / PDB).",
        }

    if record.get("pharmacology_details") or record.get("brand_names"):
        return {
            "entity_type": ENTITY_DRUG_SUBSTANCE,
            "label_ar": ENTITY_TYPE_LABELS_AR[ENTITY_DRUG_SUBSTANCE],
            "label_en": ENTITY_TYPE_LABELS_EN[ENTITY_DRUG_SUBSTANCE],
            "classification_reason_ar": "يمتلك السجل بيانات دوائية وعلاجية موثقة في ChEMBL/PubChem.",
        }

    symbol_or_formula = (record.get("symbol_or_formula") or record.get("molecular_formula") or "").strip()
    if record.get("atomic_number") or ids.get("atomic_number") or symbol_or_formula in VALID_ELEMENTS:
        parsed = parse_chemical_formula(symbol_or_formula)
        if parsed and len(parsed) == 1 and list(parsed.values())[0] == 1:
            return {
                "entity_type": ENTITY_CHEMICAL_ELEMENT,
                "label_ar": ENTITY_TYPE_LABELS_AR[ENTITY_CHEMICAL_ELEMENT],
                "label_en": ENTITY_TYPE_LABELS_EN[ENTITY_CHEMICAL_ELEMENT],
                "classification_reason_ar": f"يمثل الرمز ({symbol_or_formula}) عنصرا كيميائيا نقيا في الجدول الدوري.",
            }

    return {
        "entity_type": ENTITY_CHEMICAL_COMPOUND,
        "label_ar": ENTITY_TYPE_LABELS_AR[ENTITY_CHEMICAL_COMPOUND],
        "label_en": ENTITY_TYPE_LABELS_EN[ENTITY_CHEMICAL_COMPOUND],
        "classification_reason_ar": f"يتكون الكيان ({symbol_or_formula}) من ذرات مرتبطة كيميائيا بصيغة جزيئية محددة.",
    }


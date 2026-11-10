from dataclasses import asdict, dataclass, field
from typing import Any
from config.settings import (
    CHEMICAL_SUBCLASS_LABELS_AR,
    ENTITY_TYPE_LABELS_AR,
    ENTITY_TYPE_LABELS_EN,
    OP_SCIENTIFIC_CALCULATION,
    OP_SCIENTIFIC_PREDICTION,
    OP_SCIENTIFIC_RETRIEVAL,
    OPERATION_LABELS_AR,
)


@dataclass
class ProvenanceTag:
    operation_type: str
    operation_label_ar: str
    sources: list[str] = field(default_factory=list)
    notes_ar: str = ""

    def to_dict(self) -> dict[str, Any]:
        return asdict(self)


def build_provenance_tag(operation_type: str, sources: list[str], notes_ar: str = "") -> dict[str, Any]:
    label = OPERATION_LABELS_AR.get(operation_type, OPERATION_LABELS_AR[OP_SCIENTIFIC_RETRIEVAL])
    return ProvenanceTag(
        operation_type=operation_type,
        operation_label_ar=label,
        sources=list(dict.fromkeys(sources)),
        notes_ar=notes_ar,
    ).to_dict()


def build_unified_entity_schema(raw_data: dict[str, Any]) -> dict[str, Any]:
    entity_type = raw_data.get("entity_type", "chemical_compound")
    return {
        "key": raw_data.get("key", ""),
        "entity_type": entity_type,
        "entity_type_label_ar": ENTITY_TYPE_LABELS_AR.get(entity_type, "مادة علمية"),
        "entity_type_label_en": ENTITY_TYPE_LABELS_EN.get(entity_type, "Scientific Entity"),
        "subclass": raw_data.get("subclass", ""),
        "subclass_label_ar": CHEMICAL_SUBCLASS_LABELS_AR.get(raw_data.get("subclass", ""), ENTITY_TYPE_LABELS_AR.get(entity_type, "مادة علمية")),
        "name_ar": raw_data.get("name_ar", ""),
        "name_en": raw_data.get("name_en", ""),
        "iupac_name": raw_data.get("iupac_name", ""),
        "symbol_or_formula": raw_data.get("symbol_or_formula") or raw_data.get("molecular_formula") or raw_data.get("symbol", ""),
        "canonical_smiles": raw_data.get("canonical_smiles", ""),
        "molar_mass": raw_data.get("molar_mass"),
        "abbreviations": list(raw_data.get("abbreviations", [])),
        "brand_names": list(raw_data.get("brand_names", [])),
        "synonyms_ar": list(raw_data.get("synonyms_ar", [])),
        "synonyms_en": list(raw_data.get("synonyms_en", [])),
        "identifiers": {
            "atomic_number": raw_data.get("atomic_number"),
            "pubchem_cid": raw_data.get("pubchem_cid"),
            "chembl_id": raw_data.get("chembl_id"),
            "uniprot_id": raw_data.get("uniprot_id"),
            "pdb_ids": list(raw_data.get("pdb_ids", [])),
            "rhea_ids": list(raw_data.get("rhea_ids", [])),
        },
        "atomic_details": dict(raw_data.get("atomic_details", {})),
        "structure_details": dict(raw_data.get("structure_details", {})),
        "physical_properties": dict(raw_data.get("physical_properties", {})),
        "chemical_properties": dict(raw_data.get("chemical_properties", {})),
        "applications": dict(raw_data.get("applications", {})),
        "pharmacology_details": dict(raw_data.get("pharmacology_details", {})),
        "biology_details": dict(raw_data.get("biology_details", {})),
        "education_notes": dict(raw_data.get("education_notes", {})),
        "sources": list(raw_data.get("sources", ["Local Scientific Catalog"])),
        "matched_by": raw_data.get("matched_by", "exact"),
        "match_category": raw_data.get("match_category", "Exact Match"),
        "match_reason_ar": raw_data.get("match_reason_ar", ""),
        "match_score": float(raw_data.get("match_score", 100.0)),
    }


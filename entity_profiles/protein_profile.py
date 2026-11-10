from typing import Any
from biology.biomolecule_module import build_biology_section
from config.settings import ENTITY_BIOMOLECULE_PROTEIN, OP_SCIENTIFIC_RETRIEVAL
from core.schemas import build_provenance_tag
from education.explanation_engine import build_entity_explanation


def build_protein_scientific_profile(entity: dict[str, Any]) -> dict[str, Any]:
    struct = entity.get("structure_details", {})
    phys = entity.get("physical_properties", {})
    apps = entity.get("applications", {})
    bio_sec = build_biology_section(entity)
    explanation = build_entity_explanation(entity)
    mass = entity.get("molar_mass")

    return {
        "profile_kind": ENTITY_BIOMOLECULE_PROTEIN,
        "has_molecular_structure": True,
        "basic_info": {
            "name_en": entity.get("name_en", ""),
            "name_ar": entity.get("name_ar", ""),
            "iupac_name": entity.get("iupac_name", ""),
            "formula_or_symbol": entity.get("symbol_or_formula", ""),
            "molar_mass": f"{mass:,.1f} Da" if mass else phys.get("molar_mass_da", "غير محدد"),
            "entity_type_label_ar": entity.get("entity_type_label_ar", "بروتين أو جزيء حيوي"),
            "entity_type_label_en": entity.get("entity_type_label_en", "Protein / Biomolecule"),
        },
        "biology": bio_sec,
        "structure": {
            "molecular_geometry": struct.get("molecular_geometry", "بنية بروتينية ثلاثية الأبعاد"),
            "chemical_bonds": struct.get("chemical_bonds", "روابط ببتيدية وجسور ثنائية الكبريتيد وتفاعلات غير تساهمية"),
            "representation_3d_note": struct.get("representation_3d_note", "تتوفر البنية الفراغية ثلاثية الأبعاد في بنك معلومات البروتينات RCSB PDB."),
            "pdb_links": bio_sec.get("pdb_links", []),
        },
        "properties": {
            "physical_state": phys.get("physical_state", "جزيء حيوي في الأوساط الفسيولوجية"),
            "molar_mass_da": phys.get("molar_mass_da", f"{mass} Da" if mass else ""),
        },
        "applications": {
            "daily": apps.get("daily", ""),
            "industrial": apps.get("industrial", ""),
            "scientific": apps.get("scientific", ""),
        },
        "explanation": explanation,
        "identifiers": entity.get("identifiers", {}),
        "sources": entity.get("sources", ["UniProt", "PDB"]),
        "provenance": build_provenance_tag(
            OP_SCIENTIFIC_RETRIEVAL,
            sources=entity.get("sources", ["UniProt", "PDB", "Rhea"]),
            notes_ar="ملف جزيء حيوي/بروتين: يعرض بيانات التسلسل والوظيفة الحيوية وبنك معلومات البروتينات PDB.",
        ),
        "reaction_lab_entity": {
            "key": entity.get("key", ""),
            "name_ar": entity.get("name_ar", ""),
            "name_en": entity.get("name_en", ""),
            "symbol_or_formula": entity.get("symbol_or_formula", ""),
            "molar_mass": mass,
        },
    }


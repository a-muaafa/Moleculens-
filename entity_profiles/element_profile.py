from typing import Any
from config.settings import ENTITY_CHEMICAL_ELEMENT, OP_SCIENTIFIC_RETRIEVAL
from core.schemas import build_provenance_tag
from education.explanation_engine import build_entity_explanation


def build_element_scientific_profile(entity: dict[str, Any]) -> dict[str, Any]:
    atomic_det = entity.get("atomic_details", {})
    phys = entity.get("physical_properties", {})
    chem = entity.get("chemical_properties", {})
    apps = entity.get("applications", {})
    z = entity.get("identifiers", {}).get("atomic_number") or atomic_det.get("atomic_number")
    symbol = entity.get("symbol_or_formula", "")
    mass = entity.get("molar_mass")
    explanation = build_entity_explanation(entity)

    return {
        "profile_kind": ENTITY_CHEMICAL_ELEMENT,
        "has_molecular_structure": False,
        "basic_info": {
            "name_en": entity.get("name_en", ""),
            "name_ar": entity.get("name_ar", ""),
            "symbol": symbol,
            "atomic_number": z,
            "atomic_mass": f"{mass} u (g/mol)" if mass else atomic_det.get("atomic_mass", ""),
            "entity_type_label_ar": entity.get("entity_type_label_ar", "عنصر كيميائي"),
            "entity_type_label_en": entity.get("entity_type_label_en", "Chemical Element"),
        },
        "atomic_structure": {
            "protons_count": z,
            "electrons_count": z,
            "electron_configuration": atomic_det.get("electron_configuration", ""),
            "electron_shells": atomic_det.get("shells", ""),
            "group": atomic_det.get("group", ""),
            "period": atomic_det.get("period", ""),
            "block": atomic_det.get("block", ""),
            "oxidation_states": atomic_det.get("oxidation_states", ""),
            "crystal_structure": atomic_det.get("crystal_structure", ""),
        },
        "properties": {
            "physical_state": phys.get("physical_state", "صلب في الظروف القياسية"),
            "density": phys.get("density", "غير محدد"),
            "melting_point": phys.get("melting_point", "غير محدد"),
            "boiling_point": phys.get("boiling_point", "غير محدد"),
            "chemical_behavior": chem.get("reactivity", ""),
        },
        "applications": {
            "daily": apps.get("daily", ""),
            "industrial": apps.get("industrial", ""),
            "scientific": apps.get("scientific", ""),
        },
        "explanation": explanation,
        "identifiers": entity.get("identifiers", {}),
        "sources": entity.get("sources", ["IUPAC Periodic Table"]),
        "provenance": build_provenance_tag(
            OP_SCIENTIFIC_RETRIEVAL,
            sources=entity.get("sources", ["IUPAC Periodic Table"]),
            notes_ar="ملف عنصر كيميائي: يعرض الخصائص الذرية والتوزيع الإلكتروني ولا يتضمن بنية جزيئية تساهمية.",
        ),
        "reaction_lab_entity": {
            "key": entity.get("key", ""),
            "name_ar": entity.get("name_ar", ""),
            "name_en": entity.get("name_en", ""),
            "symbol_or_formula": symbol,
            "molar_mass": mass,
        },
    }


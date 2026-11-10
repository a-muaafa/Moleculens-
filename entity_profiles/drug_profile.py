from typing import Any
from config.settings import ENTITY_DRUG_SUBSTANCE
from entity_profiles.compound_profile import build_compound_scientific_profile
from pharmacy.pharmacology_module import build_pharmacology_section


def build_drug_scientific_profile(entity: dict[str, Any]) -> dict[str, Any]:
    base_profile = build_compound_scientific_profile(entity)
    base_profile["profile_kind"] = ENTITY_DRUG_SUBSTANCE
    base_profile["pharmacology"] = build_pharmacology_section(entity)
    base_profile["brand_names"] = entity.get("brand_names", [])
    return base_profile


from typing import Any
from config.settings import (
    ENTITY_BIOMOLECULE_PROTEIN,
    ENTITY_CHEMICAL_COMPOUND,
    ENTITY_CHEMICAL_ELEMENT,
    ENTITY_DRUG_SUBSTANCE,
)
from entity_profiles.compound_profile import build_compound_scientific_profile
from entity_profiles.drug_profile import build_drug_scientific_profile
from entity_profiles.element_profile import build_element_scientific_profile
from entity_profiles.protein_profile import build_protein_scientific_profile


def build_scientific_profile(entity: dict[str, Any]) -> dict[str, Any]:
    entity_type = entity.get("entity_type", ENTITY_CHEMICAL_COMPOUND)
    if entity_type == ENTITY_CHEMICAL_ELEMENT:
        return build_element_scientific_profile(entity)
    if entity_type == ENTITY_DRUG_SUBSTANCE:
        return build_drug_scientific_profile(entity)
    if entity_type == ENTITY_BIOMOLECULE_PROTEIN:
        return build_protein_scientific_profile(entity)
    return build_compound_scientific_profile(entity)


from typing import Any
from entity_profiles.profile_builder import build_scientific_profile
from reaction_engine.prediction import predict_reaction
from reaction_engine.reaction_service import analyze_chemical_reaction
from search.search_engine import (
    get_catalog_dropdown_options,
    get_live_suggestions,
    search_scientific_entity,
)


class MolecuLensCore:
    @staticmethod
    def explore_entity(
        query: str,
        search_type: str = "auto",
        allow_network: bool = True,
        selected_candidate: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        if selected_candidate is not None:
            profile = build_scientific_profile(selected_candidate)
            return {
                "found": True,
                "is_ambiguous": False,
                "ambiguity_message": None,
                "entity": selected_candidate,
                "profile": profile,
                "candidates": [selected_candidate],
                "did_you_mean": [],
                "message": selected_candidate.get("match_reason_ar", "تم اختيار المادة بنجاح."),
            }

        search_res = search_scientific_entity(query, search_type=search_type, allow_network=allow_network)
        if not search_res["found"] or search_res["is_ambiguous"]:
            return {
                "found": search_res["found"],
                "is_ambiguous": search_res["is_ambiguous"],
                "ambiguity_message": search_res["ambiguity_message"],
                "entity": None,
                "profile": None,
                "candidates": search_res["candidates"],
                "did_you_mean": search_res["did_you_mean"],
                "message": search_res["message"],
            }

        entity = search_res["selected"]
        profile = build_scientific_profile(entity)
        return {
            "found": True,
            "is_ambiguous": False,
            "ambiguity_message": None,
            "entity": entity,
            "profile": profile,
            "candidates": search_res["candidates"],
            "did_you_mean": [],
            "message": search_res["message"],
        }

    @staticmethod
    def run_reaction_analysis(
        reactant_inputs: list[dict[str, Any]],
        domain: str = "general_chemistry",
        custom_conditions: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        return analyze_chemical_reaction(
            reactant_inputs=reactant_inputs,
            domain=domain,
            custom_conditions=custom_conditions,
        )

    @staticmethod
    def run_reaction_prediction(
        reactant_inputs: list[dict[str, Any]],
        domain: str | None = None,
        custom_conditions: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        """Explores whether the supplied substances will react, and what the products would be."""
        return predict_reaction(
            reactant_inputs=reactant_inputs,
            domain=domain,
            custom_conditions=custom_conditions,
        )

    @staticmethod
    def fetch_suggestions(query: str, limit: int = 6) -> list[dict[str, Any]]:
        return get_live_suggestions(query, limit=limit)

    @staticmethod
    def fetch_catalog_options() -> list[dict[str, Any]]:
        return get_catalog_dropdown_options()


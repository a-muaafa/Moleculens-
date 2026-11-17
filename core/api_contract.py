"""
Unified Backend-Frontend API Contract Layer for MolecuLens.

Establishes a strict, standardized contract between the scientific backend and any
frontend client (Streamlit UI, React/Vue SPA, or REST consumer):
1. Standardized JSON response envelope (ok, status_code, ui_state, endpoint,
   message_ar, message_en, units, data, error, meta).
2. Explicit HTTP status codes:
   - 200 OK: Successful retrieval or calculation.
   - 400 Bad Request: Malformed, empty, oversized, or unsafe input.
   - 404 Not Found: Entity or record not found (ui_state = "empty").
   - 422 Unprocessable Entity: Scientifically invalid request (e.g., L_STP on solid/liquid,
     or insufficient evidence to predict a reaction without fabricating products).
   - 503 Service Unavailable: External service failure when strictly required.
3. Canonical field names and scientific units across all endpoints.
4. Explicit UI state messages (loading, empty, error, ready).
5. Zero scientific calculation logic in the frontend: all stoichiometry, QSAR,
   formula resolution, and visual experiment stages are computed in the backend.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any

from chemistry_engine.molecular_intelligence import MolecularIntelligence
from core.security import validate_safe_chemical_input
from core.validators import validate_reaction_inputs, validate_search_query
from entity_profiles.profile_builder import build_scientific_profile
from reaction_engine.prediction import (
    parse_and_resolve_reactant_text,
    predict_reaction,
)
from reaction_engine.reaction_service import analyze_chemical_reaction
from reaction_engine.scientific_audit import audit_registry_and_rules
from reaction_engine.visualizer import ReactionVisualizer
from search.search_engine import search_scientific_entity

API_VERSION = "1.0.0"

CANONICAL_UNITS: dict[str, str] = {
    "molar_mass": "g/mol",
    "amount": "mol",
    "mass": "g",
    "volume_stp": "L_STP",
    "temperature": "°C",
    "pressure": "atm",
    "enthalpy": "kJ/mol",
    "tpsa": "Å²",
    "concentration_molar": "M",
    "volume_ml": "mL",
}

UI_STATE_MESSAGES: dict[str, dict[str, str]] = {
    "loading": {
        "ar": "جارٍ معالجة البيانات العلمية وحساب النتائج في الخادم...",
        "en": "Processing scientific data and computing results in the backend...",
    },
    "empty": {
        "ar": "لا توجد نتائج مطابقة للمدخلات الحالية. جرّب البحث بالاسم العلمي أو الصيغة الكيميائية.",
        "en": "No matching records found for the current input.",
    },
    "error": {
        "ar": "تعذر إتمام العملية بسبب خطأ في المدخلات أو عدم كفاية الأدلة العلمية.",
        "en": "Operation could not be completed due to invalid input or insufficient scientific evidence.",
    },
    "ready": {
        "ar": "تم إجراء الحسابات العلمية واسترجاع البيانات الموثقة بنجاح.",
        "en": "Scientific calculations and verified data retrieved successfully.",
    },
}


def make_api_response(
    *,
    ok: bool,
    status_code: int,
    ui_state: str,
    endpoint: str,
    data: dict[str, Any] | None = None,
    message_ar: str = "",
    message_en: str = "",
    error_code: str | None = None,
    error_details_ar: str = "",
    suggestions_ar: list[str] | None = None,
    cache_hit: bool = False,
) -> dict[str, Any]:
    """Builds the unified MolecuLens API response envelope."""
    default_msgs = UI_STATE_MESSAGES.get(ui_state, UI_STATE_MESSAGES["ready"])
    final_msg_ar = message_ar or default_msgs["ar"]
    final_msg_en = message_en or default_msgs["en"]

    error_obj = None
    if not ok:
        error_obj = {
            "code": error_code or f"HTTP_{status_code}",
            "details_ar": error_details_ar or final_msg_ar,
            "suggestions_ar": list(suggestions_ar or []),
        }

    return {
        "ok": bool(ok),
        "status_code": int(status_code),
        "ui_state": ui_state,
        "endpoint": endpoint,
        "message_ar": final_msg_ar,
        "message_en": final_msg_en,
        "units": dict(CANONICAL_UNITS),
        "data": data,
        "error": error_obj,
        "meta": {
            "api_version": API_VERSION,
            "timestamp_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
            "computed_in_backend": True,
            "cache_hit": bool(cache_hit),
        },
    }


class MolecuLensAPI:
    """Backend service façade returning standardized API responses for all platform modules."""

    @staticmethod
    def explore_entity(
        query: str,
        search_type: str = "auto",
        allow_network: bool = True,
        force_network: bool = False,
    ) -> dict[str, Any]:
        endpoint = "/api/v1/explore"
        valid, cleaned, err_msg = validate_search_query(query)
        if not valid:
            return make_api_response(
                ok=False,
                status_code=400,
                ui_state="error",
                endpoint=endpoint,
                message_ar=err_msg,
                message_en="Invalid search query.",
                error_code="INVALID_SEARCH_QUERY",
                error_details_ar=err_msg,
                suggestions_ar=["أدخل اسم المادة بالعربية أو الإنجليزية، أو صيغتها الجزيئية مثل H2O أو NaCl."],
            )

        search_res = search_scientific_entity(
            cleaned,
            search_type=search_type,
            allow_network=allow_network,
            force_network=force_network,
        )
        if not search_res.get("found"):
            return make_api_response(
                ok=False,
                status_code=404,
                ui_state="empty",
                endpoint=endpoint,
                message_ar=search_res.get("message") or UI_STATE_MESSAGES["empty"]["ar"],
                message_en="Entity not found in local catalog or external databases.",
                data={
                    "query": cleaned,
                    "candidates": search_res.get("candidates", []),
                    "did_you_mean": search_res.get("did_you_mean", []),
                    "network_deferred": search_res.get("network_deferred", False),
                },
                error_code="ENTITY_NOT_FOUND",
                error_details_ar=search_res.get("message", ""),
                suggestions_ar=[
                    (c if isinstance(c, str) else str(c.get("name_ar", "")))
                    for c in search_res.get("did_you_mean", [])[:3]
                    if c
                ],
            )

        selected = search_res.get("selected")
        profile = build_scientific_profile(selected) if selected else None
        return make_api_response(
            ok=True,
            status_code=200,
            ui_state="ready",
            endpoint=endpoint,
            message_ar=search_res.get("message", UI_STATE_MESSAGES["ready"]["ar"]),
            message_en="Entity profile retrieved successfully.",
            data={
                "entity": selected,
                "profile": profile,
                "is_ambiguous": bool(search_res.get("is_ambiguous")),
                "ambiguity_message": search_res.get("ambiguity_message"),
                "candidates": search_res.get("candidates", []),
                "network_deferred": search_res.get("network_deferred", False),
            },
        )

    @staticmethod
    def predict_reaction(
        reactants_input: str | list[dict[str, Any]],
        domain: str | None = None,
        custom_conditions: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        endpoint = "/api/v1/reactions/predict"
        resolution_notes: list[str] = []

        if isinstance(reactants_input, str):
            ok_in, cleaned_text, err_in = validate_safe_chemical_input(
                reactants_input, max_len=300, field_name_ar="مدخل المتفاعلات"
            )
            if not ok_in:
                return make_api_response(
                    ok=False,
                    status_code=400,
                    ui_state="error",
                    endpoint=endpoint,
                    message_ar=err_in,
                    message_en="Invalid reactant input string.",
                    error_code="INVALID_REACTANT_INPUT",
                    error_details_ar=err_in,
                )
            resolved_items, resolution_notes = parse_and_resolve_reactant_text(cleaned_text)
        elif isinstance(reactants_input, list):
            resolved_items = reactants_input
        else:
            return make_api_response(
                ok=False,
                status_code=400,
                ui_state="error",
                endpoint=endpoint,
                message_ar="مدخل المتفاعلات غير صالح.",
                message_en="Invalid reactants payload type.",
                error_code="INVALID_PAYLOAD_TYPE",
            )

        if len(resolved_items) < 2:
            return make_api_response(
                ok=False,
                status_code=400,
                ui_state="error",
                endpoint=endpoint,
                message_ar="يلزم إدخال مادتين متفاعلتين صالحتين على الأقل لإجراء التوقع الكيميائي.",
                message_en="At least two valid reactants are required for reaction prediction.",
                data={"resolution_notes_ar": resolution_notes},
                error_code="INSUFFICIENT_REACTANTS",
                error_details_ar="تم التعرف على أقل من مادتين متفاعلتين.",
                suggestions_ar=resolution_notes or ["أدخل مادتين مفصولتين بعلامة + مثل: HCl + NaOH"],
            )

        pred = predict_reaction(
            reactant_inputs=resolved_items,
            domain=domain,
            custom_conditions=custom_conditions,
        )
        pred["resolution_notes_ar"] = resolution_notes

        if not pred.get("supported"):
            return make_api_response(
                ok=False,
                status_code=422,
                ui_state="empty",
                endpoint=endpoint,
                message_ar=pred.get("verdict_ar", "لا تتوفر أدلة علمية كافية لتنبؤ نواتج هذه المواد."),
                message_en="Insufficient scientific evidence to predict reaction products; no fabricated products returned.",
                data=pred,
                error_code="INSUFFICIENT_REACTION_EVIDENCE",
                error_details_ar=pred.get("explanation_ar", ""),
                suggestions_ar=pred.get("missing_information_ar", []),
            )

        pred["visual_experiment_stages"] = ReactionVisualizer.build_end_to_end_experiment_stages(pred)
        return make_api_response(
            ok=True,
            status_code=200,
            ui_state="ready",
            endpoint=endpoint,
            message_ar=pred.get("verdict_ar", UI_STATE_MESSAGES["ready"]["ar"]),
            message_en="Reaction prediction computed successfully.",
            data=pred,
        )

    @staticmethod
    def analyze_reaction(
        reactant_inputs: list[dict[str, Any]],
        domain: str = "general_chemistry",
        custom_conditions: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        endpoint = "/api/v1/reactions/analyze"
        valid, err_msg = validate_reaction_inputs(reactant_inputs)
        if not valid:
            status_code = 422 if "L_STP" in err_msg else 400
            err_code = "INVALID_STP_PHASE" if "L_STP" in err_msg else "INVALID_REACTION_INPUT"
            return make_api_response(
                ok=False,
                status_code=status_code,
                ui_state="error",
                endpoint=endpoint,
                message_ar=err_msg,
                message_en="Reaction input validation failed.",
                error_code=err_code,
                error_details_ar=err_msg,
            )

        rxn_res = analyze_chemical_reaction(
            reactant_inputs=reactant_inputs,
            domain=domain,
            custom_conditions=custom_conditions,
        )
        if not rxn_res.get("supported"):
            return make_api_response(
                ok=False,
                status_code=422,
                ui_state="empty",
                endpoint=endpoint,
                message_ar=rxn_res.get("rejection_reason_ar", "التفاعل غير مدعوم في السجل الموثق."),
                message_en="Reaction not supported in verified registry.",
                data=rxn_res,
                error_code="UNSUPPORTED_DOCUMENTED_REACTION",
                error_details_ar=rxn_res.get("rejection_reason_ar", ""),
            )

        rxn_res["visual_experiment_stages"] = ReactionVisualizer.build_end_to_end_experiment_stages(rxn_res)
        return make_api_response(
            ok=True,
            status_code=200,
            ui_state="ready",
            endpoint=endpoint,
            message_ar=f"تم تحليل التفاعل الموثق ({rxn_res.get('reaction_name_ar', '')}) بنجاح.",
            message_en="Verified reaction analyzed successfully.",
            data=rxn_res,
        )

    @staticmethod
    def evaluate_logp(smiles: str) -> dict[str, Any]:
        endpoint = "/api/v1/molecular-intelligence/logp"
        ok_in, clean_smi, err_in = validate_safe_chemical_input(smiles, max_len=500, field_name_ar="صيغة SMILES")
        if not ok_in:
            return make_api_response(
                ok=False,
                status_code=400,
                ui_state="error",
                endpoint=endpoint,
                message_ar=err_in,
                message_en="Invalid SMILES input.",
                error_code="INVALID_SMILES",
                error_details_ar=err_in,
            )

        res = MolecularIntelligence.evaluate_logp_prediction(clean_smi)
        if not res.get("available"):
            return make_api_response(
                ok=False,
                status_code=400,
                ui_state="error",
                endpoint=endpoint,
                message_ar=res.get("message_ar", "تعذر حساب معامل LogP لهذا المركب."),
                message_en="Could not compute LogP for the provided SMILES.",
                error_code="LOGP_COMPUTATION_FAILED",
                error_details_ar=res.get("message_ar", ""),
            )

        return make_api_response(
            ok=True,
            status_code=200,
            ui_state="ready",
            endpoint=endpoint,
            message_ar="تم حساب cLogP الحتمي وتقييم نموذج Ridge QSAR بنجاح.",
            message_en="Deterministic cLogP and Ridge QSAR evaluation completed.",
            data=res,
        )

    @staticmethod
    def modify_structure(smiles: str, modification_key: str) -> dict[str, Any]:
        endpoint = "/api/v1/molecular-intelligence/modify"
        ok_in, clean_smi, err_in = validate_safe_chemical_input(smiles, max_len=500, field_name_ar="صيغة SMILES")
        if not ok_in:
            return make_api_response(
                ok=False,
                status_code=400,
                ui_state="error",
                endpoint=endpoint,
                message_ar=err_in,
                message_en="Invalid SMILES input.",
                error_code="INVALID_SMILES",
                error_details_ar=err_in,
            )

        res = MolecularIntelligence.modify_structure(clean_smi, modification_key)
        if not res.get("success"):
            err_ar = res.get("error_ar", "تعذر تطبيق التعديل البنيوي.")
            status_code = 400 if "غير صالحة" in err_ar or "غير معروف" in err_ar else 422
            return make_api_response(
                ok=False,
                status_code=status_code,
                ui_state="error",
                endpoint=endpoint,
                message_ar=err_ar,
                message_en="Structure modification could not be applied.",
                error_code="MODIFICATION_REJECTED",
                error_details_ar=err_ar,
            )

        return make_api_response(
            ok=True,
            status_code=200,
            ui_state="ready",
            endpoint=endpoint,
            message_ar="تم تنفيذ التعديل البنيوي الافتراضي والتحقق من سلامة التكافؤ الكيميائي.",
            message_en="Virtual structural modification completed and chemically validated.",
            data=res,
        )

    @staticmethod
    def scientific_audit() -> dict[str, Any]:
        endpoint = "/api/v1/reactions/audit"
        report = audit_registry_and_rules()
        return make_api_response(
            ok=report["overall_passed"],
            status_code=200 if report["overall_passed"] else 500,
            ui_state="ready" if report["overall_passed"] else "error",
            endpoint=endpoint,
            message_ar="اكتمل التدقيق العلمي المستقل لسجل التفاعلات والقواعد الكيميائية.",
            message_en="Independent scientific audit of reactions and rules completed.",
            data=report,
        )

    @staticmethod
    def simulate_titration(
        acid_key: str,
        base_key: str,
        acid_conc_m: float = 0.1,
        acid_volume_ml: float = 25.0,
        base_conc_m: float = 0.1,
        titration_mode: str = "acid_with_base",
    ) -> dict[str, Any]:
        endpoint = "/api/v1/titration/simulate"
        from titration.titration_engine import TitrationEngine

        try:
            sim = TitrationEngine.simulate_titration(
                acid_key=acid_key,
                base_key=base_key,
                acid_conc_m=acid_conc_m,
                acid_volume_ml=acid_volume_ml,
                base_conc_m=base_conc_m,
                titration_mode=titration_mode,
            )
            return make_api_response(
                ok=True,
                status_code=200,
                ui_state="ready",
                endpoint=endpoint,
                message_ar="تمت محاكاة المعايرة وحساب نقطة التكافؤ ونقطة النهاية بنجاح.",
                message_en="Titration simulation completed successfully.",
                data=sim,
            )
        except ValueError as exc:
            return make_api_response(
                ok=False,
                status_code=400,
                ui_state="error",
                endpoint=endpoint,
                message_ar=str(exc),
                message_en=str(exc),
                error_code="INVALID_TITRATION_PARAMETERS",
                error_details_ar=str(exc),
            )

    @staticmethod
    def pharmacology_overview() -> dict[str, Any]:
        endpoint = "/api/v1/pharmacology/overview"
        from frontend.pharmacology_view import PHARMACOLOGY_PILLARS

        return make_api_response(
            ok=True,
            status_code=200,
            ui_state="ready",
            endpoint=endpoint,
            message_ar="قسم علم الأدوية (Pharmacology) في حالة Coming Soon – قيد التطوير.",
            message_en="Pharmacology module is in Coming Soon status.",
            data={
                "status": "coming_soon",
                "status_label_ar": "Coming Soon – قيد التطوير",
                "pillars": PHARMACOLOGY_PILLARS,
            },
        )

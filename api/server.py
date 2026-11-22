"""
Real HTTP REST API server for MolecuLens (FastAPI).

Enforces the layered backend architecture:
  HTTP Request -> Pydantic Validation -> MolecuLensAPI Service -> Scientific Engine -> Repository -> Response

Run locally:
    uvicorn api.server:app --host 0.0.0.0 --port 8000
"""

from __future__ import annotations

from typing import Any
import base64


def _json_safe(obj: Any) -> Any:
    if isinstance(obj, bytes):
        return base64.b64encode(obj).decode("ascii")
    if isinstance(obj, dict):
        return {str(k): _json_safe(v) for k, v in obj.items()}
    if isinstance(obj, (list, tuple)):
        return [_json_safe(v) for v in obj]
    return obj
from fastapi import FastAPI, HTTPException, Query
from pydantic import BaseModel, Field

from config.settings import ALLOW_EXTERNAL_NETWORK, SCIENTIFIC_DISCLAIMER_AR
from core.api_contract import MolecuLensAPI
from chemistry_engine.structure_analyzer import get_cheminformatics_engine_info
from titration.titration_engine import TitrationEngine, solve_equilibrium_ph_exact

app = FastAPI(
    title="MolecuLens Scientific REST API",
    version="3.3.0",
    description="Layered HTTP REST API for MolecuLens chemical exploration, reaction networks, titration equilibrium, and QSAR evaluation.",
)


class ReactantItem(BaseModel):
    formula: str = Field(..., min_length=1, max_length=120)
    amount: float | None = Field(default=None, gt=0)
    unit: str | None = Field(default=None, max_length=20)


class ReactionPredictRequest(BaseModel):
    reactants: list[ReactantItem] = Field(...)
    domain: str = Field(default="all_domains", max_length=64)


class TitrationRequest(BaseModel):
    acid_key: str = Field(..., min_length=1, max_length=40)
    base_key: str = Field(..., min_length=1, max_length=40)
    acid_conc_m: float = Field(..., gt=0, le=20.0)
    acid_volume_ml: float = Field(..., gt=0, le=10000.0)
    base_conc_m: float = Field(..., gt=0, le=20.0)
    titration_mode: str = Field(default="acid_with_base")
    temperature_c: float = Field(default=25.0, ge=0.0, le=60.0)
    use_activity_correction: bool = Field(default=False)


class LogPRequest(BaseModel):
    smiles: str = Field(..., min_length=1, max_length=500)


@app.get("/api/v1/health")
def health_check() -> dict[str, Any]:
    return {
        "ok": True,
        "service": "MolecuLens REST API",
        "version": "3.3.0",
        "offline_first_enabled": not ALLOW_EXTERNAL_NETWORK,
        "cheminformatics_engine": get_cheminformatics_engine_info(),
        "disclaimer_ar": SCIENTIFIC_DISCLAIMER_AR,
    }


@app.get("/api/v1/explore")
def explore_substance(
    q: str = Query(..., min_length=1, max_length=200),
    allow_network: bool = Query(default=False),
) -> dict[str, Any]:
    res = MolecuLensAPI.explore_entity(q, allow_network=allow_network)
    if not res.get("ok") and res.get("status_code") in {400, 404, 422}:
        raise HTTPException(status_code=int(res["status_code"]), detail=res.get("message_ar", "Invalid request"))
    return _json_safe(res)


@app.post("/api/v1/reactions/predict")
def predict_reaction_endpoint(body: ReactionPredictRequest) -> dict[str, Any]:
    if len(body.reactants) < 2 or len(body.reactants) > 6:
        raise HTTPException(status_code=422, detail='reactants must contain between 2 and 6 items')
    raw_items = [
        {"substance": {"symbol_or_formula": r.formula, "name_en": r.formula, "name_ar": r.formula}}
        for r in body.reactants
    ]
    from reaction_engine.prediction import predict_reaction

    pred = predict_reaction(raw_items, domain=body.domain)
    return {
        "ok": True,
        "status_code": 200,
        "endpoint": "/api/v1/reactions/predict",
        "domain": body.domain,
        "data": pred,
    }


@app.post("/api/v1/titration/simulate")
def simulate_titration_endpoint(body: TitrationRequest) -> dict[str, Any]:
    if body.acid_key not in TitrationEngine.ACID_DATABASE:
        raise HTTPException(status_code=400, detail=f"Unknown acid_key: {body.acid_key}")
    if body.base_key not in TitrationEngine.BASE_DATABASE:
        raise HTTPException(status_code=400, detail=f"Unknown base_key: {body.base_key}")
    sim = TitrationEngine.simulate_titration(
        acid_key=body.acid_key,
        base_key=body.base_key,
        acid_conc_m=body.acid_conc_m,
        acid_volume_ml=body.acid_volume_ml,
        base_conc_m=body.base_conc_m,
        titration_mode=body.titration_mode,
    )
    acid_spec = TitrationEngine.ACID_DATABASE[body.acid_key]
    base_spec = TitrationEngine.BASE_DATABASE[body.base_key]
    exact_initial_ph = solve_equilibrium_ph_exact(
        c_acid_m=body.acid_conc_m if body.titration_mode == "acid_with_base" else 0.0,
        c_base_m=0.0 if body.titration_mode == "acid_with_base" else body.base_conc_m,
        acid_spec=acid_spec,
        base_spec=base_spec,
        temperature_c=body.temperature_c,
        use_activity_correction=body.use_activity_correction,
    )
    return {
        "ok": True,
        "status_code": 200,
        "endpoint": "/api/v1/titration/simulate",
        "exact_initial_ph": exact_initial_ph,
        "temperature_c": body.temperature_c,
        "use_activity_correction": body.use_activity_correction,
        "data": sim,
    }


@app.post("/api/v1/molecular/logp")
def evaluate_logp_endpoint(body: LogPRequest) -> dict[str, Any]:
    res = MolecuLensAPI.evaluate_logp(body.smiles)
    if not res.get("ok"):
        raise HTTPException(status_code=int(res.get("status_code", 400)), detail=res.get("message_ar", "Invalid SMILES"))
    return _json_safe(res)


@app.get("/api/v1/audit/scientific")
def scientific_audit_endpoint() -> dict[str, Any]:
    return _json_safe(MolecuLensAPI.scientific_audit())

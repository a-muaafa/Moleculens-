import requests
from chemistry_engine.structure_analyzer import compute_inchikey, validate_and_parse_smiles
from config.settings import CHEMBL_BASE_URL, EXTERNAL_API_TIMEOUT_SECONDS
from data_integrations.cache_store import get_cached_item, set_cached_item


def search_chembl_molecules(query: str, limit: int = 3, timeout: int = EXTERNAL_API_TIMEOUT_SECONDS) -> dict:
    clean_q = (query or "").strip()
    if not clean_q:
        return {"ok": False, "records": [], "error": "Empty query"}
    cache_key = f"chembl:search:{clean_q.lower()}"
    hit, cached = get_cached_item(cache_key)
    if hit:
        return cached
    try:
        resp = requests.get(
            f"{CHEMBL_BASE_URL}/molecule/search.json",
            params={"q": clean_q, "limit": limit},
            timeout=timeout,
        )
        if resp.status_code != 200:
            return {"ok": False, "records": [], "error": f"ChEMBL HTTP {resp.status_code}"}
        records = []
        for m in resp.json().get("molecules", []):
            structs = m.get("molecule_structures") or {}
            props = m.get("molecule_properties") or {}
            smiles = structs.get("canonical_smiles")
            if not smiles:
                continue
            _, canon_smi, _ = validate_and_parse_smiles(smiles)
            final_smi = canon_smi or smiles
            inchikey = structs.get("standard_inchi_key") or compute_inchikey(final_smi) or ""
            records.append(
                {
                    "chembl_id": m.get("molecule_chembl_id"),
                    "name_en": m.get("pref_name") or clean_q,
                    "canonical_smiles": final_smi,
                    "inchikey": inchikey,
                    "symbol_or_formula": props.get("full_molformula", ""),
                    "molar_mass": float(props["full_mwt"]) if props.get("full_mwt") else None,
                    "max_phase": m.get("max_phase"),
                    "source": "ChEMBL REST API",
                }
            )
        res = {"ok": bool(records), "records": records, "error": None if records else "No ChEMBL molecules"}
        set_cached_item(cache_key, res)
        return res
    except Exception as exc:
        return {"ok": False, "records": [], "error": str(exc)}


def fetch_chembl_drug_info(chembl_id: str, timeout: int = EXTERNAL_API_TIMEOUT_SECONDS) -> dict:
    if not chembl_id:
        return {"ok": False, "data": {}, "error": "Missing ChEMBL ID"}
    clean_id = chembl_id.strip().upper()
    cache_key = f"chembl:molecule:{clean_id}"
    hit, cached = get_cached_item(cache_key)
    if hit:
        return cached

    url = f"{CHEMBL_BASE_URL}/molecule/{clean_id}.json"
    try:
        resp = requests.get(url, timeout=timeout)
        if resp.status_code != 200:
            return {"ok": False, "data": {}, "error": f"ChEMBL HTTP {resp.status_code}"}
        payload = resp.json()
        structs = payload.get("molecule_structures") or {}
        data = {
            "chembl_id": payload.get("molecule_chembl_id", clean_id),
            "pref_name": payload.get("pref_name", ""),
            "max_phase": payload.get("max_phase"),
            "molecule_type": payload.get("molecule_type", ""),
            "first_approval": payload.get("first_approval"),
            "canonical_smiles": structs.get("canonical_smiles", ""),
            "inchikey": structs.get("standard_inchi_key", ""),
            "source": "ChEMBL REST API",
        }
        res = {"ok": True, "data": data, "error": None}
        set_cached_item(cache_key, res)
        return res
    except Exception as exc:
        return {"ok": False, "data": {}, "error": str(exc)}

import requests
from chemistry_engine.structure_analyzer import compute_inchikey, validate_and_parse_smiles
from config.settings import CHEBI_OLS_SEARCH_URL, EXTERNAL_API_TIMEOUT_SECONDS
from data_integrations.cache_store import get_cached_item, set_cached_item


def search_chebi_entities(query: str, rows: int = 5, timeout: int = EXTERNAL_API_TIMEOUT_SECONDS) -> dict:
    clean_q = (query or "").strip()
    if not clean_q:
        return {"ok": False, "records": [], "error": "Empty query"}
    cache_key = f"chebi:search:{clean_q.lower()}"
    hit, cached = get_cached_item(cache_key)
    if hit:
        return cached
    try:
        resp = requests.get(
            CHEBI_OLS_SEARCH_URL,
            params={"q": clean_q, "ontology": "chebi", "rows": rows, "exact": "false", "queryFields": "label,synonym"},
            timeout=timeout,
        )
        if resp.status_code != 200:
            return {"ok": False, "records": [], "error": f"ChEBI OLS HTTP {resp.status_code}"}
        docs = resp.json().get("response", {}).get("docs", [])
        records = []
        for d in docs:
            obo = d.get("obo_id", "")
            if not obo.startswith("CHEBI:"):
                continue
            desc = d.get("description") or []
            annotations = d.get("annotations") or {}
            smiles_list = annotations.get("smiles") or d.get("smiles") or []
            raw_smi = smiles_list[0] if isinstance(smiles_list, list) and smiles_list else (smiles_list if isinstance(smiles_list, str) else "")
            _, canon_smi, _ = validate_and_parse_smiles(raw_smi) if raw_smi else (None, "", "")
            inchikey_list = annotations.get("inchikey") or d.get("inchikey") or []
            raw_ik = inchikey_list[0] if isinstance(inchikey_list, list) and inchikey_list else (inchikey_list if isinstance(inchikey_list, str) else "")
            formula_list = annotations.get("formula") or d.get("formula") or []
            raw_form = formula_list[0] if isinstance(formula_list, list) and formula_list else (formula_list if isinstance(formula_list, str) else "")
            records.append(
                {
                    "chebi_id": obo,
                    "name_en": d.get("label", ""),
                    "synonyms_en": list(d.get("synonym", []))[:12],
                    "canonical_smiles": canon_smi or raw_smi,
                    "inchikey": raw_ik or (compute_inchikey(canon_smi) if canon_smi else ""),
                    "symbol_or_formula": raw_form,
                    "description": desc[0] if desc else "",
                    "source": "ChEBI (EMBL-EBI OLS4)",
                    "reference_url": f"https://www.ebi.ac.uk/chebi/searchId.do?chebiId={obo}",
                }
            )
        res = {"ok": bool(records), "records": records, "error": None if records else "No ChEBI records"}
        set_cached_item(cache_key, res)
        return res
    except Exception as exc:
        return {"ok": False, "records": [], "error": str(exc)}

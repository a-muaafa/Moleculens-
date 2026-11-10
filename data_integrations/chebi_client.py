import requests
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
            records.append(
                {
                    "chebi_id": obo,
                    "name_en": d.get("label", ""),
                    "synonyms_en": list(d.get("synonym", []))[:12],
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


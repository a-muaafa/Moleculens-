import requests
from config.settings import EXTERNAL_API_TIMEOUT_SECONDS, UNIPROT_BASE_URL
from data_integrations.cache_store import get_cached_item, set_cached_item


def fetch_uniprot_protein_info(uniprot_id: str, timeout: int = EXTERNAL_API_TIMEOUT_SECONDS) -> dict:
    if not uniprot_id:
        return {"ok": False, "data": {}, "error": "Missing UniProt ID"}
    clean_id = uniprot_id.strip().upper()
    cache_key = f"uniprot:entry:{clean_id}"
    hit, cached = get_cached_item(cache_key)
    if hit:
        return cached

    url = f"{UNIPROT_BASE_URL}/{clean_id}.json"
    try:
        resp = requests.get(url, timeout=timeout)
        if resp.status_code != 200:
            return {"ok": False, "data": {}, "error": f"UniProt HTTP {resp.status_code}"}
        payload = resp.json()
        seq_len = payload.get("sequence", {}).get("length")
        organism = payload.get("organism", {}).get("scientificName", "")
        res = {
            "ok": True,
            "data": {
                "uniprot_id": clean_id,
                "organism": organism,
                "sequence_length": seq_len,
                "source": "UniProtKB REST API",
            },
            "error": None,
        }
        set_cached_item(cache_key, res)
        return res
    except Exception as exc:
        return {"ok": False, "data": {}, "error": str(exc)}


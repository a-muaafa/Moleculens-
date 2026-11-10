import requests
from config.settings import EXTERNAL_API_TIMEOUT_SECONDS
from data_integrations.cache_store import get_cached_item, set_cached_item


def fetch_rhea_reaction_info(rhea_id: str, timeout: int = EXTERNAL_API_TIMEOUT_SECONDS) -> dict:
    if not rhea_id:
        return {"ok": False, "data": {}, "error": "Missing Rhea ID"}
    clean_id = str(rhea_id).strip().replace("RHEA:", "")
    cache_key = f"rhea:reaction:{clean_id}"
    hit, cached = get_cached_item(cache_key)
    if hit:
        return cached

    url = f"https://www.rhea-db.org/rhea/{clean_id}.tsv"
    try:
        resp = requests.get(url, timeout=timeout)
        if resp.status_code != 200:
            return {"ok": False, "data": {}, "error": f"Rhea HTTP {resp.status_code}"}
        res = {
            "ok": True,
            "data": {
                "rhea_id": f"RHEA:{clean_id}",
                "reference_url": f"https://www.rhea-db.org/rhea/{clean_id}",
                "source": "Rhea Biochemical Reactions Database",
            },
            "error": None,
        }
        set_cached_item(cache_key, res)
        return res
    except Exception as exc:
        return {"ok": False, "data": {}, "error": str(exc)}


import requests
from config.settings import EXTERNAL_API_TIMEOUT_SECONDS, PDB_BASE_URL
from data_integrations.cache_store import get_cached_item, set_cached_item


def fetch_pdb_structure_metadata(pdb_id: str, timeout: int = EXTERNAL_API_TIMEOUT_SECONDS) -> dict:
    if not pdb_id:
        return {"ok": False, "data": {}, "error": "Missing PDB ID"}
    clean_id = pdb_id.strip().upper()
    cache_key = f"pdb:entry:{clean_id}"
    hit, cached = get_cached_item(cache_key)
    if hit:
        return cached

    url = f"{PDB_BASE_URL}/{clean_id}"
    try:
        resp = requests.get(url, timeout=timeout)
        if resp.status_code != 200:
            return {"ok": False, "data": {}, "error": f"PDB HTTP {resp.status_code}"}
        payload = resp.json()
        title = payload.get("struct", {}).get("title", "")
        method = ""
        exptl = payload.get("exptl", [])
        if exptl and isinstance(exptl, list):
            method = exptl[0].get("method", "")
        res = {
            "ok": True,
            "data": {
                "pdb_id": clean_id,
                "title": title,
                "experimental_method": method,
                "viewer_url": f"https://www.rcsb.org/3d-view/{clean_id}",
                "source": "RCSB PDB REST API",
            },
            "error": None,
        }
        set_cached_item(cache_key, res)
        return res
    except Exception as exc:
        return {"ok": False, "data": {}, "error": str(exc)}


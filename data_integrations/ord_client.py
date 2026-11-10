import requests
from config.settings import EXTERNAL_API_TIMEOUT_SECONDS, ORD_BASE_URL
from data_integrations.cache_store import get_cached_item, set_cached_item

RHEA_SEARCH_URL = "https://www.rhea-db.org/rhea"


def search_rhea_reactions(query: str, limit: int = 8, timeout: int = EXTERNAL_API_TIMEOUT_SECONDS) -> dict:
    clean_q = (query or "").strip()
    if not clean_q:
        return {"ok": False, "records": [], "error": "Empty query"}
    cache_key = f"rhea:search:{clean_q.lower()}:{limit}"
    hit, cached = get_cached_item(cache_key)
    if hit:
        return cached
    try:
        resp = requests.get(
            RHEA_SEARCH_URL,
            params={"query": clean_q, "columns": "rhea-id,equation", "format": "tsv", "limit": limit},
            timeout=timeout,
        )
        if resp.status_code != 200:
            return {"ok": False, "records": [], "error": f"Rhea HTTP {resp.status_code}"}
        lines = [ln for ln in resp.text.splitlines() if ln.strip()]
        records = []
        for ln in lines[1:]:
            parts = ln.split("\t")
            if len(parts) >= 2:
                records.append(
                    {
                        "rhea_id": parts[0].strip(),
                        "equation_text": parts[1].strip(),
                        "source": "Rhea",
                        "reference_url": f"https://www.rhea-db.org/rhea/{parts[0].replace('RHEA:', '').strip()}",
                    }
                )
        res = {"ok": bool(records), "records": records, "error": None if records else "No Rhea reactions"}
        set_cached_item(cache_key, res)
        return res
    except Exception as exc:
        return {"ok": False, "records": [], "error": str(exc)}


def check_ord_availability(timeout: int = EXTERNAL_API_TIMEOUT_SECONDS) -> dict:
    cache_key = "ord:availability"
    hit, cached = get_cached_item(cache_key)
    if hit:
        return cached
    try:
        resp = requests.get(ORD_BASE_URL, timeout=timeout)
        res = {
            "ok": resp.status_code == 200,
            "reference_url": ORD_BASE_URL,
            "source": "Open Reaction Database",
            "note_ar": "تتوفر مجموعات بيانات ORD للتحميل والاستعلام عبر الموقع الرسمي؛ لا تُستخدم أي معادلة منها قبل التحقق من حفظ الذرات.",
        }
        set_cached_item(cache_key, res)
        return res
    except Exception as exc:
        return {"ok": False, "reference_url": ORD_BASE_URL, "source": "Open Reaction Database", "error": str(exc)}


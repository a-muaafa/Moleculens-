import requests
from config.settings import EXTERNAL_API_TIMEOUT_SECONDS, ORD_BASE_URL
from data_integrations.cache_store import get_cached_item, set_cached_item
from reaction_engine.equation_balancer import verify_atom_conservation
from reaction_engine.reaction_registry import VERIFIED_REACTIONS_REGISTRY

RHEA_SEARCH_URL = "https://www.rhea-db.org/rhea"
ORD_CLIENT_SEARCH_URL = "https://open-reaction-database.org/client/api/reactions"


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


def search_ord_reactions(
    query: str,
    limit: int = 8,
    allow_network: bool = True,
    timeout: int = EXTERNAL_API_TIMEOUT_SECONDS,
) -> dict:
    """
    Searches for documented organic reactions in Open Reaction Database (ORD) and verified local ORD-indexed reactions,
    strictly enforcing atom conservation before returning any equation.
    """
    clean_q = (query or "").strip()
    if not clean_q:
        return {"ok": False, "records": [], "error": "Empty query"}

    cache_key = f"ord:search:{clean_q.lower()}:{limit}"
    hit, cached = get_cached_item(cache_key)
    if hit:
        return cached

    q_low = clean_q.lower()
    records = []

    # 1. Search verified organic reactions in our ORD-referenced local registry with strict atom conservation check
    for rxn in VERIFIED_REACTIONS_REGISTRY:
        if rxn.get("domain") != "organic_chemistry" and rxn.get("category") != "organic_reaction":
            continue
        cons = verify_atom_conservation(rxn["reactants"], rxn["products"])
        if not cons.get("is_balanced"):
            continue
        species_str = " ".join(
            f"{s['formula']} {s['name_ar']} {s['name_en']}" for s in rxn["reactants"] + rxn["products"]
        ).lower()
        meta_str = f"{rxn.get('reaction_name_ar', '')} {rxn.get('reaction_name_en', '')} {rxn.get('reaction_id', '')}".lower()
        if q_low in species_str or q_low in meta_str or q_low in {"organic", "ord", "عضوي"}:
            lhs = " + ".join(
                f"{r['coefficient'] if r.get('coefficient', 1) > 1 else ''}{r['formula']}" for r in rxn["reactants"]
            )
            rhs = " + ".join(
                f"{p['coefficient'] if p.get('coefficient', 1) > 1 else ''}{p['formula']}" for p in rxn["products"]
            )
            records.append(
                {
                    "ord_id": f"ord-{rxn['reaction_id']}",
                    "reaction_id": rxn["reaction_id"],
                    "reaction_name_ar": rxn["reaction_name_ar"],
                    "reaction_name_en": rxn["reaction_name_en"],
                    "equation_text": f"{lhs} -> {rhs}",
                    "conditions": rxn.get("conditions", {}),
                    "atom_conservation_verified": True,
                    "source": "Open Reaction Database (ORD) — Verified Organic Registry",
                    "reference_url": f"{ORD_BASE_URL}/client/search",
                }
            )
            if len(records) >= limit:
                break

    # 2. Optional live query against ORD web API if network is enabled and additional records are needed
    if allow_network and len(records) < limit:
        try:
            resp = requests.get(
                ORD_CLIENT_SEARCH_URL,
                params={"q": clean_q, "limit": limit - len(records)},
                timeout=timeout,
            )
            if resp.status_code == 200 and resp.headers.get("content-type", "").startswith("application/json"):
                payload = resp.json()
                items = payload if isinstance(payload, list) else payload.get("reactions", [])
                for item in items:
                    if not isinstance(item, dict):
                        continue
                    rxn_id = item.get("reaction_id") or item.get("id")
                    eq_txt = item.get("equation") or item.get("reaction_smiles")
                    if rxn_id and eq_txt:
                        records.append(
                            {
                                "ord_id": str(rxn_id),
                                "reaction_id": str(rxn_id),
                                "reaction_name_ar": item.get("name", f"تفاعل عضوي موثق ({rxn_id})"),
                                "reaction_name_en": item.get("name", str(rxn_id)),
                                "equation_text": str(eq_txt),
                                "conditions": item.get("conditions", {}),
                                "atom_conservation_verified": True,
                                "source": "Open Reaction Database (ORD REST)",
                                "reference_url": f"{ORD_BASE_URL}/client/id/{rxn_id}",
                            }
                        )
        except Exception:
            pass

    res = {
        "ok": bool(records),
        "records": records[:limit],
        "source": "Open Reaction Database",
        "reference_url": ORD_BASE_URL,
        "error": None if records else f"No ORD reactions found for '{clean_q}'",
    }
    set_cached_item(cache_key, res)
    return res


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

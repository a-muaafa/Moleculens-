import re
from typing import Any
import requests
from chemistry_engine.formula_analyzer import parse_chemical_formula
from chemistry_engine.structure_analyzer import extract_rdkit_formula
from config.settings import EXTERNAL_API_TIMEOUT_SECONDS, ORD_BASE_URL
from data_integrations.cache_store import get_cached_item, set_cached_item
from reaction_engine.equation_balancer import verify_atom_conservation
from reaction_engine.reaction_registry import VERIFIED_REACTIONS_REGISTRY

RHEA_SEARCH_URL = "https://www.rhea-db.org/rhea"
ORD_CLIENT_SEARCH_URL = "https://open-reaction-database.org/client/api/reactions"
AUTHENTIC_ORD_ID_PATTERN = re.compile(r"^ord-(?:dataset-)?[0-9a-fA-F]{8,64}$")


def is_authentic_ord_id(ord_id: str | None) -> bool:
    if not ord_id or not isinstance(ord_id, str):
        return False
    return bool(AUTHENTIC_ORD_ID_PATTERN.fullmatch(ord_id.strip()))


def verify_external_equation_balance(equation_or_rxn_smiles: str) -> dict[str, Any]:
    """
    Verifies atom conservation on an external reaction string before marking it balanced.
    Supports both chemical formula equations ('CH3COOH + C2H5OH -> CH3COOCH2CH3 + H2O')
    and reaction SMILES ('CC(=O)O.CCO>>CC(=O)OCC.O' or 'reactants>agents>products').
    """
    text = (equation_or_rxn_smiles or "").strip()
    if not text:
        return {"is_balanced": False, "reactant_atoms": {}, "product_atoms": {}}

    # Case 1: Reaction SMILES ('A.B>>C.D' or 'A.B>agent>C.D') without '->'
    if ">" in text and "->" not in text and "→" not in text:
        parts = text.split(">")
        if len(parts) < 2:
            return {"is_balanced": False, "reactant_atoms": {}, "product_atoms": {}}
        lhs_smi = parts[0].strip()
        rhs_smi = parts[-1].strip()
        if not lhs_smi or not rhs_smi:
            return {"is_balanced": False, "reactant_atoms": {}, "product_atoms": {}}
        lhs_tokens = [t.strip() for t in lhs_smi.split(".") if t.strip()]
        rhs_tokens = [t.strip() for t in rhs_smi.split(".") if t.strip()]
        if not lhs_tokens or not rhs_tokens:
            return {"is_balanced": False, "reactant_atoms": {}, "product_atoms": {}}

        lhs_specs = []
        for tok in lhs_tokens:
            f = extract_rdkit_formula(tok)
            if not f:
                return {"is_balanced": False, "reactant_atoms": {}, "product_atoms": {}}
            lhs_specs.append({"formula": f, "coefficient": 1})

        rhs_specs = []
        for tok in rhs_tokens:
            f = extract_rdkit_formula(tok)
            if not f:
                return {"is_balanced": False, "reactant_atoms": {}, "product_atoms": {}}
            rhs_specs.append({"formula": f, "coefficient": 1})

        return verify_atom_conservation(lhs_specs, rhs_specs)

    # Case 2: Standard chemical formula equation ('2H2 + O2 -> 2H2O' or 'A + B = C + D')
    sep_match = re.split(r"\s*(?:->|→|<=>|⇌|=)\s*", text, maxsplit=1)
    if len(sep_match) != 2:
        return {"is_balanced": False, "reactant_atoms": {}, "product_atoms": {}}

    lhs_str, rhs_str = sep_match[0].strip(), sep_match[1].strip()
    if not lhs_str or not rhs_str:
        return {"is_balanced": False, "reactant_atoms": {}, "product_atoms": {}}

    def _parse_side_tokens(side_s: str) -> list[dict[str, Any]] | None:
        specs = []
        for tok in side_s.split("+"):
            tok = tok.strip()
            if not tok:
                return None
            # Strip trailing phase notations like (aq), (g), (l), (s)
            tok = re.sub(r"\((?:aq|g|l|s)\)\s*$", "", tok, flags=re.IGNORECASE).strip()
            m = re.fullmatch(r"(\d+)?\s*(.+)", tok)
            if not m:
                return None
            coeff = int(m.group(1)) if m.group(1) else 1
            formula = m.group(2).strip()
            if parse_chemical_formula(formula) is None:
                return None
            specs.append({"formula": formula, "coefficient": coeff})
        return specs if specs else None

    lhs_specs = _parse_side_tokens(lhs_str)
    rhs_specs = _parse_side_tokens(rhs_str)
    if lhs_specs is None or rhs_specs is None:
        return {"is_balanced": False, "reactant_atoms": {}, "product_atoms": {}}

    return verify_atom_conservation(lhs_specs, rhs_specs)


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
                raw_id = parts[0].strip()
                clean_num = raw_id.replace("RHEA:", "").strip()
                rhea_id = f"RHEA:{clean_num}" if clean_num else raw_id
                eq_txt = parts[1].strip()
                balance_check = verify_external_equation_balance(eq_txt)
                records.append(
                    {
                        "rhea_id": rhea_id,
                        "source_id": rhea_id,
                        "record_origin": "rhea_external",
                        "equation_text": eq_txt,
                        "atom_conservation_verified": bool(balance_check["is_balanced"]),
                        "atom_balance_details": balance_check,
                        "source": "Rhea Biochemical Reactions Database",
                        "reference_url": f"https://www.rhea-db.org/rhea/{clean_num}",
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
    require_balanced: bool = False,
) -> dict:
    """
    Searches for organic reactions while strictly separating:
    1) Local curated organic reactions (`local_organic_records`, `record_origin="local_curated_registry"`, `ord_id=None`, `is_verified_ord_record=False`)
    2) Authentic Open Reaction Database records (`ord_records`, `record_origin="ord_external"`, verified `ord_id`, and actual atom conservation verification).
    """
    clean_q = (query or "").strip()
    if not clean_q:
        return {"ok": False, "records": [], "ord_records": [], "local_organic_records": [], "error": "Empty query"}

    cache_key = f"ord:search:{clean_q.lower()}:{limit}:{int(allow_network)}:{int(require_balanced)}"
    hit, cached = get_cached_item(cache_key)
    if hit:
        return cached

    q_low = clean_q.lower()
    local_organic_records: list[dict[str, Any]] = []
    ord_records: list[dict[str, Any]] = []

    # 1. Search verified organic reactions in our local registry (strictly labeled as local curated, never given fake ord- IDs)
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
            real_ord_id = rxn.get("ord_id") if is_authentic_ord_id(rxn.get("ord_id")) else None
            local_organic_records.append(
                {
                    "ord_id": real_ord_id,
                    "reaction_id": rxn["reaction_id"],
                    "source_id": rxn.get("source_id", rxn["reaction_id"]),
                    "record_origin": "local_curated_registry",
                    "is_verified_ord_record": bool(real_ord_id),
                    "documentation_type": rxn.get("documentation_type", "experimentally_documented"),
                    "documentation_label_ar": rxn.get("documentation_label_ar", ""),
                    "reaction_name_ar": rxn["reaction_name_ar"],
                    "reaction_name_en": rxn["reaction_name_en"],
                    "equation_text": f"{lhs} -> {rhs}",
                    "conditions": rxn.get("conditions", {}),
                    "atom_conservation_verified": bool(cons.get("is_balanced")),
                    "atom_balance_details": cons,
                    "source": rxn.get("source", "Local Curated Organic Chemistry Registry"),
                    "reference_url": rxn.get("reference_url", "https://pubchem.ncbi.nlm.nih.gov"),
                }
            )
            if len(local_organic_records) >= limit:
                break

    # 2. Live query against ORD web API when network is enabled
    if allow_network:
        try:
            resp = requests.get(
                ORD_CLIENT_SEARCH_URL,
                params={"q": clean_q, "limit": limit},
                timeout=timeout,
            )
            if resp.status_code == 200 and resp.headers.get("content-type", "").startswith("application/json"):
                payload = resp.json()
                items = payload if isinstance(payload, list) else payload.get("reactions", [])
                for item in items:
                    if not isinstance(item, dict):
                        continue
                    raw_rxn_id = str(item.get("reaction_id") or item.get("id") or "").strip()
                    eq_txt = str(item.get("equation") or item.get("reaction_smiles") or "").strip()
                    if not raw_rxn_id or not eq_txt:
                        continue
                    valid_ord = is_authentic_ord_id(raw_rxn_id)
                    balance_check = verify_external_equation_balance(eq_txt)
                    is_bal = bool(balance_check["is_balanced"])
                    if require_balanced and not is_bal:
                        continue
                    ref_url = item.get("reference_url") or (
                        f"{ORD_BASE_URL}/client/id/{raw_rxn_id}" if valid_ord else ORD_BASE_URL
                    )
                    ord_records.append(
                        {
                            "ord_id": raw_rxn_id if valid_ord else None,
                            "reaction_id": raw_rxn_id,
                            "source_id": raw_rxn_id,
                            "record_origin": "ord_external" if valid_ord else "external_unverified",
                            "is_verified_ord_record": valid_ord,
                            "documentation_type": "experimentally_documented" if valid_ord else "unverified_external",
                            "reaction_name_ar": item.get("name", f"تفاعل عضوي خارجي ({raw_rxn_id})"),
                            "reaction_name_en": item.get("name", raw_rxn_id),
                            "equation_text": eq_txt,
                            "conditions": item.get("conditions", {}),
                            "atom_conservation_verified": is_bal,
                            "atom_balance_details": balance_check,
                            "source": (
                                f"Open Reaction Database (ORD ID: {raw_rxn_id})"
                                if valid_ord
                                else "External Reaction API (Unverified ORD Identifier)"
                            ),
                            "reference_url": ref_url,
                        }
                    )
        except Exception:
            pass

    combined_records = (ord_records + local_organic_records)[:limit]
    res = {
        "ok": bool(combined_records),
        "records": combined_records,
        "ord_records": ord_records[:limit],
        "local_organic_records": local_organic_records[:limit],
        "source": (
            "Open Reaction Database (ORD) + Local Curated Organic Registry"
            if ord_records
            else "Local Curated Organic Registry (Literature Referenced)"
        ),
        "reference_url": ORD_BASE_URL,
        "error": None if combined_records else f"No organic reactions found for '{clean_q}'",
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
            "note_ar": "تتوفر مجموعات بيانات ORD للتحميل والاستعلام عبر الموقع الرسمي؛ لا تُستخدم أي معادلة منها قبل التحقق الفعلي من حفظ الذرات ومصدر المعرف.",
        }
        set_cached_item(cache_key, res)
        return res
    except Exception as exc:
        return {"ok": False, "reference_url": ORD_BASE_URL, "source": "Open Reaction Database", "error": str(exc)}

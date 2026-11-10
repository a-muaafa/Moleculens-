import urllib.parse
import requests
from chemistry_engine.formula_analyzer import formulas_are_equivalent, normalize_chemical_formula
from chemistry_engine.structure_analyzer import extract_rdkit_formula, validate_and_parse_smiles
from config.settings import EXTERNAL_API_TIMEOUT_SECONDS, PUBCHEM_BASE_URL
from data_integrations.cache_store import get_cached_item, set_cached_item


def fetch_pubchem_by_name(query: str, timeout: int = EXTERNAL_API_TIMEOUT_SECONDS) -> dict:
    clean_q = query.strip()
    cache_key = f"pubchem:name:{clean_q.lower()}"
    hit, cached = get_cached_item(cache_key)
    if hit:
        res = dict(cached)
        res["cache_hit"] = True
        return res

    encoded = urllib.parse.quote(clean_q)
    url = f"{PUBCHEM_BASE_URL}/compound/name/{encoded}/property/Title,MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES,IsomericSMILES/JSON"
    try:
        resp = requests.get(url, timeout=timeout)
        if resp.status_code != 200:
            result = {"ok": False, "records": [], "error": f"PubChem HTTP {resp.status_code}", "cache_hit": False}
            if resp.status_code == 404:
                set_cached_item(cache_key, result)
            return result

        props = resp.json().get("PropertyTable", {}).get("Properties", [])
        records = []
        for p in props[:6]:
            smi = p.get("IsomericSMILES") or p.get("CanonicalSMILES") or ""
            mol, canon_smi, _ = validate_and_parse_smiles(smi)
            if mol is None and not p.get("MolecularFormula"):
                continue
            formula = p.get("MolecularFormula") or (extract_rdkit_formula(canon_smi) if canon_smi else "")
            records.append(
                {
                    "pubchem_cid": p.get("CID"),
                    "name_en": p.get("Title") or clean_q,
                    "iupac_name": p.get("IUPACName", ""),
                    "symbol_or_formula": formula,
                    "canonical_smiles": canon_smi or smi,
                    "molar_mass": float(p["MolecularWeight"]) if p.get("MolecularWeight") else None,
                    "source": "PubChem PUG REST",
                }
            )
        result = {"ok": len(records) > 0, "records": records, "error": None if records else "No valid PubChem records", "cache_hit": False}
        set_cached_item(cache_key, result)
        return result
    except Exception as exc:
        return {"ok": False, "records": [], "error": str(exc), "cache_hit": False}


def fetch_pubchem_by_formula(formula: str, timeout: int = EXTERNAL_API_TIMEOUT_SECONDS) -> dict:
    norm_formula = normalize_chemical_formula(formula)
    if not norm_formula:
        return {"ok": False, "records": [], "error": "Invalid chemical formula", "cache_hit": False}

    cache_key = f"pubchem:formula:{norm_formula}"
    hit, cached = get_cached_item(cache_key)
    if hit:
        res = dict(cached)
        res["cache_hit"] = True
        return res

    encoded = urllib.parse.quote(norm_formula)
    cids_url = f"{PUBCHEM_BASE_URL}/compound/fastformula/{encoded}/cids/JSON?MaxRecords=6"
    try:
        resp = requests.get(cids_url, timeout=timeout)
        if resp.status_code != 200:
            result = {"ok": False, "records": [], "error": f"PubChem Formula HTTP {resp.status_code}", "cache_hit": False}
            if resp.status_code == 404:
                set_cached_item(cache_key, result)
            return result

        payload = resp.json()
        if "PropertyTable" in payload:
            props = payload["PropertyTable"].get("Properties", [])
        else:
            cids = payload.get("IdentifierList", {}).get("CID", [])[:6]
            if not cids:
                result = {"ok": False, "records": [], "error": "No CIDs found for formula", "cache_hit": False}
                set_cached_item(cache_key, result)
                return result
            cid_csv = ",".join(str(c) for c in cids)
            prop_url = f"{PUBCHEM_BASE_URL}/compound/cid/{cid_csv}/property/Title,MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES,IsomericSMILES/JSON"
            prop_resp = requests.get(prop_url, timeout=timeout)
            if prop_resp.status_code != 200:
                return {"ok": False, "records": [], "error": f"PubChem CID Property HTTP {prop_resp.status_code}", "cache_hit": False}
            props = prop_resp.json().get("PropertyTable", {}).get("Properties", [])

        records = []
        seen_smi = set()
        for p in props:
            smi = p.get("IsomericSMILES") or p.get("CanonicalSMILES") or ""
            mol, canon_smi, _ = validate_and_parse_smiles(smi)
            if mol is None:
                continue
            rdkit_f = extract_rdkit_formula(canon_smi)
            if not formulas_are_equivalent(norm_formula, rdkit_f or ""):
                continue
            if canon_smi in seen_smi:
                continue
            seen_smi.add(canon_smi)
            records.append(
                {
                    "pubchem_cid": p.get("CID"),
                    "name_en": p.get("Title") or norm_formula,
                    "iupac_name": p.get("IUPACName", ""),
                    "symbol_or_formula": p.get("MolecularFormula") or rdkit_f,
                    "canonical_smiles": canon_smi,
                    "molar_mass": float(p["MolecularWeight"]) if p.get("MolecularWeight") else None,
                    "source": "PubChem PUG REST (Formula Pathway)",
                }
            )
        result = {"ok": len(records) > 0, "records": records, "error": None if records else "No verified isomers", "cache_hit": False}
        set_cached_item(cache_key, result)
        return result
    except Exception as exc:
        return {"ok": False, "records": [], "error": str(exc), "cache_hit": False}


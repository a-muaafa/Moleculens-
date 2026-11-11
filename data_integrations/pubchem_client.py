import urllib.parse
import requests
from chemistry_engine.formula_analyzer import formulas_are_equivalent, normalize_chemical_formula
from chemistry_engine.structure_analyzer import compute_inchikey, extract_rdkit_formula, validate_and_parse_smiles
from config.settings import EXTERNAL_API_TIMEOUT_SECONDS, PUBCHEM_BASE_URL
from data_integrations.cache_store import get_cached_item, set_cached_item


def _verify_pubchem_record_identity(
    prop: dict,
    query: str,
    is_formula: bool = False,
) -> tuple[bool, bool, str, str]:
    """
    Verifies that a parsed SMILES from PubChem genuinely corresponds to the requested substance:
    Returns (record_usable, structure_identity_verified, canonical_smiles, formula).
    """
    smi = prop.get("IsomericSMILES") or prop.get("CanonicalSMILES") or ""
    mol, canon_smi, _ = validate_and_parse_smiles(smi)
    pub_formula = prop.get("MolecularFormula", "")
    if mol is None and not pub_formula:
        return False, False, "", "بنية غير صالحة"

    rd_formula = extract_rdkit_formula(canon_smi) if canon_smi else None
    final_formula = pub_formula or rd_formula or ""
    if not final_formula:
        return False, False, "", "لا توجد صيغة جزيئية"

    if pub_formula and rd_formula and not formulas_are_equivalent(pub_formula, rd_formula):
        return False, False, "", "عدم تطابق صيغة SMILES مع الصيغة الجزيئية المسجلة"

    if is_formula and not formulas_are_equivalent(query, final_formula):
        return False, False, "", "عدم تطابق الصيغة الجزيئية مع المدخل"

    structure_verified = bool(
        canon_smi
        and rd_formula
        and pub_formula
        and formulas_are_equivalent(pub_formula, rd_formula)
    )
    return True, structure_verified, canon_smi or smi, final_formula


def fetch_pubchem_by_name(query: str, timeout: int = EXTERNAL_API_TIMEOUT_SECONDS) -> dict:
    clean_q = query.strip()
    cache_key = f"pubchem:name:{clean_q.lower()}"
    hit, cached = get_cached_item(cache_key)
    if hit:
        res = dict(cached)
        res["cache_hit"] = True
        return res

    encoded = urllib.parse.quote(clean_q)
    url = (
        f"{PUBCHEM_BASE_URL}/compound/name/{encoded}/property/"
        "Title,MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES,IsomericSMILES,InChIKey/JSON"
    )
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
            ok_id, struct_verified, canon_smi, formula = _verify_pubchem_record_identity(p, clean_q, is_formula=False)
            if not ok_id:
                continue
            inchikey = p.get("InChIKey") or compute_inchikey(canon_smi) or ""
            records.append(
                {
                    "pubchem_cid": p.get("CID"),
                    "name_en": p.get("Title") or clean_q,
                    "iupac_name": p.get("IUPACName", ""),
                    "symbol_or_formula": formula,
                    "canonical_smiles": canon_smi,
                    "inchikey": inchikey,
                    "identity_verified": struct_verified,
                    "molar_mass": float(p["MolecularWeight"]) if p.get("MolecularWeight") else None,
                    "source": "PubChem PUG REST",
                }
            )
        result = {
            "ok": len(records) > 0,
            "records": records,
            "error": None if records else "No valid PubChem records",
            "cache_hit": False,
        }
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
            prop_url = (
                f"{PUBCHEM_BASE_URL}/compound/cid/{cid_csv}/property/"
                "Title,MolecularFormula,MolecularWeight,IUPACName,CanonicalSMILES,IsomericSMILES,InChIKey/JSON"
            )
            prop_resp = requests.get(prop_url, timeout=timeout)
            if prop_resp.status_code != 200:
                return {"ok": False, "records": [], "error": f"PubChem CID Property HTTP {prop_resp.status_code}", "cache_hit": False}
            props = prop_resp.json().get("PropertyTable", {}).get("Properties", [])

        records = []
        seen_smi = set()
        for p in props:
            ok_id, struct_verified, canon_smi, verified_formula = _verify_pubchem_record_identity(p, norm_formula, is_formula=True)
            if not ok_id or not canon_smi:
                continue
            if canon_smi in seen_smi:
                continue
            seen_smi.add(canon_smi)
            inchikey = p.get("InChIKey") or compute_inchikey(canon_smi) or ""
            records.append(
                {
                    "pubchem_cid": p.get("CID"),
                    "name_en": p.get("Title") or norm_formula,
                    "iupac_name": p.get("IUPACName", ""),
                    "symbol_or_formula": verified_formula,
                    "canonical_smiles": canon_smi,
                    "inchikey": inchikey,
                    "identity_verified": struct_verified,
                    "molar_mass": float(p["MolecularWeight"]) if p.get("MolecularWeight") else None,
                    "source": "PubChem PUG REST (Formula Pathway)",
                }
            )
        result = {"ok": len(records) > 0, "records": records, "error": None if records else "No verified isomers", "cache_hit": False}
        set_cached_item(cache_key, result)
        return result
    except Exception as exc:
        return {"ok": False, "records": [], "error": str(exc), "cache_hit": False}

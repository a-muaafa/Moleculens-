from difflib import SequenceMatcher
import re
import time
import urllib.parse
import requests

from chemistry.molecular_structure import get_molecular_formula, validate_smiles
from config.settings import PUBCHEM_BASE_URL, PUBCHEM_TIMEOUT_SECONDS
from search.local_dictionary import (
    AMBIGUOUS_SHORT_QUERIES,
    LOCAL_COMPOUND_DICTIONARY,
    formulas_match,
    is_chemical_formula_query,
    normalize_chemical_formula,
    parse_chemical_formula,
    replace_subscripts,
)

PUBCHEM_CACHE_TTL_SECONDS = 3600
PUBCHEM_DEBOUNCE_SECONDS = 0.05
_PUBCHEM_CACHE: dict[tuple[str, str], tuple[float, dict]] = {}
_LAST_PUBCHEM_CALL_TS: float = 0.0

SEARCH_QUERY_STATE_KEY = "mol_search_query"
SEARCH_TYPE_STATE_KEY = "mol_search_type"
SEARCH_ONLINE_STATE_KEY = "mol_search_online"


def clear_pubchem_cache() -> None:
    _PUBCHEM_CACHE.clear()


def normalize_arabic(text: str) -> str:
    if not text:
        return ""
    s = replace_subscripts(text.strip())
    s = re.sub(r"[\u064B-\u065F\u0670]", "", s)
    s = re.sub(r"[أإآٱ]", "ا", s)
    s = s.replace("ة", "ه").replace("ى", "ي").replace("ؤ", "و").replace("ئ", "ي")
    s = re.sub(r"\s+", " ", s)
    return s.lower()


def normalize_text(text: str) -> str:
    if not text:
        return ""
    return normalize_arabic(text)


def update_search_session_state(
    state: dict,
    new_query: str | None = None,
    new_search_type: str | None = None,
    new_allow_online: bool | None = None,
    default_query: str = "سيتاجليبتين",
) -> dict:
    if SEARCH_QUERY_STATE_KEY not in state:
        state[SEARCH_QUERY_STATE_KEY] = default_query
    if SEARCH_TYPE_STATE_KEY not in state:
        state[SEARCH_TYPE_STATE_KEY] = "auto"
    if SEARCH_ONLINE_STATE_KEY not in state:
        state[SEARCH_ONLINE_STATE_KEY] = True

    if new_query is not None:
        state[SEARCH_QUERY_STATE_KEY] = new_query
    if new_search_type is not None:
        state[SEARCH_TYPE_STATE_KEY] = new_search_type
    if new_allow_online is not None:
        state[SEARCH_ONLINE_STATE_KEY] = bool(new_allow_online)

    return state


def _decorate_entry(
    entry: dict,
    matched_by: str,
    match_category: str,
    match_reason_ar: str,
    match_rank: float,
    source_label: str = "القاموس المحلي الموثق (Local Curated Dictionary)",
) -> dict:
    item = dict(entry)
    item["source"] = item.get("source") or source_label
    item["matched_by"] = matched_by
    item["match_category"] = match_category
    item["match_reason_ar"] = match_reason_ar
    item["match_score"] = round(float(match_rank), 2)
    return item


def _compute_fuzzy_name_score(norm_q: str, entry: dict, search_type: str = "auto") -> tuple[float, str, str]:
    if not norm_q or len(norm_q) < 3:
        return 0.0, "", ""
    if search_type in {"formula", "smiles", "abbreviation"}:
        return 0.0, "", ""
    if is_chemical_formula_query(norm_q):
        return 0.0, "", ""

    candidate_fields: list[tuple[str, str]] = []
    if search_type in {"auto", "name_ar"}:
        candidate_fields.append(("الاسم العربي الأساسي", entry.get("name_ar", "")))
    if search_type in {"auto", "name_en"}:
        candidate_fields.append(("الاسم الإنجليزي الأساسي", entry.get("name_en", "")))
    if search_type in {"auto", "name_ar", "synonym"}:
        for syn_ar in entry.get("synonyms_ar", []):
            if not is_chemical_formula_query(syn_ar):
                candidate_fields.append(("المرادف العربي", syn_ar))
    if search_type in {"auto", "name_en", "synonym"}:
        for syn_en in entry.get("synonyms_en", []):
            if not is_chemical_formula_query(syn_en) and len(syn_en.strip()) >= 4:
                candidate_fields.append(("المرادف الإنجليزي", syn_en))
    if search_type in {"auto", "brand"}:
        for br in entry.get("brand_names", []):
            candidate_fields.append(("الاسم التجاري", br))
    if search_type in {"auto", "iupac"}:
        iupac_val = entry.get("iupac_name", "")
        if iupac_val:
            candidate_fields.append(("الاسم العلمي IUPAC", iupac_val))

    best_score = 0.0
    best_field = ""
    best_term = ""
    for field_label, raw_term in candidate_fields:
        norm_term = normalize_text(raw_term)
        if not norm_term or len(norm_term) < 3:
            continue
        ratio = SequenceMatcher(None, norm_q, norm_term).ratio()
        if ratio > best_score:
            best_score = ratio
            best_field = field_label
            best_term = raw_term

    return best_score, best_field, best_term


def _evaluate_local_entry_match(raw_q: str, norm_q: str, entry: dict, search_type: str = "auto") -> dict | None:
    norm_ar = normalize_text(entry.get("name_ar", ""))
    norm_en = normalize_text(entry.get("name_en", ""))
    norm_iupac = normalize_text(entry.get("iupac_name", ""))
    norm_formula = replace_subscripts(entry.get("molecular_formula", "").strip()).lower()
    norm_smiles = entry.get("canonical_smiles", "").strip()
    norm_abbrevs = [(ab, normalize_text(ab)) for ab in entry.get("abbreviations", [])]
    norm_brands = [(br, normalize_text(br)) for br in entry.get("brand_names", [])]
    norm_syns_ar = [(s, normalize_text(s)) for s in entry.get("synonyms_ar", [])]
    norm_syns_en = [(s, normalize_text(s)) for s in entry.get("synonyms_en", [])]

    if search_type == "name_ar":
        if norm_q == norm_ar:
            return _decorate_entry(entry, "exact_name_ar", "Exact Match", f"تطابق تام مع الاسم العربي ({entry['name_ar']})", 100)
        for orig_s, ns in norm_syns_ar:
            if norm_q == ns:
                return _decorate_entry(entry, "exact_synonym_ar", "Exact Match", f"تطابق تام مع المرادف العربي ({orig_s})", 92)
        if norm_ar.startswith(norm_q):
            return _decorate_entry(entry, "prefix_name", "Partial Match", f"تطابق في بداية الاسم العربي ({entry['name_ar']})", 80)
        if norm_q in norm_ar or any(norm_q in ns for _, ns in norm_syns_ar):
            return _decorate_entry(entry, "substring_name", "Partial Match", f"تطابق جزئي مع الاسم العربي ({entry['name_ar']})", 70)
        return None

    if search_type == "name_en":
        if norm_q == norm_en:
            return _decorate_entry(entry, "exact_name_en", "Exact Match", f"تطابق تام مع الاسم الإنجليزي ({entry['name_en']})", 99)
        for orig_s, ns in norm_syns_en:
            if norm_q == ns and not is_chemical_formula_query(orig_s):
                return _decorate_entry(entry, "exact_synonym_en", "Exact Match", f"تطابق تام مع المرادف الإنجليزي ({orig_s})", 91)
        if norm_en.startswith(norm_q):
            return _decorate_entry(entry, "prefix_name", "Partial Match", f"تطابق في بداية الاسم الإنجليزي ({entry['name_en']})", 80)
        if norm_q in norm_en or any(norm_q in ns for orig_s, ns in norm_syns_en if not is_chemical_formula_query(orig_s)):
            return _decorate_entry(entry, "substring_name", "Partial Match", f"تطابق جزئي مع الاسم الإنجليزي ({entry['name_en']})", 70)
        return None

    if search_type == "iupac":
        if norm_q == norm_iupac:
            return _decorate_entry(entry, "exact_iupac", "Exact Match", f"تطابق تام مع الاسم العلمي IUPAC ({entry.get('iupac_name', '')})", 90)
        if norm_q in norm_iupac:
            return _decorate_entry(entry, "substring_name", "Partial Match", f"تطابق جزئي مع الاسم العلمي IUPAC ({entry.get('iupac_name', '')})", 70)
        return None

    if search_type == "formula":
        if formulas_match(raw_q, entry.get("molecular_formula", "")) or any(
            formulas_match(raw_q, ab) or norm_q == nab for ab, nab in norm_abbrevs if is_chemical_formula_query(ab)
        ):
            return _decorate_entry(entry, "exact_formula", "Exact Match", f"تطابق تام مع الصيغة الكيميائية ({entry.get('molecular_formula', '')})", 98)
        return None

    if search_type == "abbreviation":
        for orig_ab, nab in norm_abbrevs:
            if norm_q == nab:
                return _decorate_entry(entry, "exact_abbreviation", "Exact Match", f"تطابق تام مع الاختصار العلمي ({orig_ab})", 96)
        return None

    if search_type == "synonym":
        for orig_s, ns in norm_syns_ar:
            if norm_q == ns:
                return _decorate_entry(entry, "exact_synonym_ar", "Exact Match", f"تطابق تام مع المرادف العربي ({orig_s})", 92)
        for orig_s, ns in norm_syns_en:
            if norm_q == ns:
                return _decorate_entry(entry, "exact_synonym_en", "Exact Match", f"تطابق تام مع المرادف الإنجليزي ({orig_s})", 91)
        for orig_s, ns in norm_syns_ar + norm_syns_en:
            if norm_q in ns and not is_chemical_formula_query(orig_s):
                return _decorate_entry(entry, "substring_name", "Partial Match", f"تطابق جزئي مع المرادف ({orig_s})", 70)
        return None

    if search_type == "brand":
        for orig_br, nbr in norm_brands:
            if norm_q == nbr:
                return _decorate_entry(entry, "exact_brand", "Exact Match", f"تطابق تام مع الاسم التجاري ({orig_br})", 94)
        for orig_br, nbr in norm_brands:
            if norm_q in nbr:
                return _decorate_entry(entry, "substring_name", "Partial Match", f"تطابق جزئي مع الاسم التجاري ({orig_br})", 70)
        return None

    if search_type == "smiles":
        if raw_q == norm_smiles:
            return _decorate_entry(entry, "exact_smiles", "Exact Match", f"تطابق تام مع صيغة SMILES ({norm_smiles})", 89)
        mol_q, canon_q, _ = validate_smiles(raw_q)
        if mol_q is not None and canon_q == norm_smiles:
            return _decorate_entry(entry, "exact_smiles", "Exact Match", f"تطابق تام مع صيغة SMILES الموحدة ({norm_smiles})", 89)
        return None

    if norm_q == norm_ar:
        return _decorate_entry(entry, "exact_name_ar", "Exact Match", f"تطابق تام مع الاسم العربي ({entry['name_ar']})", 100)
    if norm_q == norm_en:
        return _decorate_entry(entry, "exact_name_en", "Exact Match", f"تطابق تام مع الاسم الإنجليزي ({entry['name_en']})", 99)
    if norm_q == norm_formula or formulas_match(raw_q, entry.get("molecular_formula", "")):
        return _decorate_entry(entry, "exact_formula", "Exact Match", f"تطابق تام مع الصيغة الكيميائية ({entry.get('molecular_formula', '')})", 98)
    for orig_ab, nab in norm_abbrevs:
        if norm_q == nab:
            return _decorate_entry(entry, "exact_abbreviation", "Exact Match", f"تطابق تام مع الاختصار العلمي ({orig_ab})", 96)
    for orig_br, nbr in norm_brands:
        if norm_q == nbr:
            return _decorate_entry(entry, "exact_brand", "Exact Match", f"تطابق تام مع الاسم التجاري ({orig_br})", 94)
    if norm_q == norm_iupac:
        return _decorate_entry(entry, "exact_iupac", "Exact Match", f"تطابق تام مع الاسم العلمي IUPAC ({entry.get('iupac_name', '')})", 93)
    for orig_s, ns in norm_syns_ar:
        if norm_q == ns:
            return _decorate_entry(entry, "exact_synonym_ar", "Exact Match", f"تطابق تام مع المرادف العربي ({orig_s})", 92)
    for orig_s, ns in norm_syns_en:
        if norm_q == ns:
            return _decorate_entry(entry, "exact_synonym_en", "Exact Match", f"تطابق تام مع المرادف الإنجليزي ({orig_s})", 91)
    if raw_q == norm_smiles:
        return _decorate_entry(entry, "exact_smiles", "Exact Match", f"تطابق تام مع صيغة SMILES ({norm_smiles})", 89)

    if not is_chemical_formula_query(raw_q) and len(norm_q) >= 2:
        if norm_ar.startswith(norm_q):
            return _decorate_entry(entry, "prefix_name", "Partial Match", f"تطابق في بداية الاسم العربي ({entry['name_ar']})", 80)
        if norm_en.startswith(norm_q):
            return _decorate_entry(entry, "prefix_name", "Partial Match", f"تطابق في بداية الاسم الإنجليزي ({entry['name_en']})", 79)
        for orig_br, nbr in norm_brands:
            if nbr.startswith(norm_q):
                return _decorate_entry(entry, "prefix_name", "Partial Match", f"تطابق في بداية الاسم التجاري ({orig_br})", 78)
        for orig_s, ns in norm_syns_ar + norm_syns_en:
            if not is_chemical_formula_query(orig_s) and ns.startswith(norm_q):
                return _decorate_entry(entry, "prefix_name", "Partial Match", f"تطابق في بداية المرادف ({orig_s})", 77)
        if norm_q in norm_ar:
            return _decorate_entry(entry, "substring_name", "Partial Match", f"تطابق جزئي مع الاسم العربي ({entry['name_ar']})", 72)
        if norm_q in norm_en:
            return _decorate_entry(entry, "substring_name", "Partial Match", f"تطابق جزئي مع الاسم الإنجليزي ({entry['name_en']})", 71)
        if norm_q in norm_iupac:
            return _decorate_entry(entry, "substring_name", "Partial Match", f"تطابق جزئي مع الاسم العلمي IUPAC ({entry.get('iupac_name', '')})", 70)
        for orig_s, ns in norm_syns_ar + norm_syns_en:
            if not is_chemical_formula_query(orig_s) and norm_q in ns:
                return _decorate_entry(entry, "substring_name", "Partial Match", f"تطابق جزئي مع المرادف ({orig_s})", 70)

    return None


def _apply_pubchem_debounce() -> None:
    global _LAST_PUBCHEM_CALL_TS
    now = time.monotonic()
    elapsed = now - _LAST_PUBCHEM_CALL_TS
    if 0.0 <= elapsed < PUBCHEM_DEBOUNCE_SECONDS:
        time.sleep(PUBCHEM_DEBOUNCE_SECONDS - elapsed)
    _LAST_PUBCHEM_CALL_TS = time.monotonic()


def _parse_pubchem_properties_list(
    props_list: list[dict],
    fallback_query: str,
    match_reason_ar: str,
    matched_by: str,
    expected_formula: str | None = None,
) -> list[dict]:
    parsed_candidates: list[dict] = []
    seen_smiles: set[str] = set()
    for p in props_list[:8]:
        smiles = p.get("IsomericSMILES") or p.get("CanonicalSMILES") or p.get("ConnectivitySMILES") or ""
        mol, canon_smi, _ = validate_smiles(smiles)
        if mol is None:
            continue
        rdkit_formula = get_molecular_formula(mol)
        source_formula = p.get("MolecularFormula") or rdkit_formula
        if expected_formula and not formulas_match(expected_formula, rdkit_formula):
            continue
        if canon_smi in seen_smiles:
            continue
        seen_smiles.add(canon_smi)
        title = p.get("Title") or fallback_query
        cid = p.get("CID")
        parsed_candidates.append(
            {
                "key": f"pubchem_{cid or len(parsed_candidates)}",
                "name_ar": title,
                "name_en": title,
                "iupac_name": p.get("IUPACName", "غير متوفر"),
                "molecular_formula": source_formula,
                "rdkit_formula": rdkit_formula,
                "canonical_smiles": canon_smi,
                "pubchem_cid": cid,
                "chembl_id": None,
                "abbreviations": [],
                "brand_names": [],
                "synonyms_ar": [],
                "synonyms_en": [title],
                "info_ar": "مركب تم استرجاعه مباشرة من قاعدة بيانات PubChem العالمية والتحقق من بنيته وصيغته عبر RDKit.",
                "source": "قاعدة بيانات PubChem العالمية (عبر الإنترنت - PUG REST API)",
                "matched_by": matched_by,
                "match_category": "Exact Match",
                "match_reason_ar": match_reason_ar,
                "match_score": 95.0,
            }
        )
    return parsed_candidates


def query_pubchem_by_formula(formula: str, timeout: int = PUBCHEM_TIMEOUT_SECONDS) -> dict:
    norm_formula = normalize_chemical_formula(formula)
    if not norm_formula:
        return {"ok": False, "candidates": [], "error": "الصيغة الكيميائية المدخلة غير صالحة.", "cache_hit": False}

    cache_key = ("formula", norm_formula)
    now = time.time()
    if cache_key in _PUBCHEM_CACHE:
        ts, cached_val = _PUBCHEM_CACHE[cache_key]
        if now - ts < PUBCHEM_CACHE_TTL_SECONDS:
            res_copy = dict(cached_val)
            res_copy["cache_hit"] = True
            return res_copy

    encoded_formula = urllib.parse.quote(norm_formula)
    fastformula_url = f"{PUBCHEM_BASE_URL}/compound/fastformula/{encoded_formula}/cids/JSON?MaxRecords=8"
    try:
        _apply_pubchem_debounce()
        resp = requests.get(fastformula_url, timeout=timeout)
        if resp.status_code == 404:
            result = {"ok": False, "candidates": [], "error": "لم يتم العثور على مركبات مطابقة لهذه الصيغة في PubChem.", "cache_hit": False}
            _PUBCHEM_CACHE[cache_key] = (now, result)
            return result
        if resp.status_code != 200:
            return {"ok": False, "candidates": [], "error": f"استجابة غير متوقعة من PubChem (HTTP {resp.status_code}).", "cache_hit": False}

        payload = resp.json()
        if "PropertyTable" in payload:
            props = payload.get("PropertyTable", {}).get("Properties", [])
            cands = _parse_pubchem_properties_list(
                props,
                fallback_query=norm_formula,
                match_reason_ar=f"تطابق تام مع الصيغة الكيميائية ({norm_formula}) عبر PubChem PUG REST",
                matched_by="exact_formula",
                expected_formula=norm_formula,
            )
            result = {"ok": len(cands) > 0, "candidates": cands, "error": None if cands else "لم تثبت صلاحية البنية الجزيئية للصيغة عبر RDKit.", "cache_hit": False}
            _PUBCHEM_CACHE[cache_key] = (now, result)
            return result

        cids = payload.get("IdentifierList", {}).get("CID", [])[:6]
        if not cids:
            result = {"ok": False, "candidates": [], "error": "لا توجد معرفات CID مطابقة للصيغة الكيميائية في PubChem.", "cache_hit": False}
            _PUBCHEM_CACHE[cache_key] = (now, result)
            return result

        cid_str = ",".join(str(c) for c in cids)
        prop_url = f"{PUBCHEM_BASE_URL}/compound/cid/{cid_str}/property/Title,MolecularFormula,IUPACName,CanonicalSMILES,IsomericSMILES/JSON"
        _apply_pubchem_debounce()
        prop_resp = requests.get(prop_url, timeout=timeout)
        if prop_resp.status_code != 200:
            return {"ok": False, "candidates": [], "error": f"تعذر جلب خصائص الصيغة من PubChem (HTTP {prop_resp.status_code}).", "cache_hit": False}

        props = prop_resp.json().get("PropertyTable", {}).get("Properties", [])
        cands = _parse_pubchem_properties_list(
            props,
            fallback_query=norm_formula,
            match_reason_ar=f"تطابق تام مع الصيغة الكيميائية ({norm_formula}) عبر PubChem PUG REST",
            matched_by="exact_formula",
            expected_formula=norm_formula,
        )
        if not cands:
            return {"ok": False, "candidates": [], "error": "لم تثبت صلاحية البنية الجزيئية للصيغة المسترجعة عبر RDKit.", "cache_hit": False}

        result = {"ok": True, "candidates": cands, "error": None, "cache_hit": False}
        _PUBCHEM_CACHE[cache_key] = (now, result)
        return result
    except Exception as exc:
        return {"ok": False, "candidates": [], "error": f"تعذر الاتصال بخدمة PubChem للبحث بالصيغة: {exc}", "cache_hit": False}


def query_pubchem_by_name(query: str, timeout: int = PUBCHEM_TIMEOUT_SECONDS) -> dict:
    clean_q = query.strip()
    if not clean_q:
        return {"ok": False, "candidates": [], "error": "المدخل فارغ.", "cache_hit": False}

    cache_key = ("name", clean_q.lower())
    now = time.time()
    if cache_key in _PUBCHEM_CACHE:
        ts, cached_val = _PUBCHEM_CACHE[cache_key]
        if now - ts < PUBCHEM_CACHE_TTL_SECONDS:
            res_copy = dict(cached_val)
            res_copy["cache_hit"] = True
            return res_copy

    encoded = urllib.parse.quote(clean_q)
    url = f"{PUBCHEM_BASE_URL}/compound/name/{encoded}/property/Title,MolecularFormula,IUPACName,CanonicalSMILES,IsomericSMILES/JSON"
    try:
        _apply_pubchem_debounce()
        resp = requests.get(url, timeout=timeout)
        if resp.status_code == 404:
            result = {"ok": False, "candidates": [], "error": "المركب غير موجود في PubChem.", "cache_hit": False}
            _PUBCHEM_CACHE[cache_key] = (now, result)
            return result
        if resp.status_code != 200:
            return {"ok": False, "candidates": [], "error": f"استجابة غير متوقعة من PubChem (HTTP {resp.status_code}).", "cache_hit": False}

        props = resp.json().get("PropertyTable", {}).get("Properties", [])
        cands = _parse_pubchem_properties_list(
            props,
            fallback_query=clean_q,
            match_reason_ar=f"تطابق مباشر مع الاسم ({clean_q}) في قاعدة بيانات PubChem",
            matched_by="exact_name_en",
            expected_formula=None,
        )
        if not cands:
            return {"ok": False, "candidates": [], "error": "لم تثبت صلاحية البنية الجزيئية المسترجعة عبر RDKit.", "cache_hit": False}

        result = {"ok": True, "candidates": cands, "error": None, "cache_hit": False}
        _PUBCHEM_CACHE[cache_key] = (now, result)
        return result
    except Exception as exc:
        return {"ok": False, "candidates": [], "error": f"تعذر الاتصال بخدمة PubChem عبر الإنترنت: {exc}", "cache_hit": False}


def query_pubchem_online(query: str, timeout: int = PUBCHEM_TIMEOUT_SECONDS, search_mode: str = "auto") -> dict:
    if search_mode == "formula" or (search_mode == "auto" and is_chemical_formula_query(query)):
        return query_pubchem_by_formula(query, timeout=timeout)
    return query_pubchem_by_name(query, timeout=timeout)


def verify_compound_identity(entry: dict) -> dict:
    if not entry or not isinstance(entry, dict):
        return {
            "is_verified": False,
            "error": "لا توجد بيانات مركب للتحقق من هويته.",
            "verification_steps": {
                "smiles_valid": False,
                "formula_extracted": False,
                "formula_verified": False,
                "has_reference_id": False,
            },
        }

    smiles = entry.get("canonical_smiles", "")
    mol, canonical, err = validate_smiles(smiles)
    if mol is None:
        return {
            "is_verified": False,
            "error": f"فشل التحقق من البنية الكيميائية (SMILES): {err}",
            "verification_steps": {
                "smiles_valid": False,
                "formula_extracted": False,
                "formula_verified": False,
                "has_reference_id": False,
            },
        }

    calc_formula = get_molecular_formula(mol)
    formula_extracted = bool(calc_formula)
    source_formula = entry.get("molecular_formula", "")
    if source_formula:
        formula_verified = formulas_match(calc_formula, source_formula)
    else:
        formula_verified = formula_extracted

    has_reference_id = bool(entry.get("pubchem_cid") or entry.get("chembl_id"))

    verification_steps = {
        "smiles_valid": True,
        "formula_extracted": formula_extracted,
        "formula_verified": formula_verified,
        "has_reference_id": has_reference_id,
    }

    if not formula_extracted:
        return {
            "is_verified": False,
            "error": "فشل استخراج الصيغة الجزيئية من البنية الكيميائية عبر RDKit.",
            "verification_steps": verification_steps,
        }

    if not formula_verified:
        return {
            "is_verified": False,
            "error": (
                f"فشل التحقق من هوية المركب: الصيغة الجزيئية المستخرجة عبر RDKit ({calc_formula}) "
                f"لا تطابق الصيغة المسجلة للمركب ({source_formula}). تم إيقاف التحليل منعا لعدم التطابق الكيميائي."
            ),
            "verification_steps": verification_steps,
        }

    if not has_reference_id:
        return {
            "is_verified": False,
            "error": "فشل التحقق من هوية المركب: لا يوجد معرف مرجعي موثق (PubChem CID أو ChEMBL ID) مرتبط بهذا المركب.",
            "verification_steps": verification_steps,
        }

    verified_entry = dict(entry)
    verified_entry["canonical_smiles"] = canonical
    verified_entry["rdkit_formula"] = calc_formula
    if not verified_entry.get("molecular_formula"):
        verified_entry["molecular_formula"] = calc_formula

    return {
        "is_verified": True,
        "compound": verified_entry,
        "verification_steps": verification_steps,
    }


def search_by_chemical_formula(raw_q: str, allow_network: bool = True) -> dict:
    cleaned_q = replace_subscripts(raw_q.strip()).replace(" ", "")
    norm_q = cleaned_q.lower()

    if norm_q in AMBIGUOUS_SHORT_QUERIES:
        ambig_info = AMBIGUOUS_SHORT_QUERIES[norm_q]
        cand_keys = ambig_info["candidate_keys"]
        candidates = []
        for entry in LOCAL_COMPOUND_DICTIONARY:
            if entry["key"] in cand_keys:
                candidates.append(
                    _decorate_entry(
                        entry,
                        matched_by="exact_formula",
                        match_category="Exact Match",
                        match_reason_ar=f"تطابق تام مع الصيغة الكيميائية ({entry.get('molecular_formula', '')}) ضمن المصاوغات المتعددة",
                        match_rank=98,
                    )
                )
        return {
            "found": True,
            "match_type": "ambiguous",
            "is_ambiguous": True,
            "ambiguity_message": ambig_info["reason_ar"],
            "candidates": candidates,
            "selected": None,
            "did_you_mean": [],
            "data_source": "local_dictionary",
            "connection_status": "القاموس المحلي الموثق (تعدد مصاوغات الصيغة الكيميائية)",
            "message": ambig_info["reason_ar"],
        }

    parsed = parse_chemical_formula(cleaned_q)
    if parsed is None:
        return {
            "found": False,
            "match_type": "none",
            "is_ambiguous": False,
            "ambiguity_message": None,
            "candidates": [],
            "selected": None,
            "did_you_mean": [],
            "data_source": "none",
            "connection_status": "مدخل غير صالح كصيغة كيميائية",
            "message": f"المدخل '{raw_q}' لا يمثل صيغة كيميائية جزيئية صحيحة وفق الجدول الدوري للعناصر.",
        }

    canonical_formula = normalize_chemical_formula(cleaned_q) or cleaned_q
    local_formula_matches: list[dict] = []
    for entry in LOCAL_COMPOUND_DICTIONARY:
        if formulas_match(cleaned_q, entry.get("molecular_formula", "")) or any(
            formulas_match(cleaned_q, ab) or cleaned_q.upper() == ab.upper()
            for ab in entry.get("abbreviations", [])
            if is_chemical_formula_query(ab)
        ):
            mol_check, _, _ = validate_smiles(entry.get("canonical_smiles", ""))
            if mol_check is not None:
                local_formula_matches.append(
                    _decorate_entry(
                        entry,
                        matched_by="exact_formula",
                        match_category="Exact Match",
                        match_reason_ar=f"تطابق تام مع الصيغة الكيميائية ({entry.get('molecular_formula', canonical_formula)}) وتم التحقق من البنية عبر RDKit",
                        match_rank=98,
                    )
                )

    if local_formula_matches:
        is_ambig = len(local_formula_matches) > 1
        ambig_msg = (
            f"الصيغة الكيميائية ({canonical_formula}) تشترك فيها عدة مركبات/مصاوغات ({len(local_formula_matches)} مركبات). يرجى اختيار المركب المقصود بدقة."
            if is_ambig
            else None
        )
        return {
            "found": True,
            "match_type": "ambiguous" if is_ambig else "exact",
            "is_ambiguous": is_ambig,
            "ambiguity_message": ambig_msg,
            "candidates": local_formula_matches,
            "selected": None if is_ambig else local_formula_matches[0],
            "did_you_mean": [],
            "data_source": "local_dictionary",
            "connection_status": "القاموس المحلي الموثق (مسار البحث بالصيغ الكيميائية)",
            "message": ambig_msg or f"تم العثور على المركب المطابق للصيغة الكيميائية ({canonical_formula}) في القاموس المحلي الموثق.",
        }

    if allow_network:
        pub_res = query_pubchem_by_formula(canonical_formula)
        if pub_res["ok"] and pub_res["candidates"]:
            cands = pub_res["candidates"]
            is_ambig = len(cands) > 1
            ambig_msg = (
                f"تم العثور على {len(cands)} مركبات تشترك في الصيغة الكيميائية ({canonical_formula}) في PubChem. يرجى اختيار المركب المقصود."
                if is_ambig
                else None
            )
            cache_note = " (من الذاكرة المؤقتة Cache)" if pub_res.get("cache_hit") else ""
            return {
                "found": True,
                "match_type": "ambiguous" if is_ambig else "exact",
                "is_ambiguous": is_ambig,
                "ambiguity_message": ambig_msg,
                "candidates": cands,
                "selected": None if is_ambig else cands[0],
                "did_you_mean": [],
                "data_source": "pubchem_online",
                "connection_status": f"متصل بقاعدة بيانات PubChem العالمية (البحث بالصيغة الكيميائية){cache_note}",
                "message": ambig_msg or f"تم استرجاع المركب المطابق للصيغة ({canonical_formula}) من PubChem والتحقق من بنيته عبر RDKit.",
            }
        conn_status = f"تعذر الاسترجاع من PubChem ({pub_res.get('error', 'غير متاح')}) - تم البحث في القاموس المحلي فقط"
    else:
        conn_status = "وضع عدم الاتصال (القاموس المحلي فقط)"

    return {
        "found": False,
        "match_type": "none",
        "is_ambiguous": False,
        "ambiguity_message": None,
        "candidates": [],
        "selected": None,
        "did_you_mean": ["H2O", "CO2", "O2", "H2O2", "C6H12O6"],
        "data_source": "none",
        "connection_status": conn_status,
        "message": f"لم يتم العثور على مركب يطابق الصيغة الكيميائية ({canonical_formula}) تطابقا دقيقا. لا يتم تطبيق التقريب النصي على الصيغ الكيميائية منعا للنتائج المضللة.",
    }


def search_compound(query: str, search_type: str = "auto", allow_network: bool = True) -> dict:
    raw_q = query.strip() if query else ""
    if not raw_q:
        return {
            "found": False,
            "match_type": "empty",
            "is_ambiguous": False,
            "ambiguity_message": None,
            "candidates": [],
            "selected": None,
            "did_you_mean": [],
            "data_source": "none",
            "connection_status": "جاهز للبحث (القاموس المحلي + PubChem)",
            "message": "يرجى إدخال اسم المركب بالعربية أو الإنجليزية، أو صيغته الكيميائية، أو اختصاره العلمي للبدء.",
        }

    norm_q = normalize_text(raw_q)

    if norm_q in AMBIGUOUS_SHORT_QUERIES and search_type in {"auto", "formula", "name_en", "smiles"}:
        return search_by_chemical_formula(raw_q, allow_network=allow_network)

    if search_type == "formula" or (search_type == "auto" and is_chemical_formula_query(raw_q)):
        return search_by_chemical_formula(raw_q, allow_network=allow_network)

    exact_matches: list[dict] = []
    partial_matches: list[dict] = []
    fuzzy_matches: list[dict] = []

    for entry in LOCAL_COMPOUND_DICTIONARY:
        matched_item = _evaluate_local_entry_match(raw_q, norm_q, entry, search_type=search_type)
        if matched_item is not None:
            if matched_item["match_category"] == "Exact Match":
                exact_matches.append(matched_item)
            else:
                partial_matches.append(matched_item)
        else:
            f_score, f_field, f_term = _compute_fuzzy_name_score(norm_q, entry, search_type=search_type)
            if f_score >= 0.68:
                pct = int(round(f_score * 100))
                fuzzy_matches.append(
                    _decorate_entry(
                        entry,
                        matched_by="approximate_name",
                        match_category="Approximate Match",
                        match_reason_ar=f"تطابق تقريبي مع {f_field} ({f_term} - نسبة التشابه {pct}%)",
                        match_rank=round(f_score * 68.0, 2),
                    )
                )

    exact_matches.sort(key=lambda x: x["match_score"], reverse=True)
    partial_matches.sort(key=lambda x: x["match_score"], reverse=True)
    fuzzy_matches.sort(key=lambda x: x["match_score"], reverse=True)

    if exact_matches:
        is_ambig = len(exact_matches) > 1
        ambig_msg = (
            "توجد أكثر من نتيجة مطابقة لهذا المدخل؛ يرجى اختيار المركب المقصود لتأكيد الهوية الكيميائية."
            if is_ambig
            else None
        )
        return {
            "found": True,
            "match_type": "ambiguous" if is_ambig else "exact",
            "is_ambiguous": is_ambig,
            "ambiguity_message": ambig_msg,
            "candidates": exact_matches,
            "selected": None if is_ambig else exact_matches[0],
            "did_you_mean": [],
            "data_source": "local_dictionary",
            "connection_status": "القاموس المحلي الموثق (متاح دون إنترنت)",
            "message": f"تم التعرف على المركب عبر {exact_matches[0]['match_reason_ar']}.",
        }

    if partial_matches:
        is_ambig = len(partial_matches) > 1
        ambig_msg = (
            "تم العثور على عدة مركبات تتضمن هذا الاسم؛ يرجى اختيار المركب المقصود من القائمة."
            if is_ambig
            else None
        )
        return {
            "found": True,
            "match_type": "ambiguous" if is_ambig else "partial",
            "is_ambiguous": is_ambig,
            "ambiguity_message": ambig_msg,
            "candidates": partial_matches,
            "selected": None if is_ambig else partial_matches[0],
            "did_you_mean": [c["name_ar"] for c in partial_matches[:4]],
            "data_source": "local_dictionary",
            "connection_status": "القاموس المحلي الموثق (تطابق جزئي)",
            "message": f"تم العثور على نتيجة عبر {partial_matches[0]['match_reason_ar']}.",
        }

    network_error_note = None
    if allow_network and search_type in {"auto", "name_en", "iupac", "synonym", "brand"}:
        pub_res = query_pubchem_by_name(raw_q)
        if pub_res["ok"] and pub_res["candidates"]:
            cands = pub_res["candidates"]
            is_ambig = len(cands) > 1
            cache_note = " (من الذاكرة المؤقتة Cache)" if pub_res.get("cache_hit") else ""
            return {
                "found": True,
                "match_type": "ambiguous" if is_ambig else "online_exact",
                "is_ambiguous": is_ambig,
                "ambiguity_message": "توجد عدة مركبات مطابقة في PubChem؛ يرجى اختيار المركب المقصود." if is_ambig else None,
                "candidates": cands,
                "selected": None if is_ambig else cands[0],
                "did_you_mean": [],
                "data_source": "pubchem_online",
                "connection_status": f"متصل بقاعدة بيانات PubChem العالمية (عبر الإنترنت){cache_note}",
                "message": "تم استرجاع بيانات المركب من PubChem والتحقق من البنية والصيغة عبر RDKit.",
            }
        network_error_note = pub_res.get("error")

    if search_type in {"auto", "smiles"} and len(raw_q) >= 2 and not any("\u0600" <= c <= "\u06FF" for c in raw_q):
        mol_try, canon_try, _ = validate_smiles(raw_q)
        if mol_try is not None:
            for entry in LOCAL_COMPOUND_DICTIONARY:
                _, entry_canon, _ = validate_smiles(entry.get("canonical_smiles", ""))
                if entry_canon == canon_try:
                    matched_entry = _decorate_entry(
                        entry,
                        matched_by="exact_smiles",
                        match_category="Exact Match",
                        match_reason_ar=f"تطابق تام مع البنية الخطية SMILES ({canon_try})",
                        match_rank=89,
                    )
                    return {
                        "found": True,
                        "match_type": "exact",
                        "is_ambiguous": False,
                        "ambiguity_message": None,
                        "candidates": [matched_entry],
                        "selected": matched_entry,
                        "did_you_mean": [],
                        "data_source": "local_dictionary",
                        "connection_status": "القاموس المحلي الموثق (تطابق SMILES)",
                        "message": f"تم التعرف على المركب عبر صيغة SMILES الموحدة ({canon_try}).",
                    }

    if fuzzy_matches:
        did_you_mean_list = [f"{m['name_ar']} ({m['name_en']})" for m in fuzzy_matches[:4]]
        status_label = (
            "القاموس المحلي الموثق (المسار البديل عند تعذر PubChem)"
            if network_error_note
            else "القاموس المحلي الموثق (تصحيح إملائي تقريبي للأسماء)"
        )
        return {
            "found": True,
            "match_type": "approximate",
            "is_ambiguous": len(fuzzy_matches) > 1,
            "ambiguity_message": "لم يتم العثور على تطابق حرفي تام؛ هل تقصد أحد المركبات التالية بناء على التشابه الإملائي للاسم؟",
            "candidates": fuzzy_matches[:4],
            "selected": fuzzy_matches[0] if len(fuzzy_matches) == 1 else None,
            "did_you_mean": did_you_mean_list,
            "data_source": "local_fallback" if network_error_note else "local_dictionary",
            "connection_status": status_label,
            "message": f"تم اقتراح أقرب المركبات إملائيا ({fuzzy_matches[0]['match_reason_ar']}): {', '.join(did_you_mean_list)}",
        }

    default_suggestions = ["سيتاجليبتين (Sitagliptin)", "فيلداجليبتين (Vildagliptin)", "ثاني أكسيد الكربون (CO2)", "الماء (H2O)", "ATP"]
    status_msg = (
        f"غير متصل بـ PubChem ({network_error_note}) - تم البحث في القاموس المحلي فقط"
        if network_error_note
        else "القاموس المحلي الموثق"
    )
    return {
        "found": False,
        "match_type": "none",
        "is_ambiguous": False,
        "ambiguity_message": None,
        "candidates": [],
        "selected": None,
        "did_you_mean": default_suggestions,
        "data_source": "none",
        "connection_status": status_msg,
        "message": f"لم يتم العثور على مركب مطابق للمدخل '{raw_q}'. يرجى التحقق من التهجئة أو اختيار أحد المركبات المقترحة.",
    }

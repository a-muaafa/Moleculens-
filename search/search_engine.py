from difflib import SequenceMatcher
import threading
import time
from typing import Any

from chemistry_engine.formula_analyzer import (
    format_chemical_formula_display,
    normalize_chemical_formula,
    parse_chemical_formula,
)
from chemistry_engine.structure_analyzer import HAS_RDKIT, Chem, validate_and_parse_smiles
from core.validators import validate_search_query
from data_integrations.data_manager import DataManager
from search.entity_classifier import classify_scientific_entity
from search.local_catalog import (
    AMBIGUOUS_CATALOG_QUERIES,
    LOCAL_SCIENTIFIC_CATALOG,
)
from search.normalizer import (
    detect_query_input_type,
    formula_key,
    highlight_matched_text,
    is_chemical_formula_input,
    is_element_symbol_query,
    normalize_search_input,
)

AMBIGUITY_SCORE_WINDOW = 8.0
APPROX_MIN_RATIO = 0.74
DEFAULT_DID_YOU_MEAN = [
    "حمض الخليك (CH3COOH)",
    "الماء (H2O)",
    "الحديد (Fe)",
    "الأسبرين (Aspirin)",
    "الجلوكوز (Glucose)",
]

_INDEX_CACHE: list[dict[str, Any]] | None = None

# --------------------------------------------------------------------------------------
# Network debouncing — prevents a burst of external API calls while the user is typing
# --------------------------------------------------------------------------------------

NETWORK_DEBOUNCE_SECONDS = 2.0
_NETWORK_STATE: dict[str, float] = {"last_request": 0.0}
_NETWORK_LOCK = threading.RLock()


def set_network_debounce_seconds(seconds: float) -> None:
    """Sets the minimum interval (seconds) between two consecutive external lookups."""
    global NETWORK_DEBOUNCE_SECONDS
    NETWORK_DEBOUNCE_SECONDS = max(float(seconds), 0.0)


def network_debounce_status() -> dict[str, Any]:
    """Returns the current debounce state (for diagnostics and for the UI notice)."""
    with _NETWORK_LOCK:
        elapsed = time.time() - _NETWORK_STATE["last_request"]
    return {
        "interval_seconds": NETWORK_DEBOUNCE_SECONDS,
        "seconds_since_last_request": round(elapsed, 3),
        "slot_available": elapsed >= NETWORK_DEBOUNCE_SECONDS,
    }


def _claim_network_slot(force: bool = False) -> bool:
    """Claims the right to perform an external lookup, respecting the debounce interval."""
    if force:
        with _NETWORK_LOCK:
            _NETWORK_STATE["last_request"] = time.time()
        return True
    with _NETWORK_LOCK:
        now = time.time()
        if now - _NETWORK_STATE["last_request"] < NETWORK_DEBOUNCE_SECONDS:
            return False
        _NETWORK_STATE["last_request"] = now
        return True


def _canonicalize_smiles(smiles: str) -> str:
    if not smiles:
        return ""
    if HAS_RDKIT and Chem is not None:
        mol = Chem.MolFromSmiles(smiles)
        return Chem.MolToSmiles(mol) if mol is not None else ""
    _, canon, _ = validate_and_parse_smiles(smiles)
    return canon or ""


def _build_index() -> list[dict[str, Any]]:
    index = []
    for entry in LOCAL_SCIENTIFIC_CATALOG:
        formula = entry.get("symbol_or_formula", "")
        smiles = entry.get("canonical_smiles", "")
        canon = _canonicalize_smiles(smiles)
        index.append(
            {
                "entry": entry,
                "norm_ar": normalize_search_input(entry.get("name_ar", "")),
                "norm_en": normalize_search_input(entry.get("name_en", "")),
                "norm_iupac": normalize_search_input(entry.get("iupac_name", "")),
                "fkey": formula_key(formula),
                "hill": normalize_chemical_formula(formula) if entry.get("entity_type") != "biomolecule_protein" else None,
                "abbrevs": [(a, normalize_search_input(a)) for a in entry.get("abbreviations", [])],
                "brands": [(b, normalize_search_input(b)) for b in entry.get("brand_names", [])],
                "syns_ar": [(s, normalize_search_input(s)) for s in entry.get("synonyms_ar", [])],
                "syns_en": [(s, normalize_search_input(s)) for s in entry.get("synonyms_en", [])],
                "canon_smiles": canon,
                "is_ion": entry.get("subclass") == "ion",
            }
        )
    return index


def _get_index() -> list[dict[str, Any]]:
    global _INDEX_CACHE
    if _INDEX_CACHE is None:
        _INDEX_CACHE = _build_index()
    return _INDEX_CACHE


def _decorate_and_unify(
    entry: dict[str, Any],
    matched_by: str,
    match_category: str,
    match_reason_ar: str,
    match_score: float,
    allow_network: bool = False,
    query: str = "",
) -> dict[str, Any]:
    raw_copy = dict(entry)
    classification = classify_scientific_entity(raw_copy)
    raw_copy["entity_type"] = classification["entity_type"]
    raw_copy["matched_by"] = matched_by
    raw_copy["match_category"] = match_category
    raw_copy["match_reason_ar"] = match_reason_ar
    raw_copy["match_score"] = round(float(match_score), 2)
    unified = DataManager.unify_entity_record(raw_copy, allow_network_enrichment=allow_network)
    unified["classification"] = classification
    unified["formatted_formula"] = format_chemical_formula_display(unified["symbol_or_formula"])
    unified["highlighted_name_ar"] = highlight_matched_text(unified["name_ar"], query)
    unified["highlighted_name_en"] = highlight_matched_text(unified["name_en"], query)
    unified["highlighted_formula"] = highlight_matched_text(unified["symbol_or_formula"], query)
    return unified


def _match_entry(raw_q: str, norm_q: str, fkey_q: str, idx: dict[str, Any], query_is_formula: bool) -> tuple[str, str, str, float] | None:
    entry = idx["entry"]
    formula = entry.get("symbol_or_formula", "")

    if norm_q and norm_q == idx["norm_ar"]:
        return ("exact_name_ar", "Exact Match", f"تطابق تام مع الاسم العربي ({entry['name_ar']})", 100.0)
    if norm_q and norm_q == idx["norm_en"]:
        return ("exact_name_en", "Exact Match", f"تطابق تام مع الاسم الإنجليزي ({entry['name_en']})", 99.0)
    if fkey_q and fkey_q == idx["fkey"]:
        return ("exact_formula", "Exact Match", f"تطابق تام مع الصيغة أو الرمز الكيميائي ({formula})", 98.0)
    for orig, nab in idx["abbrevs"]:
        if norm_q == nab and orig.lower() != formula.lower():
            return ("exact_abbreviation", "Exact Match", f"تطابق تام مع الاختصار العلمي ({orig})", 96.0)
    if fkey_q and fkey_q.lower() == idx["fkey"].lower():
        return ("exact_formula", "Exact Match", f"تطابق مع الصيغة أو الرمز الكيميائي بصرف النظر عن حالة الأحرف ({formula})", 88.0)
    for orig, nbr in idx["brands"]:
        if norm_q == nbr:
            return ("exact_brand", "Exact Match", f"تطابق تام مع الاسم التجاري ({orig})", 94.0)
    if norm_q and norm_q == idx["norm_iupac"]:
        return ("exact_iupac", "Exact Match", f"تطابق تام مع الاسم النظامي IUPAC ({entry.get('iupac_name', '')})", 93.0)
    for orig, ns in idx["syns_ar"]:
        if norm_q == ns:
            return ("exact_synonym_ar", "Exact Match", f"تطابق تام مع المرادف العربي ({orig})", 92.0)
    for orig, ns in idx["syns_en"]:
        if norm_q == ns:
            return ("exact_synonym_en", "Exact Match", f"تطابق تام مع المرادف الإنجليزي ({orig})", 91.0)

    if query_is_formula or len(norm_q) < 2:
        return None

    candidates_ar = [(entry["name_ar"], idx["norm_ar"], "العربي", True)] + [(o, n, "العربي", False) for o, n in idx["syns_ar"]]
    candidates_en = [(entry["name_en"], idx["norm_en"], "الإنجليزي", True)] + [(o, n, "الإنجليزي", False) for o, n in idx["syns_en"] + idx["brands"]]
    for orig, n, lang, is_primary in candidates_ar + candidates_en:
        if n.startswith(norm_q):
            score = 82.0 if (is_primary and lang == "العربي") else (80.0 if is_primary else 76.0)
            return ("prefix_name", "Partial Match", f"تطابق في بداية الاسم أو المرادف {lang} ({orig})", score)
    for orig, n, lang, _ in candidates_ar + candidates_en:
        if norm_q in n and not is_chemical_formula_input(orig):
            return ("substring_name", "Partial Match", f"تطابق جزئي مع الاسم أو المرادف {lang} ({orig})", 72.0)
    return None


def _formula_prefix_match(raw_q: str, idx: dict[str, Any]) -> tuple[str, str, str, float] | None:
    q = formula_key(raw_q)
    if len(q) >= 2 and idx["fkey"].lower().startswith(q.lower()) and idx["fkey"].lower() != q.lower():
        return ("prefix_formula", "Partial Match", f"تطابق في بداية الصيغة الكيميائية ({idx['entry']['symbol_or_formula']})", 74.0)
    return None


def _compute_approximate_name_score(raw_q: str, norm_q: str, idx: dict[str, Any]) -> tuple[float, str]:
    if len(norm_q) < 3 or is_chemical_formula_input(raw_q) or is_element_symbol_query(raw_q)[0]:
        return 0.0, ""
    entry = idx["entry"]
    terms = [(entry.get("name_ar", ""), idx["norm_ar"]), (entry.get("name_en", ""), idx["norm_en"])]
    terms += [(o, n) for o, n in idx["brands"] + idx["syns_ar"] + idx["syns_en"] if len(n) >= 4 and not is_chemical_formula_input(o)]
    best_ratio, best_term = 0.0, ""
    for orig, nt in terms:
        if len(nt) < 3:
            continue
        r = SequenceMatcher(None, norm_q, nt, autojunk=False).ratio()
        if r > best_ratio:
            best_ratio, best_term = r, orig
    return best_ratio, best_term


def _smiles_match(raw_q: str, idx_list: list[dict[str, Any]]) -> list[dict[str, Any]]:
    q = raw_q.strip()
    if " " in q or len(q) < 2:
        return []
    canon = _canonicalize_smiles(q)
    if not canon:
        return []
    return [i for i in idx_list if i["canon_smiles"] and i["canon_smiles"] == canon]


def _hill_equivalents(raw_q: str, idx_list: list[dict[str, Any]]) -> list[dict[str, Any]]:
    q_hill = normalize_chemical_formula(raw_q)
    if not q_hill:
        return []
    has_charge = any(ch in formula_key(raw_q) for ch in "+-")
    out = []
    for i in idx_list:
        if i["hill"] != q_hill or i["entry"].get("entity_type") == "chemical_element":
            continue
        if i["is_ion"] and not has_charge:
            continue
        out.append(i)
    return out


def _empty_result(message: str, did_you_mean: list[str] | None = None) -> dict[str, Any]:
    return {
        "found": False,
        "is_ambiguous": False,
        "ambiguity_message": None,
        "selected": None,
        "candidates": [],
        "did_you_mean": did_you_mean if did_you_mean is not None else list(DEFAULT_DID_YOU_MEAN),
        "message": message,
        "input_type": None,
        "network_deferred": False,
    }


def _finalize_group(items: list[dict[str, Any]], ambiguity_message: str, input_type: dict[str, str], did_you_mean: list[str] | None = None) -> dict[str, Any]:
    items.sort(key=lambda x: (-x["match_score"], len(x["name_en"])))
    top = items[0]["match_score"]
    contenders = [i for i in items if top - i["match_score"] <= AMBIGUITY_SCORE_WINDOW]
    is_ambig = len(contenders) > 1
    return {
        "found": True,
        "is_ambiguous": is_ambig,
        "ambiguity_message": ambiguity_message if is_ambig else None,
        "selected": None if is_ambig else items[0],
        "candidates": contenders if is_ambig else items,
        "did_you_mean": did_you_mean or [],
        "message": items[0]["match_reason_ar"],
        "input_type": input_type,
        "network_deferred": False,
    }


def search_scientific_entity(
    query: str,
    search_type: str = "auto",
    allow_network: bool = True,
    force_network: bool = False,
) -> dict[str, Any]:
    is_valid, cleaned_q, err_msg = validate_search_query(query)
    if not is_valid:
        return _empty_result(err_msg)

    input_type = detect_query_input_type(cleaned_q)
    norm_q = normalize_search_input(cleaned_q)
    fkey_q = formula_key(cleaned_q)
    query_is_formula = is_chemical_formula_input(cleaned_q) or is_element_symbol_query(cleaned_q)[0]
    index = _get_index()

    ambig = AMBIGUOUS_CATALOG_QUERIES.get(fkey_q.lower())
    if ambig:
        cands = [
            _decorate_and_unify(
                i["entry"], "exact_formula", "Exact Match",
                f"تطابق تام مع الصيغة الكيميائية ({i['entry']['symbol_or_formula']}) ضمن مصاوغات متعددة", 98.0, query=cleaned_q,
            )
            for i in index
            if i["entry"]["key"] in ambig["candidate_keys"]
        ]
        return {
            "found": True, "is_ambiguous": True, "ambiguity_message": ambig["reason_ar"], "selected": None,
            "candidates": cands, "did_you_mean": [], "message": ambig["reason_ar"], "input_type": input_type,
            "network_deferred": False,
        }

    exact: list[dict[str, Any]] = []
    partial: list[dict[str, Any]] = []
    approx: list[dict[str, Any]] = []
    for idx in index:
        m = _match_entry(cleaned_q, norm_q, fkey_q, idx, query_is_formula)
        if m is None:
            if query_is_formula:
                continue
            ratio, term = _compute_approximate_name_score(cleaned_q, norm_q, idx)
            if ratio >= APPROX_MIN_RATIO:
                pct = int(round(ratio * 100))
                approx.append(
                    _decorate_and_unify(idx["entry"], "approximate_name", "Approximate Match",
                                        f"تطابق إملائي تقريبي مع ({term} - تشابه {pct}%)", round(ratio * 68.0, 2), query=cleaned_q)
                )
            continue
        by, cat, reason, score = m
        item = _decorate_and_unify(idx["entry"], by, cat, reason, score, query=cleaned_q)
        (exact if cat == "Exact Match" else partial).append(item)

    if exact:
        return _finalize_group(exact, "يحتمل المدخل أكثر من مادة علمية؛ يرجى اختيار المادة المقصودة من القائمة أدناه.", input_type)

    smiles_hits = _smiles_match(cleaned_q, index)
    if smiles_hits:
        items = [
            _decorate_and_unify(i["entry"], "exact_smiles", "Exact Match", "تطابق تام مع البنية الجزيئية (SMILES) بعد التوحيد القياسي عبر RDKit", 97.0, query=cleaned_q)
            for i in smiles_hits
        ]
        return _finalize_group(items, "يطابق المدخل أكثر من مادة؛ يرجى اختيار المادة المقصودة.", input_type)

    if query_is_formula:
        hill_hits = _hill_equivalents(cleaned_q, index)
        if hill_hits:
            items = [
                _decorate_and_unify(
                    i["entry"], "formula_isomer", "Exact Match",
                    f"تطابق في الصيغة الجزيئية فقط ({normalize_chemical_formula(cleaned_q)}) وقد تختلف البنية بين المصاوغات", 90.0, query=cleaned_q,
                )
                for i in hill_hits
            ]
            res = _finalize_group(items, "تشترك عدة مواد في هذه الصيغة الجزيئية وتختلف في بنيتها الفراغية/البنائية (Isomers)؛ يرجى اختيار المادة المقصودة.", input_type)
            if len(items) > 1:
                res["is_ambiguous"] = True
                res["selected"] = None
                res["candidates"] = items
                res["ambiguity_message"] = "تشترك عدة مواد في هذه الصيغة الجزيئية وتختلف في بنيتها الفراغية/البنائية (Isomers)؛ يرجى اختيار المادة المقصودة."
            return res

    if partial:
        res = _finalize_group(partial, "تم العثور على عدة نتائج متقاربة؛ يرجى اختيار المادة المقصودة بدقة.", input_type,
                              did_you_mean=[f"{c['name_ar']} ({c['symbol_or_formula']})" for c in sorted(partial, key=lambda x: -x['match_score'])[:4]])
        res["is_ambiguous"] = len(partial) > 1
        res["selected"] = None if len(partial) > 1 else partial[0]
        res["candidates"] = sorted(partial, key=lambda x: -x["match_score"])
        if res["is_ambiguous"]:
            res["ambiguity_message"] = "تم العثور على عدة نتائج متقاربة؛ يرجى اختيار المادة المقصودة بدقة."
        return res

    network_allowed = allow_network and _claim_network_slot(force=force_network)
    ext = DataManager.resolve_external_substances(cleaned_q, is_formula=query_is_formula, allow_network=network_allowed)
    if ext.get("ok") and ext.get("records"):
        origin_ar = "قاعدة البيانات المحلية الذكية" if ext["origin"] == "local_database" else "المصادر العلمية الخارجية (PubChem / ChEMBL / ChEBI)"
        ext_cands = [
            _decorate_and_unify(rec, "external_source", "Exact Match", f"تم استرجاع المادة والتحقق من هويتها وبنيتها عبر {origin_ar}", 95.0, query=cleaned_q)
            for rec in ext["records"]
        ]
        is_ambig = len(ext_cands) > 1
        return {
            "found": True, "is_ambiguous": is_ambig,
            "ambiguity_message": "تتعدد المواد المطابقة؛ يرجى اختيار المادة المقصودة." if is_ambig else None,
            "selected": None if is_ambig else ext_cands[0], "candidates": ext_cands, "did_you_mean": [],
            "message": "تم جلب بيانات المادة والتحقق من هويتها وتوحيدها عبر Data Manager وتخزينها في قاعدة البيانات المحلية.", "input_type": input_type,
            "network_deferred": False,
        }

    if approx:
        approx.sort(key=lambda x: -x["match_score"])
        return {
            "found": True, "is_ambiguous": len(approx) > 1,
            "ambiguity_message": "لم يعثر على تطابق حرفي؛ هل تقصد إحدى المواد التالية؟" if len(approx) > 1 else None,
            "selected": approx[0] if len(approx) == 1 else None, "candidates": approx[:4],
            "did_you_mean": [f"{c['name_ar']} ({c['name_en']})" for c in approx[:4]],
            "message": approx[0]["match_reason_ar"], "input_type": input_type,
            "network_deferred": False,
        }

    res = _empty_result(
        f"لم يتم العثور على مادة علمية مطابقة للمدخل '{cleaned_q}'. يرجى التأكد من صحة الإملاء أو الصيغة الكيميائية أو اختيار إحدى المقترحات أدناه."
    )
    res["input_type"] = input_type
    res["network_deferred"] = bool(allow_network and not network_allowed)
    if res["network_deferred"]:
        res["message"] = (
            "لم يتم العثور على المادة في السجل المحلي، وقد تم تأجيل الاستعلام عن المصادر الخارجية مؤقتاً "
            "لتقليل الطلبات المتكررة أثناء الكتابة. اضغط «تحديث من المصادر الخارجية» لجلب البيانات الآن."
        )
    return res


def get_live_suggestions(query: str, limit: int = 8) -> list[dict[str, Any]]:
    index = _get_index()
    if not query or not query.strip():
        picks = ["acetic_acid", "water_h2o", "sodium_hydroxide_naoh", "sodium_chloride_nacl", "benzene", "ethanol", "aspirin", "d_galactose"]
        return [
            _decorate_and_unify(i["entry"], "default", "Exact Match", "مادة مقترحة للاستكشاف", 100.0)
            for i in index if i["entry"]["key"] in picks
        ][:limit]

    cleaned = query.strip()
    norm_q = normalize_search_input(cleaned)
    fkey_q = formula_key(cleaned)
    query_is_formula = is_chemical_formula_input(cleaned) or is_element_symbol_query(cleaned)[0]
    items: dict[str, dict[str, Any]] = {}

    for idx in index:
        m = _match_entry(cleaned, norm_q, fkey_q, idx, False) or _formula_prefix_match(cleaned, idx)
        if m is None and len(norm_q) >= 3 and not query_is_formula:
            ratio, term = _compute_approximate_name_score(cleaned, norm_q, idx)
            if ratio >= APPROX_MIN_RATIO:
                m = ("approximate_name", "Approximate Match", f"قريب إملائيا من ({term})", round(ratio * 68.0, 2))
        if m is None:
            continue
        by, cat, reason, score = m
        key = idx["entry"]["key"]
        if key not in items or items[key]["match_score"] < score:
            items[key] = _decorate_and_unify(idx["entry"], by, cat, reason, score, query=cleaned)

    if len(items) < limit:
        try:
            from data_integrations.local_db import get_local_database

            for rec in get_local_database().search_prefix(cleaned, limit=limit):
                if rec.get("key") not in items and rec.get("db_source", "").startswith("External"):
                    items[rec["key"]] = _decorate_and_unify(rec, "local_database", "Partial Match", "مادة محفوظة محليا من مصدر علمي خارجي", 70.0, query=cleaned)
        except Exception:
            pass

    ranked = sorted(
        items.values(),
        key=lambda x: (
            -x["match_score"],
            0 if x["name_en"].lower().startswith(norm_q) else 1,
            len(x["name_en"].split()[0]),
            len(x["name_en"]),
        ),
    )
    return ranked[:limit]


def get_catalog_dropdown_options() -> list[dict[str, Any]]:
    items = []
    for i in _get_index():
        entry = i["entry"]
        u = _decorate_and_unify(entry, "catalog", "Exact Match", "اختيار مباشر", 100.0)
        f_disp = format_chemical_formula_display(u["symbol_or_formula"])
        u["dropdown_label"] = f"{u['name_ar']} | {u['name_en']} ({f_disp}) — [{u['entity_type_label_ar']}]"
        items.append(u)
    return items


REACTION_CATEGORY_SUBCLASSES = {
    "all": None,
    "element": {"element"},
    "acid": {"acid", "carboxylic_acid"},
    "base": {"base"},
    "salt": {"salt", "ion"},
    "inorganic_compound": {"inorganic_compound"},
    "organic_compound": {"organic_compound", "aromatic_compound", "alcohol", "aldehyde", "ketone", "ester"},
    "drug": {"drug"},
    "biomolecule": {"biomolecule"},
}


def get_reaction_substance_options(category: str = "all") -> list[dict[str, Any]]:
    allowed = REACTION_CATEGORY_SUBCLASSES.get(category)
    options = []
    for i in _get_index():
        entry = i["entry"]
        if parse_chemical_formula(entry["symbol_or_formula"]) is None:
            continue
        if allowed is not None and entry.get("subclass") not in allowed:
            continue
        f_disp = format_chemical_formula_display(entry["symbol_or_formula"])
        options.append(
            {
                "key": entry["key"],
                "formula": entry["symbol_or_formula"],
                "formatted_formula": f_disp,
                "name_ar": entry["name_ar"],
                "name_en": entry["name_en"],
                "subclass": entry.get("subclass", ""),
                "label": f"{entry['name_ar']} | {entry['name_en']} ({f_disp})",
            }
        )
    return options

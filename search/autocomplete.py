from search.compound_search import (
    _compute_fuzzy_name_score,
    _decorate_entry,
    _evaluate_local_entry_match,
    normalize_text,
)
from search.local_dictionary import (
    LOCAL_COMPOUND_DICTIONARY,
    formulas_match,
    is_chemical_formula_query,
    replace_subscripts,
)


def format_suggestion_label(entry: dict) -> str:
    formula = entry.get("molecular_formula", "")
    cid = entry.get("pubchem_cid") or "N/A"
    chembl = entry.get("chembl_id") or "N/A"
    reason = entry.get("match_reason_ar", "")
    base = f"{entry['name_ar']} | {entry['name_en']} ({formula}) [CID: {cid} | {chembl}]"
    if reason:
        return f"{base} — {reason}"
    return base


def build_instant_search_catalog() -> list[dict]:
    catalog = []
    for entry in LOCAL_COMPOUND_DICTIONARY:
        item = _decorate_entry(
            entry,
            matched_by="catalog_entry",
            match_category="Exact Match",
            match_reason_ar="اختيار مباشر من فهرس القاموس المحلي الموثق",
            match_rank=100,
        )
        abbrevs = ", ".join(entry.get("abbreviations", []))
        brands = ", ".join(entry.get("brand_names", []))
        syns_ar = "، ".join(entry.get("synonyms_ar", [])[:3])
        syns_en = ", ".join(entry.get("synonyms_en", [])[:3])
        extra_tags = " | ".join(filter(None, [abbrevs, brands, syns_ar, syns_en]))
        cid = entry.get("pubchem_cid") or "N/A"
        chembl = entry.get("chembl_id") or "N/A"
        item["instant_label"] = (
            f"{entry['name_ar']} | {entry['name_en']} | الصيغة: {entry.get('molecular_formula', '')} "
            f"| CID: {cid} | {chembl}" + (f" | ({extra_tags})" if extra_tags else "")
        )
        catalog.append(item)
    return catalog


def get_autocomplete_suggestions(query: str, limit: int = 6, search_type: str = "auto") -> list[dict]:
    raw_q = query.strip() if query else ""
    if not raw_q:
        default_items = []
        for entry in LOCAL_COMPOUND_DICTIONARY[:limit]:
            item = _decorate_entry(
                entry,
                matched_by="default_suggestion",
                match_category="Exact Match",
                match_reason_ar="اقتراح افتراضي من القاموس المحلي الموثق",
                match_rank=50,
            )
            item["display_label"] = format_suggestion_label(item)
            default_items.append(item)
        return default_items

    cleaned_q = replace_subscripts(raw_q).strip()
    norm_q = normalize_text(cleaned_q)
    is_formula_mode = search_type == "formula" or (search_type == "auto" and is_chemical_formula_query(cleaned_q))

    scored: list[tuple[float, dict]] = []
    seen_keys: set[str] = set()

    if is_formula_mode:
        clean_formula_q = cleaned_q.replace(" ", "").lower()
        for entry in LOCAL_COMPOUND_DICTIONARY:
            entry_formula = entry.get("molecular_formula", "")
            entry_formula_low = entry_formula.lower()
            abbrevs_low = [ab.lower() for ab in entry.get("abbreviations", []) if is_chemical_formula_query(ab)]
            if formulas_match(cleaned_q, entry_formula) or clean_formula_q in abbrevs_low:
                item = _decorate_entry(
                    entry,
                    matched_by="exact_formula",
                    match_category="Exact Match",
                    match_reason_ar=f"تطابق تام مع الصيغة الكيميائية ({entry_formula})",
                    match_rank=98.0,
                )
                item["display_label"] = format_suggestion_label(item)
                scored.append((98.0, item))
                seen_keys.add(entry["key"])
            elif entry_formula_low.startswith(clean_formula_q) or any(ab.startswith(clean_formula_q) for ab in abbrevs_low):
                item = _decorate_entry(
                    entry,
                    matched_by="prefix_formula",
                    match_category="Partial Match",
                    match_reason_ar=f"تطابق في بداية الصيغة الكيميائية ({entry_formula})",
                    match_rank=84.0,
                )
                item["display_label"] = format_suggestion_label(item)
                scored.append((84.0, item))
                seen_keys.add(entry["key"])

        scored.sort(key=lambda x: x[0], reverse=True)
        return [item for _, item in scored[:limit]]

    for entry in LOCAL_COMPOUND_DICTIONARY:
        matched_item = _evaluate_local_entry_match(cleaned_q, norm_q, entry, search_type=search_type)
        if matched_item is not None:
            matched_item["display_label"] = format_suggestion_label(matched_item)
            scored.append((matched_item["match_score"], matched_item))
            seen_keys.add(entry["key"])
            continue

        if search_type == "auto":
            entry_formula = entry.get("molecular_formula", "")
            if entry_formula.lower().startswith(cleaned_q.lower()):
                item = _decorate_entry(
                    entry,
                    matched_by="prefix_formula",
                    match_category="Partial Match",
                    match_reason_ar=f"تطابق في بداية الصيغة الكيميائية ({entry_formula})",
                    match_rank=82.0,
                )
                item["display_label"] = format_suggestion_label(item)
                scored.append((82.0, item))
                seen_keys.add(entry["key"])
                continue

        if len(norm_q) >= 3 and search_type not in {"formula", "smiles", "abbreviation"}:
            f_score, f_field, f_term = _compute_fuzzy_name_score(norm_q, entry, search_type=search_type)
            if f_score >= 0.65:
                pct = int(round(f_score * 100))
                rank_val = round(f_score * 65.0, 2)
                item = _decorate_entry(
                    entry,
                    matched_by="approximate_name",
                    match_category="Did You Mean",
                    match_reason_ar=f"تصحيح إملائي مقترح لـ {f_field} ({f_term} - تشابه {pct}%)",
                    match_rank=rank_val,
                )
                item["display_label"] = format_suggestion_label(item)
                scored.append((rank_val, item))
                seen_keys.add(entry["key"])

    scored.sort(key=lambda x: x[0], reverse=True)
    return [item for _, item in scored[:limit]]

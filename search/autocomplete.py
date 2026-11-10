from search.local_dictionary import LOCAL_COMPOUND_DICTIONARY, normalize_arabic_text


def get_autocomplete_suggestions(prefix: str, limit: int = 8) -> list[dict]:
    if not prefix or not isinstance(prefix, str) or not prefix.strip():
        return [dict(item) for item in LOCAL_COMPOUND_DICTIONARY[:limit]]

    raw_p = prefix.strip()
    low_p = raw_p.lower()
    norm_p = normalize_arabic_text(raw_p)

    prefix_matches = []
    substring_matches = []

    for item in LOCAL_COMPOUND_DICTIONARY:
        name_ar_norm = normalize_arabic_text(item["name_ar"])
        name_en_low = item["name_en"].lower()
        formula_low = item.get("molecular_formula", "").lower()
        abbrs_low = [a.lower() for a in item.get("abbreviations", [])]
        syns_ar_norm = [normalize_arabic_text(s) for s in item.get("synonyms_ar", [])]
        syns_en_low = [s.lower() for s in item.get("synonyms_en", [])]
        brands_low = [b.lower() for b in item.get("brand_names", [])]
        brands_ar_norm = [normalize_arabic_text(b) for b in item.get("brand_names", [])]

        is_prefix = (
            (norm_p and name_ar_norm.startswith(norm_p))
            or name_en_low.startswith(low_p)
            or formula_low.startswith(low_p)
            or any(a.startswith(low_p) for a in abbrs_low)
            or any(norm_p and s.startswith(norm_p) for s in syns_ar_norm)
            or any(s.startswith(low_p) for s in syns_en_low)
            or any(b.startswith(low_p) for b in brands_low)
            or any(norm_p and b.startswith(norm_p) for b in brands_ar_norm)
        )

        if is_prefix:
            prefix_matches.append(dict(item))
            continue

        is_sub = (
            (norm_p and norm_p in name_ar_norm)
            or (low_p in name_en_low)
            or (low_p in formula_low)
            or (low_p in item.get("iupac_name", "").lower())
            or any(low_p in a for a in abbrs_low)
            or any(norm_p and norm_p in s for s in syns_ar_norm)
            or any(low_p in s for s in syns_en_low)
            or any(low_p in b for b in brands_low)
        )
        if is_sub:
            substring_matches.append(dict(item))

    combined = prefix_matches + substring_matches
    return combined[:limit]


def format_suggestion_label(entry: dict) -> str:
    cid_part = f"CID: {entry['pubchem_cid']}" if entry.get("pubchem_cid") else "بدون CID"
    chembl_part = f" | ChEMBL: {entry['chembl_id']}" if entry.get("chembl_id") else ""
    return f"{entry['name_ar']} | {entry['name_en']} ({entry['molecular_formula']}) [{cid_part}{chembl_part}]"


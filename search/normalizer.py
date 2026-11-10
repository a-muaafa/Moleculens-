import html
import re
from difflib import SequenceMatcher
from chemistry_engine.formula_analyzer import (
    SPECIAL_FORMULA_CASE_MAP,
    SUPERSCRIPT_MAP,
    parse_chemical_formula,
    replace_subscripts,
)
from chemistry_engine.periodic_table import VALID_ELEMENTS

KNOWN_FORMULA_STRINGS = {
    "H2", "O2", "N2", "Cl2", "H2O", "CO2", "CO", "NO", "NO2", "SO2",
    "H2O2", "NH3", "CH4", "HCl", "NaCl", "NaOH", "KOH", "CaO", "MgO",
    "Fe2O3", "Fe3O4", "H2SO4", "HNO3", "H2CO3", "C6H12O6", "C2H6O",
    "C2H4O2", "C2H5OH", "CH3COOH", "CaCO3",
}


def normalize_arabic_text(text: str) -> str:
    if not text:
        return ""
    s = replace_subscripts(text.strip())
    s = re.sub(r"[\u064B-\u065F\u0670]", "", s)
    s = re.sub(r"[أإآٱ]", "ا", s)
    s = s.replace("ة", "ه").replace("ى", "ي").replace("ؤ", "و").replace("ئ", "ي")
    s = re.sub(r"\s+", " ", s)
    return s.lower()


def normalize_search_input(text: str) -> str:
    if not text:
        return ""
    return normalize_arabic_text(text)


def is_element_symbol_query(query: str) -> tuple[bool, str]:
    if not query or not isinstance(query, str):
        return False, ""
    cleaned = replace_subscripts(query.strip())
    if len(cleaned) in {1, 2} and cleaned.isalpha():
        formatted = cleaned[0].upper() + cleaned[1:].lower()
        if formatted in VALID_ELEMENTS:
            return True, formatted
    return False, ""


def formula_key(text: str) -> str:
    if not text:
        return ""
    return replace_subscripts(text).translate(SUPERSCRIPT_MAP).replace(" ", "").replace("^", "")


def _count_element_tokens(cleaned: str) -> int:
    return len(re.findall(r"[A-Z][a-z]?", cleaned))


def is_chemical_formula_input(query: str) -> bool:
    if not query or not isinstance(query, str):
        return False
    cleaned = formula_key(query.strip())
    if not cleaned:
        return False
    is_el, _ = is_element_symbol_query(cleaned)
    if is_el:
        return False
    if cleaned.upper() in {f.upper() for f in KNOWN_FORMULA_STRINGS}:
        return True
    if re.search(r"[\u0600-\u06FF]", cleaned):
        return False
    if any(ch.isdigit() for ch in cleaned):
        return parse_chemical_formula(cleaned) is not None
    if cleaned.lower() in SPECIAL_FORMULA_CASE_MAP:
        return True
    if cleaned[0].isupper() and _count_element_tokens(cleaned) >= 2 and re.fullmatch(r"[A-Za-z()+\-]+", cleaned):
        return parse_chemical_formula(cleaned) is not None
    return False


def detect_query_input_type(query: str) -> dict[str, str]:
    raw = (query or "").strip()
    if not raw:
        return {"type": "empty", "label_ar": "لا يوجد مدخل"}
    if re.search(r"[\u0600-\u06FF]", raw):
        return {"type": "arabic_name", "label_ar": "اسم عربي أو مرادف عربي"}
    is_el, _ = is_element_symbol_query(raw)
    if is_el:
        return {"type": "element_symbol", "label_ar": "رمز عنصر كيميائي"}
    key = formula_key(raw)
    if is_chemical_formula_input(raw) or re.fullmatch(r"[A-Za-z0-9()]+[+\-]?\d*[+\-]?", key) and any(c.isdigit() for c in key) and parse_chemical_formula(key):
        tokens = re.findall(r"[A-Z][a-z]?", key)
        if "(" in key or len(tokens) != len(set(tokens)):
            return {"type": "condensed_formula", "label_ar": "صيغة كيميائية بنائية مكثفة"}
        return {"type": "chemical_formula", "label_ar": "صيغة كيميائية جزيئية"}
    if re.search(r"[=#@\[\]\\/]", raw) and " " not in raw:
        return {"type": "smiles", "label_ar": "تمثيل SMILES الجزيئي"}
    if re.search(r"\d", raw) and re.search(r"[,\-]", raw) and re.search(r"[a-z]{3,}", raw):
        return {"type": "iupac_name", "label_ar": "اسم نظامي (IUPAC)"}
    if raw.replace("-", "").isalnum() and (raw.isupper() or any(c.isdigit() for c in raw)) and len(raw) <= 8:
        return {"type": "abbreviation", "label_ar": "اختصار علمي أو اسم تجاري"}
    return {"type": "english_name", "label_ar": "اسم إنجليزي أو اسم علمي أو مرادف"}


def _norm_char(ch: str) -> str:
    c = replace_subscripts(ch)
    c = c.translate(SUPERSCRIPT_MAP)
    if re.match(r"[\u064B-\u065F\u0670]", c):
        return ""
    if c in "أإآٱ":
        return "ا"
    return {"ة": "ه", "ى": "ي", "ؤ": "و", "ئ": "ي"}.get(c, c).lower()


def highlight_matched_text(text: str, query: str) -> str:
    if not text:
        return ""
    if not query or not query.strip():
        return html.escape(text)

    kept: list[tuple[int, str]] = []
    for i, ch in enumerate(text):
        n = _norm_char(ch)
        for nc in n:
            kept.append((i, nc))
    norm_text = "".join(c for _, c in kept)
    norm_q = "".join(_norm_char(c) for c in query.strip())
    ranges: list[tuple[int, int]] = []

    def add_range(start: int, length: int) -> None:
        if length <= 0 or start < 0:
            return
        s = kept[start][0]
        e = kept[start + length - 1][0] + 1
        ranges.append((s, e))

    idx = norm_text.find(norm_q) if norm_q else -1
    if idx >= 0:
        add_range(idx, len(norm_q))
    else:
        words = [w for w in re.split(r"\s+", norm_q) if len(w) >= 2]
        for w in words:
            pos = norm_text.find(w)
            if pos >= 0:
                add_range(pos, len(w))
        if not ranges and len(norm_q) >= 3:
            sm = SequenceMatcher(None, norm_q, norm_text, autojunk=False)
            for block in sm.get_matching_blocks():
                if block.size >= 2:
                    add_range(block.b, block.size)

    if not ranges:
        return html.escape(text)
    ranges.sort()
    merged: list[list[int]] = []
    for s, e in ranges:
        if merged and s <= merged[-1][1]:
            merged[-1][1] = max(merged[-1][1], e)
        else:
            merged.append([s, e])
    out = []
    cursor = 0
    for s, e in merged:
        out.append(html.escape(text[cursor:s]))
        out.append(f'<mark class="ml-highlight">{html.escape(text[s:e])}</mark>')
        cursor = e
    out.append(html.escape(text[cursor:]))
    return "".join(out)


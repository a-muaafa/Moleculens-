"""
Security, Input Validation, and External Data Sanitization utilities for MolecuLens.

Guarantees:
1. Prevention of HTML/JavaScript/XSS execution from user inputs or external API records.
2. Pre-validation of chemical strings (SMILES, formulas, names) before passing to RDKit
   or stoichiometric engines (preventing DoS via oversized inputs or control characters).
3. Strict HTTPS-only verification for external scientific API endpoints.
"""

from __future__ import annotations

import html
import re
from typing import Any

_CONTROL_CHARS_RE = re.compile(r"[\x00-\x08\x0b\x0c\x0e-\x1f\x7f]")
_SCRIPT_TAG_RE = re.compile(r"<\s*/?\s*(script|iframe|object|embed|style|link|meta)\b[^>]*>", re.IGNORECASE)
_EVENT_HANDLER_RE = re.compile(r"\bon[a-z]+\s*=", re.IGNORECASE)
_JS_PROTOCOL_RE = re.compile(r"(javascript|vbscript|data\s*:\s*text/html)\s*:", re.IGNORECASE)


def contains_unsafe_markup(text: str) -> bool:
    """Returns True if text contains HTML script/iframe tags, inline event handlers, or JS URLs."""
    if not text or not isinstance(text, str):
        return False
    return bool(
        _SCRIPT_TAG_RE.search(text)
        or _EVENT_HANDLER_RE.search(text)
        or _JS_PROTOCOL_RE.search(text)
    )


def sanitize_external_text(value: Any, max_len: int = 2000) -> str:
    """
    Sanitizes a string from an external API or user input so it cannot execute HTML or JS.
    Strips control characters, removes dangerous tags/protocols, truncates to max_len,
    and escapes HTML entities.
    """
    if value is None:
        return ""
    raw = str(value)
    cleaned = _CONTROL_CHARS_RE.sub("", raw)
    cleaned = _SCRIPT_TAG_RE.sub("", cleaned)
    cleaned = _EVENT_HANDLER_RE.sub("", cleaned)
    cleaned = _JS_PROTOCOL_RE.sub("", cleaned)
    cleaned = cleaned.strip()
    if len(cleaned) > max_len:
        cleaned = cleaned[:max_len]
    return html.escape(cleaned, quote=True)


def sanitize_external_record(obj: Any, max_len: int = 2000) -> Any:
    """
    Recursively sanitizes string fields in a dictionary or list fetched from an external API.
    """
    if isinstance(obj, str):
        return sanitize_external_text(obj, max_len=max_len)
    if isinstance(obj, dict):
        return {k: sanitize_external_record(v, max_len=max_len) for k, v in obj.items()}
    if isinstance(obj, list):
        return [sanitize_external_record(item, max_len=max_len) for item in obj]
    return obj


def validate_safe_chemical_input(
    raw: Any,
    max_len: int = 300,
    field_name_ar: str = "المدخل الكيميائي",
) -> tuple[bool, str, str]:
    """
    Validates a chemical identifier (SMILES, formula, or name) before passing it to RDKit
    or calculation engines.
    Returns (is_valid, cleaned_string, error_message_ar).
    """
    if raw is None or not isinstance(raw, str):
        return False, "", f"{field_name_ar} غير صالح: يجب أن يكون نصاً."
    cleaned = _CONTROL_CHARS_RE.sub("", raw).strip()
    if not cleaned:
        return False, "", f"{field_name_ar} فارغ."
    if len(cleaned) > max_len:
        return (
            False,
            "",
            f"{field_name_ar} يتجاوز الطول الآمن المسموح ({max_len} حرفاً).",
        )
    if contains_unsafe_markup(cleaned):
        return (
            False,
            "",
            f"رُفض {field_name_ar} لاحتوائه على وسوم HTML أو نصوص برمجية غير مسموحة.",
        )
    return True, cleaned, ""


def is_safe_https_url(url: str) -> bool:
    """Ensures an external URL uses https:// and does not point to local/private loopback."""
    if not url or not isinstance(url, str):
        return False
    u = url.strip().lower()
    if not u.startswith("https://"):
        return False
    if any(host in u for host in ("localhost", "127.0.0.1", "0.0.0.0", "[::1]")):
        return False
    return True

"""
Resilient HTTP client for MolecuLens external scientific integrations.

Provides:
1. Short configurable timeout (default 3s) and bounded retry on transient 5xx/timeout errors.
2. Fast fail on SSL certificate/handshake errors so the app falls back to local data immediately.
3. Strict HTTPS URL validation (SSRF protection).
4. Clear Arabic error classification that never leaks raw tracebacks to the UI.
"""

from __future__ import annotations

import time
from typing import Any
import requests

from config.settings import EXTERNAL_API_MAX_RETRIES, EXTERNAL_API_TIMEOUT_SECONDS
from core.security import is_safe_https_url

RETRYABLE_STATUS_CODES: set[int] = {429, 500, 502, 503, 504}


def classify_network_error_ar(exc: Exception, service_name: str = "الخدمة الخارجية") -> str:
    """Maps network/SSL/timeout exceptions to clear Arabic messages for offline-first fallback."""
    if isinstance(exc, requests.exceptions.SSLError):
        return f"تعذر التحقق من شهادة الأمان (SSL) عند الاتصال بـ {service_name}؛ تم التراجع تلقائياً إلى البيانات المحلية."
    if isinstance(exc, requests.exceptions.Timeout):
        return f"انتهت مهلة الاتصال بـ {service_name}؛ تم التراجع تلقائياً إلى البيانات المحلية دون تعطيل التطبيق."
    if isinstance(exc, requests.exceptions.ConnectionError):
        return f"لا يتوفر اتصال بالإنترنت أو تعذر الوصول إلى {service_name}؛ يعمل التطبيق الآن بالبيانات المحلية."
    return f"تعذر جلب البيانات من {service_name} مؤقتاً؛ تم الاعتماد على البيانات المحلية."


def request_with_retry(
    url: str,
    params: dict[str, Any] | None = None,
    headers: dict[str, str] | None = None,
    timeout: float = EXTERNAL_API_TIMEOUT_SECONDS,
    max_retries: int = EXTERNAL_API_MAX_RETRIES,
    backoff_seconds: float = 0.2,
) -> requests.Response:
    """
    Executes an HTTP GET request with short timeout and bounded retry on transient failures.
    Fails fast on SSLError so offline-first fallback activates immediately without hanging.
    Raises ValueError on non-HTTPS URLs and requests.RequestException if all attempts fail.
    """
    if not is_safe_https_url(url):
        raise ValueError(f"Unsafe or non-HTTPS URL rejected: {url}")

    attempts = max(1, int(max_retries) + 1)
    last_exc: Exception | None = None

    for attempt in range(attempts):
        try:
            resp = requests.get(url, params=params, headers=headers, timeout=timeout)
            if resp.status_code in RETRYABLE_STATUS_CODES and attempt < attempts - 1:
                time.sleep(backoff_seconds * (2 ** attempt))
                continue
            return resp
        except requests.exceptions.SSLError as exc:
            # Fail fast on SSL errors — retrying won't fix certificate/proxy issues during a demo.
            raise exc
        except (requests.exceptions.Timeout, requests.exceptions.ConnectionError) as exc:
            last_exc = exc
            if attempt < attempts - 1:
                time.sleep(backoff_seconds * (2 ** attempt))
                continue
            raise
        except Exception as exc:
            last_exc = exc
            raise

    if last_exc is not None:
        raise last_exc
    raise requests.exceptions.RequestException("HTTP request failed after retries")

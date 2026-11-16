"""
Resilient HTTP client for MolecuLens external scientific integrations.

Provides:
1. Configurable timeout and bounded retry on transient network/5xx errors.
2. Strict HTTPS URL validation.
3. Safe error handling that never crashes the application when offline.
"""

from __future__ import annotations

import time
from typing import Any
import requests

from config.settings import EXTERNAL_API_MAX_RETRIES, EXTERNAL_API_TIMEOUT_SECONDS
from core.security import is_safe_https_url

RETRYABLE_STATUS_CODES: set[int] = {429, 500, 502, 503, 504}


def request_with_retry(
    url: str,
    params: dict[str, Any] | None = None,
    headers: dict[str, str] | None = None,
    timeout: float = EXTERNAL_API_TIMEOUT_SECONDS,
    max_retries: int = EXTERNAL_API_MAX_RETRIES,
    backoff_seconds: float = 0.25,
) -> requests.Response:
    """
    Executes an HTTP GET request with timeout and bounded retry on transient failures.
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

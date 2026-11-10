import time
from typing import Any
from config.settings import CACHE_TTL_SECONDS

_MEMORY_CACHE: dict[str, tuple[float, Any]] = {}


def get_cached_item(cache_key: str, ttl_seconds: int = CACHE_TTL_SECONDS) -> tuple[bool, Any]:
    if cache_key in _MEMORY_CACHE:
        ts, data = _MEMORY_CACHE[cache_key]
        if time.time() - ts < ttl_seconds:
            return True, data
        del _MEMORY_CACHE[cache_key]
    return False, None


def set_cached_item(cache_key: str, value: Any) -> None:
    _MEMORY_CACHE[cache_key] = (time.time(), value)


def clear_all_caches() -> None:
    _MEMORY_CACHE.clear()


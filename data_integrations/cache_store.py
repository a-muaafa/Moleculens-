"""
Thread-safe, bounded in-memory cache for external scientific API calls.

Design goals (platform engineering):
  1. Bounded memory: an LRU eviction policy prevents unbounded growth when a user
     explores many substances in one session.
  2. Thread safety: Streamlit executes sessions in threads, so the cache is guarded
     by a re-entrant lock.
  3. Observability: `cache_statistics()` exposes hit/miss counters for diagnostics.
  4. Optional on-disk persistence: `enable_disk_persistence()` keeps hot entries in a
     small JSON file so that repeated sessions avoid re-fetching the same records.
"""

from __future__ import annotations

import json
import os
import threading
import time
from collections import OrderedDict
from typing import Any

from config.settings import CACHE_TTL_SECONDS

DEFAULT_MAX_ENTRIES = 512

_MEMORY_CACHE: "OrderedDict[str, tuple[float, Any]]" = OrderedDict()
_CACHE_LOCK = threading.RLock()
_STATS: dict[str, int] = {"hits": 0, "misses": 0, "evictions": 0, "sets": 0, "expirations": 0}
_MAX_ENTRIES = DEFAULT_MAX_ENTRIES
_DISK_CACHE_PATH: str | None = None


def configure_cache(max_entries: int = DEFAULT_MAX_ENTRIES) -> None:
    """Sets the maximum number of entries kept in memory (LRU eviction beyond it)."""
    global _MAX_ENTRIES
    with _CACHE_LOCK:
        _MAX_ENTRIES = max(int(max_entries), 1)
        _evict_if_needed()


def enable_disk_persistence(path: str) -> None:
    """Enables a small JSON file cache so cached records survive an application restart."""
    global _DISK_CACHE_PATH
    with _CACHE_LOCK:
        _DISK_CACHE_PATH = path
        _load_disk_cache()


def disable_disk_persistence() -> None:
    """Disables the disk cache (used by tests and by privacy-conscious deployments)."""
    global _DISK_CACHE_PATH
    with _CACHE_LOCK:
        _DISK_CACHE_PATH = None


def _evict_if_needed() -> None:
    while len(_MEMORY_CACHE) > _MAX_ENTRIES:
        _MEMORY_CACHE.popitem(last=False)
        _STATS["evictions"] += 1


def _load_disk_cache() -> None:
    if not _DISK_CACHE_PATH or not os.path.isfile(_DISK_CACHE_PATH):
        return
    try:
        with open(_DISK_CACHE_PATH, "r", encoding="utf-8") as handle:
            payload = json.load(handle)
        stored_at = float(payload.get("stored_at", 0.0))
        if time.time() - stored_at >= CACHE_TTL_SECONDS:
            return
        for key, value in payload.get("items", {}).items():
            _MEMORY_CACHE[key] = (stored_at, value)
        _evict_if_needed()
    except (OSError, ValueError, TypeError):
        return


def persist_to_disk() -> bool:
    """Writes the current cache to the configured disk cache file (if enabled)."""
    if not _DISK_CACHE_PATH:
        return False
    try:
        payload = {
            "stored_at": time.time(),
            "items": {key: value for key, (_, value) in _MEMORY_CACHE.items()},
        }
        os.makedirs(os.path.dirname(os.path.abspath(_DISK_CACHE_PATH)), exist_ok=True)
        with open(_DISK_CACHE_PATH, "w", encoding="utf-8") as handle:
            json.dump(payload, handle, ensure_ascii=False)
        return True
    except (OSError, TypeError, ValueError):
        return False


def get_cached_item(cache_key: str, ttl_seconds: int = CACHE_TTL_SECONDS) -> tuple[bool, Any]:
    """Returns (hit, value) for a cache key, honouring the TTL and refreshing LRU order."""
    if not cache_key:
        return False, None
    with _CACHE_LOCK:
        if cache_key in _MEMORY_CACHE:
            ts, data = _MEMORY_CACHE[cache_key]
            if time.time() - ts < ttl_seconds:
                _MEMORY_CACHE.move_to_end(cache_key)
                _STATS["hits"] += 1
                return True, data
            del _MEMORY_CACHE[cache_key]
            _STATS["expirations"] += 1
        _STATS["misses"] += 1
    return False, None


def set_cached_item(cache_key: str, value: Any) -> None:
    """Stores a value in the cache, evicting the least recently used entry if needed."""
    if not cache_key:
        return
    with _CACHE_LOCK:
        _MEMORY_CACHE[cache_key] = (time.time(), value)
        _MEMORY_CACHE.move_to_end(cache_key)
        _STATS["sets"] += 1
        _evict_if_needed()


def clear_all_caches() -> None:
    """Clears the memory cache (and the statistics counters)."""
    with _CACHE_LOCK:
        _MEMORY_CACHE.clear()
        for key in _STATS:
            _STATS[key] = 0


def cache_statistics() -> dict[str, Any]:
    """Returns cache diagnostics for monitoring and for the platform's technical views."""
    with _CACHE_LOCK:
        return {
            "entries": len(_MEMORY_CACHE),
            "max_entries": _MAX_ENTRIES,
            "hits": _STATS["hits"],
            "misses": _STATS["misses"],
            "sets": _STATS["sets"],
            "evictions": _STATS["evictions"],
            "expirations": _STATS["expirations"],
            "hit_rate_pct": round(
                (_STATS["hits"] / (_STATS["hits"] + _STATS["misses"]) * 100.0)
                if (_STATS["hits"] + _STATS["misses"])
                else 0.0,
                2,
            ),
            "disk_persistence_enabled": bool(_DISK_CACHE_PATH),
        }

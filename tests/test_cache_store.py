"""Tests for the bounded, thread-safe external-data cache."""

import os
import tempfile
import threading
import time

from data_integrations.cache_store import (
    cache_statistics,
    clear_all_caches,
    configure_cache,
    disable_disk_persistence,
    enable_disk_persistence,
    get_cached_item,
    persist_to_disk,
    set_cached_item,
)


def test_set_and_get_cached_item():
    clear_all_caches()
    set_cached_item("pubchem:313", {"name": "Hydrochloric Acid"})
    hit, value = get_cached_item("pubchem:313")
    assert hit is True
    assert value["name"] == "Hydrochloric Acid"


def test_miss_returns_false():
    clear_all_caches()
    hit, value = get_cached_item("unknown-key")
    assert hit is False
    assert value is None


def test_expired_entry_is_not_returned():
    clear_all_caches()
    set_cached_item("expired:key", 123)
    hit, _ = get_cached_item("expired:key", ttl_seconds=0)
    assert hit is False


def test_lru_eviction_keeps_memory_bounded():
    clear_all_caches()
    configure_cache(max_entries=3)
    for idx in range(5):
        set_cached_item(f"key:{idx}", idx)
    stats = cache_statistics()
    assert stats["entries"] <= 3
    assert stats["evictions"] >= 2
    # the most recently used entry survives
    assert get_cached_item("key:4")[0] is True


def test_statistics_track_hits_and_misses():
    clear_all_caches()
    set_cached_item("k", "v")
    get_cached_item("k")
    get_cached_item("absent")
    stats = cache_statistics()
    assert stats["hits"] == 1
    assert stats["misses"] == 1
    assert stats["hit_rate_pct"] == 50.0


def test_concurrent_access_is_safe():
    clear_all_caches()
    errors: list[Exception] = []

    def worker(idx: int) -> None:
        try:
            for _ in range(50):
                set_cached_item(f"thread:{idx}", idx)
                get_cached_item(f"thread:{idx}")
        except Exception as exc:  # pragma: no cover - defensive
            errors.append(exc)

    threads = [threading.Thread(target=worker, args=(i,)) for i in range(4)]
    for thread in threads:
        thread.start()
    for thread in threads:
        thread.join()
    assert not errors


def test_disk_persistence_round_trip():
    clear_all_caches()
    try:
        with tempfile.TemporaryDirectory() as tmp:
            path = os.path.join(tmp, "cache", "api_cache.json")
            enable_disk_persistence(path)
            set_cached_item("persisted:key", {"x": 1})
            assert persist_to_disk() is True

            clear_all_caches()
            enable_disk_persistence(path)
            hit, value = get_cached_item("persisted:key")
            assert hit is True
            assert value["x"] == 1
    finally:
        disable_disk_persistence()

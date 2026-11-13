"""Tests for the search network-debouncing layer (request throttling without losing results)."""

import time

from search.search_engine import (
    NETWORK_DEBOUNCE_SECONDS,
    network_debounce_status,
    search_scientific_entity,
    set_network_debounce_seconds,
)


def _reset_debounce() -> None:
    set_network_debounce_seconds(NETWORK_DEBOUNCE_SECONDS)


def test_local_results_are_not_affected_by_debouncing():
    """A material present in the curated catalogue must always resolve, network or not."""
    _reset_debounce()
    res = search_scientific_entity("CH3COOH", allow_network=True)
    assert res["found"] is True
    assert res["selected"]["symbol_or_formula"] == "CH3COOH"
    assert res.get("network_deferred") is False


def test_debounce_blocks_a_second_external_lookup():
    set_network_debounce_seconds(30.0)
    try:
        first = search_scientific_entity("zzz-nonexistent-query-1", allow_network=True)
        second = search_scientific_entity("zzz-nonexistent-query-2", allow_network=True)
        assert first.get("network_deferred") is False
        assert second.get("network_deferred") is True
        assert "تأجيل" in second["message"] or "تحديث من المصادر" in second["message"]
    finally:
        _reset_debounce()


def test_forced_network_lookup_bypasses_the_debounce():
    set_network_debounce_seconds(30.0)
    try:
        search_scientific_entity("zzz-nonexistent-query-3", allow_network=True)
        forced = search_scientific_entity("zzz-nonexistent-query-4", allow_network=True, force_network=True)
        assert forced.get("network_deferred") is False
    finally:
        _reset_debounce()


def test_debounce_status_is_reported():
    _reset_debounce()
    status = network_debounce_status()
    assert status["interval_seconds"] == NETWORK_DEBOUNCE_SECONDS
    assert "slot_available" in status
    time.sleep(0.01)
    assert isinstance(status["seconds_since_last_request"], float)

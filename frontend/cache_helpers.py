"""
Streamlit caching wrappers and cache invalidation helpers for MolecuLens.

Uses:
- @st.cache_resource for heavy reusable objects (trained QSAR model bundle, local SQLite handle)
  so models are never retrained and databases are never reloaded on every UI interaction.
- @st.cache_data for deterministic reusable computations (molecular properties, similarity
  searches, titration curves, and scientific audit reports).
- Explicit cache version counter + invalidate_all_caches() for clean invalidation whenever
  underlying data or models change.

Degrades gracefully to direct execution when imported outside Streamlit (e.g., in unit tests).
"""

from __future__ import annotations

from typing import Any, Callable, TypeVar

F = TypeVar("F", bound=Callable[..., Any])

_CACHE_VERSION: int = 1

try:
    import streamlit as st
    _HAS_STREAMLIT = hasattr(st, "cache_resource") and hasattr(st, "cache_data")
except Exception:
    st = None  # type: ignore[assignment]
    _HAS_STREAMLIT = False


def _cache_resource_deco(**kwargs: Any) -> Callable[[F], F]:
    if _HAS_STREAMLIT and st is not None:
        return st.cache_resource(**kwargs)
    return lambda fn: fn


def _cache_data_deco(**kwargs: Any) -> Callable[[F], F]:
    if _HAS_STREAMLIT and st is not None:
        return st.cache_data(**kwargs)
    return lambda fn: fn


def get_cache_version() -> int:
    """Returns the current cache schema/data version token."""
    return _CACHE_VERSION


def invalidate_all_caches() -> int:
    """
    Invalidates all Streamlit resource/data caches and increments the cache version token.
    Call this whenever underlying datasets, registries, or QSAR training data change.
    """
    global _CACHE_VERSION
    _CACHE_VERSION += 1
    if _HAS_STREAMLIT and st is not None:
        try:
            st.cache_data.clear()
        except Exception:
            pass
        try:
            st.cache_resource.clear()
        except Exception:
            pass
    return _CACHE_VERSION


@_cache_resource_deco(show_spinner=False)
def get_cached_qsar_bundle(cache_version: int = 1) -> dict[str, Any]:
    """Loads and caches the trained Ridge & Hybrid Residual QSAR model bundle once."""
    _ = cache_version
    from chemistry_engine.molecular_intelligence import MolecularIntelligence

    return MolecularIntelligence.train_and_evaluate_qsar_logp_model()


@_cache_resource_deco(show_spinner=False)
def get_cached_local_db(cache_version: int = 1) -> Any:
    """Returns the cached local SQLite database handle without reloading on every rerun."""
    _ = cache_version
    from data_integrations.local_db import get_local_database

    return get_local_database()


@_cache_data_deco(show_spinner=False, ttl=1800)
def get_cached_similarity_results(smiles: str, top_n: int = 6, cache_version: int = 1) -> list[dict[str, Any]]:
    """Caches Morgan Fingerprint / Tanimoto similarity results for a given SMILES."""
    _ = cache_version
    from chemistry_engine.molecular_intelligence import MolecularIntelligence

    return MolecularIntelligence.find_most_similar(smiles, top_n=top_n)


@_cache_data_deco(show_spinner=False, ttl=1800)
def get_cached_molecular_properties(smiles: str, cache_version: int = 1) -> dict[str, Any] | None:
    """Caches deterministic physicochemical descriptors calculated from SMILES."""
    _ = cache_version
    from chemistry_engine.structure_analyzer import calculate_rdkit_descriptors

    return calculate_rdkit_descriptors(smiles)


@_cache_data_deco(show_spinner=False, ttl=1800)
def get_cached_titration_result(
    acid_key: str,
    base_key: str,
    acid_conc_m: float,
    acid_volume_ml: float,
    base_conc_m: float,
    titration_mode: str = "acid_with_base",
    base_volume_ml: float | None = None,
    cache_version: int = 1,
) -> dict[str, Any]:
    """Caches acid-base titration curve calculations for identical input parameters."""
    _ = cache_version
    from titration.titration_engine import TitrationEngine

    kwargs: dict[str, Any] = {
        "acid_key": acid_key,
        "base_key": base_key,
        "acid_conc_m": acid_conc_m,
        "acid_volume_ml": acid_volume_ml,
        "base_conc_m": base_conc_m,
        "titration_mode": titration_mode,
    }
    if base_volume_ml is not None:
        kwargs["base_volume_ml"] = base_volume_ml
    return TitrationEngine.simulate_titration(**kwargs)


@_cache_data_deco(show_spinner=False, ttl=1800)
def get_cached_logp_evaluation(smiles: str, cache_version: int = 1) -> dict[str, Any]:
    """Caches LogP evaluation (deterministic Crippen + cached QSAR model prediction)."""
    # Warm the resource-cached QSAR model bundle so it is never retrained on rerun.
    get_cached_qsar_bundle(cache_version=cache_version)
    from chemistry_engine.molecular_intelligence import MolecularIntelligence

    return MolecularIntelligence.evaluate_logp_prediction(smiles)


@_cache_data_deco(show_spinner=False, ttl=3600)
def get_cached_reaction_audit(cache_version: int = 1) -> dict[str, Any]:
    """Caches the full scientific audit summary of the verified reaction registry."""
    _ = cache_version
    from reaction_engine.scientific_audit import audit_registry_and_rules

    return audit_registry_and_rules()

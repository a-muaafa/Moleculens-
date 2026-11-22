"""
Auditable Scientific Repository Layer for MolecuLens (v3.3.0).

Centralizes data access (local catalog + SQLite store) behind a clean Repository
boundary and attaches structured provenance/audit fields to every scientific record,
strictly separating developer-curated records from external API records.
"""

from __future__ import annotations

from datetime import datetime, timezone
from typing import Any

from data_integrations.local_db import get_local_database
from search.local_catalog import LOCAL_SCIENTIFIC_CATALOG


def attach_record_provenance(
    record: dict[str, Any],
    *,
    record_origin: str = "curated_local",
    source_version: str = "moleculens-catalog-3.3.0",
) -> dict[str, Any]:
    """
    Returns a shallow copy of `record` enriched with auditable provenance fields:
    source_id, source_url, source_version, retrieved_at, verified_at,
    verification_status, license, data_quality_status, and record_origin.
    Never mixes curated local data with external API provenance.
    """
    enriched = dict(record)
    is_external = record_origin == "external_api" or "external" in str(record.get("db_source", "")).lower()
    origin = "external_api" if is_external else "curated_local"

    cid = record.get("pubchem_cid")
    chembl_id = record.get("chembl_id")
    if origin == "external_api":
        default_src_id = f"CID:{cid}" if cid else (str(chembl_id) if chembl_id else str(record.get("key", "external")))
        default_url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}" if cid else ""
        default_license = "External public database terms (PubChem/ChEMBL/ChEBI)"
        default_verif = "external_retrieved_unreviewed"
        default_quality = "external_auto_sanitized"
    else:
        default_src_id = str(record.get("key") or record.get("symbol_or_formula") or "local_record")
        default_url = f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}" if cid else ""
        default_license = "MIT / Curated Educational Reference (IUPAC, CRC, NIST)"
        default_verif = "curated_local_verified" if record.get("canonical_smiles") or record.get("symbol_or_formula") else "curated_partial"
        default_quality = "verified_local"

    now_iso = datetime.now(timezone.utc).isoformat()
    enriched["record_origin"] = origin
    enriched["source_id"] = record.get("source_id") or default_src_id
    enriched["source_url"] = record.get("source_url") or default_url
    enriched["source_version"] = record.get("source_version") or source_version
    enriched["retrieved_at"] = record.get("retrieved_at") or now_iso
    enriched["verified_at"] = record.get("verified_at") or (now_iso if origin == "curated_local" else None)
    enriched["verification_status"] = record.get("verification_status") or default_verif
    enriched["license"] = record.get("license") or default_license
    enriched["data_quality_status"] = record.get("data_quality_status") or default_quality
    return enriched


class ScientificRepository:
    """Unified Repository facade over the local curated catalog and SQLite store."""

    @staticmethod
    def get_all_curated_substances() -> list[dict[str, Any]]:
        return [attach_record_provenance(rec, record_origin="curated_local") for rec in LOCAL_SCIENTIFIC_CATALOG]

    @staticmethod
    def find_by_term(term: str) -> list[dict[str, Any]]:
        db = get_local_database()
        rows = db.lookup_by_term(term)
        return [
            attach_record_provenance(
                r,
                record_origin="external_api" if "external" in str(r.get("db_source", "")).lower() else "curated_local",
            )
            for r in rows
        ]

    @staticmethod
    def find_by_formula(formula: str) -> list[dict[str, Any]]:
        db = get_local_database()
        rows = db.lookup_by_formula(formula)
        return [
            attach_record_provenance(
                r,
                record_origin="external_api" if "external" in str(r.get("db_source", "")).lower() else "curated_local",
            )
            for r in rows
        ]

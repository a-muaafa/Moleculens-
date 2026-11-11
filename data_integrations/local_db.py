import json
import os
import sqlite3
import threading
from contextlib import contextmanager
from datetime import datetime, timezone
from typing import Any

from chemistry_engine.formula_analyzer import normalize_chemical_formula, replace_subscripts
from config.settings import LOCAL_SQLITE_DB_PATH
from search.normalizer import normalize_search_input

_SCHEMA = """
CREATE TABLE IF NOT EXISTS substances (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    key TEXT NOT NULL UNIQUE,
    entity_type TEXT NOT NULL,
    subclass TEXT DEFAULT '',
    name_ar TEXT DEFAULT '',
    name_en TEXT DEFAULT '',
    iupac_name TEXT DEFAULT '',
    formula TEXT DEFAULT '',
    hill_formula TEXT DEFAULT '',
    smiles TEXT DEFAULT '',
    molar_mass REAL,
    pubchem_cid INTEGER,
    chebi_id TEXT,
    chembl_id TEXT,
    uniprot_id TEXT,
    pdb_ids TEXT DEFAULT '',
    source TEXT DEFAULT '',
    payload_json TEXT NOT NULL,
    fetched_at TEXT NOT NULL,
    updated_at TEXT NOT NULL
);
CREATE TABLE IF NOT EXISTS synonyms (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    substance_id INTEGER NOT NULL REFERENCES substances(id) ON DELETE CASCADE,
    term TEXT NOT NULL,
    term_norm TEXT NOT NULL,
    lang TEXT DEFAULT '',
    UNIQUE(substance_id, term_norm)
);
CREATE TABLE IF NOT EXISTS reactions (
    id INTEGER PRIMARY KEY AUTOINCREMENT,
    reaction_id TEXT NOT NULL UNIQUE,
    domain TEXT DEFAULT '',
    category TEXT DEFAULT '',
    equation_text TEXT DEFAULT '',
    reactant_formulas TEXT DEFAULT '',
    product_formulas TEXT DEFAULT '',
    source TEXT DEFAULT '',
    payload_json TEXT NOT NULL,
    fetched_at TEXT NOT NULL,
    updated_at TEXT NOT NULL
);
CREATE INDEX IF NOT EXISTS idx_syn_norm ON synonyms(term_norm);
CREATE INDEX IF NOT EXISTS idx_sub_formula ON substances(formula);
CREATE INDEX IF NOT EXISTS idx_sub_hill ON substances(hill_formula);
CREATE INDEX IF NOT EXISTS idx_sub_cid ON substances(pubchem_cid);
CREATE INDEX IF NOT EXISTS idx_sub_chebi ON substances(chebi_id);
CREATE INDEX IF NOT EXISTS idx_sub_chembl ON substances(chembl_id);
CREATE INDEX IF NOT EXISTS idx_sub_smiles ON substances(smiles);
CREATE INDEX IF NOT EXISTS idx_rxn_reactants ON reactions(reactant_formulas);
"""


def _now() -> str:
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


class SmartLocalDatabase:
    def __init__(self, db_path: str | None = None):
        self.db_path = db_path or LOCAL_SQLITE_DB_PATH
        parent = os.path.dirname(self.db_path)
        if parent:
            os.makedirs(parent, exist_ok=True)
        self._lock = threading.RLock()
        with self._session() as conn:
            conn.executescript(_SCHEMA)

    @contextmanager
    def _session(self):
        conn = sqlite3.connect(self.db_path)
        conn.row_factory = sqlite3.Row
        conn.execute("PRAGMA foreign_keys = ON")
        try:
            with conn:
                yield conn
        finally:
            conn.close()

    def _find_existing_id(self, conn: sqlite3.Connection, entry: dict[str, Any]) -> int | None:
        row = conn.execute("SELECT id FROM substances WHERE key = ?", (entry.get("key"),)).fetchone()
        if row:
            return row["id"]
        cid = entry.get("pubchem_cid")
        if cid:
            row = conn.execute("SELECT id FROM substances WHERE pubchem_cid = ?", (cid,)).fetchone()
            if row:
                return row["id"]
        chembl = entry.get("chembl_id")
        if chembl:
            row = conn.execute("SELECT id FROM substances WHERE chembl_id = ?", (chembl,)).fetchone()
            if row:
                return row["id"]
        chebi = entry.get("chebi_id")
        if chebi:
            row = conn.execute("SELECT id FROM substances WHERE chebi_id = ?", (chebi,)).fetchone()
            if row:
                return row["id"]
        return None

    def upsert_substance(self, entry: dict[str, Any], source: str = "Curated Scientific Catalog") -> int:
        formula = replace_subscripts(entry.get("symbol_or_formula", "") or "")
        hill = normalize_chemical_formula(formula) or ""
        now = _now()
        with self._lock, self._session() as conn:
            existing_id = self._find_existing_id(conn, entry)
            if existing_id is not None and not source.startswith("Curated"):
                src_row = conn.execute("SELECT source FROM substances WHERE id = ?", (existing_id,)).fetchone()
                if src_row and src_row["source"].startswith("Curated"):
                    return existing_id
            values = (
                entry.get("entity_type", "chemical_compound"),
                entry.get("subclass", ""),
                entry.get("name_ar", ""),
                entry.get("name_en", ""),
                entry.get("iupac_name", ""),
                formula,
                hill,
                entry.get("canonical_smiles", ""),
                entry.get("molar_mass"),
                entry.get("pubchem_cid"),
                entry.get("chebi_id"),
                entry.get("chembl_id"),
                entry.get("uniprot_id"),
                ",".join(entry.get("pdb_ids", []) or []),
                source,
                json.dumps(entry, ensure_ascii=False, default=str),
            )
            if existing_id is None:
                cur = conn.execute(
                    "INSERT INTO substances (entity_type, subclass, name_ar, name_en, iupac_name, formula, hill_formula, smiles, "
                    "molar_mass, pubchem_cid, chebi_id, chembl_id, uniprot_id, pdb_ids, source, payload_json, key, fetched_at, updated_at) "
                    "VALUES (?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?,?)",
                    (*values, entry.get("key"), now, now),
                )
                sub_id = cur.lastrowid
            else:
                conn.execute(
                    "UPDATE substances SET entity_type=?, subclass=?, name_ar=?, name_en=?, iupac_name=?, formula=?, hill_formula=?, "
                    "smiles=?, molar_mass=?, pubchem_cid=COALESCE(?, pubchem_cid), chebi_id=COALESCE(?, chebi_id), "
                    "chembl_id=COALESCE(?, chembl_id), uniprot_id=COALESCE(?, uniprot_id), pdb_ids=?, source=?, payload_json=?, "
                    "updated_at=? WHERE id=?",
                    (*values, now, existing_id),
                )
                sub_id = existing_id

            terms: list[tuple[str, str]] = []
            for t in entry.get("synonyms_ar", []) + [entry.get("name_ar", "")]:
                terms.append((t, "ar"))
            for t in (
                entry.get("synonyms_en", [])
                + [entry.get("name_en", ""), entry.get("iupac_name", ""), formula]
                + entry.get("abbreviations", [])
                + entry.get("brand_names", [])
            ):
                terms.append((t, "en"))
            for term, lang in terms:
                if not term:
                    continue
                norm = normalize_search_input(term)
                if norm:
                    conn.execute(
                        "INSERT OR IGNORE INTO synonyms (substance_id, term, term_norm, lang) VALUES (?,?,?,?)",
                        (sub_id, term, norm, lang),
                    )
            return sub_id

    def upsert_reaction(self, reaction: dict[str, Any], source: str = "Curated Reaction Registry") -> None:
        now = _now()
        rxn = dict(reaction)
        lhs = " + ".join(f"{r['coefficient'] if r.get('coefficient', 1) > 1 else ''}{r['formula']}" for r in rxn.get("reactants", []))
        rhs = " + ".join(f"{p['coefficient'] if p.get('coefficient', 1) > 1 else ''}{p['formula']}" for p in rxn.get("products", []))
        eq_text = rxn.get("equation_text") or f"{lhs} -> {rhs}"
        rxn["equation_text"] = eq_text
        reactants = ",".join(sorted(r["formula"] for r in rxn.get("reactants", [])))
        products = ",".join(sorted(p["formula"] for p in rxn.get("products", [])))
        with self._lock, self._session() as conn:
            row = conn.execute("SELECT id FROM reactions WHERE reaction_id = ?", (rxn["reaction_id"],)).fetchone()
            payload = json.dumps(rxn, ensure_ascii=False, default=str)
            if row:
                conn.execute(
                    "UPDATE reactions SET domain=?, category=?, equation_text=?, reactant_formulas=?, product_formulas=?, "
                    "source=?, payload_json=?, updated_at=? WHERE id=?",
                    (rxn.get("domain", ""), rxn.get("category", ""), eq_text, reactants, products, source, payload, now, row["id"]),
                )
            else:
                conn.execute(
                    "INSERT INTO reactions (reaction_id, domain, category, equation_text, reactant_formulas, product_formulas, "
                    "source, payload_json, fetched_at, updated_at) VALUES (?,?,?,?,?,?,?,?,?,?)",
                    (rxn["reaction_id"], rxn.get("domain", ""), rxn.get("category", ""), eq_text, reactants, products, source, payload, now, now),
                )

    def seed(self, catalog: list[dict[str, Any]], reactions: list[dict[str, Any]]) -> dict[str, int]:
        for entry in catalog:
            self.upsert_substance(entry, source="Curated Scientific Catalog")
        for rxn in reactions:
            self.upsert_reaction(rxn)
        return self.get_database_statistics()

    def _row_to_entry(self, row: sqlite3.Row) -> dict[str, Any]:
        entry = json.loads(row["payload_json"])
        entry["db_source"] = row["source"]
        entry["fetched_at"] = row["fetched_at"]
        entry["updated_at"] = row["updated_at"]
        return entry

    def lookup_by_term(self, query: str) -> list[dict[str, Any]]:
        norm = normalize_search_input(query)
        if not norm:
            return []
        with self._session() as conn:
            rows = conn.execute(
                "SELECT DISTINCT s.* FROM substances s JOIN synonyms y ON y.substance_id = s.id WHERE y.term_norm = ?",
                (norm,),
            ).fetchall()
        return [self._row_to_entry(r) for r in rows]

    def lookup_by_formula(self, formula: str) -> list[dict[str, Any]]:
        hill = normalize_chemical_formula(formula) or ""
        raw = replace_subscripts(formula)
        with self._session() as conn:
            rows = conn.execute(
                "SELECT * FROM substances WHERE formula = ? OR (hill_formula = ? AND hill_formula != '')",
                (raw, hill),
            ).fetchall()
        return [self._row_to_entry(r) for r in rows]

    def search_prefix(self, query: str, limit: int = 10) -> list[dict[str, Any]]:
        norm = normalize_search_input(query)
        if not norm:
            return []
        with self._session() as conn:
            rows = conn.execute(
                "SELECT DISTINCT s.* FROM substances s JOIN synonyms y ON y.substance_id = s.id "
                "WHERE y.term_norm LIKE ? ORDER BY length(y.term_norm) LIMIT ?",
                (norm + "%", limit),
            ).fetchall()
        return [self._row_to_entry(r) for r in rows]

    def list_external_substances(self) -> list[dict[str, Any]]:
        with self._session() as conn:
            rows = conn.execute("SELECT * FROM substances WHERE source NOT LIKE 'Curated%'").fetchall()
        return [self._row_to_entry(r) for r in rows]

    def search_reactions_by_formula(self, formula: str) -> list[dict[str, Any]]:
        target = normalize_chemical_formula(formula) or replace_subscripts(formula)
        results = []
        with self._session() as conn:
            rows = conn.execute("SELECT * FROM reactions").fetchall()
        for r in rows:
            formulas = [f for f in (r["reactant_formulas"] + "," + r["product_formulas"]).split(",") if f]
            if any((normalize_chemical_formula(f) or f) == target for f in formulas):
                payload = json.loads(r["payload_json"])
                if "equation_text" not in payload:
                    payload["equation_text"] = r["equation_text"]
                results.append(payload)
        return results

    def get_database_statistics(self) -> dict[str, int]:
        with self._session() as conn:
            subs = conn.execute("SELECT COUNT(*) AS n FROM substances").fetchone()["n"]
            syns = conn.execute("SELECT COUNT(*) AS n FROM synonyms").fetchone()["n"]
            rxns = conn.execute("SELECT COUNT(*) AS n FROM reactions").fetchone()["n"]
            by_type = conn.execute("SELECT entity_type, COUNT(*) AS n FROM substances GROUP BY entity_type").fetchall()
            by_sub = conn.execute("SELECT subclass, COUNT(*) AS n FROM substances GROUP BY subclass").fetchall()
            ext = conn.execute("SELECT COUNT(*) AS n FROM substances WHERE source NOT LIKE 'Curated%'").fetchone()["n"]
        return {
            "substances": subs,
            "synonyms": syns,
            "reactions": rxns,
            "external_cached_substances": ext,
            **{f"type:{r['entity_type']}": r["n"] for r in by_type},
            **{f"subclass:{r['subclass']}": r["n"] for r in by_sub},
        }


_INSTANCE: SmartLocalDatabase | None = None
_INSTANCE_LOCK = threading.Lock()


def get_local_database() -> SmartLocalDatabase:
    global _INSTANCE
    with _INSTANCE_LOCK:
        if _INSTANCE is None:
            from reaction_engine.reaction_registry import VERIFIED_REACTIONS_REGISTRY
            from search.local_catalog import LOCAL_SCIENTIFIC_CATALOG

            db = SmartLocalDatabase()
            stats = db.get_database_statistics()
            catalog_expected = len(LOCAL_SCIENTIFIC_CATALOG)
            curated = stats["substances"] - stats["external_cached_substances"]
            if curated < catalog_expected or stats["reactions"] < len(VERIFIED_REACTIONS_REGISTRY):
                db.seed(LOCAL_SCIENTIFIC_CATALOG, VERIFIED_REACTIONS_REGISTRY)
            _INSTANCE = db
        return _INSTANCE

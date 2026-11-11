from typing import Any
from chemistry_engine.formula_analyzer import formulas_are_equivalent
from chemistry_engine.structure_analyzer import (
    HAS_RDKIT,
    Chem,
    compute_inchikey,
    validate_and_parse_smiles,
    verify_substance_structure_identity,
)
from config.settings import (
    ENTITY_BIOMOLECULE_PROTEIN,
    ENTITY_CHEMICAL_COMPOUND,
    ENTITY_CHEMICAL_ELEMENT,
    ENTITY_DRUG_SUBSTANCE,
)
from core.schemas import build_unified_entity_schema
from data_integrations.chebi_client import search_chebi_entities
from data_integrations.chembl_client import fetch_chembl_drug_info, search_chembl_molecules
from data_integrations.local_db import get_local_database
from data_integrations.ord_client import search_ord_reactions, search_rhea_reactions
from data_integrations.pdb_client import fetch_pdb_structure_metadata
from data_integrations.pubchem_client import fetch_pubchem_by_formula, fetch_pubchem_by_name
from data_integrations.rhea_client import fetch_rhea_reaction_info
from data_integrations.uniprot_client import fetch_uniprot_protein_info


class DataManager:
    SOURCE_ROLES = {
        "PubChem": "الأسماء والصيغ والبنية الجزيئية والخصائص الفيزيائية المتاحة",
        "ChEMBL": "النشاط الحيوي والبيانات الدوائية الموثقة",
        "UniProt": "البروتينات والتسلسلات الأمينية والوظائف الحيوية",
        "PDB": "التراكيب الفراغية ثلاثية الأبعاد للجزيئات الحيوية",
        "Rhea": "التفاعلات الكيميائية الحيوية والإنزيمية",
        "ChEBI": "التصنيف الكيميائي الحيوي والمرادفات",
        "Open Reaction Database": "مرجع التفاعلات العضوية الموثقة",
    }

    @staticmethod
    def unify_entity_record(raw_record: dict[str, Any], allow_network_enrichment: bool = False) -> dict[str, Any]:
        unified = build_unified_entity_schema(raw_record)
        entity_type = unified["entity_type"]
        ids = unified["identifiers"]
        sources = list(unified["sources"])

        if entity_type == ENTITY_CHEMICAL_ELEMENT:
            if "IUPAC Periodic Table" not in sources:
                sources.append("IUPAC Periodic Table")
            if ids.get("pubchem_cid") and "PubChem" not in sources:
                sources.append("PubChem")

        elif entity_type == ENTITY_CHEMICAL_COMPOUND:
            if ids.get("pubchem_cid") and "PubChem" not in sources:
                sources.append("PubChem")

        elif entity_type == ENTITY_DRUG_SUBSTANCE:
            if ids.get("pubchem_cid") and "PubChem" not in sources:
                sources.append("PubChem")
            if ids.get("chembl_id") and "ChEMBL" not in sources:
                sources.append("ChEMBL")
            if allow_network_enrichment and ids.get("chembl_id"):
                ch_res = fetch_chembl_drug_info(ids["chembl_id"])
                if ch_res.get("ok"):
                    unified["pharmacology_details"]["live_chembl_metadata"] = ch_res["data"]

        elif entity_type == ENTITY_BIOMOLECULE_PROTEIN:
            if ids.get("uniprot_id") and "UniProt" not in sources:
                sources.append("UniProt")
            if ids.get("pdb_ids") and "PDB" not in sources:
                sources.append("PDB")
            if ids.get("rhea_ids") and "Rhea" not in sources:
                sources.append("Rhea")
            if allow_network_enrichment:
                if ids.get("uniprot_id"):
                    up_res = fetch_uniprot_protein_info(ids["uniprot_id"])
                    if up_res.get("ok"):
                        unified["biology_details"]["live_uniprot_metadata"] = up_res["data"]
                if ids.get("pdb_ids"):
                    pdb_res = fetch_pdb_structure_metadata(ids["pdb_ids"][0])
                    if pdb_res.get("ok"):
                        unified["biology_details"]["live_pdb_metadata"] = pdb_res["data"]
                if ids.get("rhea_ids"):
                    rhea_res = fetch_rhea_reaction_info(ids["rhea_ids"][0])
                    if rhea_res.get("ok"):
                        unified["biology_details"]["live_rhea_metadata"] = rhea_res["data"]

        unified["sources"] = list(dict.fromkeys(sources))
        unified["identity_verification"] = verify_substance_structure_identity(unified)
        return unified

    @staticmethod
    def search_external_pubchem(query: str, is_formula: bool = False) -> dict[str, Any]:
        if is_formula:
            return fetch_pubchem_by_formula(query)
        return fetch_pubchem_by_name(query)

    @staticmethod
    def infer_subclass(smiles: str, formula: str = "", is_drug: bool = False) -> str:
        if is_drug:
            return "drug"
        if HAS_RDKIT and Chem is not None:
            mol = Chem.MolFromSmiles(smiles) if smiles else None
            if mol is None:
                return "organic_compound"
            symbols = {a.GetSymbol() for a in mol.GetAtoms()}
            if "C" not in symbols:
                charged = [a for a in mol.GetAtoms() if a.GetFormalCharge() != 0]
                if len(mol.GetAtoms()) == 1 and charged:
                    return "ion"
                return "inorganic_compound"
            if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3](=O)[OX2H1]")):
                return "carboxylic_acid"
            if mol.HasSubstructMatch(Chem.MolFromSmarts("[#6][CX3](=O)[OX2][#6]")):
                return "ester"
            if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX3H1](=O)[#6,#1]")):
                return "aldehyde"
            if mol.HasSubstructMatch(Chem.MolFromSmarts("[#6][CX3](=O)[#6]")):
                return "ketone"
            if mol.GetRingInfo().NumRings() and any(a.GetIsAromatic() for a in mol.GetAtoms()):
                return "aromatic_compound"
            if mol.HasSubstructMatch(Chem.MolFromSmarts("[CX4][OX2H]")):
                return "alcohol"
            return "organic_compound"
        return "organic_compound"

    @staticmethod
    def match_records_by_structure_and_identifiers(
        primary_rec: dict[str, Any],
        candidate_records: list[dict[str, Any]],
    ) -> dict[str, Any] | None:
        """
        Matches a primary record (e.g., from PubChem) against candidate records (from ChEMBL or ChEBI)
        using scientific identifiers and chemical structure (InChIKey, Canonical SMILES, Formula + Synonyms),
        rather than relying on string name equality alone.
        """
        if not candidate_records:
            return None

        p_ik = (primary_rec.get("inchikey") or compute_inchikey(primary_rec.get("canonical_smiles", "")) or "").strip()
        _, p_smi, _ = validate_and_parse_smiles(primary_rec.get("canonical_smiles", ""))
        p_form = (primary_rec.get("symbol_or_formula") or "").strip()
        p_names = {
            s.strip().lower()
            for s in [primary_rec.get("name_en", ""), primary_rec.get("iupac_name", "")] + list(primary_rec.get("synonyms_en", []))
            if s and s.strip()
        }

        # Priority 1: Exact InChIKey match (or first 14-char connectivity block match)
        if p_ik:
            for cand in candidate_records:
                c_ik = (cand.get("inchikey") or compute_inchikey(cand.get("canonical_smiles", "")) or "").strip()
                if c_ik and (c_ik == p_ik or c_ik[:14] == p_ik[:14]):
                    return cand

        # Priority 2: Exact Canonical SMILES match
        if p_smi:
            for cand in candidate_records:
                _, c_smi, _ = validate_and_parse_smiles(cand.get("canonical_smiles", ""))
                if c_smi and c_smi == p_smi:
                    return cand

        # Priority 3: Equivalent Molecular Formula + matching name/synonym
        for cand in candidate_records:
            c_form = (cand.get("symbol_or_formula") or "").strip()
            c_names = {
                s.strip().lower()
                for s in [cand.get("name_en", "")] + list(cand.get("synonyms_en", []))
                if s and s.strip()
            }
            if p_form and c_form and formulas_are_equivalent(p_form, c_form) and (p_names & c_names):
                return cand

        # Priority 4: Synonym / Name match ONLY if candidate has no conflicting formula or SMILES
        for cand in candidate_records:
            c_form = (cand.get("symbol_or_formula") or "").strip()
            _, c_smi, _ = validate_and_parse_smiles(cand.get("canonical_smiles", ""))
            if p_smi and c_smi and p_smi != c_smi:
                continue
            if p_form and c_form and not formulas_are_equivalent(p_form, c_form):
                continue
            c_names = {
                s.strip().lower()
                for s in [cand.get("name_en", "")] + list(cand.get("synonyms_en", []))
                if s and s.strip()
            }
            if p_names & c_names:
                return cand

        return None

    @staticmethod
    def resolve_external_substances(query: str, is_formula: bool = False, allow_network: bool = True) -> dict[str, Any]:
        db = get_local_database()
        cached = db.lookup_by_formula(query) if is_formula else db.lookup_by_term(query)
        cached = [c for c in cached if c.get("db_source", "").startswith("External")]
        if cached:
            return {"ok": True, "records": cached, "origin": "local_database", "errors": []}
        if not allow_network:
            return {"ok": False, "records": [], "origin": "none", "errors": ["network disabled"]}

        errors: list[str] = []
        pc = DataManager.search_external_pubchem(query, is_formula=is_formula)
        if not pc.get("ok"):
            errors.append(f"PubChem: {pc.get('error')}")
        chembl = {"ok": False, "records": []} if is_formula else search_chembl_molecules(query)
        if not is_formula and not chembl.get("ok"):
            errors.append(f"ChEMBL: {chembl.get('error')}")
        chebi = {"ok": False, "records": []} if is_formula else search_chebi_entities(query)
        if not is_formula and not chebi.get("ok"):
            errors.append(f"ChEBI: {chebi.get('error')}")

        chembl_records = chembl.get("records", [])
        chebi_records = chebi.get("records", [])

        entries: list[dict[str, Any]] = []
        for rec in pc.get("records", []):
            ch = DataManager.match_records_by_structure_and_identifiers(rec, chembl_records)
            cb = DataManager.match_records_by_structure_and_identifiers(rec, chebi_records)
            sources = ["PubChem"]
            entry: dict[str, Any] = {
                "key": f"pubchem_{rec.get('pubchem_cid')}",
                "name_ar": rec.get("name_en", query),
                "name_en": rec.get("name_en", query),
                "iupac_name": rec.get("iupac_name", ""),
                "symbol_or_formula": rec.get("symbol_or_formula", ""),
                "canonical_smiles": rec.get("canonical_smiles", ""),
                "inchikey": rec.get("inchikey", ""),
                "molar_mass": rec.get("molar_mass"),
                "pubchem_cid": rec.get("pubchem_cid"),
                "synonyms_ar": [],
                "synonyms_en": [rec.get("name_en", "")],
            }
            is_drug = False
            if ch:
                entry["chembl_id"] = ch["chembl_id"]
                if not entry["inchikey"] and ch.get("inchikey"):
                    entry["inchikey"] = ch["inchikey"]
                sources.append("ChEMBL")
                is_drug = bool(ch.get("max_phase") and float(ch["max_phase"]) >= 1)
            if cb:
                entry["chebi_id"] = cb["chebi_id"]
                entry["synonyms_en"] = list(dict.fromkeys(entry["synonyms_en"] + cb.get("synonyms_en", [])))
                sources.append("ChEBI")
            entry["subclass"] = DataManager.infer_subclass(entry["canonical_smiles"], entry["symbol_or_formula"], is_drug)
            entry["entity_type"] = ENTITY_DRUG_SUBSTANCE if is_drug else ENTITY_CHEMICAL_COMPOUND
            entry["sources"] = sources
            entry["identity_verification"] = verify_substance_structure_identity(entry, expected_query=query)
            db.upsert_substance(entry, source="External: " + ", ".join(sources))
            entries.append(entry)
        return {"ok": bool(entries), "records": entries, "origin": "external_sources", "errors": errors}

    @staticmethod
    def search_documented_reactions_external(query: str, domain: str = "all_domains", allow_network: bool = True) -> dict[str, Any]:
        """Searches documented reactions across both ORD (for organic chemistry) and Rhea (for biochemistry)."""
        records: list[dict[str, Any]] = []
        errors: list[str] = []

        if domain in {"all_domains", "all", "organic_chemistry"}:
            ord_res = search_ord_reactions(query, allow_network=allow_network)
            if ord_res.get("ok"):
                records.extend(ord_res.get("records", []))
            elif ord_res.get("error"):
                errors.append(f"ORD: {ord_res['error']}")

        if domain in {"all_domains", "all", "biochemistry"} and allow_network:
            rhea_res = search_rhea_reactions(query)
            if rhea_res.get("ok"):
                records.extend(rhea_res.get("records", []))
            elif rhea_res.get("error"):
                errors.append(f"Rhea: {rhea_res['error']}")

        return {
            "ok": bool(records),
            "records": records,
            "errors": errors,
        }

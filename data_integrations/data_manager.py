from typing import Any
from rdkit import Chem
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
from data_integrations.ord_client import search_rhea_reactions
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

        chembl_by_name = {r["name_en"].lower(): r for r in chembl.get("records", [])}
        chebi_by_name = {r["name_en"].lower(): r for r in chebi.get("records", [])}

        entries: list[dict[str, Any]] = []
        for rec in pc.get("records", []):
            name_l = rec.get("name_en", "").lower()
            ch = chembl_by_name.get(name_l)
            cb = chebi_by_name.get(name_l)
            sources = ["PubChem"]
            entry: dict[str, Any] = {
                "key": f"pubchem_{rec.get('pubchem_cid')}",
                "name_ar": rec.get("name_en", query),
                "name_en": rec.get("name_en", query),
                "iupac_name": rec.get("iupac_name", ""),
                "symbol_or_formula": rec.get("symbol_or_formula", ""),
                "canonical_smiles": rec.get("canonical_smiles", ""),
                "molar_mass": rec.get("molar_mass"),
                "pubchem_cid": rec.get("pubchem_cid"),
                "synonyms_ar": [],
                "synonyms_en": [rec.get("name_en", "")],
            }
            is_drug = False
            if ch:
                entry["chembl_id"] = ch["chembl_id"]
                sources.append("ChEMBL")
                is_drug = bool(ch.get("max_phase") and float(ch["max_phase"]) >= 1)
            if cb:
                entry["chebi_id"] = cb["chebi_id"]
                entry["synonyms_en"] += cb.get("synonyms_en", [])
                sources.append("ChEBI")
            entry["subclass"] = DataManager.infer_subclass(entry["canonical_smiles"], entry["symbol_or_formula"], is_drug)
            entry["entity_type"] = ENTITY_DRUG_SUBSTANCE if is_drug else ENTITY_CHEMICAL_COMPOUND
            entry["sources"] = sources
            db.upsert_substance(entry, source="External: " + ", ".join(sources))
            entries.append(entry)
        return {"ok": bool(entries), "records": entries, "origin": "external_sources", "errors": errors}

    @staticmethod
    def search_documented_reactions_external(query: str) -> dict[str, Any]:
        return search_rhea_reactions(query)


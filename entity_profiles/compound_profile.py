from typing import Any
from chemistry_engine.formula_analyzer import analyze_formula_composition
from chemistry_engine.structure_analyzer import render_molecule_2d_png
from config.settings import (
    ENTITY_CHEMICAL_COMPOUND,
    OP_SCIENTIFIC_CALCULATION,
    OP_SCIENTIFIC_RETRIEVAL,
)
from core.schemas import build_provenance_tag
from education.explanation_engine import build_entity_explanation


def build_compound_scientific_profile(entity: dict[str, Any]) -> dict[str, Any]:
    formula = entity.get("symbol_or_formula", "")
    smiles = entity.get("canonical_smiles", "")
    composition = analyze_formula_composition(formula)
    molar_mass = entity.get("molar_mass") or composition.get("molar_mass")
    struct = entity.get("structure_details", {})
    phys = entity.get("physical_properties", {})
    chem = entity.get("chemical_properties", {})
    apps = entity.get("applications", {})
    explanation = build_entity_explanation(entity, composition=composition)
    image_2d = render_molecule_2d_png(smiles) if smiles else None
    cid = entity.get("identifiers", {}).get("pubchem_cid")

    return {
        "profile_kind": ENTITY_CHEMICAL_COMPOUND,
        "has_molecular_structure": True,
        "basic_info": {
            "name_en": entity.get("name_en", ""),
            "name_ar": entity.get("name_ar", ""),
            "iupac_name": entity.get("iupac_name", ""),
            "formula": formula,
            "molar_mass": f"{molar_mass} g/mol" if molar_mass else "غير متوفر",
            "molar_mass_value": molar_mass,
            "entity_type_label_ar": entity.get("entity_type_label_ar", "مركب كيميائي"),
            "entity_type_label_en": entity.get("entity_type_label_en", "Chemical Compound"),
        },
        "composition": {
            "summary_ar": composition.get("summary_ar", ""),
            "total_atoms": composition.get("total_atoms", 0),
            "components": composition.get("components", []),
            "calculation_provenance": build_provenance_tag(
                OP_SCIENTIFIC_CALCULATION,
                sources=["IUPAC Atomic Weights", "Formula Analyzer"],
                notes_ar="تم حساب المكونات الذرية والنسب المئوية الكتلية والكتلة المولية حتميا من الصيغة الجزيئية.",
            ),
        },
        "structure": {
            "molecular_geometry": struct.get("molecular_geometry", "بنية جزيئية محددة وفق تكافؤ الذرات"),
            "hybridization": struct.get("hybridization", "غير محدد"),
            "chemical_bonds": struct.get("chemical_bonds", "روابط كيميائية بين الذرات المكونة"),
            "dipole_moment": struct.get("dipole_moment", "غير محدد"),
            "representation_3d_note": struct.get(
                "representation_3d_note",
                f"يتوفر التمثيل الفراغي ثلاثي الأبعاد عبر PubChem 3D Conformer (CID: {cid})." if cid else "متاح عبر محرك البنية الجزيئية.",
            ),
            "pubchem_3d_url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=3D-Conformer" if cid else None,
            "image_2d_png": image_2d,
        },
        "properties": {
            "physical_state": phys.get("physical_state", "غير محدد"),
            "density": phys.get("density", "غير محدد"),
            "melting_point": phys.get("melting_point", "غير محدد"),
            "boiling_point": phys.get("boiling_point", "غير محدد"),
            "chemical_properties_summary": chem.get("acidity_basicity") or chem.get("reactivity", ""),
            "reactivity": chem.get("reactivity", ""),
        },
        "applications": {
            "daily": apps.get("daily", "تطبيقات عامة في الحياة اليومية."),
            "industrial": apps.get("industrial", "تطبيقات صناعية وكيميائية."),
            "scientific": apps.get("scientific", "تطبيقات بحثية ومختبرية."),
        },
        "explanation": explanation,
        "identifiers": entity.get("identifiers", {}),
        "sources": entity.get("sources", ["PubChem"]),
        "provenance": build_provenance_tag(
            OP_SCIENTIFIC_RETRIEVAL,
            sources=entity.get("sources", ["PubChem"]),
            notes_ar="المعلومات المرجعية والخصائص الفيزيائية مسترجعة من المصادر العلمية المعتمدة.",
        ),
        "reaction_lab_entity": {
            "key": entity.get("key", ""),
            "name_ar": entity.get("name_ar", ""),
            "name_en": entity.get("name_en", ""),
            "symbol_or_formula": formula,
            "molar_mass": molar_mass,
        },
    }


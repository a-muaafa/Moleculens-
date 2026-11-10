from chemistry_engine.formula_analyzer import analyze_formula_composition
from chemistry_engine.structure_analyzer import compute_rdkit_molecular_metrics


def enrich_chemical_metrics(symbol_or_formula: str, smiles: str = "") -> dict:
    composition = analyze_formula_composition(symbol_or_formula)
    rdkit_metrics = compute_rdkit_molecular_metrics(smiles) if smiles else None
    molar_mass = composition.get("molar_mass")
    if molar_mass is None and rdkit_metrics is not None:
        molar_mass = rdkit_metrics.get("molecular_weight")

    return {
        "molar_mass": molar_mass,
        "composition": composition,
        "rdkit_metrics": rdkit_metrics,
    }


"""
Independent Scientific Audit Engine for MolecuLens Reaction Registry and Rules.

Provides a deterministic, transparent, multi-criterion scientific audit across:
1. Atom & Mass Conservation (Lavoisier Stoichiometric Balance).
2. Physical State Consistency (s, l, g, aq) vs. Solubility & Reaction Medium.
3. Thermodynamic Sign & Enthalpy Consistency (Exothermic vs. Endothermic vs. Delta H).
4. Reference Operating Window Validity (Temperature & Pressure bounds vs.Optimal).
5. Bibliographic Reference Integrity (Unique Source ID, Author/Handbook Citation, URL).
6. General Chemistry Inference Rules Audit (13 rules with test reaction verification).
"""

from __future__ import annotations

import re
from typing import Any

from chemistry_engine.formula_analyzer import calculate_molar_mass, parse_chemical_formula
from reaction_engine.prediction import RULES, predict_reaction
from reaction_engine.reaction_registry import (
    INSOLUBLE_SOLIDS_AND_METALS,
    VERIFIED_REACTIONS_REGISTRY,
)

VALID_STATES: set[str] = {"s", "l", "g", "aq"}

RULE_REPRESENTATIVE_PROBES: list[tuple[str, list[str], list[str]]] = [
    ("combustion", ["C5H12", "O2"], ["CO2", "H2O"]),
    ("acid_base_neutralization", ["HBr", "LiOH"], ["LiBr", "H2O"]),
    ("acid_carbonate", ["CH3COOH", "K2CO3"], ["KCH3COO", "H2O", "CO2"]),
    ("acid_oxide", ["HBr", "MgO"], ["MgBr2", "H2O"]),
    ("acid_ammonia", ["HBr", "NH3"], ["NH4Br"]),
    ("acid_active_metal", ["H2SO4", "Mg"], ["MgSO4", "H2"]),
    ("alkali_water", ["K", "H2O"], ["KOH", "H2"]),
    ("metal_oxygen", ["Li", "O2"], ["Li2O"]),
    ("metal_halogen", ["Ca", "Cl2"], ["CaCl2"]),
    ("salt_salt_precipitation", ["Pb(NO3)2", "KI"], ["PbI2", "KNO3"]),
    ("salt_hydroxide_precipitation", ["AlCl3", "NaOH"], ["Al(OH)3", "NaCl"]),
    ("salt_metal_displacement", ["Zn", "Pb(NO3)2"], ["Zn(NO3)2", "Pb"]),
    ("esterification", ["HCOOH", "CH3OH"], ["C2H4O2", "H2O"]),
]


def _count_side_atoms(species_list: list[dict[str, Any]]) -> dict[str, int]:
    totals: dict[str, int] = {}
    for sp in species_list:
        coeff = int(sp.get("coefficient", 1) or 1)
        atoms = parse_chemical_formula(sp.get("formula", ""))
        for el, count in atoms.items():
            totals[el] = totals.get(el, 0) + coeff * count
    return totals


def _side_molar_mass_sum(species_list: list[dict[str, Any]]) -> float:
    total = 0.0
    for sp in species_list:
        coeff = float(sp.get("coefficient", 1) or 1)
        mm = calculate_molar_mass(sp.get("formula", ""))
        total += coeff * mm
    return round(total, 3)


def audit_single_reaction(rxn: dict[str, Any]) -> dict[str, Any]:
    """
    Runs an independent 5-point scientific audit on a single reaction record.
    """
    rxn_id = rxn.get("reaction_id", "unknown")
    reactants = rxn.get("reactants", [])
    products = rxn.get("products", [])
    conditions = rxn.get("conditions", {}) or {}
    spec = rxn.get("reference_conditions_spec", {}) or {}
    category = rxn.get("category", "")

    errors: list[str] = []
    warnings: list[str] = []
    checks: dict[str, bool] = {}

    # 1. Atom & Mass Conservation
    r_atoms = _count_side_atoms(reactants)
    p_atoms = _count_side_atoms(products)
    atoms_balanced = bool(r_atoms) and (r_atoms == p_atoms)
    r_mass = _side_molar_mass_sum(reactants)
    p_mass = _side_molar_mass_sum(products)
    mass_balanced = abs(r_mass - p_mass) <= 0.05
    checks["atom_and_mass_conservation"] = atoms_balanced and mass_balanced
    if not atoms_balanced:
        errors.append(f"[{rxn_id}] اختلال موازنة الذرات: المتفاعلات {r_atoms} مقابل النواتج {p_atoms}.")
    if not mass_balanced:
        errors.append(f"[{rxn_id}] اختلال حفظ الكتلة المولية: {r_mass} g مقابل {p_mass} g.")

    # 2. Physical State Validity & Phase Consistency
    states_ok = True
    cond_text = " ".join(str(v) for v in conditions.values())
    for side_name, sp_list in (("reactants", reactants), ("products", products)):
        for sp in sp_list:
            st = sp.get("state", "")
            formula = sp.get("formula", "")
            if st not in VALID_STATES:
                states_ok = False
                errors.append(f"[{rxn_id}] حالة فيزيائية غير صالحة ({st}) للمادة {formula} في {side_name}.")
            if category == "precipitation" and side_name == "products":
                if formula in {"AgCl", "BaSO4", "CaCO3", "Fe(OH)3", "Cu(OH)2"} and st != "s":
                    states_ok = False
                    errors.append(f"[{rxn_id}] الراسب {formula} يجب أن يكون في الحالة الصلبة (s) وليس ({st}).")
            if category in {"neutralization", "precipitation", "single_displacement"} and "طور غازي" not in cond_text:
                if formula in {"Mg(OH)2", "CuO", "CaCO3"} and st == "aq":
                    states_ok = False
                    errors.append(f"[{rxn_id}] المادة شحيحة الذوبان {formula} صُنفت خطأ كمحلول مائي (aq).")
    checks["physical_state_consistency"] = states_ok

    # 3. Thermodynamic Sign & Enthalpy Consistency
    thermo_ar = str(conditions.get("thermodynamics_ar", "") or "")
    thermo_ok = True
    has_quantitative_enthalpy = False
    delta_h_value: float | None = None
    if thermo_ar:
        m_num = re.search(r"([+-]?\d+(?:\.\d+)?)\s*kJ", thermo_ar, flags=re.IGNORECASE)
        if m_num:
            has_quantitative_enthalpy = True
            delta_h_value = float(m_num.group(1))
            if ("طارد" in thermo_ar or "باعث" in thermo_ar) and delta_h_value > 0:
                thermo_ok = False
                errors.append(f"[{rxn_id}] تناقض حراري: وُصف التفاعل بأنه طارد للحرارة بينما ΔH موجبة ({delta_h_value} kJ).")
            if "ماص" in thermo_ar and delta_h_value < 0:
                thermo_ok = False
                errors.append(f"[{rxn_id}] تناقض حراري: وُصف التفاعل بأنه ماص للحرارة بينما ΔH سالبة ({delta_h_value} kJ).")
    checks["thermodynamic_consistency"] = thermo_ok

    # 4. Reference Conditions Operating Window Validity
    window_ok = True
    if spec:
        t_min = float(spec.get("temp_c_min", 0.0))
        t_max = float(spec.get("temp_c_max", 0.0))
        t_opt = float(spec.get("temp_c_optimal", 0.0))
        p_min = float(spec.get("pressure_atm_min", 0.0))
        p_max = float(spec.get("pressure_atm_max", 0.0))
        p_opt = float(spec.get("pressure_atm_optimal", 0.0))
        if not (t_min <= t_opt <= t_max):
            window_ok = False
            errors.append(f"[{rxn_id}] نطاق الحرارة المرجعي غير متسق: [{t_min}, {t_opt}, {t_max}] °C.")
        if not (0.0 < p_min <= p_opt <= p_max):
            window_ok = False
            errors.append(f"[{rxn_id}] نطاق الضغط المرجعي غير متسق: [{p_min}, {p_opt}, {p_max}] atm.")
    else:
        window_ok = False
        errors.append(f"[{rxn_id}] يفتقر إلى مواصفات الظروف المرجعية (reference_conditions_spec).")
    checks["reference_window_validity"] = window_ok

    # 5. Bibliographic Reference Completeness
    source_id = str(rxn.get("source_id", "") or "").strip()
    source_text = str(rxn.get("source", "") or "").strip()
    ref_url = str(rxn.get("reference_url", "") or "").strip()
    doc_type = str(rxn.get("documentation_type", "") or "").strip()
    ref_ok = (
        len(source_id) >= 5
        and len(source_text) >= 20
        and ref_url.startswith("https://")
        and doc_type in {"experimentally_documented", "literature_referenced", "theoretical_model"}
    )
    if not ref_ok:
        errors.append(f"[{rxn_id}] نقص في بيانات التوثيق المرجعي (source_id / source / reference_url).")
    checks["bibliographic_reference_integrity"] = ref_ok

    return {
        "reaction_id": rxn_id,
        "reaction_name_ar": rxn.get("reaction_name_ar", ""),
        "category": category,
        "passed": len(errors) == 0,
        "checks": checks,
        "stoichiometric_mass_g": r_mass,
        "atom_counts": r_atoms,
        "has_quantitative_enthalpy": has_quantitative_enthalpy,
        "delta_h_kj": delta_h_value,
        "documentation_type": doc_type,
        "source_id": source_id,
        "source": source_text,
        "reference_url": ref_url,
        "errors": errors,
        "warnings": warnings,
    }


def audit_registry_and_rules() -> dict[str, Any]:
    """
    Performs a complete independent scientific audit across all 87 registry reactions
    and all 13 general chemistry inference rules.
    """
    reaction_reports: list[dict[str, Any]] = []
    all_errors: list[str] = []
    source_ids: set[str] = set()
    sources: set[str] = set()
    doc_type_counts: dict[str, int] = {}
    category_counts: dict[str, int] = {}
    quantitative_enthalpy_count = 0

    for rxn in VERIFIED_REACTIONS_REGISTRY:
        rep = audit_single_reaction(rxn)
        reaction_reports.append(rep)
        all_errors.extend(rep["errors"])

        sid = rep["source_id"]
        src = rep["source"]
        if sid in source_ids:
            all_errors.append(f"[{rep['reaction_id']}] تكرار معرّف المرجع source_id: {sid}")
        source_ids.add(sid)
        if src in sources:
            all_errors.append(f"[{rep['reaction_id']}] تكرار نص المرجع العلمي source: {src}")
        sources.add(src)

        dt = rep["documentation_type"]
        doc_type_counts[dt] = doc_type_counts.get(dt, 0) + 1
        cat = rep["category"]
        category_counts[cat] = category_counts.get(cat, 0) + 1
        if rep["has_quantitative_enthalpy"]:
            quantitative_enthalpy_count += 1

    # Audit the 13 general-chemistry rules via non-registry probe reactions
    rule_reports: list[dict[str, Any]] = []
    for rule_id, reactant_formulas, expected_products in RULE_REPRESENTATIVE_PROBES:
        rule_meta = RULES.get(rule_id, {})
        inputs = [{"substance": {"symbol_or_formula": f}} for f in reactant_formulas]
        pred = predict_reaction(inputs)
        got_products = [p["formula"] for p in pred.get("products", [])]
        atom_cons = pred.get("atom_conservation", {}) or {}
        eq_unicode = pred.get("equation_unicode", "") or pred.get("equation_text", "")
        rule_passed = (
            bool(rule_meta)
            and pred.get("supported") is True
            and pred.get("mode") in {"rule_inferred", "documented"}
            and set(got_products) == set(expected_products)
            and bool(atom_cons.get("is_balanced"))
            and bool(rule_meta.get("source"))
            and str(rule_meta.get("reference_url", "")).startswith("https://")
        )
        if not rule_passed:
            all_errors.append(
                f"[Rule:{rule_id}] فشل تدقيق القاعدة على المتفاعلات {reactant_formulas}: النواتج {got_products} مقابل المتوقع {expected_products}."
            )
        rule_reports.append(
            {
                "rule_id": rule_id,
                "rule_name_ar": rule_meta.get("rule_name_ar", rule_id),
                "probe_reactants": reactant_formulas,
                "predicted_products": got_products,
                "balanced_equation": eq_unicode,
                "source_id": rule_meta.get("source_id", ""),
                "source": rule_meta.get("source", ""),
                "passed": rule_passed,
            }
        )

    passed_reactions = sum(1 for r in reaction_reports if r["passed"])
    passed_rules = sum(1 for r in rule_reports if r["passed"])

    return {
        "overall_passed": len(all_errors) == 0,
        "total_reactions": len(reaction_reports),
        "passed_reactions": passed_reactions,
        "unique_source_ids": len(source_ids),
        "unique_bibliographic_sources": len(sources),
        "quantitative_enthalpy_reactions": quantitative_enthalpy_count,
        "documentation_type_breakdown": doc_type_counts,
        "category_breakdown": category_counts,
        "total_rules": len(RULES),
        "passed_rules": passed_rules,
        "errors": all_errors,
        "reaction_reports": reaction_reports,
        "rule_reports": rule_reports,
    }

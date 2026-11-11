from typing import Any
from chemistry_engine.formula_analyzer import calculate_molar_mass, formulas_are_equivalent, is_gaseous_at_stp
from config.settings import MOLAR_VOLUME_STP_LITERS


def convert_amount_to_moles(
    amount: float,
    unit: str,
    molar_mass: float,
    formula: str = "",
    substance: dict[str, Any] | None = None,
    rxn_state: str = "",
) -> float:
    val = float(amount)
    if unit == "mol":
        return val
    if unit == "mmol":
        return val / 1000.0
    if unit == "g":
        return val / molar_mass if molar_mass > 0 else 0.0
    if unit == "mg":
        return (val / 1000.0) / molar_mass if molar_mass > 0 else 0.0
    if unit == "L_STP":
        if (formula or substance or rxn_state) and not is_gaseous_at_stp(formula, substance=substance, rxn_state=rxn_state):
            raise ValueError(
                f"Scientific Error: STP gas volume conversion (L_STP) cannot be applied to non-gaseous substance '{formula}'."
            )
        return val / MOLAR_VOLUME_STP_LITERS
    raise ValueError(f"Unsupported unit: {unit}")


def compute_reaction_stoichiometry(
    rxn_reactants: list[dict[str, Any]],
    user_inputs: list[dict[str, Any]],
) -> dict[str, Any]:
    enriched_reactants: list[dict[str, Any]] = []
    extents: list[tuple[float, dict[str, Any]]] = []

    for r_spec in rxn_reactants:
        formula = r_spec["formula"]
        coeff = int(r_spec["coefficient"])
        mm = calculate_molar_mass(formula) or 1.0
        r_state = r_spec.get("state", "")

        matched_user = None
        for u_item in user_inputs:
            u_sub = u_item.get("substance", {})
            u_form = u_sub.get("symbol_or_formula") or u_item.get("formula", "")
            if formulas_are_equivalent(formula, u_form):
                matched_user = u_item
                break

        input_amount = float(matched_user["amount"]) if matched_user else float(coeff)
        input_unit = matched_user.get("unit", "mol") if matched_user else "mol"
        u_sub_dict = matched_user.get("substance") if matched_user else None
        initial_moles = convert_amount_to_moles(
            input_amount,
            input_unit,
            mm,
            formula=formula,
            substance=u_sub_dict,
            rxn_state=r_state,
        )
        extent_i = initial_moles / coeff

        item_record = {
            "formula": formula,
            "name_ar": r_spec["name_ar"],
            "name_en": r_spec["name_en"],
            "coefficient": coeff,
            "molar_mass": round(mm, 4),
            "input_amount": input_amount,
            "input_unit": input_unit,
            "initial_mol": round(initial_moles, 6),
            "initial_mass_g": round(initial_moles * mm, 4),
            "extent_ratio": extent_i,
        }
        enriched_reactants.append(item_record)
        extents.append((extent_i, item_record))

    min_extent = min(x[0] for x in extents)
    max_extent = max(x[0] for x in extents)
    is_stoichiometric = abs(max_extent - min_extent) <= 1e-6 or len(enriched_reactants) == 1

    limiting_reagents: list[dict[str, Any]] = []
    excess_reagents: list[dict[str, Any]] = []
    reactant_table: list[dict[str, Any]] = []

    for item in enriched_reactants:
        coeff = item["coefficient"]
        mm = item["molar_mass"]
        init_mol = item["initial_mol"]
        consumed_mol = min_extent * coeff
        remaining_mol = max(0.0, init_mol - consumed_mol)
        if remaining_mol < 1e-7:
            remaining_mol = 0.0

        is_limiting = abs(item["extent_ratio"] - min_extent) <= 1e-6
        if is_stoichiometric:
            role_ar = "متفاعل بنسبة متكافئة تماماً (Stoichiometric Ratio)"
        elif is_limiting:
            role_ar = "المتفاعل المحدد (Limiting Reagent)"
            limiting_reagents.append(item)
        else:
            role_ar = "متفاعل فائض (Excess Reagent)"
            excess_reagents.append(item)

        row = {
            **item,
            "consumed_mol": round(consumed_mol, 4),
            "consumed_mass_g": round(consumed_mol * mm, 4),
            "remaining_mol": round(remaining_mol, 4),
            "remaining_mass_g": round(remaining_mol * mm, 4),
            "role_ar": role_ar,
        }
        reactant_table.append(row)

    return {
        "reaction_extent_mol": round(min_extent, 6),
        "is_exact_stoichiometric_ratio": is_stoichiometric,
        "limiting_reagent": limiting_reagents[0] if limiting_reagents and not is_stoichiometric else None,
        "excess_reagents": excess_reagents if not is_stoichiometric else [],
        "reactants_breakdown": reactant_table,
    }

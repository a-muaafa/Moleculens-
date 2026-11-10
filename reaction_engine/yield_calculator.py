from typing import Any
from chemistry_engine.formula_analyzer import calculate_molar_mass


def compute_theoretical_yields(
    rxn_products: list[dict[str, Any]],
    reaction_extent_mol: float,
) -> list[dict[str, Any]]:
    yields_list: list[dict[str, Any]] = []
    for p_spec in rxn_products:
        formula = p_spec["formula"]
        coeff = int(p_spec["coefficient"])
        mm = calculate_molar_mass(formula) or 1.0
        produced_mol = reaction_extent_mol * coeff
        produced_mass_g = produced_mol * mm

        yields_list.append(
            {
                "formula": formula,
                "name_ar": p_spec["name_ar"],
                "name_en": p_spec["name_en"],
                "coefficient": coeff,
                "molar_mass": round(mm, 4),
                "theoretical_yield_mol": round(produced_mol, 4),
                "theoretical_yield_g": round(produced_mass_g, 4),
                "state": p_spec.get("state", ""),
            }
        )
    return yields_list


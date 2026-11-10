from fractions import Fraction
import math
from chemistry_engine.formula_analyzer import parse_chemical_formula


def verify_atom_conservation(
    reactants: list[dict],
    products: list[dict],
) -> dict:
    reactant_atoms: dict[str, int] = {}
    product_atoms: dict[str, int] = {}

    for item in reactants:
        coeff = int(item.get("coefficient", 1))
        counts = parse_chemical_formula(item["formula"]) or {}
        for el, cnt in counts.items():
            reactant_atoms[el] = reactant_atoms.get(el, 0) + coeff * cnt

    for item in products:
        coeff = int(item.get("coefficient", 1))
        counts = parse_chemical_formula(item["formula"]) or {}
        for el, cnt in counts.items():
            product_atoms[el] = product_atoms.get(el, 0) + coeff * cnt

    is_balanced = reactant_atoms == product_atoms and len(reactant_atoms) > 0
    return {
        "is_balanced": is_balanced,
        "reactant_atoms": reactant_atoms,
        "product_atoms": product_atoms,
    }


def balance_reaction_algebraically(
    reactant_formulas: list[str],
    product_formulas: list[str],
) -> tuple[list[int], list[int]] | None:
    all_formulas = reactant_formulas + product_formulas
    n_species = len(all_formulas)
    n_react = len(reactant_formulas)
    parsed_list = [parse_chemical_formula(f) for f in all_formulas]
    if any(p is None for p in parsed_list):
        return None

    elements = sorted({el for p in parsed_list if p for el in p.keys()})
    if not elements:
        return None

    matrix: list[list[Fraction]] = []
    for el in elements:
        row: list[Fraction] = []
        for idx, p in enumerate(parsed_list):
            cnt = (p or {}).get(el, 0)
            row.append(Fraction(cnt if idx < n_react else -cnt))
        matrix.append(row)

    n_rows = len(matrix)
    pivot_Row = 0
    pivot_cols: list[int] = []
    for col in range(n_species - 1):
        if pivot_Row >= n_rows:
            break
        sel = None
        for r in range(pivot_Row, n_rows):
            if matrix[r][col] != 0:
                sel = r
                break
        if sel is None:
            continue
        matrix[pivot_Row], matrix[sel] = matrix[sel], matrix[pivot_Row]
        piv_val = matrix[pivot_Row][col]
        matrix[pivot_Row] = [val / piv_val for val in matrix[pivot_Row]]
        for r in range(n_rows):
            if r != pivot_Row and matrix[r][col] != 0:
                factor = matrix[r][col]
                matrix[r] = [matrix[r][c] - factor * matrix[pivot_Row][c] for c in range(n_species)]
        pivot_cols.append(col)
        pivot_Row += 1

    if len(pivot_cols) != n_species - 1:
        return None

    free_col = n_species - 1
    sol_fracs = [Fraction(0)] * n_species
    sol_fracs[free_col] = Fraction(1)
    for idx_r, p_col in enumerate(pivot_cols):
        sol_fracs[p_col] = -matrix[idx_r][free_col]

    if any(val <= 0 for val in sol_fracs):
        return None

    lcm_denom = 1
    for val in sol_fracs:
        lcm_denom = (lcm_denom * val.denominator) // math.gcd(lcm_denom, val.denominator)

    int_coeffs = [int(val * lcm_denom) for val in sol_fracs]
    common_gcd = int_coeffs[0]
    for c in int_coeffs[1:]:
        common_gcd = math.gcd(common_gcd, c)
    int_coeffs = [c // common_gcd for c in int_coeffs]

    return int_coeffs[:n_react], int_coeffs[n_react:]


def _formula_to_latex(formula: str) -> str:
    import re
    parts = re.findall(r"([A-Za-z]+|\d+)", formula)
    out = []
    for p in parts:
        if p.isdigit():
            out.append(f"_{{{p}}}" if len(p) > 1 else f"_{p}")
        else:
            out.append(f"\\text{{{p}}}")
    return "".join(out)


def format_balanced_equation(reactants: list[dict], products: list[dict]) -> dict[str, str]:
    r_text_parts = []
    r_latex_parts = []
    for r in reactants:
        c = int(r.get("coefficient", 1))
        f = r["formula"]
        r_text_parts.append(f"{c if c > 1 else ''}{f}")
        r_latex_parts.append(f"{c if c > 1 else ''}{_formula_to_latex(f)}")

    p_text_parts = []
    p_latex_parts = []
    for p in products:
        c = int(p.get("coefficient", 1))
        f = p["formula"]
        p_text_parts.append(f"{c if c > 1 else ''}{f}")
        p_latex_parts.append(f"{c if c > 1 else ''}{_formula_to_latex(f)}")

    return {
        "text": " + ".join(r_text_parts) + " -> " + " + ".join(p_text_parts),
        "latex": " + ".join(r_latex_parts) + " \\rightarrow " + " + ".join(p_latex_parts),
    }


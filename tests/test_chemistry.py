from chemistry_engine.formula_analyzer import (
    analyze_formula_composition,
    calculate_molar_mass,
    format_chemical_equation_display,
    format_chemical_formula_display,
    formulas_are_equivalent,
    is_gaseous_at_stp,
    normalize_chemical_formula,
    parse_chemical_formula,
)
from chemistry_engine.periodic_table import VALID_ELEMENTS, get_element_basic_info
from chemistry_engine.structure_analyzer import verify_substance_structure_identity
from search.local_catalog import get_catalog_entity_by_key


def test_periodic_table_elements_and_formula_parser():
    assert len(VALID_ELEMENTS) == 118
    fe_info = get_element_basic_info("Fe")
    assert fe_info is not None
    assert fe_info["atomic_number"] == 26
    assert parse_chemical_formula("CH3COOH") == {"C": 2, "H": 4, "O": 2}
    assert normalize_chemical_formula("CH3COOH") == "C2H4O2"
    assert formulas_are_equivalent("CH3COOH", "C2H4O2") is True
    assert abs((calculate_molar_mass("H2O") or 0.0) - 18.015) < 0.01


def test_subscript_and_equation_display_formatting():
    assert format_chemical_formula_display("H2O") == "H₂O"
    assert format_chemical_formula_display("CH3COOH") == "CH₃COOH"
    assert format_chemical_formula_display("Ca(OH)2") == "Ca(OH)₂"
    assert format_chemical_formula_display("SO4^2-") == "SO₄²⁻"
    eq = format_chemical_equation_display("2H2 + O2 -> 2H2O")
    assert eq == "2H₂ + O₂ → 2H₂O"


def test_stp_gas_phase_discrimination():
    assert is_gaseous_at_stp("H2") is True
    assert is_gaseous_at_stp("O2") is True
    assert is_gaseous_at_stp("CO2") is True
    assert is_gaseous_at_stp("CH4") is True
    assert is_gaseous_at_stp("H2O") is False
    assert is_gaseous_at_stp("NaCl") is False
    assert is_gaseous_at_stp("Fe") is False
    assert is_gaseous_at_stp("CH3COOH") is False


def test_compound_identity_verification():
    acetic = get_catalog_entity_by_key("acetic_acid")
    ver = verify_substance_structure_identity(acetic)
    assert ver["is_verified"] is True
    assert ver["formula_matches"] is True

    # Tampered formula should fail identity verification
    tampered = dict(acetic, symbol_or_formula="C6H12O6")
    ver_bad = verify_substance_structure_identity(tampered)
    assert ver_bad["is_verified"] is False
    assert ver_bad["formula_matches"] is False

    comp = analyze_formula_composition("H2SO4")
    assert comp["is_valid"] is True
    assert comp["formatted_formula"] == "H₂SO₄"

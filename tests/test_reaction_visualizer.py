"""
Unit tests for Reaction Visualizer module in Chemical Playground.
Verifies SVG generation for chemical species, atom conservation auditing,
bond transformation analysis, and Plotly reaction energy coordinate diagrams.
"""

import plotly.graph_objects as go
from reaction_engine.visualizer import ReactionVisualizer


def test_get_species_svg_known_and_fallback():
    # Known compound with SMILES (e.g. H2O or CH3COOH)
    svg_water = ReactionVisualizer.get_species_svg("H2O")
    assert "<svg" in svg_water
    assert "</svg>" in svg_water

    # Formula without RDKit representation falls back cleanly to atomic badge SVG
    svg_fallback = ReactionVisualizer.get_species_svg("CustomFormulaX")
    assert "<svg" in svg_fallback
    assert "CustomFormulaX" in svg_fallback


def test_audit_atom_conservation_balanced_and_unbalanced():
    # Balanced combustion of hydrogen: 2H2 + O2 -> 2H2O
    reactants_balanced = [
        {"formula": "H2", "coefficient": 2},
        {"formula": "O2", "coefficient": 1},
    ]
    products_balanced = [
        {"formula": "H2O", "coefficient": 2},
    ]
    audit_bal = ReactionVisualizer.audit_atom_conservation(reactants_balanced, products_balanced)
    assert len(audit_bal) == 2  # H and O
    for row in audit_bal:
        assert row["is_conserved"] is True
        assert "محفوظ 100%" in row["status_ar"]

    # Unbalanced: H2 + O2 -> H2O
    reactants_unbalanced = [
        {"formula": "H2", "coefficient": 1},
        {"formula": "O2", "coefficient": 1},
    ]
    products_unbalanced = [
        {"formula": "H2O", "coefficient": 1},
    ]
    audit_unbal = ReactionVisualizer.audit_atom_conservation(reactants_unbalanced, products_unbalanced)
    o_row = next((r for r in audit_unbal if r["element_symbol"] == "O"), None)
    assert o_row is not None
    assert o_row["is_conserved"] is False


def test_analyze_bond_transformations_hydrogen_combustion():
    reactants = [{"formula": "H2"}, {"formula": "O2"}]
    products = [{"formula": "H2O"}]
    analysis = ReactionVisualizer.analyze_bond_transformations(reactants, products)

    assert "bonds_broken" in analysis
    assert "bonds_formed" in analysis
    assert len(analysis["bonds_broken"]) > 0
    assert len(analysis["bonds_formed"]) > 0
    assert analysis["is_exothermic"] is True
    assert analysis["estimated_delta_h_kj"] < 0
    assert "طارد للحرارة" in analysis["thermodynamic_type_ar"]


def test_create_reaction_energy_diagram():
    fig = ReactionVisualizer.create_reaction_energy_diagram(
        reactants_str="2H2 + O2",
        products_str="2H2O",
        is_exothermic=True,
        delta_h_kj=-483.6,
        activation_energy_kj=65.0,
    )
    assert isinstance(fig, go.Figure)
    assert len(fig.data) > 0
    # Check annotations exist for Delta H
    layout_dict = fig.to_dict()["layout"]
    assert "annotations" in layout_dict
    assert len(layout_dict["annotations"]) >= 1

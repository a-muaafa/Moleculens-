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
    assert analysis["is_illustrative"] is False
    assert analysis["calculation_source"] == "bond_dissociation_table"
    assert analysis["is_exothermic"] is True
    assert analysis["estimated_delta_h_kj"] < 0
    assert "طارد للحرارة" in analysis["thermodynamic_type_ar"]


def test_analyze_bond_transformations_illustrative_when_incomplete():
    reactants = [{"formula": "AgNO3"}, {"formula": "NaCl"}]
    products = [{"formula": "AgCl"}, {"formula": "NaNO3"}]
    analysis = ReactionVisualizer.analyze_bond_transformations(reactants, products)

    assert analysis["is_illustrative"] is True
    assert analysis["calculation_source"] == "illustrative_qualitative"
    assert analysis["estimated_delta_h_kj"] is None
    assert analysis["total_broken_energy_kj"] is None
    assert analysis["total_formed_energy_kj"] is None


def test_create_reaction_energy_diagram():
    fig_no_ea = ReactionVisualizer.create_reaction_energy_diagram(
        reactants_str="2H2 + O2",
        products_str="2H2O",
        is_exothermic=True,
        delta_h_kj=-483.6,
        activation_energy_kj=None,
        is_illustrative=False,
    )
    assert isinstance(fig_no_ea, go.Figure)
    assert len(fig_no_ea.data) > 0
    layout_dict = fig_no_ea.to_dict()["layout"]
    assert "annotations" in layout_dict
    assert len(layout_dict["annotations"]) >= 1
    assert layout_dict["yaxis"]["showticklabels"] is False
    all_text = " ".join(
        str(ann.get("text", "")) for ann in layout_dict["annotations"]
    ) + " " + " ".join(str(getattr(tr, "text", "")) for tr in fig_no_ea.data)
    assert "65.0" not in all_text
    assert "-483.6" in all_text

    # When both Delta H and Ea are provided, Y axis displays quantitative kJ ticks
    fig_full = ReactionVisualizer.create_reaction_energy_diagram(
        reactants_str="2H2 + O2",
        products_str="2H2O",
        is_exothermic=True,
        delta_h_kj=-483.6,
        activation_energy_kj=85.0,
        is_illustrative=False,
    )
    assert fig_full.to_dict()["layout"]["yaxis"]["showticklabels"] is True


def test_create_illustrative_energy_diagram_hides_fake_numbers():
    fig = ReactionVisualizer.create_reaction_energy_diagram(
        reactants_str="A + B",
        products_str="C + D",
        is_exothermic=None,
        delta_h_kj=None,
        activation_energy_kj=None,
        is_illustrative=True,
    )
    assert isinstance(fig, go.Figure)
    layout_dict = fig.to_dict()["layout"]
    assert layout_dict["yaxis"]["showticklabels"] is False
    annotations_text = " ".join(str(ann.get("text", "")) for ann in layout_dict.get("annotations", []))
    assert "kJ" not in annotations_text
    assert "-150" not in annotations_text


def test_end_to_end_experiment_stages_and_mass_conservation():
    from reaction_engine.prediction import predict_reaction

    pred = predict_reaction([
        {"substance": {"symbol_or_formula": "HCl"}, "amount": 2.0, "unit": "mol"},
        {"substance": {"symbol_or_formula": "NaOH"}, "amount": 1.0, "unit": "mol"},
    ])
    assert pred["supported"] is True

    stages = ReactionVisualizer.build_end_to_end_experiment_stages(pred)
    assert "stage1_preparation" in stages
    assert "stage2_active_vessel" in stages
    assert "stage3_mass_balance" in stages
    assert "stage4_final_inventory" in stages

    s1 = stages["stage1_preparation"]
    assert len(s1["reactant_cards"]) == 2

    s2 = stages["stage2_active_vessel"]
    assert "<svg" in s2["vessel_svg"]
    assert "HCl" in s2["vessel_svg"]

    s3 = stages["stage3_mass_balance"]
    assert s3["mass_conserved"] is True
    assert abs(s3["total_initial_mass_g"] - s3["total_final_vessel_mass_g"]) <= 0.05
    # Since HCl had 2 mol and NaOH had 1 mol, 1 mol of HCl remains in excess
    assert s3["total_remaining_mass_g"] > 30.0

    s4 = stages["stage4_final_inventory"]
    # 2 products (NaCl, H2O) + 1 excess reactant (HCl) = 3 items in final vessel
    assert len(s4["final_inventory"]) == 3

    fig_mass = ReactionVisualizer.create_stoichiometry_mass_chart(s3)
    assert isinstance(fig_mass, go.Figure)
    assert len(fig_mass.data) >= 3


def test_independent_scientific_audit_registry_and_rules():
    from reaction_engine.scientific_audit import audit_registry_and_rules

    report = audit_registry_and_rules()
    assert report["overall_passed"] is True, f"Audit errors: {report['errors']}"
    assert report["total_reactions"] == 87
    assert report["passed_reactions"] == 87
    assert report["unique_source_ids"] == 87
    assert report["unique_bibliographic_sources"] == 87
    assert report["total_rules"] == 13
    assert report["passed_rules"] == 13
    assert report["quantitative_enthalpy_reactions"] >= 60


def test_preset8_ch3cooh_k2co3_and_resolve_token():
    from frontend.prediction_view import _resolve_token
    from reaction_engine.prediction import predict_reaction

    ent, note = _resolve_token("K2CO3")
    assert ent is not None
    assert ent["symbol_or_formula"] == "K2CO3"
    assert ent["molar_mass"] > 130.0

    pred = predict_reaction([
        {"substance": {"symbol_or_formula": "CH3COOH"}},
        {"substance": ent},
    ])
    assert pred["supported"] is True
    assert pred["mode"] == "rule_inferred"
    prod_formulas = {p["formula"] for p in pred["products"]}
    assert prod_formulas == {"KCH3COO", "H2O", "CO2"}
    assert pred["atom_conservation"]["is_balanced"] is True


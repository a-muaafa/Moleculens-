"""
Unit tests for the TitrationEngine module in MolecuLens.
Verifies thermodynamic acid-base equilibria, equivalence point calculations,
buffer half-equivalence relations, and indicator matching.
"""

import pytest
from titration.titration_engine import TitrationEngine


def test_strong_acid_strong_base_titration():
    res = TitrationEngine.simulate_titration(
        acid_key="HCl",
        base_key="NaOH",
        acid_conc_m=0.100,
        acid_volume_ml=25.0,
        base_conc_m=0.100,
    )
    assert res["success"] is True
    assert res["v_eq_ml"] == pytest.approx(25.0, abs=0.01)
    assert res["ph_eq"] == pytest.approx(7.00, abs=0.05)
    assert res["initial_ph"] == pytest.approx(1.00, abs=0.05)
    assert res["half_v_eq_ml"] is None


def test_weak_acid_strong_base_titration():
    res = TitrationEngine.simulate_titration(
        acid_key="CH3COOH",
        base_key="NaOH",
        acid_conc_m=0.100,
        acid_volume_ml=25.0,
        base_conc_m=0.100,
    )
    assert res["success"] is True
    assert res["v_eq_ml"] == pytest.approx(25.0, abs=0.01)
    assert res["ph_eq"] > 8.0
    assert res["half_v_eq_ml"] == pytest.approx(12.5, abs=0.01)
    assert res["half_ph"] == pytest.approx(4.76, abs=0.05)
    assert res["initial_ph"] > 2.0


def test_weak_base_strong_acid_titration():
    res = TitrationEngine.simulate_titration(
        acid_key="HCl",
        base_key="NH3",
        acid_conc_m=0.100,
        acid_volume_ml=25.0,
        base_conc_m=0.100,
    )
    assert res["success"] is True
    assert res["v_eq_ml"] == pytest.approx(25.0, abs=0.01)
    assert res["ph_eq"] < 6.0


def test_plotly_curve_generation():
    sim = TitrationEngine.simulate_titration("HCl", "NaOH", 0.1, 25.0, 0.1)
    fig = TitrationEngine.create_plotly_curve(sim, current_v=15.0, indicator_key="phenolphthalein")
    assert fig is not None
    assert len(fig.data) >= 2


def test_invalid_parameters():
    with pytest.raises(ValueError):
        TitrationEngine.simulate_titration("UnknownAcid", "NaOH", 0.1, 25.0, 0.1)

    with pytest.raises(ValueError):
        TitrationEngine.simulate_titration("HCl", "NaOH", -0.1, 25.0, 0.1)


def test_expanded_acids_bases_indicators_and_polyprotic_simulation():
    required_acids = {
        "HCl", "H2SO4", "HNO3", "H3PO4", "CH3COOH", "HCOOH",
        "C6H8O7", "H2C2O4", "C3H6O3", "C7H6O2", "H2CO3", "H3BO3",
    }
    required_bases = {
        "NaOH", "KOH", "Ca(OH)2", "Ba(OH)2", "NH3",
        "Mg(OH)2", "CH3NH2", "C2H5NH2", "NaHCO3", "Na2CO3",
    }
    required_indicators = {
        "phenolphthalein", "methyl_orange", "methyl_red", "bromothymol_blue",
        "bromocresol_green", "bromocresol_purple", "thymol_blue", "litmus",
    }

    assert required_acids <= set(TitrationEngine.ACID_DATABASE.keys())
    assert required_bases <= set(TitrationEngine.BASE_DATABASE.keys())
    assert required_indicators <= set(TitrationEngine.INDICATORS.keys())

    # Polyprotic acid H3PO4 (3 protons) titrated with NaOH: 25 mL of 0.1 M H3PO4 requires 75 mL of 0.1 M NaOH for total equivalence
    sim_poly = TitrationEngine.simulate_titration("H3PO4", "NaOH", 0.1, 25.0, 0.1)
    assert sim_poly["success"] is True
    assert sim_poly["is_polyprotic"] is True
    assert sim_poly["v_eq_ml"] == pytest.approx(75.0, abs=0.1)
    assert sim_poly["step_equivalence_volumes_ml"] == [25.0, 50.0, 75.0]
    assert "endpoint_volume_ml" in sim_poly

    # Flask color SVG simulation (at pH=10.5 above phenolphthalein range 8.2-10.0, basic pink #EC4899 is shown)
    svg = TitrationEngine.create_flask_svg(current_ph=10.5, indicator_key="phenolphthalein", current_vol_ml=26.0, v_eq_ml=25.0)
    assert "<svg" in svg
    assert "#EC4899" in svg


def test_comprehensive_acid_base_and_indicator_catalogue():
    acids = TitrationEngine.ACID_DATABASE
    bases = TitrationEngine.BASE_DATABASE
    indicators = TitrationEngine.INDICATORS

    assert len(acids) >= 40
    assert len(bases) >= 30
    assert len(indicators) >= 23

    # Every entry has mandatory metadata
    for key, spec in acids.items():
        assert "name_ar" in spec and "name_en" in spec
        assert "pKa" in spec and "pKa_list" in spec
        assert "molar_mass" in spec and spec["molar_mass"] > 0
        assert "category_ar" in spec
        assert "source" in spec and len(spec["source"]) >= 10
        assert "solubility_ar" in spec

    for key, spec in bases.items():
        assert "name_ar" in spec and "name_en" in spec
        assert "pKb" in spec and "pKb_list" in spec
        assert "molar_mass" in spec and spec["molar_mass"] > 0
        assert "category_ar" in spec
        assert "source" in spec and len(spec["source"]) >= 10
        assert "solubility_ar" in spec

    for key, spec in indicators.items():
        assert "name_ar" in spec and "name_en" in spec
        assert "pH_range" in spec and len(spec["pH_range"]) == 2
        assert "color_acid" in spec and "color_base" in spec
        assert "source" in spec and len(spec["source"]) >= 10

    # Simulation with new acid and base
    res1 = TitrationEngine.simulate_titration("HBr", "LiOH", 0.1, 20.0, 0.1)
    assert res1["success"] is True
    assert res1["v_eq_ml"] == pytest.approx(20.0, abs=0.01)

    res2 = TitrationEngine.simulate_titration("C3H6O2", "NaOH", 0.1, 25.0, 0.1)
    assert res2["success"] is True
    assert res2["ph_eq"] > 8.0


def test_polyprotic_step_equivalence_buffer_capacity_dilution_and_limitations():
    # 1. Polyprotic H3PO4 stepwise equivalence points
    sim_p = TitrationEngine.simulate_titration("H3PO4", "NaOH", 0.1, 25.0, 0.1)
    steps = sim_p["step_equivalence_points"]
    assert len(steps) == 3
    assert steps[0]["volume_ml"] == pytest.approx(25.0, abs=0.01)
    assert steps[0]["ph"] == pytest.approx(0.5 * (2.15 + 7.20), abs=0.05)
    assert steps[1]["volume_ml"] == pytest.approx(50.0, abs=0.01)
    assert steps[1]["ph"] == pytest.approx(0.5 * (7.20 + 12.35), abs=0.05)

    # 2. Buffer capacity & endpoint error analysis for weak acid
    sim_w = TitrationEngine.simulate_titration("CH3COOH", "NaOH", 0.1, 25.0, 0.1)
    assert sim_w["buffer_analysis"]["has_buffer_region"] is True
    assert sim_w["buffer_analysis"]["max_buffer_capacity_mol_l_ph"] > 0.02
    assert "relative_error_pct" in sim_w["titration_error_analysis"]
    assert len(sim_w["step_by_step_calculation"]) >= 5
    assert len(sim_w["model_limitations_ar"]) >= 3

    # 3. Dilution effect on initial pH (25 mL of 0.1 M HCl diluted with 50 mL water -> 0.0333 M -> pH ~ 1.48)
    sim_undil = TitrationEngine.simulate_titration("HCl", "NaOH", 0.1, 25.0, 0.1, dilution_water_ml=0.0)
    sim_dil = TitrationEngine.simulate_titration("HCl", "NaOH", 0.1, 25.0, 0.1, dilution_water_ml=50.0)
    assert sim_dil["v_eq_ml"] == pytest.approx(sim_undil["v_eq_ml"], abs=0.01)
    assert sim_dil["initial_ph"] > sim_undil["initial_ph"] + 0.4


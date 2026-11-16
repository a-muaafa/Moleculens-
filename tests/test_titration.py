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


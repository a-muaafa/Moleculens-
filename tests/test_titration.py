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


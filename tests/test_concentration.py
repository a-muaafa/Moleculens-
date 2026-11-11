"""
Unit tests for the ConcentrationCalculator module in MolecuLens.
Verifies Molarity, Normality, ppm/ppb, %, and Dilution calculations.
"""

import pytest
from calculators.concentration_calculator import ConcentrationCalculator


def test_calculate_molarity_solve_mass():
    res = ConcentrationCalculator.calculate_molarity(
        molar_mass_g_mol=58.44, volume_l=0.5, molarity_m=0.1
    )
    assert res["success"] is True
    assert res["solved_for"] == "mass_g"
    assert res["mass_g"] == pytest.approx(2.922, abs=0.001)
    assert "recipe_ar" in res and "recipe_en" in res


def test_calculate_molarity_solve_molarity():
    res = ConcentrationCalculator.calculate_molarity(
        mass_g=5.844, molar_mass_g_mol=58.44, volume_l=1.0
    )
    assert res["success"] is True
    assert res["solved_for"] == "molarity_m"
    assert res["molarity_m"] == pytest.approx(0.1, abs=0.0001)


def test_calculate_molarity_solve_volume():
    res = ConcentrationCalculator.calculate_molarity(
        mass_g=5.844, molar_mass_g_mol=58.44, molarity_m=0.1
    )
    assert res["success"] is True
    assert res["solved_for"] == "volume_l"
    assert res["volume_l"] == pytest.approx(1.0, abs=0.001)


def test_calculate_normality():
    res = ConcentrationCalculator.calculate_normality(molarity_m=0.1, equivalents_per_mole=2.0)
    assert res["success"] is True
    assert res["normality_n"] == pytest.approx(0.2, abs=0.0001)


def test_calculate_ppm_ppb():
    res = ConcentrationCalculator.calculate_ppm_ppb(solute_mg=15.0, solution_volume_l=1.0)
    assert res["success"] is True
    assert res["ppm"] == pytest.approx(15.0, abs=0.001)
    assert res["ppb"] == pytest.approx(15000.0, abs=0.1)


def test_calculate_percentage():
    res = ConcentrationCalculator.calculate_percentage(
        percent_type="w/v", solute_amount=5.0, solution_amount=100.0
    )
    assert res["success"] is True
    assert res["percentage"] == pytest.approx(5.0, abs=0.01)


def test_calculate_dilution_solve_v1():
    res = ConcentrationCalculator.calculate_dilution(c1=5.0, c2=0.5, v2=250.0)
    assert res["success"] is True
    assert res["solved_for"] == "v1"
    assert res["v1"] == pytest.approx(25.0, abs=0.01)
    assert res["solvent_to_add"] == pytest.approx(225.0, abs=0.01)


def test_calculate_dilution_invalid_concentrations():
    res = ConcentrationCalculator.calculate_dilution(c1=0.2, c2=0.5, v2=250.0)
    assert res["success"] is False
    assert "greater" in res["error"]


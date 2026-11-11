from reaction_engine.equation_balancer import (
    balance_reaction_algebraically,
    verify_atom_conservation,
)
from reaction_engine.reaction_service import analyze_chemical_reaction


def test_algebraic_equation_balancer_and_atom_conservation():
    r_coeffs, p_coeffs = balance_reaction_algebraically(["H2", "O2"], ["H2O"])
    assert r_coeffs == [2, 1]
    assert p_coeffs == [2]

    r_fe, p_fe = balance_reaction_algebraically(["Fe", "O2"], ["Fe2O3"])
    assert r_fe == [4, 3]
    assert p_fe == [2]

    cons = verify_atom_conservation(
        [{"formula": "H2", "coefficient": 2}, {"formula": "O2", "coefficient": 1}],
        [{"formula": "H2O", "coefficient": 2}],
    )
    assert cons["is_balanced"] is True
    assert cons["reactant_atoms"] == {"H": 4, "O": 2}
    assert cons["product_atoms"] == {"H": 4, "O": 2}


def test_water_synthesis_exact_stoichiometric_ratio():
    res = analyze_chemical_reaction(
        [
            {"formula": "H2", "amount": 2.0, "unit": "mol"},
            {"formula": "O2", "amount": 1.0, "unit": "mol"},
        ]
    )
    assert res["supported"] is True
    assert res["equation_text"] == "2H2 + O2 -> 2H2O"
    assert res["equation_unicode"] == "2H₂ + O₂ → 2H₂O"
    assert res["is_exact_stoichiometric_ratio"] is True
    assert res["limiting_reagent"] is None
    assert res["theoretical_yields"][0]["formula"] == "H2O"
    assert res["theoretical_yields"][0]["theoretical_yield_mol"] == 2.0
    assert abs(res["theoretical_yields"][0]["theoretical_yield_g"] - 36.03) < 0.05
    for row in res["reactants_table"]:
        assert row["remaining_mol"] == 0.0


def test_limiting_and_excess_reagent_calculation():
    res = analyze_chemical_reaction(
        [
            {"formula": "H2", "amount": 5.0, "unit": "mol"},
            {"formula": "O2", "amount": 1.0, "unit": "mol"},
        ]
    )
    assert res["supported"] is True
    assert res["is_exact_stoichiometric_ratio"] is False
    assert res["limiting_reagent"]["formula"] == "O2"
    assert len(res["excess_reagents"]) == 1
    assert res["excess_reagents"][0]["formula"] == "H2"
    assert res["theoretical_yields"][0]["theoretical_yield_mol"] == 2.0

    h2_row = next(r for r in res["reactants_table"] if r["formula"] == "H2")
    assert h2_row["consumed_mol"] == 2.0
    assert h2_row["remaining_mol"] == 3.0


def test_unsupported_and_inert_reactions_refuse_to_guess():
    res_inert = analyze_chemical_reaction(
        [
            {"formula": "Fe", "amount": 1.0, "unit": "mol"},
            {"formula": "He", "amount": 1.0, "unit": "mol"},
        ]
    )
    assert res_inert["supported"] is False
    assert res_inert["error_type"] == "chemically_inert"

    res_unsupported = analyze_chemical_reaction(
        [
            {"formula": "NaCl", "amount": 1.0, "unit": "mol"},
            {"formula": "O2", "amount": 1.0, "unit": "mol"},
        ]
    )
    assert res_unsupported["supported"] is False
    assert res_unsupported["error_type"] == "unsupported_reaction"
    assert "ترفض تخمين" in res_unsupported["rejection_reason_ar"]


def test_domain_filtering_and_condition_compatibility():
    # HCl + NaOH is acid_base_chemistry / general_chemistry; should fail if filtered by organic_chemistry
    res_org_mismatch = analyze_chemical_reaction(
        [
            {"formula": "HCl", "amount": 1.0, "unit": "mol"},
            {"formula": "NaOH", "amount": 1.0, "unit": "mol"},
        ],
        domain="organic_chemistry",
    )
    assert res_org_mismatch["supported"] is False
    assert res_org_mismatch["error_type"] == "reaction_domain_mismatch"

    # Esterification CH3COOH + C2H5OH succeeds under organic_chemistry
    res_org_ok = analyze_chemical_reaction(
        [
            {"formula": "CH3COOH", "amount": 1.0, "unit": "mol"},
            {"formula": "C2H5OH", "amount": 1.0, "unit": "mol"},
        ],
        domain="organic_chemistry",
        custom_conditions={"catalyst": "بدون عامل حفاز إضافي (No Catalyst)", "temperature": "25 °C"},
    )
    assert res_org_ok["supported"] is True
    assert res_org_ok["conditions_compatibility"]["is_compatible"] is False
    assert len(res_org_ok["conditions_compatibility"]["warnings_ar"]) >= 1


def test_stp_volume_rejected_for_solids_and_liquids_allowed_for_gases():
    res_bad_stp = analyze_chemical_reaction(
        [
            {"formula": "NaCl", "amount": 22.414, "unit": "L_STP"},
            {"formula": "AgNO3", "amount": 1.0, "unit": "mol"},
        ]
    )
    assert res_bad_stp["supported"] is False
    assert res_bad_stp["error_type"] == "validation_error"
    assert "L_STP" in res_bad_stp["rejection_reason_ar"]

    res_gas_stp = analyze_chemical_reaction(
        [
            {"formula": "H2", "amount": 44.828, "unit": "L_STP"},
            {"formula": "O2", "amount": 22.414, "unit": "L_STP"},
        ]
    )
    assert res_gas_stp["supported"] is True
    assert abs(res_gas_stp["theoretical_yields"][0]["theoretical_yield_mol"] - 2.0) < 1e-3

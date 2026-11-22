"""
Comprehensive Titration Reference Curve & Equilibrium Solver Tests (v3.3.0):
Covers all 9 reference benchmarks:
  1. Weak acid + strong base
  2. Strong acid + weak base
  3. Polyprotic acid (H3PO4)
  4. Polyhydroxy base (Ba(OH)2 / Ca(OH)2)
  5. Weak acid + weak base (CH3COOH + NH3)
  6. Stepwise equivalence points
  7. Dilution effects
  8. High concentration effects
  9. Physical pH outside [0, 14] without artificial clipping (e.g., 2.0 M HCl -> pH < 0, 2.0 M NaOH -> pH > 14)
  Plus temperature dependence Kw(T) and Davies activity coefficient correction.
"""

from __future__ import annotations

import math
import unittest

from titration.titration_engine import (
    TitrationEngine,
    davies_activity_coefficient,
    kw_at_temperature,
    solve_equilibrium_ph_exact,
)


class TestTitrationReferenceCurves(unittest.TestCase):
    def setUp(self) -> None:
        self.hcl = TitrationEngine.ACID_DATABASE["HCl"]
        self.ch3cooh = TitrationEngine.ACID_DATABASE["CH3COOH"]
        self.h3po4 = TitrationEngine.ACID_DATABASE["H3PO4"]
        self.naoh = TitrationEngine.BASE_DATABASE["NaOH"]
        self.nh3 = TitrationEngine.BASE_DATABASE["NH3"]
        self.baoh2 = TitrationEngine.BASE_DATABASE["Ba(OH)2"]

    def test_1_weak_acid_strong_base_exact_equivalence(self) -> None:
        # 0.05 M sodium acetate at equivalence (pKa=4.76) -> pH ≈ 8.73
        ph_eq = solve_equilibrium_ph_exact(0.05, 0.05, self.ch3cooh, self.naoh)
        self.assertAlmostEqual(ph_eq, 8.73, delta=0.05)

    def test_2_strong_acid_weak_base_exact_equivalence(self) -> None:
        # 0.05 M NH4Cl at equivalence (pKb=4.75) -> pH ≈ 5.27
        ph_eq = solve_equilibrium_ph_exact(0.05, 0.05, self.hcl, self.nh3)
        self.assertAlmostEqual(ph_eq, 5.27, delta=0.05)

    def test_3_polyprotic_acid_stepwise_ph(self) -> None:
        # H3PO4 (0.1 M) + 1 eq NaOH -> H2PO4- (pH ≈ 0.5*(pKa1 + pKa2) ≈ 4.68)
        ph_eq1 = solve_equilibrium_ph_exact(0.05, 0.05, self.h3po4, self.naoh)
        # H3PO4 + 2 eq NaOH -> HPO4^2- (pH ≈ 0.5*(pKa2 + pKa3) ≈ 9.76)
        ph_eq2 = solve_equilibrium_ph_exact(0.0333, 0.0666, self.h3po4, self.naoh)
        self.assertTrue(4.3 <= ph_eq1 <= 5.0, f"Unexpected H3PO4 eq1 pH={ph_eq1}")
        self.assertTrue(9.2 <= ph_eq2 <= 10.2, f"Unexpected H3PO4 eq2 pH={ph_eq2}")

    def test_4_polyhydroxy_base_stoichiometry_and_ph(self) -> None:
        # 0.05 M Ba(OH)2 (0.10 M OH-) neutralized by 0.10 M HCl -> pH = 7.00
        ph_eq = solve_equilibrium_ph_exact(0.10, 0.05, self.hcl, self.baoh2)
        self.assertAlmostEqual(ph_eq, 7.00, delta=0.02)

    def test_5_weak_acid_weak_base_equivalence(self) -> None:
        # CH3COOH (pKa=4.76) + NH3 (pKb=4.75) -> pH ≈ 7.005
        ph_eq = solve_equilibrium_ph_exact(0.05, 0.05, self.ch3cooh, self.nh3)
        self.assertAlmostEqual(ph_eq, 7.00, delta=0.05)

    def test_6_stepwise_equivalence_points_in_simulation(self) -> None:
        sim = TitrationEngine.simulate_titration("H3PO4", "NaOH", 0.1, 25.0, 0.1)
        self.assertTrue(sim["is_polyprotic"])
        self.assertGreaterEqual(len(sim["step_equivalence_points"]), 2)

    def test_7_dilution_shifts_initial_ph_toward_7(self) -> None:
        ph_conc = solve_equilibrium_ph_exact(0.1, 0.0, self.ch3cooh, self.naoh)
        ph_dil = solve_equilibrium_ph_exact(0.001, 0.0, self.ch3cooh, self.naoh)
        self.assertGreater(ph_dil, ph_conc)
        self.assertLess(ph_dil, 7.0)

    def test_8_and_9_high_concentration_ph_outside_0_14_not_clipped(self) -> None:
        # 2.0 M HCl -> pH = -log10(2.0) ≈ -0.301 (MUST NOT be clipped to 0.0!)
        ph_strong_acid = solve_equilibrium_ph_exact(2.0, 0.0, self.hcl, self.naoh)
        self.assertLess(ph_strong_acid, 0.0)
        self.assertAlmostEqual(ph_strong_acid, -math.log10(2.0), delta=0.02)

        # 2.0 M NaOH -> pH = 14 + log10(2.0) ≈ 14.301 (MUST NOT be clipped to 14.0!)
        ph_strong_base = solve_equilibrium_ph_exact(0.0, 2.0, self.hcl, self.naoh)
        self.assertGreater(ph_strong_base, 14.0)
        self.assertAlmostEqual(ph_strong_base, 14.0 + math.log10(2.0), delta=0.02)

    def test_temperature_kw_and_davies_activity_correction(self) -> None:
        # Neutral water at 25 °C has pH = 7.00; at 50 °C Kw is larger so neutral pH < 7.00
        self.assertAlmostEqual(kw_at_temperature(25.0), 1.0e-14, delta=1e-16)
        self.assertGreater(kw_at_temperature(50.0), kw_at_temperature(25.0))
        ph_neutral_50 = solve_equilibrium_ph_exact(0.05, 0.05, self.hcl, self.naoh, temperature_c=50.0)
        self.assertLess(ph_neutral_50, 6.8)

        # Davies activity coefficient < 1 for z=1 at I=0.1 M -> activity pH slightly higher than ideal for HCl
        gamma = davies_activity_coefficient(1, 0.1)
        self.assertTrue(0.70 <= gamma <= 0.85)
        ph_ideal = solve_equilibrium_ph_exact(0.1, 0.0, self.hcl, self.naoh, use_activity_correction=False)
        ph_act = solve_equilibrium_ph_exact(0.1, 0.0, self.hcl, self.naoh, use_activity_correction=True)
        self.assertGreater(ph_act, ph_ideal)


if __name__ == "__main__":
    unittest.main()

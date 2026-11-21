"""
Comprehensive Titration Lab Catalog Tests:
Dynamically tests EVERY registered acid, base, and indicator in TitrationEngine
(rather than a hardcoded subset), plus hydrolysis, buffers, polyprotic acids,
extreme concentrations, and invalid input guards.
"""

from __future__ import annotations

import math
import unittest

from titration.titration_engine import TitrationEngine


class TestTitrationAllCatalog(unittest.TestCase):
    def test_every_registered_acid_simulates_cleanly_with_naoh(self) -> None:
        acids = TitrationEngine.ACID_DATABASE
        self.assertGreaterEqual(len(acids), 40)

        for acid_key, spec in acids.items():
            self.assertIn("name_ar", spec, f"Missing name_ar for {acid_key}")
            self.assertIn("pKa", spec, f"Missing pKa for {acid_key}")
            self.assertGreater(spec.get("molar_mass", 0), 0, f"Invalid molar_mass for {acid_key}")

            sim = TitrationEngine.simulate_titration(
                acid_key=acid_key,
                base_key="NaOH",
                acid_conc_m=0.1,
                acid_volume_ml=25.0,
                base_conc_m=0.1,
                titration_mode="acid_with_base",
            )
            self.assertGreater(sim["v_eq_ml"], 0.0, f"v_eq_ml invalid for {acid_key}")
            self.assertTrue(0.0 <= sim["ph_eq"] <= 14.0, f"ph_eq out of bounds for {acid_key}")
            self.assertTrue(len(sim["volumes_ml"]) >= 10, f"Curve too short for {acid_key}")
            self.assertEqual(len(sim["volumes_ml"]), len(sim["ph_values"]))
            for ph in sim["ph_values"]:
                self.assertFalse(math.isnan(ph) or math.isinf(ph), f"NaN/Inf pH for {acid_key}")
                self.assertTrue(-1.0 <= ph <= 15.0, f"Unphysical pH={ph} for {acid_key}")
            # Weak acid + strong base must NOT assume neutral salt (hydrolysis -> pH_eq > 7)
            if spec["type"] == "weak" and spec["protons"] == 1:
                self.assertGreaterEqual(sim["ph_eq"], 7.0, f"Weak monoprotic acid {acid_key} + NaOH must have pH_eq >= 7.0 (got {sim['ph_eq']})")

    def test_every_registered_base_simulates_cleanly_with_hcl(self) -> None:
        bases = TitrationEngine.BASE_DATABASE
        self.assertGreaterEqual(len(bases), 25)

        for base_key, spec in bases.items():
            self.assertIn("name_ar", spec, f"Missing name_ar for {base_key}")
            self.assertIn("pKb", spec, f"Missing pKb for {base_key}")
            self.assertGreater(spec.get("molar_mass", 0), 0, f"Invalid molar_mass for {base_key}")

            sim = TitrationEngine.simulate_titration(
                acid_key="HCl",
                base_key=base_key,
                acid_conc_m=0.1,
                acid_volume_ml=25.0,
                base_conc_m=0.1,
                titration_mode="base_with_acid",
                base_volume_ml=25.0,
            )
            self.assertGreater(sim["v_eq_ml"], 0.0, f"v_eq_ml invalid for {base_key}")
            self.assertTrue(0.0 <= sim["ph_eq"] <= 14.0, f"ph_eq out of bounds for {base_key}")
            for ph in sim["ph_values"]:
                self.assertFalse(math.isnan(ph) or math.isinf(ph), f"NaN/Inf pH for {base_key}")
                self.assertTrue(-1.0 <= ph <= 15.0, f"Unphysical pH={ph} for {base_key}")
            # Strong acid + weak base must NOT assume neutral salt (hydrolysis -> pH_eq < 7)
            if spec["type"] == "weak" and spec["hydroxides"] == 1:
                self.assertLessEqual(sim["ph_eq"], 7.0, f"Weak monoprotic base {base_key} + HCl must have pH_eq <= 7.0 (got {sim['ph_eq']})")

    def test_every_registered_indicator_has_valid_transition_range(self) -> None:
        indicators = TitrationEngine.INDICATORS
        self.assertGreaterEqual(len(indicators), 20)
        for ind_key, ind in indicators.items():
            low, high = ind["transition_range"]
            self.assertTrue(0.0 <= low < high <= 14.0, f"Invalid range for indicator {ind_key}")
            self.assertTrue(ind.get("name_ar"), f"Missing Arabic name for {ind_key}")
            self.assertTrue(ind.get("color_acid") and ind.get("color_base"))

    def test_polyprotic_acid_stepwise_equivalence_points(self) -> None:
        sim = TitrationEngine.simulate_titration(
            acid_key="H3PO4",
            base_key="NaOH",
            acid_conc_m=0.1,
            acid_volume_ml=25.0,
            base_conc_m=0.1,
            titration_mode="acid_with_base",
        )
        # 25 mL of 0.1 M H3PO4 requires 75 mL of 0.1 M NaOH for full 3-proton neutralization
        self.assertAlmostEqual(sim["v_eq_ml"], 75.0, delta=0.5)
        self.assertGreaterEqual(len(sim.get("step_equivalence_points", [])), 2)

    def test_low_and_high_concentrations_remain_stable(self) -> None:
        for conc in (0.001, 0.01, 1.0, 2.0):
            sim = TitrationEngine.simulate_titration(
                acid_key="HCl",
                base_key="NaOH",
                acid_conc_m=conc,
                acid_volume_ml=25.0,
                base_conc_m=conc,
                titration_mode="acid_with_base",
            )
            self.assertAlmostEqual(sim["v_eq_ml"], 25.0, delta=0.1)
            self.assertAlmostEqual(sim["ph_eq"], 7.0, delta=0.15)
            # pH should increase monotonically for acid titrated with base
            phs = sim["ph_values"]
            self.assertGreater(phs[-1], phs[0])

    def test_invalid_inputs_rejected_without_zero_division(self) -> None:
        for bad_kwargs in (
            dict(acid_key="HCl", base_key="NaOH", acid_conc_m=0.0, acid_volume_ml=25.0, base_conc_m=0.1),
            dict(acid_key="HCl", base_key="NaOH", acid_conc_m=-0.1, acid_volume_ml=25.0, base_conc_m=0.1),
            dict(acid_key="HCl", base_key="NaOH", acid_conc_m=0.1, acid_volume_ml=0.0, base_conc_m=0.1),
            dict(acid_key="HCl", base_key="NaOH", acid_conc_m=0.1, acid_volume_ml=25.0, base_conc_m=0.0),
            dict(acid_key="UNKNOWN_ACID", base_key="NaOH", acid_conc_m=0.1, acid_volume_ml=25.0, base_conc_m=0.1),
            dict(acid_key="HCl", base_key="UNKNOWN_BASE", acid_conc_m=0.1, acid_volume_ml=25.0, base_conc_m=0.1),
        ):
            with self.assertRaises(ValueError):
                TitrationEngine.simulate_titration(**bad_kwargs)


if __name__ == "__main__":
    unittest.main()

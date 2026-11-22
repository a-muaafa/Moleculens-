"""
MolecuLens v3.3.0 Scientific Hardening & Integration Test Suite:
Covers:
  1. RDKit Fallback strict-mode blocking & schematic warning in SVG drawings
  2. Lipinski missing-values -> unavailable (never pass) & decoupling from Drug-likeness
  3. 3-Layer Grounded Scientific Assistant actually using compound_data & no overstated claims
  4. QSAR Hybrid Scaffold Split, Hybrid LOCO, Multi-Seed, External Test Set, Learning Curves,
     95% PI Empirical Coverage, Per-Compound Provenance, and Deterministic Reproducibility
  5. Reaction Network (competing pathways, shared reactants, unreacted species, no double consumption)
  6. Domain filtering actually filtering documented & rule-inferred reactions + explicit no-reaction priority
  7. Scientific Audit separation of citation_exists vs source_manually_verified & theoretical model guard
  8. Auditable Repository provenance separation (curated_local vs external_api)
  9. Real HTTP REST API integration tests via FastAPI TestClient
"""

from __future__ import annotations

import unittest

from fastapi.testclient import TestClient

from api.server import app as fastapi_app
from chemistry_engine.molecular_intelligence import MolecularIntelligence
from chemistry_engine.structure_analyzer import (
    USING_RDKIT_FALLBACK,
    assert_native_rdkit_if_strict,
    render_molecule_2d_svg,
)
from data_integrations.repository import ScientificRepository, attach_record_provenance
from education.ai_assistant import AIScientificAssistant
from reaction_engine.prediction import predict_reaction
from reaction_engine.scientific_audit import audit_registry_and_rules


class TestScientificHardeningV33(unittest.TestCase):
    def test_1_fallback_strict_guard_and_schematic_svg_warning(self) -> None:
        if USING_RDKIT_FALLBACK:
            with self.assertRaises(RuntimeError):
                assert_native_rdkit_if_strict("qsar_scientific_evaluation", strict=True)
            svg = render_molecule_2d_svg("CCO")
            self.assertIsNotNone(svg)
            self.assertIn("Schematic Only", svg)
        else:
            # Native RDKit allows strict scientific operations
            assert_native_rdkit_if_strict("qsar_scientific_evaluation", strict=True)

    def test_2_lipinski_missing_values_never_pass_and_decoupled_from_druglikeness(self) -> None:
        # Missing logp, hbd, hba MUST return unavailable, NEVER pass!
        partial = AIScientificAssistant.explain_lipinski_evaluation({"molecular_weight": 180.0})
        self.assertFalse(partial["available"])
        self.assertEqual(partial["status"], "unavailable")
        self.assertIsNone(partial["is_compliant"])
        self.assertIn("LogP", partial["missing_parameters"])
        self.assertTrue(any(d["status"] == "unavailable" for d in partial["details"]))

        # Complete valid descriptors -> evaluated, but is_druglike explicitly decoupled (None)
        full = AIScientificAssistant.explain_lipinski_evaluation(
            {"molecular_weight": 180.16, "logp": 1.19, "hbd": 1, "hba": 4, "tpsa": 63.6, "rotatable_bonds": 3}
        )
        self.assertTrue(full["available"])
        self.assertTrue(full["lipinski_compliant"])
        self.assertIsNone(full["is_druglike"])
        self.assertIn("مرشّح إرشادي", full["drug_likeness_note_ar"])

    def test_3_three_layer_grounded_assistant_uses_compound_data(self) -> None:
        comp_data = {
            "name_ar": "الأسبرين",
            "name_en": "Aspirin",
            "symbol_or_formula": "C9H8O4",
            "sources": ["PubChem", "ChEMBL"],
            "properties": {"molecular_weight": 180.16, "logp": 1.19, "tpsa": 63.6, "hbd": 1, "hba": 4},
        }
        res = AIScientificAssistant.answer_student_question_structured("ما قيمة LogP لهذا المركب؟", compound_data=comp_data)
        self.assertEqual(res["assistant_architecture"], "rule_based_retrieval_grounded")
        self.assertTrue(res["grounded_in_compound_data"])
        self.assertIn("1.19", res["answer_text"])
        self.assertIn("الأسبرين", res["answer_text"])
        self.assertTrue(len(res["layer2_sources_used"]) >= 1)
        # Must NOT claim LogP 0-5 is universally ideal
        self.assertNotIn("يُعد النطاق بين 0 و 5 مثالياً", res["answer_text"])

        # Unknown clinical question -> honest unavailability
        unk = AIScientificAssistant.answer_student_question_structured("ما الجرعة السريرية لعلاج الصداع؟", compound_data=comp_data)
        self.assertIn("لا تتوفر", unk["answer_text"])

    def test_4_qsar_generalization_provenance_and_reproducibility(self) -> None:
        prov = MolecularIntelligence.get_logp_dataset_with_provenance()
        self.assertEqual(len(prov), 120)
        for row in prov:
            self.assertEqual(row["measurement_type"], "LogP")
            self.assertEqual(row["temperature_c"], 25.0)
            self.assertTrue(row["source_reference"])
            self.assertTrue(row["ionization_state"])

        # Determinism check: clear cache and train twice -> identical metrics
        MolecularIntelligence._CACHED_QSAR_BUNDLE = None
        b1 = MolecularIntelligence.train_and_evaluate_qsar_logp_model()
        MolecularIntelligence._CACHED_QSAR_BUNDLE = None
        b2 = MolecularIntelligence.train_and_evaluate_qsar_logp_model()

        self.assertEqual(b1["test_rmse"], b2["test_rmse"])
        self.assertEqual(b1["hybrid_residual_qsar"]["hybrid_scaffold_rmse"], b2["hybrid_residual_qsar"]["hybrid_scaffold_rmse"])
        self.assertEqual(b1["hybrid_residual_qsar"]["hybrid_loco_rmse"], b2["hybrid_residual_qsar"]["hybrid_loco_rmse"])
        self.assertEqual(b1["external_test_evaluation"]["n_external"], 15)
        self.assertEqual(len(b1["multi_seed_evaluation"]["seeds"]), 5)
        self.assertEqual(len(b1["learning_curves"]), 4)
        self.assertIn("empirical_coverage_pct", b1["prediction_interval_calibration"])

    def test_5_reaction_network_competing_pathways_and_unreacted_species(self) -> None:
        # HCl + NaOH + KOH + NaCl: HCl can react with BOTH NaOH and KOH (competing pathways!), while NaCl is unreacted!
        items = [
            {"substance": {"symbol_or_formula": f, "name_en": f, "name_ar": f}}
            for f in ("HCl", "NaOH", "KOH", "NaCl")
        ]
        res = predict_reaction(items, domain="all_domains")
        self.assertIn("reaction_network", res)
        net = res["reaction_network"]
        self.assertGreaterEqual(net["pathway_count"], 2)
        self.assertIn("HCl", net["shared_reactants"])
        self.assertIn("NaCl", net["unreacted_formulas"])
        self.assertTrue(net["has_unresolved_competition"])
        self.assertTrue(net["double_consumption_prevented"])
        self.assertTrue(res.get("single_outcome_withheld"))

    def test_6_domain_filtering_changes_results_and_no_reaction_priority(self) -> None:
        # HCl + NaOH is acid-base; filtering to "organic_chemistry" MUST NOT return it!
        acid_base = [
            {"substance": {"symbol_or_formula": "HCl", "name_en": "HCl", "name_ar": "HCl"}},
            {"substance": {"symbol_or_formula": "NaOH", "name_en": "NaOH", "name_ar": "NaOH"}},
        ]
        res_all = predict_reaction(acid_base, domain="all_domains")
        res_org = predict_reaction(acid_base, domain="organic_chemistry")
        self.assertTrue(res_all["supported"])
        self.assertFalse(res_org["supported"])

        # Explicit no-reaction: Cu + HCl must be rejected BEFORE general acid+metal rule
        cu_hcl = [
            {"substance": {"symbol_or_formula": "Cu", "name_en": "Copper", "name_ar": "نحاس"}},
            {"substance": {"symbol_or_formula": "HCl", "name_en": "HCl", "name_ar": "حمض الهيدروكلوريك"}},
        ]
        res_cu = predict_reaction(cu_hcl)
        self.assertFalse(res_cu["supported"])
        self.assertTrue(res_cu.get("is_no_reaction") or res_cu.get("mode") in {"no_reaction", "insufficient_evidence"})

    def test_7_scientific_audit_separates_citation_exists_from_manually_verified(self) -> None:
        audit = audit_registry_and_rules()
        for rep in audit["reaction_reports"]:
            self.assertIn("citation_exists", rep)
            self.assertIn("source_manually_verified", rep)
            self.assertIn("is_theoretical_model", rep)
            # Theoretical model must NEVER be reported as experimentally_documented
            if rep["is_theoretical_model"]:
                self.assertNotEqual(rep["documentation_type"], "experimentally_documented")

    def test_8_repository_provenance_fields_separate_curated_vs_external(self) -> None:
        curated = ScientificRepository.get_all_curated_substances()
        self.assertGreaterEqual(len(curated), 240)
        sample = curated[0]
        for field in ("source_id", "source_url", "source_version", "retrieved_at", "verification_status", "license", "data_quality_status", "record_origin"):
            self.assertIn(field, sample)
        self.assertEqual(sample["record_origin"], "curated_local")

        ext_rec = attach_record_provenance({"key": "pubchem_2244", "pubchem_cid": 2244, "db_source": "External: PubChem"}, record_origin="external_api")
        self.assertEqual(ext_rec["record_origin"], "external_api")
        self.assertIsNone(ext_rec["verified_at"])

    def test_9_real_fastapi_http_rest_api_integration(self) -> None:
        client = TestClient(fastapi_app)

        r_health = client.get("/api/v1/health")
        self.assertEqual(r_health.status_code, 200)
        self.assertTrue(r_health.json()["ok"])

        r_exp = client.get("/api/v1/explore", params={"q": "CH3COOH", "allow_network": "false"})
        self.assertEqual(r_exp.status_code, 200)
        self.assertTrue(r_exp.json()["ok"])

        r_rxn = client.post(
            "/api/v1/reactions/predict",
            json={"reactants": [{"formula": "HCl"}, {"formula": "NaOH"}, {"formula": "KOH"}], "domain": "all_domains"},
        )
        self.assertEqual(r_rxn.status_code, 200)
        self.assertIn("reaction_network", r_rxn.json()["data"])

        r_tit = client.post(
            "/api/v1/titration/simulate",
            json={"acid_key": "HCl", "base_key": "NaOH", "acid_conc_m": 2.0, "acid_volume_ml": 25.0, "base_conc_m": 0.1},
        )
        self.assertEqual(r_tit.status_code, 200)
        # 2.0 M HCl exact initial pH < 0 (not clipped!)
        self.assertLess(r_tit.json()["exact_initial_ph"], 0.0)

        r_bad = client.post("/api/v1/molecular/logp", json={"smiles": "INVALID_SMILES((("})
        self.assertIn(r_bad.status_code, (400, 422))

        r_audit = client.get("/api/v1/audit/scientific")
        self.assertEqual(r_audit.status_code, 200)
        self.assertTrue(r_audit.json()["ok"])


if __name__ == "__main__":
    unittest.main()

"""
Phase 1 Acceptance & Regression Tests:
- Offline-First fallback (3s timeout, SSL/Timeout/ConnectionError handling without Tracebacks)
- Offline Presets for Reaction Lab, Prediction Engine, and Titration Lab
- Streamlit cache wrappers & explicit cache invalidation
- Incomplete UI hiding (Pharmacology hidden by default) & exact scientific disclaimer
"""

from __future__ import annotations

import unittest
from unittest.mock import patch
import requests

from config.settings import (
    EXTERNAL_API_TIMEOUT_SECONDS,
    SCIENTIFIC_DISCLAIMER_AR,
    SHOW_EXPERIMENTAL_PHARMACOLOGY,
)
from data_integrations.chembl_client import fetch_chembl_drug_info, search_chembl_molecules
from data_integrations.data_manager import DataManager
from data_integrations.http_client import classify_network_error_ar, request_with_retry
from data_integrations.pubchem_client import fetch_pubchem_by_formula, fetch_pubchem_by_name
from frontend.cache_helpers import (
    get_cache_version,
    get_cached_logp_evaluation,
    get_cached_qsar_bundle,
    get_cached_reaction_audit,
    get_cached_similarity_results,
    get_cached_titration_result,
    invalidate_all_caches,
)
from reaction_engine.prediction import predict_reaction
from reaction_engine.reaction_service import analyze_chemical_reaction
from search.search_engine import search_scientific_entity
from titration.titration_engine import TitrationEngine


class TestOfflineFirstAndCaching(unittest.TestCase):
    def test_default_timeout_is_short_and_pharmacology_hidden(self) -> None:
        self.assertLessEqual(EXTERNAL_API_TIMEOUT_SECONDS, 3)
        self.assertFalse(SHOW_EXPERIMENTAL_PHARMACOLOGY)
        self.assertIn("هذه المنصة مخصصة للمحاكاة التعليمية والبحث الأولي فقط", SCIENTIFIC_DISCLAIMER_AR)
        self.assertIn("الاستشارة الطبية أو الصيدلانية المعتمدة", SCIENTIFIC_DISCLAIMER_AR)

    def test_ssl_and_timeout_error_classification_in_arabic(self) -> None:
        ssl_msg = classify_network_error_ar(requests.exceptions.SSLError("cert failed"), "PubChem")
        self.assertIn("SSL", ssl_msg)
        self.assertIn("البيانات المحلية", ssl_msg)

        timeout_msg = classify_network_error_ar(requests.exceptions.Timeout("timed out"), "ChEMBL")
        self.assertIn("انتهت مهلة الاتصال", timeout_msg)

        conn_msg = classify_network_error_ar(requests.exceptions.ConnectionError("offline"), "PubChem")
        self.assertIn("البيانات المحلية", conn_msg)

    def test_ssl_error_fails_fast_without_retries(self) -> None:
        call_count = 0

        def _fake_get(*args, **kwargs):
            nonlocal call_count
            call_count += 1
            raise requests.exceptions.SSLError("SSL handshake failed")

        with patch("requests.get", side_effect=_fake_get):
            with self.assertRaises(requests.exceptions.SSLError):
                request_with_retry("https://pubchem.ncbi.nlm.nih.gov/rest/pug/compound/name/water/JSON", max_retries=2)
        self.assertEqual(call_count, 1)

    def test_external_clients_handle_offline_errors_gracefully(self) -> None:
        for exc in (
            requests.exceptions.SSLError("SSL error"),
            requests.exceptions.Timeout("Read timeout"),
            requests.exceptions.ConnectionError("DNS failure"),
        ):
            with patch("requests.get", side_effect=exc):
                res_pc = fetch_pubchem_by_name("nonexistent_offline_probe_xyz", timeout=1)
                self.assertFalse(res_pc["ok"])
                self.assertEqual(res_pc["records"], [])

                res_pcf = fetch_pubchem_by_formula("C99H198O12", timeout=1)
                self.assertFalse(res_pcf["ok"])

                res_ch = search_chembl_molecules("nonexistent_offline_probe_xyz", timeout=1)
                self.assertFalse(res_ch["ok"])

                res_chd = fetch_chembl_drug_info("CHEMBL999999999", timeout=1)
                self.assertFalse(res_chd["ok"])

    def test_data_manager_offline_mode_does_not_claim_external_origin(self) -> None:
        res = DataManager.resolve_external_substances("totally_unknown_compound_999", allow_network=False)
        self.assertFalse(res["ok"])
        self.assertNotEqual(res["origin"], "external_sources")

    def test_five_offline_presets_work_completely_without_network(self) -> None:
        # 1. Local search presets (Arabic, English, Formula)
        for q in ("H2O", "CH3COOH", "حمض الخليك", "Aspirin", "NaOH"):
            sr = search_scientific_entity(q, allow_network=False)
            self.assertTrue(sr["found"], f"Offline search failed for {q}")

        # 2. Reaction prediction presets (Element+Element, Element+Compound, Compound+Compound, No-reaction)
        for reactants in (
            ["H2", "O2"],
            ["Mg", "HCl"],
            ["HCl", "NaOH"],
            ["CH3COOH", "K2CO3"],
            ["Cu", "HCl"],
        ):
            items = [{"substance": {"symbol_or_formula": r, "name_en": r, "name_ar": r}} for r in reactants]
            pr = predict_reaction(items)
            self.assertIn("mode", pr)
            self.assertTrue(pr.get("mode") in {"documented", "rule_inferred", "no_reaction", "insufficient_evidence"})

        # 3. Verified reaction stoichiometry preset
        rxn_res = analyze_chemical_reaction(
            [{"formula": "HCl", "amount": 1.0, "unit": "mol"}, {"formula": "NaOH", "amount": 1.0, "unit": "mol"}],
        )
        self.assertTrue(rxn_res.get("supported") or rxn_res.get("found") or rxn_res.get("matched"))

        # 4. Five offline titration presets
        titration_presets = [
            ("HCl", "NaOH", 25.0, 0.10, 0.10, "acid_with_base"),
            ("CH3COOH", "NaOH", 25.0, 0.10, 0.10, "acid_with_base"),
            ("H3PO4", "NaOH", 25.0, 0.10, 0.10, "acid_with_base"),
            ("HCl", "NH3", 25.0, 0.10, 0.10, "base_with_acid"),
            ("H2SO4", "KOH", 25.0, 0.10, 0.10, "acid_with_base"),
        ]
        for acid_k, base_k, vol, c1, c2, mode in titration_presets:
            sim = TitrationEngine.simulate_titration(
                acid_key=acid_k,
                base_key=base_k,
                acid_conc_m=c1 if mode == "acid_with_base" else c2,
                acid_volume_ml=vol,
                base_conc_m=c2 if mode == "acid_with_base" else c1,
                titration_mode=mode,
                base_volume_ml=vol if mode == "base_with_acid" else None,
            )
            self.assertGreater(sim["v_eq_ml"], 0.0)
            self.assertTrue(0.0 <= sim["ph_eq"] <= 14.0)

    def test_cache_helpers_and_invalidation(self) -> None:
        v1 = get_cache_version()
        bundle = get_cached_qsar_bundle(cache_version=v1)
        self.assertTrue(bundle.get("available"))

        sims = get_cached_similarity_results("CC(=O)Oc1ccccc1C(=O)O", top_n=3, cache_version=v1)
        self.assertIsInstance(sims, list)

        logp = get_cached_logp_evaluation("CCO", cache_version=v1)
        self.assertTrue(logp.get("available"))

        tit = get_cached_titration_result("HCl", "NaOH", 0.1, 25.0, 0.1, cache_version=v1)
        self.assertAlmostEqual(tit["ph_eq"], 7.0, delta=0.1)

        audit = get_cached_reaction_audit(cache_version=v1)
        self.assertGreaterEqual(audit.get("summary", {}).get("reactions_total", audit.get("total_reactions", 0)), 150)

        v2 = invalidate_all_caches()
        self.assertEqual(v2, v1 + 1)


if __name__ == "__main__":
    unittest.main()

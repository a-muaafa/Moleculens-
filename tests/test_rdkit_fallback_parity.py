"""
RDKit & Pure-Python Fallback Parity and Transparency Tests:
Verifies SMILES parsing/validation, canonicalization, molecular weight & formula,
Morgan Fingerprints & Tanimoto similarity, SMARTS substructure matching,
Bemis-Murcko scaffolds, and transparent fallback version reporting.
"""

from __future__ import annotations

import unittest

from chemistry_engine.molecular_intelligence import MolecularIntelligence
from chemistry_engine.structure_analyzer import (
    Chem,
    Descriptors,
    USING_RDKIT_FALLBACK,
    get_cheminformatics_engine_info,
    validate_and_parse_smiles,
)


class TestRDKitFallbackParity(unittest.TestCase):
    def test_engine_info_transparency_never_masquerades_as_native_rdkit(self) -> None:
        info = get_cheminformatics_engine_info()
        self.assertEqual(info["using_fallback"], USING_RDKIT_FALLBACK)
        self.assertTrue(info["engine_label_ar"])
        self.assertTrue(info["equivalence_note_ar"])
        if USING_RDKIT_FALLBACK:
            self.assertFalse(info["native_rdkit_available"])
            self.assertIn("fallback", info["engine_version"].lower())
            # Must not look like a native RDKit calendar release version (e.g. 2026.09.1)
            self.assertFalse(info["engine_version"].startswith("202"))

    def test_smiles_validation_accepts_valid_and_rejects_invalid_with_arabic_error(self) -> None:
        for valid_smi in ("O", "CCO", "CC(=O)O", "c1ccccc1", "CC(=O)Oc1ccccc1C(=O)O"):
            mol, canon, err = validate_and_parse_smiles(valid_smi)
            self.assertIsNotNone(mol, f"Failed to parse valid SMILES: {valid_smi}")
            self.assertTrue(canon)
            self.assertEqual(err, "")

        for bad_smi in ("", "   ", "C(C", "invalid_smiles_123!", "C1CC"):
            mol, canon, err = validate_and_parse_smiles(bad_smi)
            self.assertIsNone(mol, f"Should reject invalid SMILES: {bad_smi}")
            self.assertFalse(canon)
            self.assertTrue(err, f"Missing Arabic error message for invalid SMILES: {bad_smi}")

    def test_molecular_weight_parity_against_iupac_benchmarks(self) -> None:
        benchmarks = [
            ("O", 18.015),
            ("CCO", 46.069),
            ("CC(=O)O", 60.052),
            ("c1ccccc1", 78.114),
            ("CC(=O)Oc1ccccc1C(=O)O", 180.159),  # Aspirin
            ("Cn1cnc2c1c(=O)n(C)c(=O)n2C", 194.194),  # Caffeine
        ]
        for smi, expected_mw in benchmarks:
            mol = Chem.MolFromSmiles(smi)
            self.assertIsNotNone(mol)
            mw = Descriptors.MolWt(mol)
            self.assertAlmostEqual(mw, expected_mw, delta=0.25, msg=f"MW mismatch for {smi}")

    def test_morgan_fingerprint_and_tanimoto_properties(self) -> None:
        aspirin = "CC(=O)Oc1ccccc1C(=O)O"
        salicylic = "O=C(O)c1ccccc1O"
        methane = "C"

        sim_self = MolecularIntelligence.calculate_similarity(aspirin, aspirin)
        self.assertAlmostEqual(sim_self, 1.0, places=4)

        sim_ab = MolecularIntelligence.calculate_similarity(aspirin, salicylic)
        sim_ba = MolecularIntelligence.calculate_similarity(salicylic, aspirin)
        self.assertAlmostEqual(sim_ab, sim_ba, places=6)

        sim_unrelated = MolecularIntelligence.calculate_similarity(aspirin, methane)
        self.assertGreater(sim_ab, sim_unrelated)

    def test_smarts_substructure_matching_parity(self) -> None:
        acid_mol = Chem.MolFromSmiles("CC(=O)O")
        ethanol_mol = Chem.MolFromSmiles("CCO")
        benzene_mol = Chem.MolFromSmiles("c1ccccc1")

        carboxyl = Chem.MolFromSmarts("[CX3](=O)[OX2H1]")
        hydroxyl = Chem.MolFromSmarts("[OX2H]")

        self.assertTrue(acid_mol.HasSubstructMatch(carboxyl))
        self.assertFalse(ethanol_mol.HasSubstructMatch(carboxyl))
        self.assertTrue(ethanol_mol.HasSubstructMatch(hydroxyl))
        self.assertFalse(benzene_mol.HasSubstructMatch(hydroxyl))

    def test_bemis_murcko_scaffold_extraction(self) -> None:
        # Aromatic ring compounds should yield a non-empty ring scaffold; acyclic should be empty or acyclic
        s_toluene = MolecularIntelligence._extract_murcko_scaffold(Chem.MolFromSmiles("Cc1ccccc1"))
        s_aspirin = MolecularIntelligence._extract_murcko_scaffold(Chem.MolFromSmiles("CC(=O)Oc1ccccc1C(=O)O"))
        self.assertTrue(s_toluene)
        self.assertTrue(s_aspirin)
        self.assertEqual(s_toluene, s_aspirin)


if __name__ == "__main__":
    unittest.main()

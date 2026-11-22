"""
150-Molecule Diverse Benchmark: Native C++ RDKit vs. Pure-Python RDKit Fallback.

Compares Native RDKit against chemistry_engine._rdkit_fallback across 8 cheminformatics
dimensions:
  1. Molecular Formula
  2. Molecular Weight (within ±0.1 Da)
  3. Canonical SMILES
  4. Morgan Fingerprints (radius=2, 2048 bits)
  5. Tanimoto Similarity (50 diverse pairs, within ±0.05)
  6. SMARTS Substructure Matching (8 functional groups x 150 molecules)
  7. Bemis-Murcko Scaffolds
  8. Molecular Descriptors (MolLogP, TPSA, HBD, HBA, RotatableBonds, RingCount)

Usage:
    python scripts/benchmark_rdkit_parity.py
    python scripts/benchmark_rdkit_parity.py --require-native-rdkit
"""

from __future__ import annotations

import argparse
import json
import sys
from datetime import datetime, timezone
from pathlib import Path
from typing import Any

ROOT = Path(__file__).resolve().parents[1]
sys.path.insert(0, str(ROOT))

# 150 diverse organic, inorganic-organic, heterocyclic, and pharmaceutical molecules
BENCHMARK_150_SMILES: list[tuple[str, str]] = [
    # Aliphatic hydrocarbons & alcohols (1-20)
    ("Methane", "C"), ("Ethane", "CC"), ("Propane", "CCC"), ("Butane", "CCCC"),
    ("Pentane", "CCCCC"), ("Hexane", "CCCCCC"), ("Isobutane", "CC(C)C"), ("Neopentane", "CC(C)(C)C"),
    ("Cyclopentane", "C1CCCC1"), ("Cyclohexane", "C1CCCCC1"), ("Methanol", "CO"), ("Ethanol", "CCO"),
    ("1-Propanol", "CCCO"), ("Isopropanol", "CC(C)O"), ("1-Butanol", "CCCCO"), ("tert-Butanol", "CC(C)(C)O"),
    ("Ethylene glycol", "OCCO"), ("Glycerol", "OCC(O)CO"), ("Cyclohexanol", "OC1CCCCC1"), ("Allyl alcohol", "C=CCO"),
    # Aldehydes, ketones, ethers, esters (21-40)
    ("Formaldehyde", "C=O"), ("Acetaldehyde", "CC=O"), ("Propionaldehyde", "CCC=O"), ("Acetone", "CC(=O)C"),
    ("Butanone", "CCC(=O)C"), ("Cyclohexanone", "O=C1CCCCC1"), ("Dimethyl ether", "COC"), ("Diethyl ether", "CCOCC"),
    ("THF", "C1CCOC1"), ("1,4-Dioxane", "C1COCCO1"), ("Methyl formate", "COC=O"), ("Methyl acetate", "CC(=O)OC"),
    ("Ethyl acetate", "CC(=O)OCC"), ("Butyl acetate", "CCCCOC(=O)C"), ("Methyl propionate", "CCC(=O)OC"),
    ("Acrolein", "C=CC=O"), ("Diisopropyl ether", "CC(C)OC(C)C"), ("MTBE", "COC(C)(C)C"),
    ("Acetylacetone", "CC(=O)CC(=O)C"), ("Gamma-butyrolactone", "O=C1CCCO1"),
    # Carboxylic acids & amines/amides (41-60)
    ("Formic acid", "O=CO"), ("Acetic acid", "CC(=O)O"), ("Propionic acid", "CCC(=O)O"), ("Butyric acid", "CCCC(=O)O"),
    ("Valeric acid", "CCCCC(=O)O"), ("Oxalic acid", "O=C(O)C(=O)O"), ("Malonic acid", "O=C(O)CC(=O)O"),
    ("Succinic acid", "O=C(O)CCC(=O)O"), ("Lactic acid", "CC(O)C(=O)O"), ("Citric acid", "O=C(O)CC(O)(C(=O)O)CC(=O)O"),
    ("Methylamine", "CN"), ("Ethylamine", "CCN"), ("Dimethylamine", "CNC"), ("Triethylamine", "CCN(CC)CC"),
    ("Ethylenediamine", "NCCN"), ("Urea", "NC(=O)N"), ("Acetamide", "CC(=O)N"), ("DMF", "CN(C)C=O"),
    ("Glycine", "NCC(=O)O"), ("Alanine", "CC(N)C(=O)O"),
    # Halogenated & sulfur/nitro/nitrile aliphatics (61-75)
    ("Chloromethane", "CCl"), ("Dichloromethane", "ClCCl"), ("Chloroform", "ClC(Cl)Cl"), ("Carbon tetrachloride", "ClC(Cl)(Cl)Cl"),
    ("Bromoethane", "CCBr"), ("Iodomethane", "CI"), ("Trifluoroacetic acid", "O=C(O)C(F)(F)F"),
    ("Acetonitrile", "CC#N"), ("Acrylonitrile", "C=CC#N"), ("Nitromethane", "C[N+](=O)[O-]"),
    ("DMSO", "CS(=O)C"), ("Methanethiol", "CS"), ("Ethanethiol", "CCS"), ("2-Fluoroethanol", "FCCO"), ("Chloroacetic acid", "O=C(O)CCl"),
    # Aromatics & substituted benzenes (76-105)
    ("Benzene", "c1ccccc1"), ("Toluene", "Cc1ccccc1"), ("Ethylbenzene", "CCc1ccccc1"), ("o-Xylene", "Cc1ccccc1C"),
    ("p-Xylene", "Cc1ccc(C)cc1"), ("Styrene", "C=Cc1ccccc1"), ("Phenol", "Oc1ccccc1"), ("Catechol", "Oc1ccccc1O"),
    ("Resorcinol", "Oc1cccc(O)c1"), ("Hydroquinone", "Oc1ccc(O)cc1"), ("Anisole", "COc1ccccc1"), ("Benzyl alcohol", "OCc1ccccc1"),
    ("Benzaldehyde", "O=Cc1ccccc1"), ("Acetophenone", "CC(=O)c1ccccc1"), ("Benzoic acid", "O=C(O)c1ccccc1"),
    ("Salicylic acid", "O=C(O)c1ccccc1O"), ("Aniline", "Nc1ccccc1"), ("Nitrobenzene", "O=[N+]([O-])c1ccccc1"),
    ("Chlorobenzene", "Clc1ccccc1"), ("Bromobenzene", "Brc1ccccc1"), ("Fluorobenzene", "Fc1ccccc1"),
    ("p-Chlorophenol", "Oc1ccc(Cl)cc1"), ("p-Nitrophenol", "O=[N+]([O-])c1ccc(O)cc1"), ("Benzonitrile", "N#Cc1ccccc1"),
    ("Phenylacetic acid", "O=C(O)Cc1ccccc1"), ("Cinnamic acid", "O=C(O)C=Cc1ccccc1"), ("Vanillin", "COc1cc(C=O)ccc1O"),
    ("Eugenol", "C=CCc1ccc(O)c(OC)c1"), ("Thymol", "Cc1ccc(C(C)C)c(O)c1"), ("Biphenyl", "c1ccc(-c2ccccc2)cc1"),
    # Polycyclics & Heterocycles (106-125)
    ("Naphthalene", "c1ccc2ccccc2c1"), ("Anthracene", "c1ccc2cc3ccccc3cc2c1"), ("Pyridine", "c1ccncc1"),
    ("Pyrimidine", "c1cncnc1"), ("Pyrazine", "c1cnccn1"), ("Pyrrole", "c1cc[nH]c1"), ("Furan", "c1ccoc1"),
    ("Thiophene", "c1ccsc1"), ("Imidazole", "c1c[nH]cn1"), ("Oxazole", "c1cocn1"), ("Thiazole", "c1cscn1"),
    ("Indole", "c1ccc2[nH]ccc2c1"), ("Quinoline", "c1ccc2ncccc2c1"), ("Isoquinoline", "c1ccc2cnccc2c1"),
    ("Benzimidazole", "c1ccc2[nH]cnc2c1"), ("Piperidine", "C1CCNCC1"), ("Morpholine", "C1COCCN1"),
    ("Piperazine", "C1CNCCN1"), ("Pyrrolidine", "C1CCNC1"), ("Nicotinic acid", "O=C(O)c1cccnc1"),
    # Drugs & Bioactive Compounds (126-150)
    ("Aspirin", "CC(=O)Oc1ccccc1C(=O)O"),
    ("Paracetamol", "CC(=O)Nc1ccc(O)cc1"),
    ("Ibuprofen", "CC(C)Cc1ccc(C(C)C(=O)O)cc1"),
    ("Naproxen", "COc1ccc2cc(C(C)C(=O)O)ccc2c1"),
    ("Caffeine", "Cn1cnc2c1c(=O)n(C)c(=O)n2C"),
    ("Theophylline", "Cn1c(=O)c2[nH]cnc2n(C)c1=O"),
    ("Theobromine", "Cn1cnc2c1c(=O)[nH]c(=O)n2C"),
    ("Metformin", "CN(C)C(=N)NC(=N)N"),
    ("Sulfanilamide", "Nc1ccc(S(N)(=O)=O)cc1"),
    ("Isoniazid", "NNC(=O)c1ccncc1"),
    ("Lidocaine", "CCN(CC)CC(=O)Nc1c(C)cccc1C"),
    ("Procaine", "CCN(CC)CCOC(=O)c1ccc(N)cc1"),
    ("Phenacetin", "CCOc1ccc(NC(C)=O)cc1"),
    ("Diclofenac", "O=C(O)Cc1ccccc1Nc1c(Cl)cccc1Cl"),
    ("Flurbiprofen", "CC(C(=O)O)c1ccc(-c2ccccc2)c(F)c1"),
    ("Ketoprofen", "CC(C(=O)O)c1cccc(C(=O)c2ccccc2)c1"),
    ("Indomethacin", "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c1ccc(Cl)cc1"),
    ("Dopamine", "NCCc1ccc(O)c(O)c1"),
    ("Epinephrine", "CNC[C@H](O)c1ccc(O)c(O)c1"),
    ("Serotonin", "NCCc1c[nH]c2ccc(O)cc12"),
    ("Histamine", "NCCc1c[nH]cn1"),
    ("Nicotine", "CN1CCC[C@H]1c1cccnc1"),
    ("Warfarin", "CC(=O)CC(c1ccccc1)c1c(O)c2ccccc2oc1=O"),
    ("Omeprazole", "COc1ccc2[nH]c(S(=O)Cc3ncc(C)c(OC)c3C)nc2c1"),
    ("Diazepam", "CN1C(=O)CN=C(c2ccccc2)c2cc(Cl)ccc21"),
]

SMARTS_PROBES: list[tuple[str, str]] = [
    ("carboxylic_acid", "[CX3](=O)[OX2H1]"),
    ("hydroxyl", "[OX2H]"),
    ("carbonyl", "[CX3]=[OX1]"),
    ("amine_nh", "[NX3;H2,H1]"),
    ("aromatic_carbon", "[c]"),
    ("ester", "[#6][CX3](=O)[OX2][#6]"),
    ("halogen", "[F,Cl,Br,I]"),
    ("nitro", "[$([NX3](=O)=O),$([NX3+](=O)[O-])]"),
]


def _import_native_rdkit() -> dict[str, Any] | None:
    """Attempts to import genuine native C++ RDKit (rejecting our pure-Python fallback)."""
    existing = sys.modules.get("rdkit")
    if existing is not None and getattr(existing, "__is_fallback__", False):
        return None
    try:
        import rdkit
        if getattr(rdkit, "__is_fallback__", False):
            return None
        from rdkit import Chem, DataStructs, RDLogger
        RDLogger.DisableLog("rdApp.warning")
        from rdkit.Chem import Crippen, Descriptors, Lipinski, rdMolDescriptors
        from rdkit.Chem.Scaffolds import MurckoScaffold
        return {
            "version": getattr(rdkit, "__version__", "unknown"),
            "Chem": Chem,
            "DataStructs": DataStructs,
            "Crippen": Crippen,
            "Descriptors": Descriptors,
            "Lipinski": Lipinski,
            "rdMolDescriptors": rdMolDescriptors,
            "MurckoScaffold": MurckoScaffold,
        }
    except ImportError:
        return None


def run_parity_benchmark(require_native_rdkit: bool = False) -> dict[str, Any]:
    """
    Runs the 150-molecule benchmark comparing Native C++ RDKit against Pure-Python Fallback.
    Never fabricates parity percentages when native RDKit is absent.
    """
    native = _import_native_rdkit()

    from chemistry_engine import _rdkit_fallback as fb

    total = len(BENCHMARK_150_SMILES)
    fb_parsed = 0
    for _, smi in BENCHMARK_150_SMILES:
        if fb.Chem.MolFromSmiles(smi) is not None:
            fb_parsed += 1

    if native is None:
        report: dict[str, Any] = {
            "generated_at_utc": datetime.now(timezone.utc).isoformat(),
            "benchmark_size": total,
            "native_rdkit_available": False,
            "native_rdkit_version": None,
            "fallback_version": "pure-python-fallback-1.0",
            "parity_verified": False,
            "fallback_parse_success_count": fb_parsed,
            "fallback_parse_success_pct": round(100.0 * fb_parsed / total, 2),
            "note_ar": (
                "لم يتم العثور على مكتبة RDKit الأصلية في هذه البيئة؛ لذا لم تُحسب نسب المطابقة الفعلية بين RDKit والبديل. "
                "شغّل هذا السكربت في بيئة تحتوي على RDKit الأصلي (أو عبر وظيفة native-rdkit في GitHub Actions) "
                "لقياس نسب المطابقة الحقيقية عبر الـ 150 جزيئاً."
            ),
            "property_match_rates_pct": None,
        }
        if require_native_rdkit:
            report["error"] = "Native C++ RDKit is required for --require-native-rdkit parity verification."
        return report

    n_Chem = native["Chem"]
    n_DS = native["DataStructs"]
    n_Desc = native["Descriptors"]
    n_Crip = native["Crippen"]
    n_Lip = native["Lipinski"]
    n_rdMD = native["rdMolDescriptors"]
    n_Scaf = native["MurckoScaffold"]

    counts = {
        "parsed_both": 0,
        "formula_exact": 0,
        "mw_within_0_1_da": 0,
        "canonical_smiles_exact": 0,
        "morgan_fp_exact": 0,
        "scaffold_exact": 0,
        "logp_within_0_5": 0,
        "tpsa_within_5_0": 0,
        "hbd_exact": 0,
        "hba_exact": 0,
        "rotatable_exact": 0,
        "ring_count_exact": 0,
    }
    smarts_total = 0
    smarts_match = 0
    native_fps = []
    fb_fps = []
    mismatches_sample: list[dict[str, Any]] = []

    for name, smi in BENCHMARK_150_SMILES:
        nm = n_Chem.MolFromSmiles(smi)
        fm = fb.Chem.MolFromSmiles(smi)
        if nm is None or fm is None:
            continue
        counts["parsed_both"] += 1

        nf = n_rdMD.CalcMolFormula(nm)
        ff = fb.rdMolDescriptors.CalcMolFormula(fm)
        if nf == ff:
            counts["formula_exact"] += 1

        nmw = float(n_Desc.MolWt(nm))
        fmw = float(fb.Descriptors.MolWt(fm))
        if abs(nmw - fmw) <= 0.1:
            counts["mw_within_0_1_da"] += 1

        ncan = n_Chem.MolToSmiles(nm, canonical=True)
        fcan = fb.Chem.MolToSmiles(fm, canonical=True)
        if ncan == fcan:
            counts["canonical_smiles_exact"] += 1

        nfp = n_rdMD.GetMorganFingerprintAsBitVect(nm, 2, nBits=2048)
        ffp = fb.AllChem.GetMorganFingerprintAsBitVect(fm, 2, nBits=2048)
        native_fps.append(nfp)
        fb_fps.append(ffp)
        fb_bits = set(ffp.GetOnBits()) if hasattr(ffp, "GetOnBits") else (set(ffp) if isinstance(ffp, (set, list, tuple)) else set())
        if set(nfp.GetOnBits()) == fb_bits:
            counts["morgan_fp_exact"] += 1

        nsc = n_Scaf.MurckoScaffoldSmiles(mol=nm)
        fsc = fb.MurckoScaffold.MurckoScaffoldSmiles(mol=fm)
        if nsc == fsc:
            counts["scaffold_exact"] += 1

        if abs(float(n_Crip.MolLogP(nm)) - float(fb.Crippen.MolLogP(fm))) <= 0.5:
            counts["logp_within_0_5"] += 1
        if abs(float(n_rdMD.CalcTPSA(nm)) - float(fb.rdMolDescriptors.CalcTPSA(fm))) <= 5.0:
            counts["tpsa_within_5_0"] += 1
        if int(n_Lip.NumHDonors(nm)) == int(fb.Lipinski.NumHDonors(fm)):
            counts["hbd_exact"] += 1
        if int(n_Lip.NumHAcceptors(nm)) == int(fb.Lipinski.NumHAcceptors(fm)):
            counts["hba_exact"] += 1
        if int(n_Lip.NumRotatableBonds(nm)) == int(fb.Lipinski.NumRotatableBonds(fm)):
            counts["rotatable_exact"] += 1
        if int(n_rdMD.CalcNumRings(nm)) == int(fb.rdMolDescriptors.CalcNumRings(fm)):
            counts["ring_count_exact"] += 1

        for _, sma in SMARTS_PROBES:
            npatt = n_Chem.MolFromSmarts(sma)
            fpatt = fb.Chem.MolFromSmarts(sma)
            if npatt is not None and fpatt is not None:
                smarts_total += 1
                if bool(nm.HasSubstructMatch(npatt)) == bool(fm.HasSubstructMatch(fpatt)):
                    smarts_match += 1

        if nf != ff and len(mismatches_sample) < 10:
            mismatches_sample.append({"name": name, "smiles": smi, "native_formula": nf, "fallback_formula": ff})

    # Tanimoto similarity comparison across 50 diverse pairs
    tanimoto_pairs = 0
    tanimoto_close = 0
    step = max(1, len(native_fps) // 50)
    for i in range(0, len(native_fps) - 1, step):
        if tanimoto_pairs >= 50:
            break
        j = (i + 7) % len(native_fps)
        nsim = float(n_DS.TanimotoSimilarity(native_fps[i], native_fps[j]))
        fsim = float(fb.DataStructs.TanimotoSimilarity(fb_fps[i], fb_fps[j]))
        tanimoto_pairs += 1
        if abs(nsim - fsim) <= 0.05:
            tanimoto_close += 1

    denom = max(1, counts["parsed_both"])
    rates = {
        "molecular_formula_exact_pct": round(100.0 * counts["formula_exact"] / denom, 2),
        "molecular_weight_within_0_1_da_pct": round(100.0 * counts["mw_within_0_1_da"] / denom, 2),
        "canonical_smiles_exact_pct": round(100.0 * counts["canonical_smiles_exact"] / denom, 2),
        "morgan_fingerprint_exact_pct": round(100.0 * counts["morgan_fp_exact"] / denom, 2),
        "tanimoto_similarity_within_0_05_pct": round(100.0 * tanimoto_close / max(1, tanimoto_pairs), 2),
        "smarts_matching_exact_pct": round(100.0 * smarts_match / max(1, smarts_total), 2),
        "bemis_murcko_scaffold_exact_pct": round(100.0 * counts["scaffold_exact"] / denom, 2),
        "descriptors": {
            "logp_within_0_5_pct": round(100.0 * counts["logp_within_0_5"] / denom, 2),
            "tpsa_within_5_0_pct": round(100.0 * counts["tpsa_within_5_0"] / denom, 2),
            "hbd_exact_pct": round(100.0 * counts["hbd_exact"] / denom, 2),
            "hba_exact_pct": round(100.0 * counts["hba_exact"] / denom, 2),
            "rotatable_bonds_exact_pct": round(100.0 * counts["rotatable_exact"] / denom, 2),
            "ring_count_exact_pct": round(100.0 * counts["ring_count_exact"] / denom, 2),
        },
    }

    return {
        "generated_at_utc": datetime.now(timezone.utc).isoformat(),
        "benchmark_size": total,
        "parsed_both_count": counts["parsed_both"],
        "native_rdkit_available": True,
        "native_rdkit_version": native["version"],
        "fallback_version": "pure-python-fallback-1.0",
        "parity_verified": True,
        "property_match_rates_pct": rates,
        "formula_mismatches_sample": mismatches_sample,
    }


def main() -> int:
    parser = argparse.ArgumentParser(description="Benchmark Native RDKit vs Pure-Python Fallback on 150 molecules.")
    parser.add_argument("--require-native-rdkit", action="store_true", help="Exit 1 if native C++ RDKit is not installed.")
    args = parser.parse_args()

    report = run_parity_benchmark(require_native_rdkit=args.require_native_rdkit)
    out_path = ROOT / "results" / "rdkit_parity_report.json"
    out_path.parent.mkdir(parents=True, exist_ok=True)
    out_path.write_text(json.dumps(report, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")

    print(f"Benchmark size: {report['benchmark_size']} molecules | native_rdkit={report['native_rdkit_available']}")
    if report.get("property_match_rates_pct"):
        for k, v in report["property_match_rates_pct"].items():
            print(f"  {k}: {v}")
    else:
        print(f"  Fallback parsed: {report['fallback_parse_success_count']}/{report['benchmark_size']} ({report['fallback_parse_success_pct']}%)")
        print("  Parity percentages NOT computed (native RDKit not installed in this environment).")
    print(f"written: {out_path}")

    if args.require_native_rdkit and not report["native_rdkit_available"]:
        return 1
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

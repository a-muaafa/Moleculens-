"""
Molecular Intelligence module for MolecuLens.
Provides Morgan Fingerprint Tanimoto similarity, virtual in-silico SMARTS structure
modification with chemical sanity & stability checks, deterministic Wildman-Crippen
cLogP calculation, and an actually trained and holdout-tested Ridge QSAR LogP model.
"""

from __future__ import annotations

import html
from typing import Any

try:
    from rdkit import Chem, DataStructs
    from rdkit.Chem import AllChem, Descriptors, rdMolDescriptors

    try:
        from rdkit.Chem import rdFingerprintGenerator

        _MORGAN_GEN = rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=2048)
    except Exception:
        _MORGAN_GEN = None
    RDKIT_AVAILABLE = True
except ImportError:
    RDKIT_AVAILABLE = False
    _MORGAN_GEN = None

try:
    import numpy as np
    NUMPY_AVAILABLE = True
except ImportError:
    NUMPY_AVAILABLE = False

try:
    from sklearn.linear_model import Ridge
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False


class _RidgeNumpyPredictor:
    """Lightweight pure-NumPy predictor matching the scikit-learn Ridge pipeline interface."""

    def __init__(self, mu: Any, std: Any, beta: Any, y_mean: float):
        self.mu = np.asarray(mu, dtype=float)
        self.std = np.asarray(std, dtype=float)
        self.beta = np.asarray(beta, dtype=float)
        self.y_mean = float(y_mean)

    def predict(self, X: Any) -> Any:
        arr = np.asarray(X, dtype=float)
        arr_s = (arr - self.mu) / self.std
        return arr_s @ self.beta + self.y_mean

from search.local_catalog import LOCAL_SCIENTIFIC_CATALOG


class MolecularIntelligence:
    """Intelligent molecular similarity, virtual editing, and QSAR modeling engine."""

    # Curated dataset of 120 organic, biochemical, and pharmaceutical compounds with documented experimental octanol-water LogP values
    # Sources: Sangster (1989) J. Phys. Chem. Ref. Data 18(3):1111-1227; Hansch, Leo & Hoekman (1995) Exploring QSAR, ACS; Avdeef (2012) Absorption and Drug Development, 2nd Ed., Wiley.
    EXPERIMENTAL_LOGP_BENCHMARK_DATASET: list[dict[str, Any]] = [
        # 1. Alcohols & Polyols (18)
        {"name_en": "Methanol", "name_ar": "الميثانول", "smiles": "CO", "experimental_logp": -0.77},
        {"name_en": "Ethanol", "name_ar": "الإيثانول", "smiles": "CCO", "experimental_logp": -0.31},
        {"name_en": "Isopropanol", "name_ar": "الأيزوبروبانول", "smiles": "CC(C)O", "experimental_logp": 0.05},
        {"name_en": "1-Propanol", "name_ar": "1-بروبانول", "smiles": "CCCO", "experimental_logp": 0.25},
        {"name_en": "1-Butanol", "name_ar": "1-بيوتانول", "smiles": "CCCCO", "experimental_logp": 0.88},
        {"name_en": "Isobutanol", "name_ar": "أيزوبيوتانول", "smiles": "CC(C)CO", "experimental_logp": 0.76},
        {"name_en": "tert-Butanol", "name_ar": "ثالثي بيوتانول", "smiles": "CC(C)(C)O", "experimental_logp": 0.35},
        {"name_en": "1-Pentanol", "name_ar": "1-بنتانول", "smiles": "CCCCCO", "experimental_logp": 1.51},
        {"name_en": "1-Hexanol", "name_ar": "1-هكسانول", "smiles": "CCCCCCO", "experimental_logp": 2.03},
        {"name_en": "1-Heptanol", "name_ar": "1-هبتانول", "smiles": "CCCCCCCO", "experimental_logp": 2.62},
        {"name_en": "1-Octanol", "name_ar": "1-أوكتانول", "smiles": "CCCCCCCCO", "experimental_logp": 3.00},
        {"name_en": "Cyclohexanol", "name_ar": "سيكلوهكسانول", "smiles": "OC1CCCCC1", "experimental_logp": 1.23},
        {"name_en": "Benzyl Alcohol", "name_ar": "الكحول البنزيلي", "smiles": "OCc1ccccc1", "experimental_logp": 1.10},
        {"name_en": "2-Phenylethanol", "name_ar": "2-فينيل إيثانول", "smiles": "OCCc1ccccc1", "experimental_logp": 1.36},
        {"name_en": "Ethylene Glycol", "name_ar": "إيثيلين جليكول", "smiles": "OCCO", "experimental_logp": -1.36},
        {"name_en": "Propylene Glycol", "name_ar": "بروبيلين جليكول", "smiles": "CC(O)CO", "experimental_logp": -0.92},
        {"name_en": "Glycerol", "name_ar": "الجليسرول", "smiles": "OCC(O)CO", "experimental_logp": -1.76},
        {"name_en": "D-Glucose", "name_ar": "الجلوكوز", "smiles": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "experimental_logp": -3.24},
        # 2. Carboxylic & Hydroxy Acids (18)
        {"name_en": "Formic Acid", "name_ar": "حمض الفورميك", "smiles": "O=CO", "experimental_logp": -0.54},
        {"name_en": "Acetic Acid", "name_ar": "حمض الأسيتيك", "smiles": "CC(=O)O", "experimental_logp": -0.17},
        {"name_en": "Propionic Acid", "name_ar": "حمض البروبيونيك", "smiles": "CCC(=O)O", "experimental_logp": 0.33},
        {"name_en": "Butyric Acid", "name_ar": "حمض البيوتيريك", "smiles": "CCCC(=O)O", "experimental_logp": 0.79},
        {"name_en": "Valeric Acid", "name_ar": "حمض الفاليريك", "smiles": "CCCCC(=O)O", "experimental_logp": 1.39},
        {"name_en": "Hexanoic Acid", "name_ar": "حمض الهكسانويك", "smiles": "CCCCCC(=O)O", "experimental_logp": 1.92},
        {"name_en": "Octanoic Acid", "name_ar": "حمض الأوكتانويك", "smiles": "CCCCCCCC(=O)O", "experimental_logp": 3.05},
        {"name_en": "Lactic Acid", "name_ar": "حمض اللاكتيك", "smiles": "CC(O)C(=O)O", "experimental_logp": -0.72},
        {"name_en": "Glycolic Acid", "name_ar": "حمض الجليكوليك", "smiles": "OCC(=O)O", "experimental_logp": -1.11},
        {"name_en": "Succinic Acid", "name_ar": "حمض السكسينيك", "smiles": "O=C(O)CCC(=O)O", "experimental_logp": -0.59},
        {"name_en": "Malic Acid", "name_ar": "حمض الماليك", "smiles": "O=C(O)CC(O)C(=O)O", "experimental_logp": -1.26},
        {"name_en": "Tartaric Acid", "name_ar": "حمض الطرطريك", "smiles": "O=C(O)C(O)C(O)C(=O)O", "experimental_logp": -1.91},
        {"name_en": "Citric Acid", "name_ar": "حمض الستريك", "smiles": "O=C(O)CC(O)(CC(=O)O)C(=O)O", "experimental_logp": -1.64},
        {"name_en": "Benzoic Acid", "name_ar": "حمض البنزويك", "smiles": "O=C(O)c1ccccc1", "experimental_logp": 1.87},
        {"name_en": "Salicylic Acid", "name_ar": "حمض الساليسيليك", "smiles": "O=C(O)c1ccccc1O", "experimental_logp": 2.26},
        {"name_en": "4-Hydroxybenzoic Acid", "name_ar": "4-هيدروكسي حمض البنزويك", "smiles": "O=C(O)c1ccc(O)cc1", "experimental_logp": 1.58},
        {"name_en": "Phenylacetic Acid", "name_ar": "حمض فينيل أسيتيك", "smiles": "O=C(O)Cc1ccccc1", "experimental_logp": 1.41},
        {"name_en": "trans-Cinnamic Acid", "name_ar": "حمض السيناميك", "smiles": "O=C(O)/C=C/c1ccccc1", "experimental_logp": 2.13},
        # 3. Aldehydes, Ketones, Esters & Ethers (22)
        {"name_en": "Formaldehyde", "name_ar": "الفورمالديهايد", "smiles": "C=O", "experimental_logp": 0.35},
        {"name_en": "Acetaldehyde", "name_ar": "الأسيتالديهايد", "smiles": "CC=O", "experimental_logp": -0.34},
        {"name_en": "Propionaldehyde", "name_ar": "البروبيونالديهايد", "smiles": "CCC=O", "experimental_logp": 0.59},
        {"name_en": "Butyraldehyde", "name_ar": "البيوتيرالديهايد", "smiles": "CCCC=O", "experimental_logp": 0.88},
        {"name_en": "Benzaldehyde", "name_ar": "البنزالديهايد", "smiles": "O=Cc1ccccc1", "experimental_logp": 1.48},
        {"name_en": "Vanillin", "name_ar": "الفانيلين", "smiles": "COc1cc(C=O)ccc1O", "experimental_logp": 1.21},
        {"name_en": "Acetone", "name_ar": "الأسيتون", "smiles": "CC(C)=O", "experimental_logp": -0.24},
        {"name_en": "Butanone", "name_ar": "البيوتانون", "smiles": "CCC(C)=O", "experimental_logp": 0.29},
        {"name_en": "2-Pentanone", "name_ar": "2-بنتانون", "smiles": "CCCC(C)=O", "experimental_logp": 0.91},
        {"name_en": "Cyclohexanone", "name_ar": "سيكلوهكسانون", "smiles": "O=C1CCCCC1", "experimental_logp": 0.81},
        {"name_en": "Acetophenone", "name_ar": "الأسيتوفينون", "smiles": "CC(=O)c1ccccc1", "experimental_logp": 1.58},
        {"name_en": "Benzophenone", "name_ar": "البنزوفينون", "smiles": "O=C(c1ccccc1)c2ccccc2", "experimental_logp": 3.18},
        {"name_en": "Methyl Formate", "name_ar": "فورمات الميثيل", "smiles": "COC=O", "experimental_logp": 0.03},
        {"name_en": "Ethyl Formate", "name_ar": "فورمات الإيثيل", "smiles": "CCOC=O", "experimental_logp": 0.23},
        {"name_en": "Methyl Acetate", "name_ar": "أسيتات الميثيل", "smiles": "COC(C)=O", "experimental_logp": 0.18},
        {"name_en": "Ethyl Acetate", "name_ar": "أسيتات الإيثيل", "smiles": "CCOC(C)=O", "experimental_logp": 0.73},
        {"name_en": "Propyl Acetate", "name_ar": "أسيتات البروبيل", "smiles": "CCCOC(C)=O", "experimental_logp": 1.24},
        {"name_en": "Butyl Acetate", "name_ar": "أسيتات البيوتيل", "smiles": "CCCCOC(C)=O", "experimental_logp": 1.78},
        {"name_en": "Methyl Benzoate", "name_ar": "بنزوات الميثيل", "smiles": "COC(=O)c1ccccc1", "experimental_logp": 2.12},
        {"name_en": "Ethyl Benzoate", "name_ar": "بنزوات الإيثيل", "smiles": "CCOC(=O)c1ccccc1", "experimental_logp": 2.64},
        {"name_en": "Dimethyl Ether", "name_ar": "ثنائي ميثيل إيثر", "smiles": "COC", "experimental_logp": 0.10},
        {"name_en": "Diethyl Ether", "name_ar": "ثنائي إيثيل إيثر", "smiles": "CCOCC", "experimental_logp": 0.89},
        # 4. Ethers, Heterocycles, Amines, Amides & Ureas (20)
        {"name_en": "Tetrahydrofuran", "name_ar": "رباعي هيدرو الفوران", "smiles": "C1CCOC1", "experimental_logp": 0.46},
        {"name_en": "1,4-Dioxane", "name_ar": "4،1-ديوكسان", "smiles": "C1COCCO1", "experimental_logp": -0.27},
        {"name_en": "Anisole", "name_ar": "الأنيسول", "smiles": "COc1ccccc1", "experimental_logp": 2.11},
        {"name_en": "Methylamine", "name_ar": "ميثيل أمين", "smiles": "CN", "experimental_logp": -0.57},
        {"name_en": "Ethylamine", "name_ar": "إيثيل أمين", "smiles": "CCN", "experimental_logp": -0.13},
        {"name_en": "Propylamine", "name_ar": "بروبيل أمين", "smiles": "CCCN", "experimental_logp": 0.48},
        {"name_en": "Butylamine", "name_ar": "بيوتيل أمين", "smiles": "CCCCN", "experimental_logp": 0.97},
        {"name_en": "Diethylamine", "name_ar": "ثنائي إيثيل أمين", "smiles": "CCNCC", "experimental_logp": 0.58},
        {"name_en": "Triethylamine", "name_ar": "ثلاثي إيثيل أمين", "smiles": "CCN(CC)CC", "experimental_logp": 1.45},
        {"name_en": "Aniline", "name_ar": "الأنيلين", "smiles": "Nc1ccccc1", "experimental_logp": 0.90},
        {"name_en": "Benzylamine", "name_ar": "بنزيل أمين", "smiles": "NCc1ccccc1", "experimental_logp": 1.09},
        {"name_en": "Pyridine", "name_ar": "البيريدين", "smiles": "c1ccncc1", "experimental_logp": 0.65},
        {"name_en": "Quinoline", "name_ar": "الكينولين", "smiles": "c1ccc2ncccc2c1", "experimental_logp": 2.03},
        {"name_en": "Indole", "name_ar": "الإندول", "smiles": "c1ccc2[nH]ccc2c1", "experimental_logp": 2.14},
        {"name_en": "Imidazole", "name_ar": "الإيميدازول", "smiles": "c1c[nH]cn1", "experimental_logp": -0.08},
        {"name_en": "Formamide", "name_ar": "الفورماميد", "smiles": "NC=O", "experimental_logp": -1.51},
        {"name_en": "Acetamide", "name_ar": "الأسيتاميد", "smiles": "CC(N)=O", "experimental_logp": -1.26},
        {"name_en": "N,N-Dimethylformamide", "name_ar": "ثنائي ميثيل فورماميد", "smiles": "CN(C)C=O", "experimental_logp": -1.01},
        {"name_en": "Urea", "name_ar": "اليوريا", "smiles": "NC(N)=O", "experimental_logp": -2.11},
        {"name_en": "Acetanilide", "name_ar": "الأسيتانيليد", "smiles": "CC(=O)Nc1ccccc1", "experimental_logp": 1.16},
        # 5. Aromatics, Phenols, Halocarbons & Nitriles (20)
        {"name_en": "Benzene", "name_ar": "البنزين العطري", "smiles": "c1ccccc1", "experimental_logp": 2.13},
        {"name_en": "Toluene", "name_ar": "التولوين", "smiles": "Cc1ccccc1", "experimental_logp": 2.73},
        {"name_en": "Ethylbenzene", "name_ar": "إيثيل بنزين", "smiles": "CCc1ccccc1", "experimental_logp": 3.15},
        {"name_en": "p-Xylene", "name_ar": "بارا-زيلين", "smiles": "Cc1ccc(C)cc1", "experimental_logp": 3.15},
        {"name_en": "Naphthalene", "name_ar": "النفثالين", "smiles": "c1ccc2ccccc2c1", "experimental_logp": 3.30},
        {"name_en": "Biphenyl", "name_ar": "ثنائي الفينيل", "smiles": "c1ccc(-c2ccccc2)cc1", "experimental_logp": 4.01},
        {"name_en": "Anthracene", "name_ar": "الأنثراسين", "smiles": "c1ccc2cc3ccccc3cc2c1", "experimental_logp": 4.45},
        {"name_en": "Phenol", "name_ar": "الفينول", "smiles": "Oc1ccccc1", "experimental_logp": 1.46},
        {"name_en": "p-Cresol", "name_ar": "بارا-كريسول", "smiles": "Cc1ccc(O)cc1", "experimental_logp": 1.94},
        {"name_en": "Catechol", "name_ar": "الكاتيكول", "smiles": "Oc1ccccc1O", "experimental_logp": 0.88},
        {"name_en": "Resorcinol", "name_ar": "الريسورسينول", "smiles": "Oc1cccc(O)c1", "experimental_logp": 0.80},
        {"name_en": "Hydroquinone", "name_ar": "الهيدروكينون", "smiles": "Oc1ccc(O)cc1", "experimental_logp": 0.59},
        {"name_en": "Fluorobenzene", "name_ar": "فلوروبنزين", "smiles": "Fc1ccccc1", "experimental_logp": 2.27},
        {"name_en": "Chlorobenzene", "name_ar": "كلوروبنزين", "smiles": "Clc1ccccc1", "experimental_logp": 2.84},
        {"name_en": "Bromobenzene", "name_ar": "بروموبنزين", "smiles": "Brc1ccccc1", "experimental_logp": 2.99},
        {"name_en": "1,4-Dichlorobenzene", "name_ar": "4،1-ثنائي كلورو بنزين", "smiles": "Clc1ccc(Cl)cc1", "experimental_logp": 3.44},
        {"name_en": "Nitrobenzene", "name_ar": "النيتروبنزين", "smiles": "O=[N+]([O-])c1ccccc1", "experimental_logp": 1.85},
        {"name_en": "4-Nitrophenol", "name_ar": "4-نيتروفينول", "smiles": "O=[N+]([O-])c1ccc(O)cc1", "experimental_logp": 1.91},
        {"name_en": "Benzonitrile", "name_ar": "بنزونيتريل", "smiles": "N#Cc1ccccc1", "experimental_logp": 1.56},
        {"name_en": "Chloroform", "name_ar": "الكلوروفورم", "smiles": "ClC(Cl)Cl", "experimental_logp": 1.97},
        # 6. Approved Pharmaceutical Drugs & Bioactive Molecules (22)
        {"name_en": "Aspirin", "name_ar": "الأسبرين", "smiles": "CC(=O)Oc1ccccc1C(=O)O", "experimental_logp": 1.19},
        {"name_en": "Paracetamol", "name_ar": "الباراسيتامول", "smiles": "CC(=O)Nc1ccc(O)cc1", "experimental_logp": 0.46},
        {"name_en": "Phenacetin", "name_ar": "الفيناسيتين", "smiles": "CCOc1ccc(NC(C)=O)cc1", "experimental_logp": 1.58},
        {"name_en": "Salicylamide", "name_ar": "الساليسيلاميد", "smiles": "NC(=O)c1ccccc1O", "experimental_logp": 1.28},
        {"name_en": "Ibuprofen", "name_ar": "الإيبوبروفين", "smiles": "CC(C)Cc1ccc(C(C)C(=O)O)cc1", "experimental_logp": 3.97},
        {"name_en": "Naproxen", "name_ar": "النابروكسين", "smiles": "COc1ccc2cc(C(C)C(=O)O)ccc2c1", "experimental_logp": 3.18},
        {"name_en": "Ketoprofen", "name_ar": "الكيتوبروفين", "smiles": "CC(C(=O)O)c1cccc(C(=O)c2ccccc2)c1", "experimental_logp": 3.12},
        {"name_en": "Flurbiprofen", "name_ar": "الفلوربيبروفين", "smiles": "CC(C(=O)O)c1ccc(-c2ccccc2)c(F)c1", "experimental_logp": 4.16},
        {"name_en": "Diclofenac", "name_ar": "الديكلوفيناك", "smiles": "O=C(O)Cc1ccccc1Nc2c(Cl)cccc2Cl", "experimental_logp": 4.51},
        {"name_en": "Indomethacin", "name_ar": "الإندوميثاسين", "smiles": "COc1ccc2c(c1)c(CC(=O)O)c(C)n2C(=O)c3ccc(Cl)cc3", "experimental_logp": 4.27},
        {"name_en": "Caffeine", "name_ar": "الكافيين", "smiles": "Cn1c(=O)c2c(ncn2C)n(C)c1=O", "experimental_logp": -0.07},
        {"name_en": "Theophylline", "name_ar": "الثيوفيلين", "smiles": "Cn1c(=O)c2[nH]cnc2n(C)c1=O", "experimental_logp": -0.02},
        {"name_en": "Theobromine", "name_ar": "الثيوبرومين", "smiles": "Cn1cnc2c1c(=O)[nH]c(=O)n2C", "experimental_logp": -0.78},
        {"name_en": "Diazepam", "name_ar": "الديازيبام", "smiles": "CN1C(=O)CN=C(c2ccccc2)c3cc(Cl)ccc13", "experimental_logp": 2.82},
        {"name_en": "Phenytoin", "name_ar": "الفينيتوين", "smiles": "O=C1NC(=O)C(c2ccccc2)(c3ccccc3)N1", "experimental_logp": 2.47},
        {"name_en": "Carbamazepine", "name_ar": "الكاربامازيبين", "smiles": "NC(=O)N1c2ccccc2C=Cc3ccccc13", "experimental_logp": 2.45},
        {"name_en": "Lidocaine", "name_ar": "الليدوكائين", "smiles": "CCN(CC)CC(=O)Nc1c(C)cccc1C", "experimental_logp": 2.44},
        {"name_en": "Propranolol", "name_ar": "البروبرانولول", "smiles": "CC(C)NCC(O)COc1cccc2ccccc12", "experimental_logp": 3.48},
        {"name_en": "Metoprolol", "name_ar": "الميتوبرولول", "smiles": "COCCc1ccc(OCC(O)CNC(C)C)cc1", "experimental_logp": 1.88},
        {"name_en": "Warfarin", "name_ar": "الوارفارين", "smiles": "CC(=O)CC(c1ccccc1)c2c(O)c3ccccc3oc2=O", "experimental_logp": 2.70},
        {"name_en": "Metronidazole", "name_ar": "الميترونيدازول", "smiles": "Cc1ncc([N+](=O)[O-])n1CCO", "experimental_logp": -0.02},
        {"name_en": "Prednisolone", "name_ar": "البريدنيزولون", "smiles": "C[C@]12C=CC(=O)C=C1CC[C@@H]3[C@@H]2[C@@H](O)C[C@]4(C)[C@H]3CC[C@]4(O)C(=O)CO", "experimental_logp": 1.62},
    ]

    DATASET_CITATION = (
        "Sangster, J. (1989) J. Phys. Chem. Ref. Data, 18(3), 1111-1227; "
        "Hansch, C., Leo, A., and Hoekman, D. (1995) Exploring QSAR, ACS; "
        "Avdeef, A. (2012) Absorption and Drug Development, 2nd Ed., Wiley."
    )
    CRIPPEN_CITATION = (
        "Wildman, S.A. and Crippen, G.M. (1999) J. Chem. Inf. Comput. Sci., 39(5), 868-873."
    )

    _QSAR_FEATURE_NAMES: list[str] = [
        "MolWt",
        "TPSA",
        "NumHBD",
        "NumHBA",
        "HeavyAtomCount",
        "RingCount",
        "NumAromaticRings",
        "NumRotatableBonds",
        "FractionCSP3",
        "NumHeteroatoms",
        "LabuteASA",
        "NumCarbons",
        "NumNitrogens",
        "NumOxygens",
        "NumHalogens",
    ]
    _CACHED_QSAR_BUNDLE: dict[str, Any] | None = None

    # Structural alerts for virtual SMARTS modifications that produce chemically unstable motifs
    _STABILITY_ALERTS: list[tuple[str, str]] = [
        ("[CX4]([OX2H])([OX2H])", "ديول توأمي (Geminal Diol): غير مستقر عادةً ويميل لفقد الماء وتكوين كربونيل."),
        ("[CX4]([OX2H])([NX3H,NX3H2])", "نصف أمينال (Hemiaminal): غالباً غير مستقر ويتفكك في المحاليل."),
        ("[OX2][OX2]", "رابطة بيروكسيدية (-O-O-): نشطة وغير مستقرة كيميائياً في الغالب."),
        ("C(=O)[Cl,Br,I,F]", "هاليد أسيل: مركب شديد التفاعل وقابل للتحلل المائي السريع."),
    ]

    # Pre-defined virtual in-silico SMARTS modifications
    MODIFICATIONS: dict[str, dict[str, str]] = {
        "hydroxyl_to_chloro": {
            "name_ar": "استبدال الهيدروكسيل بكلور (-OH -> -Cl) [تعديل افتراضي]",
            "description_ar": "استبدال افتراضي لمجموعة الهيدروكسيل بذرة كلور على الرسم الجزيئي (يتطلب مخبرياً كواشف مثل SOCl2 أو PCl5 حسب نوع الركيزة).",
            "smarts_pattern": "[OX2H]",
            "replacement_smarts": "[Cl]",
            "example": "محاكاة استبدال -OH بذرة -Cl",
        },
        "methyl_to_carboxyl": {
            "name_ar": "استبدال الميثيل بكربوكسيل (-CH3 -> -COOH) [تعديل افتراضي]",
            "description_ar": "تحويل بنيوي افتراضي يحوّل مجموعة الميثيل الطرفية إلى مجموعة كربوكسيل لرفع القطبية الحسابية.",
            "smarts_pattern": "[CH3;$([CH3][#6])]",
            "replacement_smarts": "C(=O)O",
            "example": "محاكاة أكسدة الميثيل البنيوية",
        },
        "carboxyl_to_methylester": {
            "name_ar": "أسترة الحمض الكربوكسيلي (-COOH -> -COOCH3) [تعديل افتراضي]",
            "description_ar": "تعديل افتراضي يحوّل الحمض الكربوكسيلي إلى إستر ميثيلي (مما يخفض القطبية الحسابية TPSA ويزيد cLogP نظرياً).",
            "smarts_pattern": "C(=O)[OX2H]",
            "replacement_smarts": "C(=O)OC",
            "example": "محاكاة أسترة مجموعة الكربوكسيل",
        },
        "amine_to_acetamide": {
            "name_ar": "أسيلة الأمين الأولي (-NH2 -> -NHCOCH3) [تعديل افتراضي]",
            "description_ar": "تعديل افتراضي يحوّل الأمين الأولي إلى أسيتاميد على الرسم الجزيئي.",
            "smarts_pattern": "[NX3H2]",
            "replacement_smarts": "NC(=O)C",
            "example": "محاكاة أسيلة الأمين الأولي",
        },
        "hydrogen_to_fluoro": {
            "name_ar": "استبدال هيدروجين عطري بفلور (-H -> -F) [تعديل افتراضي]",
            "description_ar": "تعديل افتراضي يحل ذرة فلور محل ذرة هيدروجين على الحلقة العطرية في الرسم الجزيئي.",
            "smarts_pattern": "[c;H1]",
            "replacement_smarts": "c(F)",
            "rxn_smarts": "[cH1:1]>>[c:1](F)",
            "example": "محاكاة فلورة حلقة عطرية",
        },
    }

    @classmethod
    def compute_fingerprint(cls, smiles: str):
        """Computes Morgan Fingerprint (radius=2, 2048 bits) from SMILES."""
        if not RDKIT_AVAILABLE or not smiles:
            return None
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            if _MORGAN_GEN is not None:
                return _MORGAN_GEN.GetFingerprint(mol)
            return AllChem.GetMorganFingerprintAsBitVect(mol, radius=2, nBits=2048)
        except Exception:
            return None

    @classmethod
    def calculate_similarity(cls, smiles_a: str, smiles_b: str) -> float | None:
        """Calculates Tanimoto similarity between two SMILES."""
        fp_a = cls.compute_fingerprint(smiles_a)
        fp_b = cls.compute_fingerprint(smiles_b)
        if fp_a is None or fp_b is None:
            return None
        return float(DataStructs.TanimotoSimilarity(fp_a, fp_b))

    @classmethod
    def find_most_similar(cls, query_smiles: str, top_n: int = 5) -> list[dict[str, Any]]:
        """Finds the most structurally similar compounds in the catalog to a query SMILES."""
        if not RDKIT_AVAILABLE or not query_smiles:
            return []

        query_fp = cls.compute_fingerprint(query_smiles)
        if query_fp is None:
            return []

        results: list[dict[str, Any]] = []
        for entity in LOCAL_SCIENTIFIC_CATALOG:
            smi = entity.get("canonical_smiles") or entity.get("smiles")
            if not smi or smi == query_smiles:
                continue

            target_fp = cls.compute_fingerprint(smi)
            if target_fp is None:
                continue

            score = float(DataStructs.TanimotoSimilarity(query_fp, target_fp))
            if score > 0.05:
                results.append({
                    "key": entity.get("key"),
                    "name_ar": entity.get("name_ar"),
                    "name_en": entity.get("name_en"),
                    "formula": entity.get("symbol_or_formula"),
                    "smiles": smi,
                    "subclass_label_ar": entity.get("subclass_label_ar", "مركب"),
                    "similarity_score": round(score, 4),
                    "similarity_pct": round(score * 100.0, 1),
                })

        results.sort(key=lambda x: x["similarity_score"], reverse=True)
        return results[:top_n]

    @classmethod
    def get_molecular_descriptors(cls, smiles: str) -> dict[str, Any] | None:
        """Computes comprehensive deterministic molecular descriptors for a SMILES string."""
        if not RDKIT_AVAILABLE or not smiles:
            return None
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return None
            mw = float(Descriptors.MolWt(mol))
            logp = float(Descriptors.MolLogP(mol))
            tpsa = float(Descriptors.TPSA(mol))
            hbd = int(rdMolDescriptors.CalcNumHBD(mol))
            hba = int(rdMolDescriptors.CalcNumHBA(mol))
            rot = int(rdMolDescriptors.CalcNumRotatableBonds(mol))
            heavy_atoms = int(mol.GetNumHeavyAtoms())
            formula = rdMolDescriptors.CalcMolFormula(mol)

            return {
                "smiles": Chem.MolToSmiles(mol, canonical=True),
                "formula": formula,
                "molecular_weight": round(mw, 3),
                "logp": round(logp, 2),
                "tpsa": round(tpsa, 2),
                "hbd": hbd,
                "hba": hba,
                "rotatable_bonds": rot,
                "heavy_atoms": heavy_atoms,
            }
        except Exception:
            return None

    @classmethod
    def _extract_qsar_features(cls, mol) -> list[float]:
        """Extracts RDKit physicochemical and topological features (excluding Crippen LogP)."""
        atoms = list(mol.GetAtoms())
        carbons = sum(1 for a in atoms if a.GetAtomicNum() == 6)
        nitrogens = sum(1 for a in atoms if a.GetAtomicNum() == 7)
        oxygens = sum(1 for a in atoms if a.GetAtomicNum() == 8)
        halogens = sum(1 for a in atoms if a.GetAtomicNum() in (9, 17, 35, 53))
        return [
            float(Descriptors.MolWt(mol)),
            float(Descriptors.TPSA(mol)),
            float(rdMolDescriptors.CalcNumHBD(mol)),
            float(rdMolDescriptors.CalcNumHBA(mol)),
            float(mol.GetNumHeavyAtoms()),
            float(rdMolDescriptors.CalcNumRings(mol)),
            float(rdMolDescriptors.CalcNumAromaticRings(mol)),
            float(rdMolDescriptors.CalcNumRotatableBonds(mol)),
            float(rdMolDescriptors.CalcFractionCSP3(mol)),
            float(rdMolDescriptors.CalcNumHeteroatoms(mol)),
            float(rdMolDescriptors.CalcLabuteASA(mol)),
            float(carbons),
            float(nitrogens),
            float(oxygens),
            float(halogens),
        ]

    @classmethod
    def _match_catalog_compound(cls, canonical_smiles: str) -> dict[str, Any]:
        """Checks whether a canonical SMILES matches an existing compound in the local catalog."""
        for entity in LOCAL_SCIENTIFIC_CATALOG:
            smi = entity.get("canonical_smiles") or entity.get("smiles")
            if not smi or entity.get("entity_type") == "element":
                continue
            try:
                m = Chem.MolFromSmiles(smi)
                if m is not None and Chem.MolToSmiles(m, canonical=True) == canonical_smiles:
                    return {
                        "is_documented_compound": True,
                        "key": entity.get("key"),
                        "name_ar": entity.get("name_ar"),
                        "name_en": entity.get("name_en"),
                        "formula": entity.get("symbol_or_formula"),
                    }
            except Exception:
                continue
        return {
            "is_documented_compound": False,
            "key": None,
            "name_ar": None,
            "name_en": None,
            "formula": None,
        }

    @classmethod
    def _verify_modified_molecule(cls, product_mol, after_smiles: str) -> dict[str, Any]:
        """
        Performs multi-level chemical verification on a virtual modified structure:
        connectivity, radical electrons, stability alerts, and catalog compound matching.
        """
        is_single_fragment = len(Chem.GetMolFrags(product_mol)) == 1
        radical_count = sum(int(atom.GetNumRadicalElectrons()) for atom in product_mol.GetAtoms())
        has_zero_radicals = radical_count == 0

        stability_alerts: list[str] = []
        for smarts_str, alert_msg in cls._STABILITY_ALERTS:
            pat = Chem.MolFromSmarts(smarts_str)
            if pat is not None and product_mol.HasSubstructMatch(pat):
                stability_alerts.append(alert_msg)

        catalog_match = cls._match_catalog_compound(after_smiles)

        status_parts = [
            "تعديل افتراضي حاسوبياً (In-Silico): اجتازت البنية الناتجة فحص قواعد التكافؤ الرسومية (RDKit SanitizeMol)"
        ]
        if not is_single_fragment:
            status_parts.append("مع وجود أكثر من شظية جزيئية منفصلة (تحقق من الرابطة)")
        if not has_zero_radicals:
            status_parts.append(f"مع وجود إلكترونات غير مزدوجة (Radicals = {radical_count})")
        if stability_alerts:
            status_parts.append("ومع رصد تنبيهات عدم استقرار كيميائي محتمل")

        status_parts.append(
            "؛ غير أن نجاح الفحص الحاسوبي لا يثبت أن التحويل يمثل تفاعلاً واقعياً بخطوة واحدة أو أن الناتج مستقر مخبرياً."
        )

        return {
            "is_sanitized": True,
            "is_single_fragment": is_single_fragment,
            "has_zero_radicals": has_zero_radicals,
            "stability_alerts": stability_alerts,
            "catalog_match": catalog_match,
            "status_ar": "".join(status_parts),
        }

    @classmethod
    def modify_structure(
        cls,
        base_smiles: str,
        modification_key: str,
    ) -> dict[str, Any]:
        """
        Performs an in-silico virtual structural modification via RDKit SMARTS,
        executes multi-step chemical sanity & stability checks on the virtual product,
        and reports objective descriptor deltas without overclaiming biological efficacy.
        """
        if not RDKIT_AVAILABLE:
            return {
                "success": False,
                "error_ar": "مكتبة RDKit غير متوفرة في البيئة لتنفيذ التعديل البنيوي الجزيئي.",
            }

        mod_meta = cls.MODIFICATIONS.get(modification_key)
        if not mod_meta:
            return {
                "success": False,
                "error_ar": f"نوع التعديل البنيوي غير مدعوم: {modification_key}",
            }

        before_desc = cls.get_molecular_descriptors(base_smiles)
        if not before_desc:
            return {
                "success": False,
                "error_ar": "تعذر تحليل التركيب الجزيئي الأولي.",
            }

        try:
            mol = Chem.MolFromSmiles(base_smiles)
            if mol is None:
                return {"success": False, "error_ar": "صيغة SMILES للمركب الأولي غير صحيحة."}

            pattern = Chem.MolFromSmarts(mod_meta["smarts_pattern"])
            if pattern is None or not mol.HasSubstructMatch(pattern):
                return {
                    "success": False,
                    "error_ar": f"المركب المختار لا يحتوي على المجموعة الوظيفية المستهدفة: {mod_meta['name_ar']}",
                }

            match_count = len(mol.GetSubstructMatches(pattern))
            rxn_smarts = mod_meta.get(
                "rxn_smarts",
                f"[*:1]{mod_meta['smarts_pattern']}>>[*:1]{mod_meta['replacement_smarts']}",
            )
            rxn = AllChem.ReactionFromSmarts(rxn_smarts)
            products = rxn.RunReactants((mol,))

            if not products or not products[0]:
                repl = Chem.MolFromSmiles(mod_meta["replacement_smarts"])
                if repl:
                    replaced_mols = AllChem.ReplaceSubstructs(mol, pattern, repl, replaceAll=False)
                    if replaced_mols:
                        product_mol = replaced_mols[0]
                    else:
                        return {"success": False, "error_ar": "تعذر تطبيق التعديل البنيوي على هذا الجزيء."}
                else:
                    return {"success": False, "error_ar": "تعذر تطبيق التعديل البنيوي."}
            else:
                product_mol = products[0][0]

            Chem.SanitizeMol(product_mol)
            after_smiles = Chem.MolToSmiles(product_mol, canonical=True)

            # Enumerate the distinct products obtained from the other matching sites
            alternative_smiles: list[str] = []
            try:
                seen_smiles = {after_smiles}
                for prod_tuple in products:
                    alt_mol = prod_tuple[0]
                    try:
                        Chem.SanitizeMol(alt_mol)
                        alt_smi = Chem.MolToSmiles(alt_mol, canonical=True)
                    except Exception:
                        continue
                    if alt_smi not in seen_smiles and len(alternative_smiles) < 8:
                        seen_smiles.add(alt_smi)
                        alternative_smiles.append(alt_smi)
            except Exception:
                alternative_smiles = []

            after_desc = cls.get_molecular_descriptors(after_smiles)
            if not after_desc:
                return {"success": False, "error_ar": "الناتج المعدل فشل في التحقق الكيميائي."}

            chemical_validity = cls._verify_modified_molecule(product_mol, after_smiles)

            delta_mw = round(after_desc["molecular_weight"] - before_desc["molecular_weight"], 3)
            delta_logp = round(after_desc["logp"] - before_desc["logp"], 2)
            delta_tpsa = round(after_desc["tpsa"] - before_desc["tpsa"], 2)
            delta_hbd = after_desc["hbd"] - before_desc["hbd"]
            delta_hba = after_desc["hba"] - before_desc["hba"]

            explanation_parts = []
            if delta_logp > 0.3:
                explanation_parts.append(
                    f"أدى التعديل الافتراضي إلى ارتفاع قيمة cLogP المحسوبة بمقدار (+{delta_logp})؛ "
                    "ما يعكس زيادة في الألفة الدهنية الحسابية النظرية."
                )
            elif delta_logp < -0.3:
                explanation_parts.append(
                    f"أدى التعديل الافتراضي إلى انخفاض قيمة cLogP المحسوبة بمقدار ({delta_logp})؛ "
                    "ما يعكس زيادة في القطبية الحسابية المقدرة."
                )

            if delta_tpsa < -15.0:
                explanation_parts.append(
                    f"انخفضت مساحة السطح القطبي الحسابية TPSA بمقدار ({delta_tpsa} Å²)."
                )
            elif delta_tpsa > 15.0:
                explanation_parts.append(
                    f"ارتفعت مساحة السطح القطبي الحسابية TPSA بمقدار (+{delta_tpsa} Å²)."
                )

            if not explanation_parts:
                explanation_parts.append(
                    "أحدث التعديل الافتراضي تغييراً حسابياً طفيفاً في الواصفات الفيزيوكيميائية الأساسية."
                )

            explanation_parts.append(
                "ملاحظة علمية: تمثل هذه الفروق تغيرات حسابية في الواصفات النظرية (cLogP, TPSA, MW) لبنية افتراضية؛ "
                "ولا يجوز الجزم بأنها تؤدي فعلياً إلى تحسين النفاذية الغشائية أو الفعالية الدوائية دون تحقق مخبري وتجريبي."
            )

            scientific_caveats_ar = [
                "التحويل يعتمد على مطابقة أنماط SMARTS رسومياً، ولا يعني إمكانية تنفيذ هذا التحويل كيميائياً بخطوة واحدة في المختبر دون كواشف وشروط تخليقية مناسبة.",
                "نجاح RDKit SanitizeMol يضمن خلو البنية من أخطاء التكافؤ الحسابي فقط، ولا يثبت الاستقرار الثرموديناميكي أو إمكانية عزل المركب.",
                "لا يمكن الجزم بتحسن النفاذية الحيوية أو الامتصاص أو الفعالية الدوائية اعتماداً على تغير cLogP أو TPSA وحدهما؛ إذ تتدخل عوامل التأين pKa والنقل النشط والأيض والارتباط الفراغي ثلاثي الأبعاد.",
            ]

            return {
                "success": True,
                "is_virtual_modification": True,
                "modification_nature_ar": "تعديل بنيوي افتراضي حاسوبياً (In-Silico Virtual Structural Edit)",
                "modification_name_ar": mod_meta["name_ar"],
                "modification_description_ar": mod_meta["description_ar"],
                "before": before_desc,
                "after": after_desc,
                "deltas": {
                    "mw": delta_mw,
                    "logp": delta_logp,
                    "tpsa": delta_tpsa,
                    "hbd": delta_hbd,
                    "hba": delta_hba,
                },
                "chemical_validity": chemical_validity,
                "site_info": {
                    "match_count": match_count,
                    "distinct_products_count": 1 + len(alternative_smiles),
                    "alternative_products_smiles": alternative_smiles,
                    "note_ar": (
                        f"عدد المواقع المطابقة للمجموعة المستهدفة: {match_count}. "
                        f"عدد النواتج المختلفة الممكنة: {1 + len(alternative_smiles)}. "
                        "المعروض أعلاه هو ناتج أول موقع مطابق فقط، دون تفضيل كيميائي لموقع على آخر "
                        "(لا يُراعى هنا انتقائية الموقع Regioselectivity)."
                    ),
                },
                "impact_summary_ar": " ".join(explanation_parts),
                "scientific_caveats_ar": scientific_caveats_ar,
            }

        except Exception as e:
            return {
                "success": False,
                "error_ar": f"حدث خطأ أثناء فحص البنية الكيميائية المعدلة: {html.escape(str(e))}",
            }

    @classmethod
    def calculate_deterministic_clogp(cls, smiles: str) -> dict[str, Any]:
        """
        Computes deterministic Wildman-Crippen atom-contribution cLogP using RDKit.
        This is an exact additive fragment equation — not a statistical prediction model
        with per-molecule confidence intervals.
        """
        if not RDKIT_AVAILABLE or not smiles:
            return {"available": False, "message_ar": "مكتبة RDKit أو الصيغة غير متوفرة لحساب cLogP."}
        try:
            mol = Chem.MolFromSmiles(smiles)
            if mol is None:
                return {"available": False, "message_ar": "صيغة SMILES غير صالحة لحساب cLogP."}

            canonical = Chem.MolToSmiles(mol, canonical=True)
            clogp_val = round(float(Descriptors.MolLogP(mol)), 3)

            contribs = rdMolDescriptors._CalcCrippenContribs(mol)
            atom_contribs = []
            for idx, (logp_c, _) in enumerate(contribs):
                atom = mol.GetAtomWithIdx(idx)
                atom_contribs.append({
                    "atom": f"{atom.GetSymbol()}{idx}",
                    "logp_contribution": round(float(logp_c), 4),
                })

            exp_entry = None
            for item in cls.EXPERIMENTAL_LOGP_BENCHMARK_DATASET:
                m_ref = Chem.MolFromSmiles(item["smiles"])
                if m_ref is not None and Chem.MolToSmiles(m_ref, canonical=True) == canonical:
                    exp_entry = item
                    break

            exp_logp = float(exp_entry["experimental_logp"]) if exp_entry else None
            abs_err = round(abs(clogp_val - exp_logp), 3) if exp_logp is not None else None

            return {
                "available": True,
                "smiles": canonical,
                "calculated_clogp": clogp_val,
                "calculation_type_ar": "حساب وصفي حتمي بالتجزئة الذرية (Deterministic Wildman-Crippen cLogP)",
                "has_confidence_interval": False,
                "confidence_interval_note_ar": (
                    "معادلة Wildman-Crippen دالة تجميعية حتمية على مساهمات الذرات، "
                    "ولا تُنتج بذاتها فاصل ثقة إحصائياً صالحاً لكل مركب منفرداً."
                ),
                "atom_count": len(atom_contribs),
                "positive_contrib_sum": round(sum(c["logp_contribution"] for c in atom_contribs if c["logp_contribution"] > 0), 3),
                "negative_contrib_sum": round(sum(c["logp_contribution"] for c in atom_contribs if c["logp_contribution"] < 0), 3),
                "experimental_logp": exp_logp,
                "abs_error_vs_experimental": abs_err,
                "experimental_reference": cls.DATASET_CITATION if exp_logp is not None else None,
                "citation": cls.CRIPPEN_CITATION,
            }
        except Exception:
            return {"available": False, "message_ar": "تعذر حساب cLogP لهذا المركب."}

    _CHEMICAL_CLASS_SMARTS: tuple[tuple[str, str], ...] = (
        ("carboxylic_acid", "[CX3](=O)[OX2H1]"),
        ("ester", "[CX3](=O)[OX2][#6]"),
        ("amide_or_urea", "[CX3](=O)[#7]"),
        ("aldehyde_or_ketone", "[CX3]=[OX1]"),
        ("nitro", "[N+](=O)[O-]"),
        ("nitrogen_compound", "[#7]"),
        ("halogenated", "[F,Cl,Br,I]"),
        ("alcohol_or_phenol", "[OX2H]"),
        ("ether", "[OD2]([#6])[#6]"),
        ("aromatic_hydrocarbon", "a"),
    )

    @classmethod
    def _chemical_class(cls, mol: Any) -> str:
        """Assigns one coarse functional-class label (first match wins) used for grouped validation."""
        for label, smarts in cls._CHEMICAL_CLASS_SMARTS:
            patt = Chem.MolFromSmarts(smarts)
            if patt is not None and mol.HasSubstructMatch(patt):
                return label
        return "other"

    @classmethod
    def train_and_evaluate_qsar_logp_model(cls) -> dict[str, Any]:
        """
        Trains and evaluates a Ridge QSAR LogP regression model on an explicit train/test
        split of the experimental LogP reference dataset, returning genuine measured
        holdout test metrics (R2, RMSE, MAE) alongside Crippen performance on the same test set.
        """
        if cls._CACHED_QSAR_BUNDLE is not None:
            return cls._CACHED_QSAR_BUNDLE

        if not RDKIT_AVAILABLE or not (SKLEARN_AVAILABLE or NUMPY_AVAILABLE):
            return {
                "available": False,
                "message_ar": "مكتبات RDKit و scikit-learn أو NumPy مطلوبة لتدريب وتقييم نموذج QSAR.",
            }

        X_rows: list[list[float]] = []
        y_vals: list[float] = []
        crippen_vals: list[float] = []
        valid_records: list[dict[str, Any]] = []

        for entry in cls.EXPERIMENTAL_LOGP_BENCHMARK_DATASET:
            mol = Chem.MolFromSmiles(entry["smiles"])
            if mol is None:
                continue
            X_rows.append(cls._extract_qsar_features(mol))
            y_vals.append(float(entry["experimental_logp"]))
            crippen_vals.append(float(Descriptors.MolLogP(mol)))
            valid_records.append(entry)

        X_arr = np.asarray(X_rows, dtype=float)
        y_arr = np.asarray(y_vals, dtype=float)
        crippen_arr = np.asarray(crippen_vals, dtype=float)

        def _fit_predict_numpy(X_train: np.ndarray, y_train: np.ndarray) -> _RidgeNumpyPredictor:
            mu_v = np.mean(X_train, axis=0)
            std_v = np.std(X_train, axis=0, ddof=0)
            std_v[std_v == 0] = 1.0
            X_tr_s = (X_train - mu_v) / std_v
            y_mean_v = float(np.mean(y_train))
            y_c_v = y_train - y_mean_v
            alpha_v = 1.0
            beta_v = np.linalg.solve(
                X_tr_s.T @ X_tr_s + alpha_v * np.eye(X_tr_s.shape[1]),
                X_tr_s.T @ y_c_v,
            )
            return _RidgeNumpyPredictor(mu_v, std_v, beta_v, y_mean_v)

        def _r2_np(y_t: np.ndarray, y_p: np.ndarray) -> float:
            ss_res = float(np.sum((y_t - y_p) ** 2))
            ss_tot = float(np.sum((y_t - np.mean(y_t)) ** 2))
            return 1.0 - (ss_res / ss_tot) if ss_tot != 0 else 0.0

        n_samples = len(valid_records)
        perm = np.random.RandomState(42).permutation(n_samples)
        n_test_count = int(round(n_samples * 0.25))
        test_idx = perm[:n_test_count]
        train_idx = perm[n_test_count:]

        if SKLEARN_AVAILABLE:
            pipeline = Pipeline([
                ("scaler", StandardScaler()),
                ("ridge", Ridge(alpha=1.0)),
            ])
            pipeline.fit(X_arr[train_idx], y_arr[train_idx])

            y_pred_test = pipeline.predict(X_arr[test_idx])
            y_true_test = y_arr[test_idx]
            crippen_test = crippen_arr[test_idx]

            test_r2 = round(float(r2_score(y_true_test, y_pred_test)), 3)
            test_rmse = round(float(np.sqrt(mean_squared_error(y_true_test, y_pred_test))), 3)
            test_mae = round(float(mean_absolute_error(y_true_test, y_pred_test)), 3)

            crippen_test_r2 = round(float(r2_score(y_true_test, crippen_test)), 3)
            crippen_test_rmse = round(float(np.sqrt(mean_squared_error(y_true_test, crippen_test))), 3)
            crippen_test_mae = round(float(mean_absolute_error(y_true_test, crippen_test)), 3)
            model_predictor = pipeline
        else:
            model_predictor = _fit_predict_numpy(X_arr[train_idx], y_arr[train_idx])
            y_pred_test = model_predictor.predict(X_arr[test_idx])
            y_true_test = y_arr[test_idx]
            crippen_test = crippen_arr[test_idx]

            test_r2 = round(_r2_np(y_true_test, y_pred_test), 3)
            test_rmse = round(float(np.sqrt(np.mean((y_true_test - y_pred_test) ** 2))), 3)
            test_mae = round(float(np.mean(np.abs(y_true_test - y_pred_test))), 3)

            crippen_test_r2 = round(_r2_np(y_true_test, crippen_test), 3)
            crippen_test_rmse = round(float(np.sqrt(np.mean((y_true_test - crippen_test) ** 2))), 3)
            crippen_test_mae = round(float(np.mean(np.abs(y_true_test - crippen_test))), 3)

        # 5-Fold Cross-Validation across all 120 compounds
        cv_folds = 5
        fold_slices = np.array_split(perm, cv_folds)
        cv_preds = np.zeros_like(y_arr)
        all_idx_set = set(range(n_samples))
        for val_idx_arr in fold_slices:
            val_set = set(int(x) for x in val_idx_arr)
            tr_idx_arr = np.array([i for i in range(n_samples) if i not in val_set], dtype=int)
            fold_model = _fit_predict_numpy(X_arr[tr_idx_arr], y_arr[tr_idx_arr])
            cv_preds[val_idx_arr] = fold_model.predict(X_arr[val_idx_arr])

        cv_r2 = round(_r2_np(y_arr, cv_preds), 3)
        cv_rmse = round(float(np.sqrt(np.mean((y_arr - cv_preds) ** 2))), 3)
        cv_mae = round(float(np.mean(np.abs(y_arr - cv_preds))), 3)

        # Leave-One-Chemical-Class-Out validation: the model never sees the class it predicts.
        # This removes near-duplicate leakage between homologous series (e.g. 1-butanol vs 1-pentanol).
        class_labels = np.array(
            [cls._chemical_class(Chem.MolFromSmiles(rec["smiles"])) for rec in valid_records]
        )
        loco_preds = np.zeros_like(y_arr)
        group_cv_by_class: list[dict[str, Any]] = []
        for grp in sorted(set(class_labels.tolist())):
            val_mask = class_labels == grp
            tr_ids = np.where(~val_mask)[0]
            va_ids = np.where(val_mask)[0]
            if len(tr_ids) < 10:
                loco_preds[va_ids] = float(np.mean(y_arr))
            else:
                grp_model = _fit_predict_numpy(X_arr[tr_ids], y_arr[tr_ids])
                loco_preds[va_ids] = grp_model.predict(X_arr[va_ids])
            group_cv_by_class.append({
                "class": grp,
                "n": int(len(va_ids)),
                "qsar_rmse": round(float(np.sqrt(np.mean((y_arr[va_ids] - loco_preds[va_ids]) ** 2))), 3),
                "crippen_rmse": round(float(np.sqrt(np.mean((y_arr[va_ids] - crippen_arr[va_ids]) ** 2))), 3),
            })

        group_cv_r2 = round(_r2_np(y_arr, loco_preds), 3)
        group_cv_rmse = round(float(np.sqrt(np.mean((y_arr - loco_preds) ** 2))), 3)
        group_cv_mae = round(float(np.mean(np.abs(y_arr - loco_preds))), 3)
        crippen_all_r2 = round(_r2_np(y_arr, crippen_arr), 3)
        crippen_all_rmse = round(float(np.sqrt(np.mean((y_arr - crippen_arr) ** 2))), 3)
        crippen_all_mae = round(float(np.mean(np.abs(y_arr - crippen_arr))), 3)

        test_set_table = []
        for idx_i, pred_val in zip(test_idx, y_pred_test):
            rec = valid_records[int(idx_i)]
            exp_v = float(rec["experimental_logp"])
            crip_v = round(float(crippen_arr[int(idx_i)]), 2)
            pred_r = round(float(pred_val), 2)
            test_set_table.append({
                "المركب (Compound)": f"{rec['name_ar']} ({rec['name_en']})",
                "LogP التجريبي الموثق": exp_v,
                "تنبؤ نموذج QSAR": pred_r,
                "حساب Crippen cLogP": crip_v,
                "خطأ QSAR المطلق": round(abs(pred_r - exp_v), 2),
            })

        mw_train = X_arr[train_idx, 0]
        heavy_train = X_arr[train_idx, 4]

        bundle = {
            "available": True,
            "pipeline": model_predictor,
            "model_name": "Ridge Regression QSAR Model (StandardScaler + Ridge alpha=1.0, 15 Descriptors)",
            "descriptor_names": list(cls._QSAR_FEATURE_NAMES),
            "n_total": int(len(valid_records)),
            "n_train": int(len(train_idx)),
            "n_test": int(len(test_idx)),
            "split_method": (
                f"{len(train_idx)} Train / {len(test_idx)} Holdout Test (random permutation, seed=42) "
                "+ 5-Fold CV + Leave-One-Chemical-Class-Out"
            ),
            "test_r2": test_r2,
            "test_rmse": test_rmse,
            "test_mae": test_mae,
            "cv_folds": cv_folds,
            "cv_r2": cv_r2,
            "cv_rmse": cv_rmse,
            "cv_mae": cv_mae,
            "group_cv_r2": group_cv_r2,
            "group_cv_rmse": group_cv_rmse,
            "group_cv_mae": group_cv_mae,
            "group_cv_by_class": group_cv_by_class,
            "class_counts": {k: int(v) for k, v in zip(*np.unique(class_labels, return_counts=True))},
            "crippen_all_r2": crippen_all_r2,
            "crippen_all_rmse": crippen_all_rmse,
            "crippen_all_mae": crippen_all_mae,
            "crippen_test_r2": crippen_test_r2,
            "crippen_test_rmse": crippen_test_rmse,
            "crippen_test_mae": crippen_test_mae,
            "test_set_table": test_set_table,
            "applicability_domain": {
                "mw_min": round(float(mw_train.min()), 2),
                "mw_max": round(float(mw_train.max()), 2),
                "heavy_min": int(heavy_train.min()),
                "heavy_max": int(heavy_train.max()),
            },
            "dataset_citation": cls.DATASET_CITATION,
        }
        cls._CACHED_QSAR_BUNDLE = bundle
        return bundle

    @classmethod
    def evaluate_logp_prediction(cls, smiles: str) -> dict[str, Any]:
        """
        Explicitly separates deterministic Wildman-Crippen cLogP calculation from an
        actually trained and holdout-evaluated QSAR LogP regression model.
        """
        det = cls.calculate_deterministic_clogp(smiles)
        if not det.get("available"):
            return {
                "available": False,
                "message_ar": det.get("message_ar", "تعذر حساب معامل LogP لهذا المركب."),
            }

        qsar_bundle = cls.train_and_evaluate_qsar_logp_model()
        qsar_eval: dict[str, Any] = {"model_available": False}

        if qsar_bundle.get("available"):
            mol = Chem.MolFromSmiles(smiles)
            feats = np.asarray([cls._extract_qsar_features(mol)], dtype=float)
            qsar_pred = round(float(qsar_bundle["pipeline"].predict(feats)[0]), 3)

            mw = float(feats[0, 0])
            heavy = int(feats[0, 4])
            ad = qsar_bundle["applicability_domain"]
            in_ad = (ad["mw_min"] <= mw <= ad["mw_max"]) and (ad["heavy_min"] <= heavy <= ad["heavy_max"])

            ad_note = (
                f"المركب ضمن مجال انطباق مجموعة التدريب ({qsar_bundle['n_train']} مركباً، الكتلة المولية MW بين {ad['mw_min']} و {ad['mw_max']} g/mol، وعدد الذرات الثقيلة بين {ad['heavy_min']} و {ad['heavy_max']})."
                if in_ad
                else f"تنبيه مجال الانطباق: المركب خارج نطاق تدريب النموذج (MW = {mw:.1f} g/mol، الذرات الثقيلة = {heavy} مقارنة بنطاق التدريب [{ad['mw_min']}-{ad['mw_max']}] g/mol)، لذا قد تنخفض موثوقية التنبؤ بالاستقراء."
            )

            qsar_eval = {
                "model_available": True,
                "qsar_predicted_logp": qsar_pred,
                "in_applicability_domain": in_ad,
                "applicability_note_ar": ad_note,
                "model_name": qsar_bundle["model_name"],
                "descriptor_names": qsar_bundle["descriptor_names"],
                "n_total": qsar_bundle["n_total"],
                "n_train": qsar_bundle["n_train"],
                "n_test": qsar_bundle["n_test"],
                "split_method": qsar_bundle["split_method"],
                "test_r2": qsar_bundle["test_r2"],
                "test_rmse": qsar_bundle["test_rmse"],
                "test_mae": qsar_bundle["test_mae"],
                "cv_folds": qsar_bundle["cv_folds"],
                "cv_r2": qsar_bundle["cv_r2"],
                "cv_rmse": qsar_bundle["cv_rmse"],
                "cv_mae": qsar_bundle["cv_mae"],
                "group_cv_r2": qsar_bundle["group_cv_r2"],
                "group_cv_rmse": qsar_bundle["group_cv_rmse"],
                "group_cv_mae": qsar_bundle["group_cv_mae"],
                "group_cv_by_class": qsar_bundle["group_cv_by_class"],
                "crippen_all_r2": qsar_bundle["crippen_all_r2"],
                "crippen_all_rmse": qsar_bundle["crippen_all_rmse"],
                "crippen_all_mae": qsar_bundle["crippen_all_mae"],
                "crippen_test_r2": qsar_bundle["crippen_test_r2"],
                "crippen_test_rmse": qsar_bundle["crippen_test_rmse"],
                "crippen_test_mae": qsar_bundle["crippen_test_mae"],
                "test_set_table": qsar_bundle["test_set_table"],
                "dataset_citation": qsar_bundle["dataset_citation"],
            }

        clogp_val = det["calculated_clogp"]
        exp_val = det.get("experimental_logp")
        exp_sentence = (
            f" القيمة التجريبية المرجعية المسجلة لهذا المركب هي {exp_val} (بفارق مطلق {det['abs_error_vs_experimental']} عن cLogP)."
            if exp_val is not None
            else " لا تتوفر قيمة تجريبية مسجلة لهذا المركب في مجموعة البيانات المرجعية المحلية البالغ حجمها " + str(len(cls.EXPERIMENTAL_LOGP_BENCHMARK_DATASET)) + " مركباً."
        )

        if qsar_eval.get("model_available"):
            pedagogical_note = (
                f"1) الحساب الحتمي (RDKit Wildman-Crippen cLogP): القيمة = {clogp_val} وهي ناتجة عن جمع مساهمات الذرات بدون تدريب إحصائي أو فاصل ثقة لكل مركب. "
                f"2) نموذج التنبؤ المختبر فعلياً (Ridge QSAR على {qsar_eval['n_total']} مركباً): التنبؤ = {qsar_eval['qsar_predicted_logp']}؛ "
                f"جرى تدريبه على {qsar_eval['n_train']} مركباً واختباره على {qsar_eval['n_test']} مركباً مستقلاً "
                f"فحقق على مجموعة الاختبار المستقلة R² = {qsar_eval['test_r2']} و RMSE = {qsar_eval['test_rmse']} "
                f"(وفي التحقق المتقاطع الخماسي العشوائي: R² = {qsar_eval['cv_r2']} و RMSE = {qsar_eval['cv_rmse']}). "
                f"أما عند حجب فئة كيميائية كاملة عن التدريب (Leave-One-Class-Out، وهو الاختبار الأصعب والأقرب للواقع) فينخفض الأداء إلى "
                f"R² = {qsar_eval['group_cv_r2']} و RMSE = {qsar_eval['group_cv_rmse']}، "
                f"مقابل Wildman-Crippen على كل المركبات: R² = {qsar_eval['crippen_all_r2']} و RMSE = {qsar_eval['crippen_all_rmse']}.{exp_sentence}"
            )
        else:
            pedagogical_note = (
                f"الحساب الحتمي (RDKit Wildman-Crippen cLogP) = {clogp_val}.{exp_sentence}"
            )

        return {
            "available": True,
            "calculated_clogp": clogp_val,
            "qsar_predicted_logp": qsar_eval.get("qsar_predicted_logp"),
            "experimental_logp": exp_val,
            "deterministic_clogp": det,
            "qsar_model_evaluation": qsar_eval,
            "methodology_ar": (
                "فصل منهجي بين: (أ) حساب cLogP الحتمي بطريقة التجزئة الذرية Wildman-Crippen، "
                "و (ب) نموذج Ridge QSAR المدرَّب والمقيَّم فعلياً "
                + (
                    f"على {qsar_eval['n_total']} مركباً مرجعياً ({qsar_eval['n_train']} للتدريب / {qsar_eval['n_test']} للاختبار المحجوز) "
                    "مع تحقق متقاطع خماسي وتحقق مُجمَّع حسب الفئة الكيميائية (Leave-One-Class-Out)."
                    if qsar_eval.get("model_available")
                    else "(غير متاح في هذه البيئة)."
                )
            ),
            "limitations_ar": (
                "حساب Wildman-Crippen cLogP قد يختلف عن القياس التجريبي عند وجود روابط هيدروجينية داخلية أو تأين يعتمد على pH (مثل الأحماض والقواعد والأيونات الثنائية). "
                + (
                    f"أما نموذج QSAR المحلي فمدرَّب على {qsar_eval['n_train']} مركباً ومختبر على {qsar_eval['n_test']} مركباً، "
                    "ومجموعة البيانات تضم سلاسل متجانسة قريبة بنيوياً، لذا فالتقسيم العشوائي متفائل؛ "
                    "والمقياس الأصدق هو التحقق بحجب الفئة الكيميائية. وقد لا يتفوق النموذج على Crippen، "
                    "ويجب مراعاة مجال الانطباق (Applicability Domain) قبل أي تعميم."
                    if qsar_eval.get("model_available")
                    else "نموذج QSAR المحلي غير متاح في هذه البيئة."
                )
            ),
            "pedagogical_note_ar": pedagogical_note,
        }

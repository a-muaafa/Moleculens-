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
    from sklearn.linear_model import Ridge
    from sklearn.metrics import mean_absolute_error, mean_squared_error, r2_score
    from sklearn.model_selection import train_test_split
    from sklearn.pipeline import Pipeline
    from sklearn.preprocessing import StandardScaler

    SKLEARN_AVAILABLE = True
except ImportError:
    SKLEARN_AVAILABLE = False

from search.local_catalog import LOCAL_SCIENTIFIC_CATALOG


class MolecularIntelligence:
    """Intelligent molecular similarity, virtual editing, and QSAR modeling engine."""

    # Curated dataset of 32 catalog compounds with documented experimental octanol-water LogP values
    # Sources: Sangster (1989) J. Phys. Chem. Ref. Data 18(3):1111-1227; Hansch, Leo & Hoekman (1995) ACS.
    EXPERIMENTAL_LOGP_BENCHMARK_DATASET: list[dict[str, Any]] = [
        {"name_en": "Methanol", "name_ar": "الميثانول", "smiles": "CO", "experimental_logp": -0.77},
        {"name_en": "Ethanol", "name_ar": "الإيثانول", "smiles": "CCO", "experimental_logp": -0.31},
        {"name_en": "Isopropanol", "name_ar": "الأيزوبروبانول", "smiles": "CC(C)O", "experimental_logp": 0.05},
        {"name_en": "1-Propanol", "name_ar": "البروبانول", "smiles": "CCCO", "experimental_logp": 0.25},
        {"name_en": "Glycerol", "name_ar": "الجليسرول", "smiles": "OCC(O)CO", "experimental_logp": -1.76},
        {"name_en": "Formic Acid", "name_ar": "حمض الفورميك", "smiles": "O=CO", "experimental_logp": -0.54},
        {"name_en": "Acetic Acid", "name_ar": "حمض الأسيتيك", "smiles": "CC(=O)O", "experimental_logp": -0.17},
        {"name_en": "Lactic Acid", "name_ar": "حمض اللاكتيك", "smiles": "CC(O)C(=O)O", "experimental_logp": -0.72},
        {"name_en": "Citric Acid", "name_ar": "حمض الستريك", "smiles": "O=C(O)CC(O)(CC(=O)O)C(=O)O", "experimental_logp": -1.64},
        {"name_en": "Formaldehyde", "name_ar": "الفورمالديهايد", "smiles": "C=O", "experimental_logp": 0.35},
        {"name_en": "Acetaldehyde", "name_ar": "الأسيتالديهايد", "smiles": "CC=O", "experimental_logp": -0.34},
        {"name_en": "Acetone", "name_ar": "الأسيتون", "smiles": "CC(C)=O", "experimental_logp": -0.24},
        {"name_en": "Butanone", "name_ar": "البيوتانون", "smiles": "CCC(C)=O", "experimental_logp": 0.29},
        {"name_en": "Methyl Formate", "name_ar": "فورمات الميثيل", "smiles": "COC=O", "experimental_logp": 0.03},
        {"name_en": "Ethyl Acetate", "name_ar": "أسيتات الإيثيل", "smiles": "CCOC(C)=O", "experimental_logp": 0.73},
        {"name_en": "Methyl Acetate", "name_ar": "أسيتات الميثيل", "smiles": "COC(C)=O", "experimental_logp": 0.18},
        {"name_en": "Dimethyl Ether", "name_ar": "ثنائي ميثيل إيثر", "smiles": "COC", "experimental_logp": 0.10},
        {"name_en": "Urea", "name_ar": "اليوريا", "smiles": "NC(N)=O", "experimental_logp": -2.11},
        {"name_en": "D-Glucose", "name_ar": "الجلوكوز", "smiles": "OC[C@H]1OC(O)[C@H](O)[C@@H](O)[C@@H]1O", "experimental_logp": -3.24},
        {"name_en": "Benzene", "name_ar": "البنزين العطري", "smiles": "c1ccccc1", "experimental_logp": 2.13},
        {"name_en": "Toluene", "name_ar": "التولوين", "smiles": "Cc1ccccc1", "experimental_logp": 2.73},
        {"name_en": "Phenol", "name_ar": "الفينول", "smiles": "Oc1ccccc1", "experimental_logp": 1.46},
        {"name_en": "Aniline", "name_ar": "الأنيلين", "smiles": "Nc1ccccc1", "experimental_logp": 0.90},
        {"name_en": "Nitrobenzene", "name_ar": "النيتروبنزين", "smiles": "O=[N+]([O-])c1ccccc1", "experimental_logp": 1.85},
        {"name_en": "Benzoic Acid", "name_ar": "حمض البنزويك", "smiles": "O=C(O)c1ccccc1", "experimental_logp": 1.87},
        {"name_en": "Salicylic Acid", "name_ar": "حمض الساليسيليك", "smiles": "O=C(O)c1ccccc1O", "experimental_logp": 2.26},
        {"name_en": "Benzaldehyde", "name_ar": "البنزالديهايد", "smiles": "O=Cc1ccccc1", "experimental_logp": 1.48},
        {"name_en": "Naphthalene", "name_ar": "النفثالين", "smiles": "c1ccc2ccccc2c1", "experimental_logp": 3.30},
        {"name_en": "Aspirin", "name_ar": "الأسبرين", "smiles": "CC(=O)Oc1ccccc1C(=O)O", "experimental_logp": 1.19},
        {"name_en": "Paracetamol", "name_ar": "الباراسيتامول", "smiles": "CC(=O)Nc1ccc(O)cc1", "experimental_logp": 0.46},
        {"name_en": "Ibuprofen", "name_ar": "الإيبوبروفين", "smiles": "CC(C)Cc1ccc(C(C)C(=O)O)cc1", "experimental_logp": 3.97},
        {"name_en": "Caffeine", "name_ar": "الكافيين", "smiles": "Cn1c(=O)c2c(ncn2C)n(C)c1=O", "experimental_logp": -0.07},
    ]

    DATASET_CITATION = (
        "Sangster, J. (1989) J. Phys. Chem. Ref. Data, 18(3), 1111-1227; "
        "Hansch, C., Leo, A., and Hoekman, D. (1995) Exploring QSAR, ACS."
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
            "smarts_pattern": "[CH3]",
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

    @classmethod
    def train_and_evaluate_qsar_logp_model(cls) -> dict[str, Any]:
        """
        Trains and evaluates a Ridge QSAR LogP regression model on an explicit train/test
        split of the 32-compound experimental LogP dataset, returning genuine measured
        holdout test metrics (R2, RMSE, MAE) alongside Crippen performance on the same test set.
        """
        if cls._CACHED_QSAR_BUNDLE is not None:
            return cls._CACHED_QSAR_BUNDLE

        if not RDKIT_AVAILABLE or not SKLEARN_AVAILABLE:
            return {
                "available": False,
                "message_ar": "مكتبات RDKit و scikit-learn مطلوبة لتدريب وتقييم نموذج QSAR.",
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

        indices = np.arange(len(valid_records))
        train_idx, test_idx = train_test_split(indices, test_size=0.25, random_state=42)

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
            "pipeline": pipeline,
            "model_name": "Ridge Regression QSAR Model (StandardScaler + Ridge alpha=1.0)",
            "descriptor_names": list(cls._QSAR_FEATURE_NAMES),
            "n_total": int(len(valid_records)),
            "n_train": int(len(train_idx)),
            "n_test": int(len(test_idx)),
            "split_method": "75% Train / 25% Holdout Test Split (random_state=42)",
            "test_r2": test_r2,
            "test_rmse": test_rmse,
            "test_mae": test_mae,
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
                f"المركب ضمن نطاق الكتلة والذرات لمجموعة التدريب (MW بين {ad['mw_min']} و {ad['mw_max']} g/mol)."
                if in_ad
                else f"تنبيه مجال الانطباق: المركب خارج نطاق تدريب النموذج الصغير (MW = {mw:.1f} مقارنة بنطاق التدريب [{ad['mw_min']}-{ad['mw_max']}])، لذا قد تنخفض موثوقية التنبؤ بالاستقراء."
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
            else " لا تتوفر قيمة تجريبية مسجلة لهذا المركب في مجموعة الاختبار المحلية البالغ حجمها 32 مركباً."
        )

        if qsar_eval.get("model_available"):
            pedagogical_note = (
                f"1) الحساب الحتمي (RDKit Wildman-Crippen cLogP): القيمة = {clogp_val} وهي ناتجة عن جمع مساهمات الذرات بدون تدريب إحصائي أو فاصل ثقة لكل مركب. "
                f"2) نموذج التنبؤ المختبر فعلياً (Ridge QSAR): التنبؤ = {qsar_eval['qsar_predicted_logp']}؛ "
                f"جرى تدريبه على {qsar_eval['n_train']} مركباً واختباره على {qsar_eval['n_test']} مركبات مستقلة "
                f"فحقق على مجموعة الاختبار R² = {qsar_eval['test_r2']} و RMSE = {qsar_eval['test_rmse']}.{exp_sentence}"
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
                "و (ب) نموذج Ridge QSAR المدرَّب والمقيَّم فعلياً على مجموعة اختبار مستقلة."
            ),
            "limitations_ar": (
                "حساب Wildman-Crippen cLogP قد يختلف عن القياس التجريبي عند وجود روابط هيدروجينية داخلية أو تأين يعتمد على pH (مثل الأحماض والقواعد والأيونات الثنائية). "
                "أما نموذج QSAR المحلي فمدرَّب على 24 مركباً ومختبر على 8 مركبات مرجعية فقط لتوضيح منهجية التقييم الحقيقي، ولا يُعمَّم كأداة حسم دوائي على مركبات خارج مجال انطباقه."
            ),
            "pedagogical_note_ar": pedagogical_note,
        }

from typing import Any
from chemistry_engine.formula_analyzer import analyze_formula_composition, format_chemical_formula_display
from chemistry_engine.structure_analyzer import render_molecule_2d_png, verify_substance_structure_identity
from config.settings import (
    ENTITY_CHEMICAL_COMPOUND,
    OP_SCIENTIFIC_CALCULATION,
    OP_SCIENTIFIC_RETRIEVAL,
)
from core.schemas import build_provenance_tag
from education.explanation_engine import build_entity_explanation

SUBCLASS_SAFETY_NOTES: dict[str, dict[str, str]] = {
    "acid": {
        "hazard_ar": "مادة حمضية أكالة أو مهيجة (Corrosive / Irritant — GHS05/GHS07): تسبب حروقاً كيميائية للجلد والعينين عند التلامس المباشر.",
        "ppe_ar": "معطف المختبر، قفازات مقاومة للأحماض، نظارات وقاية محكمة، والعمل داخل خزانة سحب الغازات (Fume Hood) للأحماض المركزة.",
        "handling_ar": "قاعدة ذهبية في المختبر: يُضاف الحمض ببطء إلى الماء (Add Acid to Water) مع التحريك المستمر لتجنب التناثر الحراري.",
    },
    "carboxylic_acid": {
        "hazard_ar": "حمض عضوي (مهيج في التراكيز العالية مثل حمض الخليك الثلجي Glacial Acetic Acid — GHS05/GHS07).",
        "ppe_ar": "معطف المختبر، قفازات النتريل، ونظارات السلامة.",
        "handling_ar": "يُتداول في مكان جيد التهوية بعيداً عن القواعد القوية والمؤكسدات.",
    },
    "base": {
        "hazard_ar": "مادة قاعدية كاوية أو مهيجة (Caustic / Alkaline — GHS05): تسبب تآكل الأنسجة وتصبن الدهون الجلدية.",
        "ppe_ar": "معطف المختبر، قفازات النتريل، ونظارات الوقاية الكيميائية.",
        "handling_ar": "تذوب القواعد القوية (مثل NaOH و KOH) في الماء بتفاعل طارد للحرارة؛ تُضاف تدريجياً مع التبريد والتحريك.",
    },
    "alcohol": {
        "hazard_ar": "سائل عضوي قابل للاشتعال (Flammable Liquid — GHS02) وقد يسبب تهيج العينين.",
        "ppe_ar": "معطف المختبر، قفازات النتريل، ونظارات السلامة.",
        "handling_ar": "يُبعد تماماً عن مصادر اللهب المكشوف والأسطح الساخنة؛ يُسخن فقط عبر حمام مائي مضبوط.",
    },
    "aromatic_compound": {
        "hazard_ar": "مركب عطري متطاير أو مهيج (GHS02/GHS07/GHS08): يلزم الحذر من استنشاق الأبخرة العضوية.",
        "ppe_ar": "معطف المختبر، قفازات مناسبة للمذيبات العضوية، والعمل حصراً تحت ساحبة الغازات (Fume Hood).",
        "handling_ar": "يُحفظ في عبوات زجاجية داكنة محكمة الإغلاق ويُتخلص من بقاياه في حاويات النفايات العضوية المخصصة.",
    },
    "salt": {
        "hazard_ar": "ملح كيميائي (تختلف الخطورة بحسب الأيونات؛ الأملاح المؤكسدة أو أملاح المعادن الثقيلة تتطلب حذراً خاصاً).",
        "ppe_ar": "معطف المختبر، قفازات النتريل، ونظارات السلامة.",
        "handling_ar": "تُوزن الأملاح الصلبة باستخدام ملعقة نظيفة وجافة (Spatula) وتُحفظ العبوات مغلقة لمنع التميع (Hygroscopicity).",
    },
}


def build_compound_scientific_profile(entity: dict[str, Any]) -> dict[str, Any]:
    formula = entity.get("symbol_or_formula", "")
    formatted_formula = format_chemical_formula_display(formula)
    smiles = entity.get("canonical_smiles", "")
    composition = analyze_formula_composition(formula)
    molar_mass = entity.get("molar_mass") or composition.get("molar_mass")
    struct = entity.get("structure_details", {})
    phys = entity.get("physical_properties", {})
    chem = entity.get("chemical_properties", {})
    apps = entity.get("applications", {})
    explanation = build_entity_explanation(entity, composition=composition)
    image_2d = render_molecule_2d_png(smiles) if smiles else None
    cid = entity.get("identifiers", {}).get("pubchem_cid") or entity.get("pubchem_cid")
    identity_check = entity.get("identity_verification") or verify_substance_structure_identity(entity)
    subclass = entity.get("subclass", "")
    lab_safety = SUBCLASS_SAFETY_NOTES.get(
        subclass,
        {
            "hazard_ar": "مادة كيميائية مخبرية: يُلتزم بقواعد السلامة القياسية وتجنب الملامسة المباشرة أو الاستنشاق.",
            "ppe_ar": "معطف المختبر (Lab Coat)، قفازات واقية (Gloves)، ونظارات السلامة (Safety Goggles).",
            "handling_ar": "راجع بطاقة بيانات السلامة المادية (MSDS / SDS) قبل الاستخدام في المختبر.",
        },
    )

    return {
        "profile_kind": ENTITY_CHEMICAL_COMPOUND,
        "has_molecular_structure": True,
        "identity_verification": identity_check,
        "lab_safety": lab_safety,
        "basic_info": {
            "name_en": entity.get("name_en", ""),
            "name_ar": entity.get("name_ar", ""),
            "iupac_name": entity.get("iupac_name", ""),
            "formula": formula,
            "formatted_formula": formatted_formula,
            "molar_mass": f"{molar_mass} g/mol" if molar_mass else "غير متوفر",
            "molar_mass_value": molar_mass,
            "entity_type_label_ar": entity.get("entity_type_label_ar", "مركب كيميائي"),
            "entity_type_label_en": entity.get("entity_type_label_en", "Chemical Compound"),
        },
        "composition": {
            "summary_ar": composition.get("summary_ar", ""),
            "total_atoms": composition.get("total_atoms", 0),
            "components": composition.get("components", []),
            "calculation_provenance": build_provenance_tag(
                OP_SCIENTIFIC_CALCULATION,
                sources=["IUPAC Atomic Weights", "Formula Analyzer"],
                notes_ar="تم حساب المكونات الذرية والنسب المئوية الكتلية والكتلة المولية (Molar Mass) حتمياً من الصيغة الجزيئية.",
            ),
        },
        "structure": {
            "molecular_geometry": struct.get("molecular_geometry", "بنية جزيئية محددة وفق تكافؤ الذرات"),
            "hybridization": struct.get("hybridization", "غير محدد"),
            "chemical_bonds": struct.get("chemical_bonds", "روابط كيميائية بين الذرات المكونة"),
            "dipole_moment": struct.get("dipole_moment", "غير محدد"),
            "representation_3d_note": struct.get(
                "representation_3d_note",
                f"يتوفر التمثيل الفراغي ثلاثي الأبعاد عبر PubChem 3D Conformer (CID: {cid})." if cid else "متاح عبر محرك البنية الجزيئية.",
            ),
            "pubchem_3d_url": f"https://pubchem.ncbi.nlm.nih.gov/compound/{cid}#section=3D-Conformer" if cid else None,
            "image_2d_png": image_2d,
        },
        "properties": {
            "physical_state": phys.get("physical_state", "غير محدد"),
            "density": phys.get("density", "غير محدد"),
            "melting_point": phys.get("melting_point", "غير محدد"),
            "boiling_point": phys.get("boiling_point", "غير محدد"),
            "chemical_properties_summary": chem.get("acidity_basicity") or chem.get("reactivity", ""),
            "reactivity": chem.get("reactivity", ""),
        },
        "applications": {
            "daily": apps.get("daily", "تطبيقات عامة في الحياة اليومية."),
            "industrial": apps.get("industrial", "تطبيقات صناعية وكيميائية."),
            "scientific": apps.get("scientific", "تطبيقات بحثية ومختبرية."),
        },
        "explanation": explanation,
        "identifiers": entity.get("identifiers", {}),
        "sources": entity.get("sources", ["PubChem"]),
        "provenance": build_provenance_tag(
            OP_SCIENTIFIC_RETRIEVAL,
            sources=entity.get("sources", ["PubChem"]),
            notes_ar="المعلومات المرجعية والخصائص الفيزيائية مسترجعة من المصادر العلمية المعتمدة.",
        ),
        "reaction_lab_entity": {
            "key": entity.get("key", ""),
            "name_ar": entity.get("name_ar", ""),
            "name_en": entity.get("name_en", ""),
            "symbol_or_formula": formula,
            "formatted_formula": formatted_formula,
            "molar_mass": molar_mass,
        },
    }

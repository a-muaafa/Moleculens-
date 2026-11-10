from typing import Any
from config.settings import (
    ENTITY_BIOMOLECULE_PROTEIN,
    ENTITY_CHEMICAL_COMPOUND,
    ENTITY_CHEMICAL_ELEMENT,
    ENTITY_DRUG_SUBSTANCE,
)


def build_entity_explanation(entity: dict[str, Any], composition: dict[str, Any] | None = None) -> dict[str, str]:
    notes = entity.get("education_notes", {})
    student_ar = notes.get("student_level_ar", "")
    advanced_ar = notes.get("advanced_level_ar", "")

    name_ar = entity.get("name_ar", "")
    name_en = entity.get("name_en", "")
    sym = entity.get("symbol_or_formula", "")
    entity_type = entity.get("entity_type", ENTITY_CHEMICAL_COMPOUND)

    if not student_ar:
        if entity_type == ENTITY_CHEMICAL_ELEMENT:
            z = entity.get("identifiers", {}).get("atomic_number") or entity.get("atomic_details", {}).get("atomic_number", "")
            student_ar = f"يمثل {name_ar} ({sym}) عنصرا كيميائيا نقيا في الجدول الدوري عدده الذري {z}، وتحدد خصائصه الذرية سلوكه الفيزيائي والكيميائي."
        elif entity_type == ENTITY_DRUG_SUBSTANCE:
            student_ar = f"يمثل {name_ar} ({name_en}) مادة دوائية ذات صيغة جزيئية ({sym}) صممت للارتباط بهدف حيوي محدد في الجسم لإحداث تأثير علاجي مدروس."
        elif entity_type == ENTITY_BIOMOLECULE_PROTEIN:
            student_ar = f"يمثل {name_ar} ({name_en}) جزيئا حيويا بروتينيا يؤدي وظيفة تنظيمية أو إنزيمية متخصصة داخل الكائن الحي."
        else:
            comp_str = composition.get("summary_ar", "") if composition else ""
            student_ar = (
                f"يتكون المركب الكيميائي {name_ar} ({sym}) من ({comp_str}) مرتبطة بروابط كيميائية محددة تمنحه خصائصه الفيزيائية والكيميائية المميزة."
            )

    if not advanced_ar:
        molar_mass = entity.get("molar_mass") or (composition.get("molar_mass") if composition else None)
        advanced_ar = (
            f"من المنظور الكيميائي المتقدم، يخضع الكيان {name_ar} ({sym}، الكتلة المولية: {molar_mass} g/mol) "
            "لقواعد الترابط الإلكتروني والديناميكا الحرارية التي تحدد استقراره وتفاعليته في الأوساط المختلفة."
        )

    return {
        "student_level_ar": student_ar,
        "advanced_level_ar": advanced_ar,
        "methodology_distinction_ar": (
            "تلتزم منصة MolecuLens بالفصل المنهجي بين ثلاث عمليات علمية: "
            "(1) البحث العلمي: استرجاع الحقائق والبيانات المسجلة من قواعد البيانات المعتمدة، "
            "(2) الحساب العلمي: تطبيق القوانين والمعادلات الرياضية والكيميائية الحتمية مثل حفظ الكتلة وحساب المولات، "
            "(3) التنبؤ العلمي: استخدام النماذج الحاسوبية لتقدير خاصية غير مقيسة مع توضيح حدود النموذج."
        ),
    }


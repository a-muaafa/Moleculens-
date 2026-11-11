from chemistry_engine.formula_analyzer import is_gaseous_at_stp
from config.settings import SUPPORTED_AMOUNT_UNITS


def validate_search_query(query: str) -> tuple[bool, str, str]:
    if query is None or not isinstance(query, str):
        return (
            False,
            "",
            "مدخل البحث غير صالح: يرجى إدخال اسم المادة (عربي/إنجليزي)، أو صيغتها الجزيئية (Molecular Formula مثل H₂O)، أو رمزها الكيميائي.",
        )
    cleaned = query.strip()
    if not cleaned:
        return (
            False,
            "",
            "حقل البحث فارغ: يرجى كتابة اسم العنصر أو المركب أو المادة الدوائية أو البروتين للبدء في الاستكشاف.",
        )
    if len(cleaned) > 250:
        return (
            False,
            "",
            "النص المدخل يتجاوز الطول المسموح (250 حرفاً): يرجى إدخال اسم علمي محدد أو صيغة كيميائية مختصرة.",
        )
    return True, cleaned, ""


def validate_reaction_inputs(reactants: list[dict]) -> tuple[bool, str]:
    if not reactants or not isinstance(reactants, list):
        return (
            False,
            "بيانات التفاعل ناقصة: يرجى اختيار المواد المتفاعلة (Reactants) وتحديد كمياتها من القوائم المغلقة قبل بدء التحليل.",
        )
    valid_items = [r for r in reactants if r and (r.get("substance") or r.get("formula"))]
    if len(valid_items) < 1:
        return (
            False,
            "لم يتم تحديد أي مادة متفاعلة: يرجى اختيار مادة متفاعلة واحدة على الأقل لإجراء الحسابات.",
        )
    for idx, item in enumerate(valid_items, start=1):
        sub = item.get("substance") or {}
        formula = sub.get("symbol_or_formula") or item.get("formula", f"#{idx}")
        name_ar = sub.get("name_ar") or formula
        amount = item.get("amount")
        unit = item.get("unit", "mol")
        if amount is None:
            return (
                False,
                f"الكمية مفقودة للمادة رقم {idx} ({name_ar} — {formula}): يرجى إدخال قيمة رقمية موجبة للكمية.",
            )
        try:
            val = float(amount)
        except (TypeError, ValueError):
            return (
                False,
                f"الكمية المدخلة للمادة رقم {idx} ({name_ar} — {formula}) غير رقمية: يرجى إدخال رقم صحيح أو عشري موجب.",
            )
        if val <= 0:
            return (
                False,
                f"قيمة غير صالحة للمادة رقم {idx} ({name_ar} — {formula}): يجب أن تكون الكمية أكبر من الصفر (> 0).",
            )
        if unit not in SUPPORTED_AMOUNT_UNITS:
            return (
                False,
                f"الوحدة '{unit}' للمادة رقم {idx} ({name_ar} — {formula}) غير مدعومة. الوحدات المعتمدة هي: {', '.join(SUPPORTED_AMOUNT_UNITS)}.",
            )
        if unit == "L_STP":
            rxn_state = item.get("state", "")
            if not is_gaseous_at_stp(formula, substance=sub, rxn_state=rxn_state):
                return (
                    False,
                    f"خطأ علمي في وحدة الحجم للمادة ({name_ar} — {formula}): وحدة الحجم الغازي في الظروف القياسية "
                    "(L_STP = 22.414 L/mol) تنطبق حصرياً على الغازات المثالية والحقيقية في الطور الغازي (Gas Phase)، "
                    "ولا يجوز تطبيقها لتحويل حجوم المواد الصلبة أو السائلة أو المحاليل المائية إلى مولات. "
                    "الحل المقترح: يرجى تغيير نوع الكمية إلى المولات (mol / mmol) أو الكتلة الوزنية (g / mg).",
                )
    return True, ""

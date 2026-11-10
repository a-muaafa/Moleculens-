from config.settings import SUPPORTED_AMOUNT_UNITS


def validate_search_query(query: str) -> tuple[bool, str, str]:
    if query is None or not isinstance(query, str):
        return False, "", "يرجى إدخال اسم المادة أو صيغتها الكيميائية أو رمزها العلمي."
    cleaned = query.strip()
    if not cleaned:
        return False, "", "حقل البحث فارغ؛ يرجى كتابة اسم العنصر أو المركب أو الدواء أو البروتين."
    if len(cleaned) > 250:
        return False, "", "النص المدخل طويل جدا؛ يرجى إدخال اسم أو صيغة علمية صالحة."
    return True, cleaned, ""


def validate_reaction_inputs(reactants: list[dict]) -> tuple[bool, str]:
    if not reactants or not isinstance(reactants, list):
        return False, "يرجى اختيار المواد المتفاعلة وتحديد كمياتها قبل تحليل التفاعل."
    valid_items = [r for r in reactants if r and (r.get("substance") or r.get("formula"))]
    if len(valid_items) < 1:
        return False, "يرجى تحديد مادة متفاعلة واحدة على الأقل."
    for idx, item in enumerate(valid_items, start=1):
        amount = item.get("amount")
        unit = item.get("unit", "mol")
        if amount is None:
            return False, f"يرجى تحديد الكمية للمادة المتفاعلة رقم {idx}."
        try:
            val = float(amount)
        except (TypeError, ValueError):
            return False, f"الكمية المدخلة للمادة رقم {idx} غير رقمية."
        if val <= 0:
            return False, f"يجب أن تكون كمية المادة رقم {idx} أكبر من الصفر."
        if unit not in SUPPORTED_AMOUNT_UNITS:
            return False, f"الوحدة '{unit}' للمادة رقم {idx} غير مدعومة. الوحدات المدعومة: {', '.join(SUPPORTED_AMOUNT_UNITS)}."
    return True, ""

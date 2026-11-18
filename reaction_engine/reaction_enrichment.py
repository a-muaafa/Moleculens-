"""
Scientific Enrichment Module for MolecuLens Reaction Lab.

Closes key scientific gaps in chemical reaction analysis and prediction:
1. Complete Ionic Equation, Net Ionic Equation, and Aqueous Spectator Ions (`build_ionic_equation_analysis`).
2. Chemical Equilibrium & Reversibility Analysis (`evaluate_reaction_equilibrium`).
3. Competing Reactions & Multi-Pathway Detection (`detect_competing_pathways`).
4. Laboratory Hazard & Safety Warnings + Anti-Operational-Recipe Guard (`evaluate_reaction_safety_and_hazards`).
5. Explicit Reaction Lab Coverage Boundaries (`get_reaction_lab_coverage_boundaries`).
"""

from __future__ import annotations

from typing import Any

from chemistry_engine.formula_analyzer import format_chemical_formula_display

# Strong electrolytes in aqueous solution that dissociate completely into ions
STRONG_ACIDS_AQUEOUS: set[str] = {"HCl", "HBr", "HI", "HNO3", "H2SO4", "HClO4"}
STRONG_BASES_AQUEOUS: set[str] = {"NaOH", "KOH", "LiOH", "RbOH", "CsOH", "Ca(OH)2", "Ba(OH)2", "Sr(OH)2"}
WEAK_ELECTROLYTES_AQUEOUS: set[str] = {
    "H2O", "CH3COOH", "HCOOH", "HF", "HCN", "H2S", "H2CO3", "H3PO4", "H2SO3",
    "H2C2O4", "HNO2", "H3BO3", "NH3", "NH4OH", "C6H8O7", "C3H6O3", "C7H6O2",
}
GASEOUS_SPECIES: set[str] = {
    "H2", "O2", "N2", "CO2", "CO", "SO2", "SO3", "NH3", "Cl2", "F2", "H2S", "NO", "NO2", "CH4", "C2H6", "C3H8", "C4H10", "C2H4", "C2H2",
}

REVERSIBLE_REACTION_IDS: set[str] = {
    "rxn_007",  # N2 + 3H2 <-> 2NH3 (Haber-Bosch)
    "rxn_013",  # 2SO2 + O2 <-> 2SO3 (Contact process)
    "rxn_016",  # CO2 + H2O <-> H2CO3
    "rxn_023",  # NH3 + H2O <-> NH4OH
    "rxn_071",  # CH3COOH + C2H5OH <-> CH3COOCH2CH3 + H2O
    "rxn_072",  # HCOOH + CH3OH <-> HCOOCH3 + H2O
    "rxn_073",  # CH3COOH + CH3OH <-> CH3COOCH3 + H2O
    "rxn_088",  # H2 + I2 <-> 2HI
}

TOXIC_OR_HAZARDOUS_GASES: dict[str, str] = {
    "Cl2": "غاز الكلور (Cl₂) سام وشديد التخريش للجهاز التنفسي والأغشية المخاطية.",
    "Br2": "أبخرة البروم (Br₂) سامة وأكّالة وتسبب حروقاً كيميائية شديدة.",
    "CO": "غاز أول أكسيد الكربون (CO) سام وقاتل بلا لون أو رائحة (يرتبط بالهيموجلوبين).",
    "SO2": "غاز ثاني أكسيد الكبريت (SO₂) خانق ومهيّج شديد للمسالك التنفسية.",
    "SO3": "ثالث أكسيد الكبريت (SO₃) يتفاعل بعنف مع الرطوبة مكوناً ضباب حمض الكبريتيك الأكّال.",
    "NO2": "غاز ثاني أكسيد النيتروجين (NO₂) بني محمر سام جداً للرئتين.",
    "NO": "أكسيد النيتريك (NO) غاز سام يتأكسد سريعاً في الهواء إلى NO₂.",
    "H2S": "غاز كبريتيد الهيدروجين (H₂S) شديد السمية وقابل للاشتعال.",
    "HCN": "سيانيد الهيدروجين (HCN) غاز قاتل فائق السمية يثبط التنفس الخلوي.",
    "NH3": "غاز الأمونيا (NH₃) كاوٍ ومهيّج للعينين والجهاز التنفسي.",
    "H2": "غاز الهيدروجين (H₂) شديد القابلية للاشتعال والانفجار عند اختلاطه بالهواء.",
    "CH4": "غاز الميثان (CH₄) سريع الاشتعال.",
}

CORROSIVE_SPECIES: dict[str, str] = {
    "H2SO4": "حمض الكبريتيك (H₂SO₄) أكّال جداً وتخفيفه بالماء طارد للحرارة بشدة.",
    "HNO3": "حمض النيتريك (HNO₃) حمض معدني أكّال ومؤكسد قوي.",
    "HCl": "حمض الهيدروكلوريك (HCl) حمض أكّال وأبخرته المركزّة مهيّجة.",
    "HF": "حمض الهيدروفلوريك (HF) شديد الخطورة يخترق الأنسجة ويسبب تسمماً جهازياً بالكالسيوم.",
    "HClO4": "حمض البيركلوريك (HClO₄) حمض فائق ومؤكسد قوي قد يشكل مخاليط متفجرة مع المواد العضوية.",
    "NaOH": "هيدروكسيد الصوديوم (NaOH) قاعدة كاوية تسبب حروقاً كيميائية شديدة للجلد والعيون.",
    "KOH": "هيدروكسيد البوتاسيوم (KOH) قاعدة كاوية شديدة القلوية.",
}

REACTIVE_WATER_METALS: set[str] = {"Li", "Na", "K", "Rb", "Cs"}
HEAVY_METAL_ELEMENTS: set[str] = {"Pb", "Hg", "Cd", "Ba", "Cr"}


def _charge_sup(charge: int) -> str:
    if charge == 0:
        return ""
    mag = abs(charge)
    sup = {1: "", 2: "²", 3: "³", 4: "⁴"}.get(mag, str(mag))
    return f"{sup}⁺" if charge > 0 else f"{sup}⁻"


def _dissociate_aqueous_species(
    formula: str,
    coeff: int,
    state: str,
    classify_fn: Any,
    is_insoluble_fn: Any,
) -> list[dict[str, Any]]:
    """
    Returns dissociated ion tokens if the species is a strong soluble electrolyte in aqueous phase;
    otherwise returns the undissociated molecular/solid/gas species token.
    """
    clean_state = (state or "").strip().lower().strip("()")
    if clean_state in {"s", "l", "g"} or formula in GASEOUS_SPECIES or formula in WEAK_ELECTROLYTES_AQUEOUS:
        phase = clean_state or ("l" if formula == "H2O" else ("g" if formula in GASEOUS_SPECIES else "s"))
        disp = f"{format_chemical_formula_display(formula)}({phase})"
        return [{"key": disp, "display": disp, "count": coeff, "is_ion": False}]

    info = classify_fn(formula)
    kind = info.get("kind", "")
    cation = info.get("cation", "")
    anion = info.get("anion", "")
    cat_chg = int(info.get("cation_charge", 0) or 0)
    an_chg = int(info.get("anion_charge", 0) or 0)
    cat_cnt = int(info.get("cation_count", 1) or 1)
    an_cnt = int(info.get("anion_count", 1) or 1)

    if not cation or not anion:
        phase = clean_state or "aq"
        disp = f"{format_chemical_formula_display(formula)}({phase})"
        return [{"key": disp, "display": disp, "count": coeff, "is_ion": False}]

    # Check if insoluble precipitate
    if cation != "H" and is_insoluble_fn(cation, anion):
        disp = f"{format_chemical_formula_display(formula)}(s)↓"
        return [{"key": disp, "display": disp, "count": coeff, "is_ion": False}]

    # Check if strong acid, strong base, or soluble salt
    is_strong_electrolyte = (
        formula in STRONG_ACIDS_AQUEOUS
        or formula in STRONG_BASES_AQUEOUS
        or kind in {"strong_acid", "strong_base", "salt"}
    ) and (formula not in WEAK_ELECTROLYTES_AQUEOUS)

    if not is_strong_electrolyte:
        phase = clean_state or "aq"
        disp = f"{format_chemical_formula_display(formula)}({phase})"
        return [{"key": disp, "display": disp, "count": coeff, "is_ion": False}]

    cat_disp = f"{format_chemical_formula_display(cation)}{_charge_sup(abs(cat_chg))}(aq)"
    an_disp = f"{format_chemical_formula_display(anion)}{_charge_sup(-abs(an_chg))}(aq)"
    return [
        {"key": cat_disp, "display": cat_disp, "count": coeff * cat_cnt, "is_ion": True},
        {"key": an_disp, "display": an_disp, "count": coeff * an_cnt, "is_ion": True},
    ]


def _format_ion_side(tokens: list[dict[str, Any]]) -> str:
    merged: list[dict[str, Any]] = []
    by_key: dict[str, dict[str, Any]] = {}
    for tok in tokens:
        k = tok["key"]
        if k in by_key:
            by_key[k]["count"] += int(tok["count"])
        else:
            entry = dict(tok)
            by_key[k] = entry
            merged.append(entry)
    parts = []
    for item in merged:
        c = int(item["count"])
        if c <= 0:
            continue
        pref = f"{c} " if c > 1 else ""
        parts.append(f"{pref}{item['display']}")
    return " + ".join(parts)


def build_ionic_equation_analysis(
    reactants: list[dict[str, Any]],
    products: list[dict[str, Any]],
    category: str = "",
    existing_net_ionic: str = "",
) -> dict[str, Any]:
    """
    Constructs the Complete Ionic Equation, Net Ionic Equation, and explicit list of
    aqueous Spectator Ions for ionic aqueous reactions.
    """
    from reaction_engine.prediction import classify_species, is_salt_insoluble

    ionic_categories = {"neutralization", "precipitation", "single_displacement"}
    has_aq = any((x.get("state", "") or "").lower().strip("()") == "aq" for x in list(reactants) + list(products))

    if category not in ionic_categories and not has_aq:
        msg = existing_net_ionic or "لا ينطبق (تفاعل جزيئي غير أيوني أو في طور غير مائي) — Non-ionic / molecular reaction."
        return {
            "is_ionic_aqueous": False,
            "complete_ionic_equation": msg,
            "net_ionic_equation": msg,
            "spectator_ions": [],
            "spectator_ions_summary_ar": "لا توجد أيونات متفرجة مائية (التفاعل جزيئي غير أيوني).",
        }

    lhs_tokens: list[dict[str, Any]] = []
    rhs_tokens: list[dict[str, Any]] = []

    for r in reactants:
        c = int(r.get("coefficient", 1) or 1)
        lhs_tokens.extend(_dissociate_aqueous_species(r["formula"], c, r.get("state", ""), classify_species, is_salt_insoluble))

    for p in products:
        c = int(p.get("coefficient", 1) or 1)
        rhs_tokens.extend(_dissociate_aqueous_species(p["formula"], c, p.get("state", ""), classify_species, is_salt_insoluble))

    has_any_ion = any(t["is_ion"] for t in lhs_tokens + rhs_tokens)
    if not has_any_ion:
        msg = existing_net_ionic or "لا ينطبق (معادلة جزيئية عامة) — Molecular equation."
        return {
            "is_ionic_aqueous": False,
            "complete_ionic_equation": msg,
            "net_ionic_equation": msg,
            "spectator_ions": [],
            "spectator_ions_summary_ar": "لا توجد أيونات متفرجة مائية.",
        }

    complete_eq = f"{_format_ion_side(lhs_tokens)} → {_format_ion_side(rhs_tokens)}"

    # Aggregate counts by key on LHS and RHS to cancel spectator ions
    lhs_counts: dict[str, int] = {}
    for t in lhs_tokens:
        lhs_counts[t["key"]] = lhs_counts.get(t["key"], 0) + int(t["count"])

    rhs_counts: dict[str, int] = {}
    for t in rhs_tokens:
        rhs_counts[t["key"]] = rhs_counts.get(t["key"], 0) + int(t["count"])

    spectator_ions: list[str] = []
    common_keys = [k for k in lhs_counts if k in rhs_counts]
    cancel_map: dict[str, int] = {}
    for k in common_keys:
        c_min = min(lhs_counts[k], rhs_counts[k])
        if c_min > 0:
            cancel_map[k] = c_min
            spectator_ions.append(k)

    net_lhs: list[dict[str, Any]] = []
    rem_cancel_l = dict(cancel_map)
    for t in lhs_tokens:
        k = t["key"]
        c_cancel = min(int(t["count"]), rem_cancel_l.get(k, 0))
        rem = int(t["count"]) - c_cancel
        rem_cancel_l[k] = rem_cancel_l.get(k, 0) - c_cancel
        if rem > 0:
            net_lhs.append({**t, "count": rem})

    net_rhs: list[dict[str, Any]] = []
    rem_cancel_r = dict(cancel_map)
    for t in rhs_tokens:
        k = t["key"]
        c_cancel = min(int(t["count"]), rem_cancel_r.get(k, 0))
        rem = int(t["count"]) - c_cancel
        rem_cancel_r[k] = rem_cancel_r.get(k, 0) - c_cancel
        if rem > 0:
            net_rhs.append({**t, "count": rem})

    if net_lhs and net_rhs:
        net_eq = f"{_format_ion_side(net_lhs)} → {_format_ion_side(net_rhs)}"
    else:
        net_eq = existing_net_ionic or complete_eq

    if spectator_ions:
        spec_summary = (
            "الأيونات المتفرجة في المحلول المائي (Spectator Ions): "
            + " ، ".join(spectator_ions)
            + " — بقيت ذائبة في الطور المائي قبل التفاعل وبعده دون تغير في حالتها الأيونية."
        )
    else:
        spec_summary = "جميع الأيونات المدخلة شاركت فعلياً في التفاعل الصافي (لا توجد أيونات متفرجة)."

    return {
        "is_ionic_aqueous": True,
        "complete_ionic_equation": complete_eq,
        "net_ionic_equation": net_eq,
        "spectator_ions": spectator_ions,
        "spectator_ions_summary_ar": spec_summary,
    }


def evaluate_reaction_equilibrium(
    reactants: list[dict[str, Any]],
    products: list[dict[str, Any]],
    category: str = "",
    equation_unicode: str = "",
    rxn_id: str = "",
    rule_id: str = "",
) -> dict[str, Any]:
    """
    Evaluates chemical reversibility and equilibrium behavior so that reversible reactions
    are never misrepresented as 100%-complete irreversible transformations.
    """
    r_formulas = {r.get("formula", "") for r in reactants}
    p_formulas = {p.get("formula", "") for p in products}

    is_reversible = False
    driving_force_ar = ""
    le_chatelier_ar = ""

    if (
        rxn_id in REVERSIBLE_REACTION_IDS
        or rule_id == "esterification"
        or ({("N2"), ("H2")} <= r_formulas and "NH3" in p_formulas)
        or ({"SO2", "O2"} <= r_formulas and "SO3" in p_formulas)
        or ({"CO2", "H2O"} == r_formulas and "H2CO3" in p_formulas)
        or ({"H2", "I2"} <= r_formulas and "HI" in p_formulas)
        or (category == "organic_reaction" and "H2O" in p_formulas and any("COO" in f for f in p_formulas))
    ):
        is_reversible = True
        if rule_id == "esterification" or any("COO" in f for f in p_formulas):
            driving_force_ar = (
                "تفاعل أسترة عكوس (Reversible Equilibrium, K_eq ≈ 4): يتوازن تكثيف الحمض الكربوكسيلي والكحول "
                "مع التحلل المائي العكسي للإستر."
            )
            le_chatelier_ar = (
                "وفق مبدأ لوشاتيليه، يزاح الاتزان نحو تكوين الإستر بإزالة الماء الناتج باستمرار أو استخدام فائض من الكحول."
            )
        elif "NH3" in p_formulas:
            driving_force_ar = (
                "تفاعل غازي عكوس طارد للحرارة (N₂ + 3H₂ ⇌ 2NH₃) محكوم باتزان ثرموديناميكي وحركي."
            )
            le_chatelier_ar = (
                "يزاح الاتزان نحو الأمونيا برفع الضغط (انخفاض مولات الغاز من 4 إلى 2) وخفض الحرارة المعتدل مع حفاز الحديد."
            )
        elif "SO3" in p_formulas:
            driving_force_ar = "أكسدة غازية عكوسة طاردة للحرارة (2SO₂ + O₂ ⇌ 2SO₃)."
            le_chatelier_ar = "يزاح الاتزان نحو SO₃ بزيادة الأكسجين والضغط وضبط الحرارة حول 450 °C بوجود حفاز V₂O₅."
        else:
            driving_force_ar = "تفاعل عكوس يصل إلى حالة اتزان كيميائي ديناميكي بين المتفاعلات والنواتج."
            le_chatelier_ar = "تتأثر نسبة التحويل الفعلية عند الاتزان بتركيز المتفاعلات ودرجة الحرارة والضغط."
    else:
        if category == "precipitation":
            driving_force_ar = "تفاعل تام عملياً (Irreversible Completion): انفصال الراسب شحيح الذوبان (Ksp صغير جداً) يسحب الأيونات من المحلول."
        elif category == "combustion":
            driving_force_ar = "تفاعل احتراق تام غير عكوس ثرموديناميكياً (ΔG° سالبة وكبيرة جداً)."
        elif any(f in {"CO2", "H2", "O2"} for f in p_formulas) and category in {"neutralization", "single_displacement", "decomposition"}:
            driving_force_ar = "تفاعل يتجه نحو الاكتمال بسبب تصاعد الناتج الغازي وخروجه من حيز التفاعل."
        elif category == "neutralization":
            driving_force_ar = "تفاعل تعادل يتجه نحو الاكتمال لتكوّن جزيئات الماء ضعيفة التأين (K ≈ 10¹⁴ للأحماض والقواعد القوية)."
        else:
            driving_force_ar = "تفاعل يتجه في الاتجاه الأمامي تحت الظروف المرجعية المحددة."
        le_chatelier_ar = "التفاعل محكوم بالمتفاعل المحدد ويتجه عملياً نحو الاكتمال في الاتجاه الأمامي."

    arrow = "⇌" if is_reversible else "→"
    eq_equil = equation_unicode.replace("→", arrow).replace("->", arrow) if equation_unicode else ""

    return {
        "is_reversible": is_reversible,
        "reaction_arrow": arrow,
        "equilibrium_equation_unicode": eq_equil,
        "equilibrium_type_ar": (
            "تفاعل عكوس يصل إلى اتزان كيميائي ديناميكي (Reversible Chemical Equilibrium ⇌)"
            if is_reversible
            else "تفاعل يتجه نحو الاكتمال في الاتجاه الأمامي (Completion-Driven Reaction →)"
        ),
        "driving_force_ar": driving_force_ar,
        "le_chatelier_guidance_ar": le_chatelier_ar,
        "yield_equilibrium_note_ar": (
            "تنبيه المردود والاتزان: يمثل المردود النظري المحسوب (Theoretical Yield) الحد الأقصى عند تحويل 100% من المتفاعل المحدد؛ "
            + (
                "ولأن هذا التفاعل عكوس (⇌)، فإن التحويل الفعلي يتوقف عند نقطة الاتزان الكيميائي (< 100%) ما لم يُسحب أحد النواتج، "
                "ولا يُعرض أي مردود فعلي (Actual Yield) بدون قياس تجريبي مخبري."
                if is_reversible
                else "ولا يُعرض أي مردود فعلي (Actual Yield) أو نسبة مئوية للمردود بدون بيانات قياس تجريبية من المختبر."
            )
        ),
    }


def detect_competing_pathways(
    species: list[dict[str, Any]],
    reactants: list[dict[str, Any]],
    products: list[dict[str, Any]],
    category: str = "",
    alternatives: list[dict[str, Any]] | None = None,
) -> dict[str, Any]:
    """
    Detects and documents competing or parallel chemical pathways depending on molar ratios
    and reaction conditions (e.g., complete vs incomplete combustion, stepwise polyprotic
    neutralization, dilute vs oxidizing acid with metals, and multi-reactant mixtures).
    """
    pathways: list[dict[str, str]] = []
    r_formulas = {r.get("formula", "") for r in reactants}

    # 1. Registry alternatives for the same reactants
    for alt in alternatives or []:
        pathways.append({
            "pathway_name_ar": alt.get("reaction_name_ar", "مسار تفاعلي بديل موثق"),
            "condition_trigger_ar": f"تغيير مجال/نوع التفاعل أو الظروف ({alt.get('category', '')})",
            "products_summary_ar": f"راجع التفاعل الموثق البديل رقم ({alt.get('reaction_id', '')})",
            "source_basis_ar": "سجل التفاعلات الموثقة (VERIFIED_REACTIONS_REGISTRY)",
        })

    # 2. Complete vs Incomplete Combustion of hydrocarbons/alcohols
    if category == "combustion" or ("O2" in r_formulas and any("C" in f for f in r_formulas if f != "O2")):
        pathways.append({
            "pathway_name_ar": "الاحتراق غير الكامل عند نقص الأكسجين (Incomplete Combustion)",
            "condition_trigger_ar": "نقص نسبة الأكسجين (O₂ متفاعل محدد بنسبة أقل من النسبة الستوكيومترية الكاملة)",
            "products_summary_ar": "أول أكسيد الكربون السام (CO) و/أو السخام الكربوني (C) مع الماء (H₂O) بدلاً من CO₂ الكامل",
            "source_basis_ar": "CRC Handbook / Brown LeMay Ch. 3 & Ch. 5 — Hydrocarbon Combustion Stoichiometry",
        })

    # 3. Polyprotic Acid Stepwise Partial vs Complete Neutralization
    polyprotic_acids = {"H2SO4": "NaHSO₄ / KHSO₄", "H3PO4": "H₂PO₄⁻ / HPO₄²⁻", "H2CO3": "HCO₃⁻", "H2C2O4": "HC₂O₄⁻"}
    for pa, inter_salt in polyprotic_acids.items():
        if pa in r_formulas and category == "neutralization":
            pathways.append({
                "pathway_name_ar": f"التعادل المرحلي الجزئي للحمض متعدد البروتونات ({pa})",
                "condition_trigger_ar": "إضافة القاعدة بنسبة مولية (1:1) غير كافية لمعادلة جميع بروتونات الحمض",
                "products_summary_ar": f"تكوّن الملح الحمضي الوسطي ({inter_salt}) بدلاً من الملح المتعادل الكامل",
                "source_basis_ar": "Harris Quantitative Chemical Analysis — Polyprotic Acid-Base Equilibria",
            })

    # 4. Metal + Acid: Dilute Non-Oxidizing vs Concentrated Oxidizing Acid
    if category == "single_displacement" and any(f in {"HCl", "H2SO4"} for f in r_formulas):
        pathways.append({
            "pathway_name_ar": "مسار الأكسدة بالأحماض المؤكسدة أو المركزة الساخنة",
            "condition_trigger_ar": "استخدام حمض الكبريتيك المركز الساخن أو حمض النيتريك (HNO₃) بدلاً من الحمض المخفف",
            "products_summary_ar": "تصاعد غازات SO₂ أو NO₂ / NO مع الماء بدلاً من تصاعد غاز الهيدروجين (H₂)",
            "source_basis_ar": "Greenwood & Earnshaw, Chemistry of the Elements — Metal-Acid Redox Pathways",
        })

    # 5. Multi-reactant mixture competing pairs
    if len(species) >= 3:
        sp_names = " + ".join(s.get("formula", "") for s in species)
        pathways.append({
            "pathway_name_ar": f"التنافس في مخلوط متعدد المواد ({sp_names})",
            "condition_trigger_ar": "وجود أكثر من مادتين في نفس وعاء التفاعل قد يتيح مسارات تفاعل متوازية بحسب سرعة التفاعل وثوابت الاتزان",
            "products_summary_ar": "الأولوية الثرموديناميكية لتفاعلات التعادل السريعة والترسيب ذي Ksp الأدنى، بينما تعزل المواد الفائضة أو الخاملة كمواد متفرجة",
            "source_basis_ar": "Vogel's Qualitative Inorganic Analysis — Selective Precipitation & Competing Equilibria",
        })

    return {
        "has_competing_pathways": len(pathways) > 0,
        "pathways": pathways,
        "summary_ar": (
            f"تم رصد ({len(pathways)}) مسار/مسارات تفاعلية متنافسة أو بديلة تعتمد على نسبة المولات وظروف الوسط."
            if pathways
            else "مسار تفاعلي وحيد مهيمن تحت الظروف القياسية المدخلة."
        ),
    }


def evaluate_reaction_safety_and_hazards(
    reactants: list[dict[str, Any]],
    products: list[dict[str, Any]],
    category: str = "",
    conditions: dict[str, Any] | None = None,
) -> dict[str, Any]:
    """
    Generates laboratory hazard warnings for dangerous reactants/products/conditions,
    enforces the virtual-lab simulation disclaimer, and withholds operational recipes
    for hazardous reactions.
    """
    all_species = list(reactants) + list(products)
    all_formulas = {s.get("formula", "") for s in all_species}
    warnings: list[str] = []
    hazard_level = "low"

    # 1. Alkali metal + water violent reaction
    if "H2O" in all_formulas and (all_formulas & REACTIVE_WATER_METALS):
        metals_found = "، ".join(sorted(all_formulas & REACTIVE_WATER_METALS))
        warnings.append(
            f"تحذير خطورة عالية (Pyrophoric Alkali Metal + Water): تفاعل الفلز القلوي ({metals_found}) مع الماء "
            "عنيف وشديد الطرد للحرارة وقد يؤدي إلى اشتعال غاز الهيدروجين المتصاعد وانفجاره."
        )
        hazard_level = "high"

    # 2. Toxic / flammable gases
    for f, msg in TOXIC_OR_HAZARDOUS_GASES.items():
        if f in all_formulas:
            warnings.append(f"تحذير غاز خطر ({format_chemical_formula_display(f)}): {msg}")
            if f in {"Cl2", "Br2", "CO", "HCN", "H2S", "NO2", "SO3", "H2"}:
                hazard_level = "high"
            elif hazard_level == "low":
                hazard_level = "moderate"

    # 3. Corrosive strong acids / bases
    for f, msg in CORROSIVE_SPECIES.items():
        if f in all_formulas:
            warnings.append(f"تحذير مادة كاوية/أكّالة ({format_chemical_formula_display(f)}): {msg}")
            if f in {"HF", "HClO4", "H2SO4", "HNO3"}:
                hazard_level = "high"
            elif hazard_level == "low":
                hazard_level = "moderate"

    # 4. Heavy metals & oxidizing salts
    for f in all_formulas:
        if any(hm in f for hm in HEAVY_METAL_ELEMENTS):
            warnings.append(
                f"تحذير سمية بيئية وصحية ({format_chemical_formula_display(f)}): يحتوي على معادن ثقيلة أو مؤكسدات سامة تتطلب بروتوكول تخلص خاص من النفايات الكيميائية."
            )
            if hazard_level == "low":
                hazard_level = "moderate"
            break

    # 5. Combustion / high temperature & pressure
    if category == "combustion":
        warnings.append("تحذير اشتعال وحرارة (Combustion Hazard): تفاعل احتراق طارد للحرارة بشدة يتطلب عزلاً عن الأبخرة القابلة للاشتعال.")
        if hazard_level == "low":
            hazard_level = "moderate"

    if not warnings:
        warnings.append("تنبيه سلامة مخبرية عام: يلزم ارتداء معدات الوقاية الشخصية القياسية (المعطف، النظارات الواقية، والقفازات) عند التعامل مع أي مواد كيميائية.")

    return {
        "has_hazards": hazard_level in {"high", "moderate"},
        "hazard_level": hazard_level,
        "hazard_level_ar": {
            "high": "خطورة مخبرية عالية (High Chemical Hazard)",
            "moderate": "خطورة مخبرية متوسطة (Moderate Chemical Hazard)",
            "low": "خطورة مخبرية منخفضة / قياسية (Standard Lab Precaution)",
        }[hazard_level],
        "warnings_ar": warnings,
        "operational_instructions_withheld": True,
        "operational_safety_policy_ar": (
            "سياسة الأمان الأكاديمي: تكتفي المنصة بالتحليل العلمي والمعادلات الستوكيومترية النظرية، "
            "وتمتنع صراحةً عن تقديم وصفات تشغيلية خطوة بخطوة لتجارب أو مخاليط خطرة قد تسبب إصابات أو أضراراً."
        ),
        "virtual_lab_disclaimer_ar": (
            "إخلاء مسؤولية علمي: جميع نتائج هذا المختبر محاكاة افتراضية تعليمية (In-Silico Virtual Simulation) "
            "وليست بديلاً عن التجارب الواقعية تحت إشراف مخبري معتمد أو أوراق بيانات السلامة الرسمية (SDS)."
        ),
    }


def get_reaction_lab_coverage_boundaries() -> dict[str, Any]:
    """Explicitly documents what the Reaction Lab engine covers and what it cannot analyze."""
    return {
        "documented_reactions_count": 152,
        "experimentally_documented_count": 102,
        "literature_referenced_count": 41,
        "theoretical_model_count": 9,
        "manually_verified_primary_sources_count": 36,
        "general_inference_rules_count": 13,
        "explicit_no_reaction_rules_count": 6,
        "local_catalog_substances_count": 241,
        "max_simultaneous_reactants": 4,
        "covered_capabilities_ar": [
            "152 تفاعلاً موثقاً بمرجع فردي ومصنفاً حسب قوة الدليل (102 تجريبي، 41 مرجعي في الأدبيات، 9 نماذج نظرية).",
            "13 قاعدة كيميائية عامة معلنة المصدر والحدود (التعادل، الكربونات، الترسيب، الإحلال الفلزي والهالوجيني، الأكاسيد، الاحتراق، والأسترة).",
            "6 قواعد صريحة لكشف عدم حدوث التفاعل (الغازات النبيلة، الفلزات أسفل الهيدروجين مع الأحماض غير المؤكسدة، الإحلال العكسي، والأيونات المتفرجة الذائبة).",
            "دعم تفاعلات مادتين أو 3 أو 4 مواد مع عزل المواد المتفرجة وحساب المعادلات الأيونية الكاملة والصافية والمتفاعل المحدد والمردود النظري.",
        ],
        "unsupported_boundaries_ar": [
            "التخليق العضوي متعدد الخطوات والتفاعلات الفراغية المعقدة غير المدرجة في السجل.",
            "حساب المردود الفعلي (Actual Yield) أو ثوابت سرعة التفاعل الحركية (Kinetics) دون معطيات تجريبية.",
            "معقدات العناصر الانتقالية التساندية غير القياسية وتفاعلات البلمرة الجذرية.",
        ],
    }

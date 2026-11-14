"""
Unit tests for What-If scenario analysis engine in Chemical Playground.
Verifies baseline vs modified scenario evaluations, limiting reagent shifts,
yield variations, and condition impact analysis.
"""

import pytest
from reaction_engine.what_if import WhatIfAnalyzer


def test_what_if_stoichiometry_and_yield_increase():
    # Baseline: 1.0 mol CH3COOH + 1.0 mol NaOH
    base_inputs = [
        {"formula": "CH3COOH", "amount": 1.0, "unit": "mol", "substance": {"name_ar": "حمض الأسيتيك"}},
        {"formula": "NaOH", "amount": 1.0, "unit": "mol", "substance": {"name_ar": "هيدروكسيد الصوديوم"}},
    ]
    # Modified: 2.0 mol CH3COOH + 2.0 mol NaOH (doubled yield)
    mod_inputs = [
        {"formula": "CH3COOH", "amount": 2.0, "unit": "mol", "substance": {"name_ar": "حمض الأسيتيك"}},
        {"formula": "NaOH", "amount": 2.0, "unit": "mol", "substance": {"name_ar": "هيدروكسيد الصوديوم"}},
    ]

    result = WhatIfAnalyzer.evaluate_what_if(base_inputs, mod_inputs)
    assert result["success"] is True
    assert "yield_comparisons" in result
    assert len(result["yield_comparisons"]) > 0

    # Yield of products should have doubled (+100%)
    for yc in result["yield_comparisons"]:
        assert yc["pct_change"] == 100.0
        assert yc["modified_moles"] == pytest.approx(yc["baseline_moles"] * 2.0, 0.001)

    assert "table_rows" in result
    assert len(result["table_rows"]) >= 3
    assert result["insights_ar"] != ""


def test_what_if_limiting_reagent_shift():
    # Baseline: 1.0 mol H2SO4 + 1.0 mol NaOH -> NaOH is limiting (ratio requires 2 NaOH per H2SO4)
    base_inputs = [
        {"formula": "H2SO4", "amount": 1.0, "unit": "mol", "substance": {"name_ar": "حمض الكبريتيك"}},
        {"formula": "NaOH", "amount": 1.0, "unit": "mol", "substance": {"name_ar": "هيدروكسيد الصوديوم"}},
    ]
    # Modified: 1.0 mol H2SO4 + 3.0 mol NaOH -> H2SO4 is now limiting
    mod_inputs = [
        {"formula": "H2SO4", "amount": 1.0, "unit": "mol", "substance": {"name_ar": "حمض الكبريتيك"}},
        {"formula": "NaOH", "amount": 3.0, "unit": "mol", "substance": {"name_ar": "هيدروكسيد الصوديوم"}},
    ]

    result = WhatIfAnalyzer.evaluate_what_if(base_inputs, mod_inputs)
    assert result["success"] is True
    assert result["limiting_reagent_shifted"] is True
    assert "تحول" in result["insights_ar"] or "المتفاعل المحدد" in result["insights_ar"]


def test_what_if_temperature_and_pressure_impact():
    base_inputs = [
        {"formula": "HCl", "amount": 1.0, "unit": "mol", "substance": {"name_ar": "حمض الهيدروكلوريك"}},
        {"formula": "NaOH", "amount": 1.0, "unit": "mol", "substance": {"name_ar": "هيدروكسيد الصوديوم"}},
    ]
    mod_inputs = [
        {"formula": "HCl", "amount": 1.0, "unit": "mol", "substance": {"name_ar": "حمض الهيدروكلوريك"}},
        {"formula": "NaOH", "amount": 1.0, "unit": "mol", "substance": {"name_ar": "هيدروكسيد الصوديوم"}},
    ]

    base_cond = {"temperature_c": 25.0, "pressure_atm": 1.0}
    mod_cond = {"temperature_c": 80.0, "pressure_atm": 5.0}

    result = WhatIfAnalyzer.evaluate_what_if(
        base_inputs,
        mod_inputs,
        baseline_conditions=base_cond,
        modified_conditions=mod_cond,
    )
    assert result["success"] is True
    assert "درجة الحرارة" in result["insights_ar"]
    assert "الضغط" in result["insights_ar"]


def test_what_if_unreactive_baseline_handled_safely():
    # Baseline with inert / unreactive mix
    base_inputs = [
        {"formula": "Fe", "amount": 1.0, "unit": "mol", "substance": {"name_ar": "الحديد"}},
        {"formula": "He", "amount": 1.0, "unit": "mol", "substance": {"name_ar": "الهيليوم"}},
    ]
    mod_inputs = [
        {"formula": "Fe", "amount": 2.0, "unit": "mol", "substance": {"name_ar": "الحديد"}},
        {"formula": "He", "amount": 2.0, "unit": "mol", "substance": {"name_ar": "الهيليوم"}},
    ]

    result = WhatIfAnalyzer.evaluate_what_if(base_inputs, mod_inputs)
    assert result["success"] is False
    assert "error_ar" in result

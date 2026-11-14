"""
What-If Scenario Analysis engine for MolecuLens Chemical Playground.
Compares baseline vs modified reaction conditions: reactant amounts, temperatures,
pressures, limiting reagent shifts, and theoretical yield variations with Le Chatelier insights.
"""

from __future__ import annotations

from typing import Any
import html

from reaction_engine.prediction import (
    predict_reaction,
    compute_reaction_stoichiometry,
    compute_theoretical_yields,
)


class WhatIfAnalyzer:
    """Analyzes comparative 'what-if' scenarios in chemical reactions."""

    @classmethod
    def evaluate_what_if(
        cls,
        baseline_inputs: list[dict[str, Any]],
        modified_inputs: list[dict[str, Any]],
        baseline_conditions: dict[str, Any] | None = None,
        modified_conditions: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        """
        Executes both baseline and modified scenarios and generates a comprehensive
        side-by-side scientific comparison and thermodynamic analysis.
        """
        base_pred = predict_reaction(baseline_inputs, custom_conditions=baseline_conditions)
        mod_pred = predict_reaction(modified_inputs, custom_conditions=modified_conditions)

        if base_pred.get("mode") not in ["documented", "rule_inferred"]:
            return {
                "success": False,
                "error_ar": "التجربة المرجعية الأصلية لا تؤدي إلى تفاعل كيميائي، لذا يتعذر إجراء مقارنة السيناريوهات البديلة.",
            }

        # Compute stoichiometry and yields for baseline
        base_stoich = compute_reaction_stoichiometry(base_pred["reactants"], baseline_inputs)
        base_yields = compute_theoretical_yields(base_pred["products"], base_stoich["reaction_extent_mol"])

        # Compute stoichiometry and yields for modified
        if mod_pred.get("mode") in ["documented", "rule_inferred"]:
            mod_stoich = compute_reaction_stoichiometry(mod_pred["reactants"], modified_inputs)
            mod_yields = compute_theoretical_yields(mod_pred["products"], mod_stoich["reaction_extent_mol"])
        else:
            mod_stoich = {}
            mod_yields = []

        base_lr = base_stoich.get("limiting_reagent") or {}
        mod_lr = mod_stoich.get("limiting_reagent") or {}

        base_lr_formula = base_lr.get("formula", "متكافئ تماماً")
        mod_lr_formula = mod_lr.get("formula", "متكافئ تماماً")

        lr_shifted = (base_lr_formula != mod_lr_formula)

        yield_comparisons = []
        for b_y in base_yields:
            p_form = b_y.get("formula")
            m_y = next((y for y in mod_yields if y.get("formula") == p_form), None)

            b_mol = float(b_y.get("theoretical_yield_mol") or 0.0)
            b_mass = float(b_y.get("theoretical_yield_g") or 0.0)
            m_mol = float(m_y.get("theoretical_yield_mol") or 0.0) if m_y else 0.0
            m_mass = float(m_y.get("theoretical_yield_g") or 0.0) if m_y else 0.0

            delta_mass = round(m_mass - b_mass, 3)
            delta_mol = round(m_mol - b_mol, 4)
            pct_change = round(((m_mass - b_mass) / b_mass * 100.0), 1) if b_mass > 0 else 0.0

            yield_comparisons.append({
                "product_formula": p_form,
                "product_name_ar": b_y.get("name_ar", p_form),
                "baseline_moles": b_mol,
                "baseline_mass_g": b_mass,
                "modified_moles": m_mol,
                "modified_mass_g": m_mass,
                "delta_mass_g": delta_mass,
                "delta_moles": delta_mol,
                "pct_change": pct_change,
            })

        # Conditions impact analysis (Temperature and Pressure)
        b_cond = baseline_conditions or {}
        m_cond = modified_conditions or {}

        b_temp = float(b_cond.get("temperature_c", 25.0))
        m_temp = float(m_cond.get("temperature_c", 25.0))
        delta_temp = m_temp - b_temp

        b_press = float(b_cond.get("pressure_atm", 1.0))
        m_press = float(m_cond.get("pressure_atm", 1.0))
        delta_press = m_press - b_press

        temp_impact_ar = []
        if abs(delta_temp) >= 5.0:
            if delta_temp > 0:
                temp_impact_ar.append(
                    f"رفع درجة الحرارة بمقدار (+{delta_temp:.1f} °C) يرفع الطاقة الحركية وسرعة تصادم الجزيئات لتجاوز حاجز طاقة التنشيط."
                )
            else:
                temp_impact_ar.append(
                    f"خفض درجة الحرارة بمقدار ({delta_temp:.1f} °C) يقلل من معدل التصادمات والسرعة الحركية."
                )

        press_impact_ar = []
        if abs(delta_press) >= 0.5:
            if delta_press > 0:
                press_impact_ar.append(
                    f"زيادة الضغط بمقدار (+{delta_press:.1f} atm) توجه الاتزان نحو الحجم والطرف الأقل مولات غازية (مبدأ لوشاتيليه)."
                )
            else:
                press_impact_ar.append(
                    f"خفض الضغط بمقدار ({delta_press:.1f} atm) يوجه التفاعل لصالح الطرف ذي الحجم الغازي الأكبر."
                )

        # Synthesize analytical explanation
        insights = []
        if lr_shifted:
            insights.append(
                f"حدث تحوّل في هوية المتفاعل المحدد من ({base_lr.get('name_ar', base_lr_formula)}) "
                f"إلى ({mod_lr.get('name_ar', mod_lr_formula)})؛ "
                "بسبب تعديل النسب المولية بحيث أصبح المتفاعل الآخر هو الذي ينفد أولاً ويحدد سقف الناتج."
            )
        elif yield_comparisons:
            p1 = yield_comparisons[0]
            if p1["pct_change"] > 0:
                insights.append(
                    f"ارتفعت كمية الناتج النظري ({p1['product_name_ar']}) بنسبة (+{p1['pct_change']}%) "
                    f"نتيجة زيادة كمية المتفاعل المحدد المغذي للتفاعل."
                )
            elif p1["pct_change"] < 0:
                insights.append(
                    f"انخفضت كمية الناتج النظري ({p1['product_name_ar']}) بنسبة ({p1['pct_change']}%) "
                    "بسبب خفض كمية المتفاعل المحدد المستهلك أولاً."
                )
            else:
                insights.append("لم تتغير كمية الناتج النظري لأن التعديل طال مادة فائضة دون المساس بكمية المتفاعل المحدد.")

        if temp_impact_ar:
            insights.extend(temp_impact_ar)
        if press_impact_ar:
            insights.extend(press_impact_ar)

        # Build comparison summary rows for table display
        table_rows = []
        table_rows.append({
            "المتغير (Parameter)": "المتفاعل المحدد (Limiting Reagent)",
            "التجربة الأصلية (Baseline)": f"{base_lr.get('name_ar', '')} ({base_lr_formula})" if base_lr else "نسبة متكافئة تماماً",
            "التجربة المعدلة (Modified Scenario)": f"{mod_lr.get('name_ar', '')} ({mod_lr_formula})" if mod_lr else "نسبة متكافئة تماماً",
            "التغير والملاحظة": "تحول المتفاعل المحدد" if lr_shifted else "لا تغيير في الهوية",
        })

        for yc in yield_comparisons:
            table_rows.append({
                "المتغير (Parameter)": f"الناتج النظري ({yc['product_formula']})",
                "التجربة الأصلية (Baseline)": f"{yc['baseline_mass_g']} g ({yc['baseline_moles']} mol)",
                "التجربة المعدلة (Modified Scenario)": f"{yc['modified_mass_g']} g ({yc['modified_moles']} mol)",
                "التغير والملاحظة": f"{'+' if yc['pct_change'] > 0 else ''}{yc['pct_change']}% كتلياً",
            })

        table_rows.append({
            "المتغير (Parameter)": "درجة الحرارة (Temperature)",
            "التجربة الأصلية (Baseline)": f"{b_temp} °C",
            "التجربة المعدلة (Modified Scenario)": f"{m_temp} °C",
            "التغير والملاحظة": f"{'+' if delta_temp > 0 else ''}{delta_temp:.1f} °C",
        })

        table_rows.append({
            "المتغير (Parameter)": "الضغط (Pressure)",
            "التجربة الأصلية (Baseline)": f"{b_press} atm",
            "التجربة المعدلة (Modified Scenario)": f"{m_press} atm",
            "التغير والملاحظة": f"{'+' if delta_press > 0 else ''}{delta_press:.1f} atm",
        })

        return {
            "success": True,
            "baseline_prediction": base_pred,
            "modified_prediction": mod_pred,
            "baseline_stoichiometry": base_stoich,
            "modified_stoichiometry": mod_stoich,
            "limiting_reagent_shifted": lr_shifted,
            "yield_comparisons": yield_comparisons,
            "table_rows": table_rows,
            "insights_ar": " ".join(insights),
        }

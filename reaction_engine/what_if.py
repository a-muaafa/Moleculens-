"""
What-If Scenario Analysis engine for MolecuLens Chemical Playground.
Compares baseline vs modified reaction conditions: reactant amounts, temperatures,
pressures, limiting reagent shifts, and theoretical yield variations strictly grounded
in actual gas stoichiometry, enthalpy direction, and documented reference condition limits.
"""

from __future__ import annotations

import re
from typing import Any

from reaction_engine.prediction import (
    compute_reaction_stoichiometry,
    compute_theoretical_yields,
    predict_reaction,
)


class WhatIfAnalyzer:
    """Analyzes comparative 'what-if' scenarios in chemical reactions."""

    _VALID_STATES = {"g", "l", "s", "aq"}

    @classmethod
    def _analyze_gas_moles(
        cls,
        reactants: list[dict[str, Any]],
        products: list[dict[str, Any]],
    ) -> dict[str, Any]:
        """Counts gaseous stoichiometric moles on each side when physical states are known."""
        if not reactants or not products:
            return {"known_states": False, "r_gas": None, "p_gas": None, "delta_n_gas": None}

        r_gas = 0
        p_gas = 0
        for item in reactants:
            st_code = str(item.get("state") or "").strip("() ").lower()
            if st_code not in cls._VALID_STATES:
                return {"known_states": False, "r_gas": None, "p_gas": None, "delta_n_gas": None}
            if st_code == "g":
                r_gas += int(item.get("coefficient", 1) or 1)

        for item in products:
            st_code = str(item.get("state") or "").strip("() ").lower()
            if st_code not in cls._VALID_STATES:
                return {"known_states": False, "r_gas": None, "p_gas": None, "delta_n_gas": None}
            if st_code == "g":
                p_gas += int(item.get("coefficient", 1) or 1)

        return {
            "known_states": True,
            "r_gas": r_gas,
            "p_gas": p_gas,
            "delta_n_gas": p_gas - r_gas,
        }

    @classmethod
    def _extract_thermal_nature(cls, prediction: dict[str, Any]) -> tuple[str, str]:
        """
        Determines whether the reaction is documented/inferred as exothermic, endothermic,
        or thermodynamically unknown based solely on actual reaction data.
        """
        conds = prediction.get("conditions") or {}
        thermo_text = str(conds.get("thermodynamics_ar", "")) if isinstance(conds, dict) else ""
        obs_texts = [
            str(o.get("text_ar", ""))
            for o in (prediction.get("observations") or [])
            if isinstance(o, dict)
        ]
        temp_text = str(conds.get("temperature", "")) if isinstance(conds, dict) else ""
        combined = " ".join([thermo_text, temp_text] + obs_texts)

        m = re.search(r"ΔH°?[a-zA-Z]*\s*[=≈]\s*([+\-−]?\d+(?:\.\d+)?)\s*kJ", combined)
        if m:
            val = float(m.group(1).replace("−", "-"))
            if val < 0:
                return "exothermic", "موثق كيميائياً كتفاعل طارد للحرارة"
            if val > 0:
                return "endothermic", "موثق كيميائياً كتفاعل ماص للحرارة"

        if "طارد للحرارة" in combined or prediction.get("category") == "combustion":
            return "exothermic", "موثق كيميائياً كتفاعل طارد للحرارة"
        if "ماص للحرارة" in combined:
            return "endothermic", "موثق كيميائياً كتفاعل ماص للحرارة"

        return "unknown", "غير محدد ثرموديناميكياً في بيانات التفاعل"

    @classmethod
    def _analyze_temperature_change(
        cls,
        delta_temp: float,
        m_temp: float,
        thermal_nature: str,
        ref_spec: dict[str, Any],
    ) -> tuple[list[str], bool]:
        """Generates evidence-grounded temperature analysis without unsupported claims."""
        notes: list[str] = []
        supported = False
        if abs(delta_temp) < 5.0:
            return notes, True

        if thermal_nature == "exothermic":
            supported = True
            if delta_temp > 0:
                notes.append(
                    f"بما أن التفاعل طارد للحرارة وفق بياناته، فإن رفع درجة الحرارة بمقدار (+{delta_temp:.1f} °C) "
                    "يزيد سرعة الوصول إلى الاتزان حركياً، لكنه يزيح موضع الاتزان الثرموديناميكي عكسياً نحو المتفاعلات وفق مبدأ لوشاتيليه."
                )
            else:
                notes.append(
                    f"بما أن التفاعل طارد للحرارة وفق بياناته، فإن خفض درجة الحرارة بمقدار ({delta_temp:.1f} °C) "
                    "يدعم الاتجاه الأمامي (نحو النواتج) ثرموديناميكياً، مع تباطؤ السرعة الحركية للتصادمات."
                )
        elif thermal_nature == "endothermic":
            supported = True
            if delta_temp > 0:
                notes.append(
                    f"بما أن التفاعل ماص للحرارة وفق بياناته، فإن رفع درجة الحرارة بمقدار (+{delta_temp:.1f} °C) "
                    "يزيح موضع الاتزان أمامياً نحو النواتج ويزيد من سرعة التفاعل."
                )
            else:
                notes.append(
                    f"بما أن التفاعل ماص للحرارة وفق بياناته، فإن خفض درجة الحرارة بمقدار ({delta_temp:.1f} °C) "
                    "يزيح موضع الاتزان عكسياً نحو المتفاعلات ويقلل من معدل التحول."
                )
        else:
            supported = False
            notes.append(
                f"تغيرت درجة الحرارة بمقدار ({'+' if delta_temp > 0 else ''}{delta_temp:.1f} °C)، "
                "ولكن لا تتوفر بيانات إنثالبي (ΔH) موثقة تحدد ما إذا كان التفاعل طارداً أم ماصاً للحرارة؛ "
                "لذا لا يمكن الجزم باتجاه انزياح الاتزان الحراري دون بيانات ثرموديناميكية."
            )

        t_min = ref_spec.get("temp_c_min")
        t_max = ref_spec.get("temp_c_max")
        if t_min is not None and t_max is not None and (m_temp < float(t_min) or m_temp > float(t_max)):
            notes.append(
                f"تنبيه مرجعي: درجة الحرارة المعدلة ({m_temp:.1f} °C) خارج النطاق الموثق للتفاعل "
                f"[{float(t_min):g} إلى {float(t_max):g} °C]."
            )

        return notes, supported

    @classmethod
    def _analyze_pressure_change(
        cls,
        delta_press: float,
        m_press: float,
        gas_info: dict[str, Any],
        ref_spec: dict[str, Any],
    ) -> tuple[list[str], bool]:
        """Generates evidence-grounded pressure analysis based strictly on Δn_gas."""
        notes: list[str] = []
        supported = False
        if abs(delta_press) < 0.5:
            return notes, True

        if not gas_info["known_states"]:
            supported = False
            notes.append(
                f"تغير الضغط بمقدار ({'+' if delta_press > 0 else ''}{delta_press:.1f} atm)، "
                "ولكن لا يمكن استنتاج أثر الضغط على الاتزان لعدم تحديد الحالات الفيزيائية (الأطوار الغازية) لجميع مواد التفاعل."
            )
        else:
            supported = True
            r_gas = int(gas_info["r_gas"])
            p_gas = int(gas_info["p_gas"])
            delta_n_gas = int(gas_info["delta_n_gas"])

            if r_gas == 0 and p_gas == 0:
                notes.append(
                    f"تغير الضغط بمقدار ({'+' if delta_press > 0 else ''}{delta_press:.1f} atm): "
                    "بما أن التفاعل يجري في الطور المكثف/المائي ولا يتضمن أي مواد غازية (مولات الغاز = 0 في الطرفين)، "
                    "فإن تغير الضغط لا يؤثر عملياً على موضع الاتزان."
                )
            elif delta_n_gas == 0:
                notes.append(
                    f"تغير الضغط بمقدار ({'+' if delta_press > 0 else ''}{delta_press:.1f} atm): "
                    f"يتساوى عدد مولات الغاز في طرفي المعادلة ({r_gas} مول غاز في المتفاعلات مقابل {p_gas} مول غاز في النواتج، Δn_gas = 0)؛ "
                    "لذا لا يزيح الضغط موضع الاتزان."
                )
            elif delta_n_gas < 0:
                if delta_press > 0:
                    notes.append(
                        f"زيادة الضغط بمقدار (+{delta_press:.1f} atm) تزيح الاتزان أمامياً نحو النواتج، "
                        f"لأن عدد مولات الغاز في النواتج ({p_gas}) أقل منه في المتفاعلات ({r_gas}؛ Δn_gas = {delta_n_gas})."
                    )
                else:
                    notes.append(
                        f"خفض الضغط بمقدار ({delta_press:.1f} atm) يزيح الاتزان عكسياً نحو المتفاعلات "
                        f"ذات العدد الأكبر من مولات الغاز ({r_gas} مقابل {p_gas}؛ Δn_gas = {delta_n_gas})."
                    )
            else:  # delta_n_gas > 0
                if delta_press > 0:
                    notes.append(
                        f"زيادة الضغط بمقدار (+{delta_press:.1f} atm) تزيح الاتزان عكسياً نحو المتفاعلات، "
                        f"لأن عدد مولات الغاز في النواتج ({p_gas}) أكبر منه في المتفاعلات ({r_gas}؛ Δn_gas = +{delta_n_gas})."
                    )
                else:
                    notes.append(
                        f"خفض الضغط بمقدار ({delta_press:.1f} atm) يزيح الاتزان أمامياً نحو النواتج "
                        f"ذات العدد الأكبر من مولات الغاز ({p_gas} مقابل {r_gas}؛ Δn_gas = +{delta_n_gas})."
                    )

        p_min = ref_spec.get("pressure_atm_min")
        p_max = ref_spec.get("pressure_atm_max")
        if p_min is not None and p_max is not None and (m_press < float(p_min) or m_press > float(p_max)):
            notes.append(
                f"تنبيه مرجعي: الضغط المعدل ({m_press:.1f} atm) خارج النطاق الموثق للتفاعل "
                f"[{float(p_min):g} إلى {float(p_max):g} atm]."
            )

        return notes, supported

    @classmethod
    def evaluate_what_if(
        cls,
        baseline_inputs: list[dict[str, Any]],
        modified_inputs: list[dict[str, Any]],
        baseline_conditions: dict[str, Any] | None = None,
        modified_conditions: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        """
        Executes both baseline and modified scenarios and generates a side-by-side
        stoichiometric comparison and reaction-grounded thermodynamic analysis.
        """
        base_pred = predict_reaction(baseline_inputs, custom_conditions=baseline_conditions)
        mod_pred = predict_reaction(modified_inputs, custom_conditions=modified_conditions)

        if base_pred.get("mode") not in ["documented", "rule_inferred"]:
            return {
                "success": False,
                "error_ar": "التجربة المرجعية الأصلية لا تؤدي إلى تفاعل كيميائي، لذا يتعذر إجراء مقارنة السيناريوهات البديلة.",
            }

        base_stoich = compute_reaction_stoichiometry(base_pred["reactants"], baseline_inputs)
        base_yields = compute_theoretical_yields(base_pred["products"], base_stoich["reaction_extent_mol"])

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
        lr_shifted = base_lr_formula != mod_lr_formula

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

        b_cond = baseline_conditions or {}
        m_cond = modified_conditions or {}

        b_temp = float(b_cond.get("temperature_c", 25.0))
        m_temp = float(m_cond.get("temperature_c", 25.0))
        delta_temp = m_temp - b_temp

        b_press = float(b_cond.get("pressure_atm", 1.0))
        m_press = float(m_cond.get("pressure_atm", 1.0))
        delta_press = m_press - b_press

        thermal_nature, thermal_basis_ar = cls._extract_thermal_nature(base_pred)
        gas_info = cls._analyze_gas_moles(base_pred["reactants"], base_pred["products"])
        ref_spec = base_pred.get("reference_conditions_spec") or {}

        temp_impact_ar, temp_supported = cls._analyze_temperature_change(
            delta_temp, m_temp, thermal_nature, ref_spec
        )
        press_impact_ar, press_supported = cls._analyze_pressure_change(
            delta_press, m_press, gas_info, ref_spec
        )

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
                    "نتيجة زيادة كمية المتفاعل المحدد المغذي للتفاعل."
                )
            elif p1["pct_change"] < 0:
                insights.append(
                    f"انخفضت كمية الناتج النظري ({p1['product_name_ar']}) بنسبة ({p1['pct_change']}%) "
                    "بسبب خفض كمية المتفاعل المحدد المستهلك أولاً."
                )
            else:
                insights.append("لم تتغير كمية الناتج النظري لأن التعديل طال مادة فائضة دون المساس بكمية المتفاعل المحدد.")

        insights.extend(temp_impact_ar)
        insights.extend(press_impact_ar)

        table_rows = [
            {
                "المتغير (Parameter)": "المتفاعل المحدد (Limiting Reagent)",
                "التجربة الأصلية (Baseline)": f"{base_lr.get('name_ar', '')} ({base_lr_formula})" if base_lr else "نسبة متكافئة تماماً",
                "التجربة المعدلة (Modified Scenario)": f"{mod_lr.get('name_ar', '')} ({mod_lr_formula})" if mod_lr else "نسبة متكافئة تماماً",
                "التغير والملاحظة": "تحول المتفاعل المحدد" if lr_shifted else "لا تغيير في الهوية",
            }
        ]

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
            "التغير والملاحظة": f"{'+' if delta_temp > 0 else ''}{delta_temp:.1f} °C ({thermal_basis_ar})",
        })

        gas_note = (
            f"Δn_gas = {gas_info['delta_n_gas']:+d}"
            if gas_info["known_states"]
            else "أطوار غير محددة"
        )
        table_rows.append({
            "المتغير (Parameter)": "الضغط (Pressure)",
            "التجربة الأصلية (Baseline)": f"{b_press} atm",
            "التجربة المعدلة (Modified Scenario)": f"{m_press} atm",
            "التغير والملاحظة": f"{'+' if delta_press > 0 else ''}{delta_press:.1f} atm ({gas_note})",
        })

        return {
            "success": True,
            "baseline_prediction": base_pred,
            "modified_prediction": mod_pred,
            "baseline_stoichiometry": base_stoich,
            "modified_stoichiometry": mod_stoich,
            "limiting_reagent_shifted": lr_shifted,
            "yield_comparisons": yield_comparisons,
            "condition_analysis": {
                "thermal_nature": thermal_nature,
                "thermal_basis_ar": thermal_basis_ar,
                "known_states": gas_info["known_states"],
                "reactant_gas_moles": gas_info["r_gas"],
                "product_gas_moles": gas_info["p_gas"],
                "delta_n_gas": gas_info["delta_n_gas"],
                "temperature_effect_supported": temp_supported,
                "pressure_effect_supported": press_supported,
            },
            "table_rows": table_rows,
            "insights_ar": " ".join(insights),
        }

"""
Reaction Visualizer module for MolecuLens Chemical Playground.
Generates 2D molecular cards for reactants and products, analyzes bond transformations,
audits atom conservation, and generates interactive Plotly reaction energy coordinate diagrams.
"""

from __future__ import annotations

import re
from typing import Any
import html
import plotly.graph_objects as go

from chemistry_engine.formula_analyzer import (
    format_chemical_formula_display,
    parse_chemical_formula,
)
from chemistry_engine.periodic_table import ELEMENT_NAMES_AR, ELEMENT_NAMES_EN
from search.local_catalog import LOCAL_SCIENTIFIC_CATALOG

try:
    from chemistry_engine.structure_analyzer import render_molecule_2d_svg
    RDKIT_AVAILABLE = True
except Exception:
    render_molecule_2d_svg = None
    RDKIT_AVAILABLE = False


class ReactionVisualizer:
    """Provides rich visual representations for chemical reactions."""

    BOND_ENERGY_DATABASE = {
        "H-H": 436,
        "C-H": 413,
        "C-C": 348,
        "C=C": 614,
        "C#C": 839,
        "O=O": 498,
        "O-H": 463,
        "C-O": 358,
        "C=O": 799,
        "N#N": 945,
        "N-H": 391,
        "H-Cl": 431,
        "H-Br": 366,
        "Cl-Cl": 242,
        "Br-Br": 193,
    }

    @classmethod
    def get_species_svg(cls, formula: str) -> str:
        """Finds or renders an SVG representation for a chemical species."""
        clean_formula = formula.strip()
        # Look up in local catalog for canonical smiles
        smiles = None
        for entity in LOCAL_SCIENTIFIC_CATALOG:
            if entity.get("symbol_or_formula") == clean_formula:
                smiles = entity.get("canonical_smiles") or entity.get("smiles")
                break

        if smiles and RDKIT_AVAILABLE and render_molecule_2d_svg:
            try:
                svg = render_molecule_2d_svg(smiles, width=180, height=140)
                if svg:
                    return svg
            except Exception:
                pass

        # Clean fallback atomic badge SVG
        f_disp = format_chemical_formula_display(clean_formula)
        return f"""
        <svg width="180" height="140" viewBox="0 0 180 140" xmlns="http://www.w3.org/2000/svg">
            <rect width="180" height="140" rx="12" fill="#1e293b" stroke="#334155" stroke-width="1.5"/>
            <circle cx="90" cy="60" r="34" fill="#0f4c5c" opacity="0.65"/>
            <text x="90" y="66" font-family="'Cairo', sans-serif" font-size="20" font-weight="bold" fill="#38bdf8" text-anchor="middle">{clean_formula}</text>
            <text x="90" y="112" font-family="'Cairo', sans-serif" font-size="12" fill="#94a3b8" text-anchor="middle">تمثيل صيغة جزيئية</text>
        </svg>
        """

    @classmethod
    def audit_atom_conservation(
        cls,
        reactants: list[dict[str, Any]],
        products: list[dict[str, Any]],
    ) -> list[dict[str, Any]]:
        """
        Audits atom conservation by counting each element on the reactant vs product side.
        Returns a detailed balance audit table.
        """
        reactant_counts: dict[str, int] = {}
        for r in reactants:
            coeff = int(r.get("coefficient", 1))
            parsed = parse_chemical_formula(r.get("formula", ""))
            for elem, cnt in parsed.items():
                reactant_counts[elem] = reactant_counts.get(elem, 0) + (cnt * coeff)

        product_counts: dict[str, int] = {}
        for p in products:
            coeff = int(p.get("coefficient", 1))
            parsed = parse_chemical_formula(p.get("formula", ""))
            for elem, cnt in parsed.items():
                product_counts[elem] = product_counts.get(elem, 0) + (cnt * coeff)

        all_elements = sorted(list(set(reactant_counts.keys()) | set(product_counts.keys())))
        audit_rows = []

        for elem in all_elements:
            r_cnt = reactant_counts.get(elem, 0)
            p_cnt = product_counts.get(elem, 0)
            is_conserved = (r_cnt == p_cnt)
            audit_rows.append({
                "element_symbol": elem,
                "element_name_ar": ELEMENT_NAMES_AR.get(elem, elem),
                "reactant_atoms": r_cnt,
                "product_atoms": p_cnt,
                "is_conserved": is_conserved,
                "status_ar": "محفوظ 100% (قانون حفظ الكتلة)" if is_conserved else "غير متطابق",
            })

        return audit_rows

    @classmethod
    def _extract_documented_enthalpy_and_direction(
        cls,
        reaction_data: dict[str, Any] | None,
        category: str = "",
    ) -> tuple[float | None, bool | None]:
        """
        Extracts documented Delta H (kJ) if explicitly present in verified reaction conditions,
        plus qualitative exothermic/endothermic nature when supported by reaction evidence.
        Returns (documented_delta_h_kj, is_exothermic_or_none).
        """
        documented_delta_h: float | None = None
        is_exothermic: bool | None = None

        text_blobs: list[str] = []
        if reaction_data:
            conds = reaction_data.get("conditions") or {}
            if isinstance(conds, dict):
                for key in ("thermodynamics_ar", "temperature", "notes"):
                    val = conds.get(key)
                    if val:
                        text_blobs.append(str(val))
            for obs in reaction_data.get("observations") or []:
                if isinstance(obs, dict) and obs.get("text_ar"):
                    text_blobs.append(str(obs["text_ar"]))
            if reaction_data.get("explanation_ar"):
                text_blobs.append(str(reaction_data["explanation_ar"]))

            thermo_text = str(conds.get("thermodynamics_ar", "")) if isinstance(conds, dict) else ""
            if thermo_text:
                m = re.search(r"ΔH°?[a-zA-Z]*\s*=\s*([+\-−]?\d+(?:\.\d+)?)\s*kJ", thermo_text)
                if m:
                    documented_delta_h = float(m.group(1).replace("−", "-"))
                    is_exothermic = documented_delta_h < 0

        joined = " ".join(text_blobs)
        if is_exothermic is None:
            if "طارد للحرارة" in joined or category == "combustion":
                is_exothermic = True
            elif "ماص للحرارة" in joined:
                is_exothermic = False

        return documented_delta_h, is_exothermic

    @classmethod
    def analyze_bond_transformations(
        cls,
        reactants: list[dict[str, Any]],
        products: list[dict[str, Any]],
        category: str = "",
        reaction_data: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        """
        Identifies bond transformations during a reaction and strictly separates
        genuinely calculated/documented thermodynamic values from qualitative
        illustrative descriptions. Never fabricates numeric fallback enthalpies.
        """
        r_formulas = [r.get("formula", "") for r in reactants]
        p_formulas = [p.get("formula", "") for p in products]

        documented_delta_h, qualitative_exothermic = cls._extract_documented_enthalpy_and_direction(
            reaction_data, category=category
        )

        bonds_broken: list[dict[str, Any]] = []
        bonds_formed: list[dict[str, Any]] = []
        complete_bde_coverage = False

        # Case 1: Complete covalent bond dissociation table coverage
        if set(r_formulas) == {"H2", "O2"} and set(p_formulas) == {"H2O"}:
            complete_bde_coverage = True
            bonds_broken = [
                {"bond": "H-H (تساهمية أحادية)", "count": 2, "energy_kj": 436, "total_kj": 872},
                {"bond": "O=O (تساهمية ثنائية)", "count": 1, "energy_kj": 498, "total_kj": 498},
            ]
            bonds_formed = [
                {"bond": "O-H (تساهمية قطبية في H2O)", "count": 4, "energy_kj": 463, "total_kj": 1852},
            ]
        elif set(r_formulas) == {"H2", "Cl2"} and set(p_formulas) == {"HCl"}:
            complete_bde_coverage = True
            bonds_broken = [
                {"bond": "H-H (تساهمية أحادية)", "count": 1, "energy_kj": 436, "total_kj": 436},
                {"bond": "Cl-Cl (تساهمية أحادية)", "count": 1, "energy_kj": 242, "total_kj": 242},
            ]
            bonds_formed = [
                {"bond": "H-Cl (تساهمية قطبية)", "count": 2, "energy_kj": 431, "total_kj": 862},
            ]
        elif set(r_formulas) == {"H2", "Br2"} and set(p_formulas) == {"HBr"}:
            complete_bde_coverage = True
            bonds_broken = [
                {"bond": "H-H (تساهمية أحادية)", "count": 1, "energy_kj": 436, "total_kj": 436},
                {"bond": "Br-Br (تساهمية أحادية)", "count": 1, "energy_kj": 193, "total_kj": 193},
            ]
            bonds_formed = [
                {"bond": "H-Br (تساهمية قطبية)", "count": 2, "energy_kj": 366, "total_kj": 732},
            ]
        elif set(r_formulas) == {"CH4", "O2"} and set(p_formulas) == {"CO2", "H2O"}:
            complete_bde_coverage = True
            bonds_broken = [
                {"bond": "C-H (تساهمية في CH4)", "count": 4, "energy_kj": 413, "total_kj": 1652},
                {"bond": "O=O (تساهمية ثنائية)", "count": 2, "energy_kj": 498, "total_kj": 996},
            ]
            bonds_formed = [
                {"bond": "C=O (تساهمية ثنائية في CO2)", "count": 2, "energy_kj": 799, "total_kj": 1598},
                {"bond": "O-H (تساهمية قطبية في H2O)", "count": 4, "energy_kj": 463, "total_kj": 1852},
            ]
        # Case 2: Qualitative ionic/aqueous/displacement mechanisms (no complete BDE balance)
        elif ("HCl" in r_formulas or "HNO3" in r_formulas or "CH3COOH" in r_formulas) and ("NaOH" in r_formulas or "KOH" in r_formulas):
            if qualitative_exothermic is None:
                qualitative_exothermic = True
            bonds_broken = [
                {"bond": "رابطة المانح البروتوني H-A", "count": 1, "description": "تحرر أيون الهيدروجين المميأ في الوسط المائي"},
                {"bond": "تآثر القاعدة مع OH-", "count": 1, "description": "توافر أيون الهيدروكسيد في المحلول المائي"},
            ]
            bonds_formed = [
                {"bond": "O-H (تكوين H2O)", "count": 1, "description": "اتحاد H+ و OH- لتكوين جزيء الماء المستقر (مع طاقات تميؤ الأيونات)"},
                {"bond": "أيونات الملح المميأة", "count": 1, "description": "بقاء الأيونات المتفرجة مميأة في المحلول المائي"},
            ]
        elif "AgCl" in p_formulas or "BaSO4" in p_formulas:
            bonds_broken = [
                {"bond": "تآثرات الإماهة المائية للأيونات", "count": 2, "description": "تحرر الأيونات من أغلفة التميؤ المائية"},
            ]
            bonds_formed = [
                {"bond": "شبكة أيونية بلورية في الراسب", "count": 1, "description": "تكوّن شبيكة بلورية صلبة شحيحة الذوبان"},
            ]
        elif "H2" in p_formulas and any(r in r_formulas for r in ["Zn", "Fe", "Mg"]):
            if qualitative_exothermic is None:
                qualitative_exothermic = True
            bonds_broken = [
                {"bond": "H-A (في الحمض) + الرابطة الفلزية للمعدن", "count": 2, "description": "أكسدة ذرات الفلز واختزال البروتونات"},
            ]
            bonds_formed = [
                {"bond": "H-H (تساهمية في H2) + كاتيون مميأ", "count": 2, "description": "تصاعد غاز الهيدروجين وتكوّن أيون الفلز المميأ في المحلول"},
            ]
        else:
            bonds_broken = [
                {"bond": "الروابط الأولية في المتفاعلات", "count": len(reactants), "description": "تفكك أو إعادة ترتيب الروابط في المواد المتفاعلة"},
            ]
            bonds_formed = [
                {"bond": "الروابط المتكونة في النواتج", "count": len(products), "description": "تشكل روابط المواد الناتجة"},
            ]

        total_broken_energy: int | None = None
        total_formed_energy: int | None = None
        est_delta_h: float | None = None
        is_exothermic: bool | None = qualitative_exothermic
        is_illustrative = True
        calculation_source = "illustrative_qualitative"

        if documented_delta_h is not None:
            is_illustrative = False
            calculation_source = "documented_reference"
            est_delta_h = round(documented_delta_h, 2)
            is_exothermic = est_delta_h < 0
            if complete_bde_coverage:
                total_broken_energy = sum(int(b.get("total_kj", 0)) for b in bonds_broken)
                total_formed_energy = sum(int(b.get("total_kj", 0)) for b in bonds_formed)
        elif complete_bde_coverage:
            is_illustrative = False
            calculation_source = "bond_dissociation_table"
            total_broken_energy = sum(int(b.get("total_kj", 0)) for b in bonds_broken)
            total_formed_energy = sum(int(b.get("total_kj", 0)) for b in bonds_formed)
            est_delta_h = float(total_broken_energy - total_formed_energy)
            is_exothermic = est_delta_h < 0

        if not is_illustrative:
            if calculation_source == "documented_reference":
                data_classification_ar = "بيانات طاقة موثقة مرجعياً (Documented Reference Enthalpy — ΔH° من السجل المعتمد)"
                summary_ar = (
                    f"قيمة التغير في المحتوى الحراري الموثقة مرجعياً: ΔH° = {est_delta_h} kJ "
                    f"({'تفاعل طارد للحرارة' if is_exothermic else 'تفاعل ماص للحرارة'}). "
                    "تفصيل الروابط أدناه يوضح التحولات البنيوية الأساسية؛ أما طاقة التنشيط Ea فتحتاج قياساً حركياً مستقلاً."
                )
            else:
                data_classification_ar = "بيانات طاقة محسوبة من جدول طاقات تفكك الروابط التساهمية القياسي (Calculated from Bond Dissociation Energies)"
                summary_ar = (
                    f"تقدير محسوب من متوسط طاقات الروابط القياسية: "
                    f"طاقة الروابط المتكسرة ({total_broken_energy} kJ/mol) مقابل طاقة الروابط المتكونة ({total_formed_energy} kJ/mol)، "
                    f"ما يعطي ΔH ≈ {est_delta_h:.1f} kJ/mol ({'طارد للحرارة' if is_exothermic else 'ماص للحرارة'} في الطور الغازي القياسي)."
                )
            thermodynamic_type_ar = (
                "طارد للحرارة (Exothermic — محسوب/موثق)"
                if is_exothermic
                else "ماص للحرارة (Endothermic — محسوب/موثق)"
            )
        else:
            data_classification_ar = (
                "تمثيل توضيحي كيفي فقط (Illustrative Qualitative Only) — "
                "لا تتوفر طاقات روابط مكتملة (مثل طاقات الشبيكة البلورية أو الإماهة) ولا ΔH° موثق لحساب الطاقة رقمياً"
            )
            if is_exothermic is True:
                thermodynamic_type_ar = "طارد للحرارة كيفياً (Qualitative Exothermic — دون قيمة ΔH رقمية)"
                summary_ar = (
                    "وصف توضيحي كيفي للروابط فقط: التفاعل معروف كيفياً بأنه طارد للحرارة، "
                    "لكن لا تتوفر بيانات كمية مكتملة (كطاقات الشبيكة الأيونية أو الإماهة) لحساب قيمة ΔH رقمياً هنا."
                )
            elif is_exothermic is False:
                thermodynamic_type_ar = "ماص للحرارة كيفياً (Qualitative Endothermic — دون قيمة ΔH رقمية)"
                summary_ar = (
                    "وصف توضيحي كيفي للروابط فقط: التفاعل معروف كيفياً بأنه ماص للحرارة، "
                    "لكن لا تتوفر بيانات كمية مكتملة لحساب قيمة ΔH رقمياً."
                )
            else:
                thermodynamic_type_ar = "غير محدد حرارياً (Thermodynamic Direction Unknown)"
                summary_ar = (
                    "تمثيل توضيحي كيفي لتغير الروابط فقط: لا تتوفر في قاعدة البيانات الحالية قيم إنثالبي موثقة "
                    "أو طاقات روابط مكتملة لهذا التفاعل لتحديد ΔH أو اتجاه التغير الحراري رقمياً."
                )

        return {
            "bonds_broken": bonds_broken,
            "bonds_formed": bonds_formed,
            "total_broken_energy_kj": total_broken_energy,
            "total_formed_energy_kj": total_formed_energy,
            "estimated_delta_h_kj": est_delta_h,
            "is_exothermic": is_exothermic,
            "is_illustrative": is_illustrative,
            "calculation_source": calculation_source,
            "data_classification_ar": data_classification_ar,
            "thermodynamic_type_ar": thermodynamic_type_ar,
            "summary_ar": summary_ar,
        }

    @classmethod
    def create_reaction_energy_diagram(
        cls,
        reactants_str: str,
        products_str: str,
        is_exothermic: bool | None = True,
        delta_h_kj: float | None = None,
        activation_energy_kj: float | None = None,
        is_illustrative: bool = False,
    ) -> go.Figure:
        """
        Builds an interactive Plotly reaction coordinate energy profile.
        Strictly separates quantified ΔH/Ea values from illustrative schematic diagrams
        and never fabricates numeric fallback enthalpies or activation barriers.
        """
        has_quant_dh = (delta_h_kj is not None) and (not is_illustrative)
        has_quant_ea = activation_energy_kj is not None

        if has_quant_dh and delta_h_kj is not None:
            eff_exothermic: bool | None = delta_h_kj < 0
            e_reactants = 0.0
            e_products = float(delta_h_kj)
            if has_quant_ea and activation_energy_kj is not None:
                e_transition = e_reactants + float(activation_energy_kj)
            else:
                # Visual barrier placement relative to max level without fabricating numeric Ea
                e_transition = max(e_reactants, e_products) + max(40.0, abs(delta_h_kj) * 0.25)
        else:
            eff_exothermic = is_exothermic
            e_reactants = 50.0
            if eff_exothermic is True:
                e_products = 20.0
                e_transition = 85.0
            elif eff_exothermic is False:
                e_products = 75.0
                e_transition = 95.0
            else:
                e_products = 50.0
                e_transition = 85.0

        span_r_to_ts = e_transition - e_reactants
        span_ts_to_p = e_transition - e_products
        xs = [0.0, 0.8, 1.2, 1.8, 2.3, 2.5, 2.7, 3.2, 3.8, 4.2, 5.0]
        ys = [
            e_reactants,
            e_reactants,
            e_reactants + span_r_to_ts * 0.15,
            e_reactants + span_r_to_ts * 0.72,
            e_transition - span_r_to_ts * 0.05,
            e_transition,
            e_transition - span_ts_to_p * 0.05,
            e_products + span_ts_to_p * 0.45,
            e_products + span_ts_to_p * 0.08,
            e_products,
            e_products,
        ]

        fig = go.Figure()

        fig.add_trace(
            go.Scatter(
                x=xs,
                y=ys,
                mode="lines",
                line=dict(color="#0f4c5c", width=4, shape="spline"),
                name="مسار طاقة الوضع (Potential Energy)",
                hoverinfo="skip",
            )
        )

        fig.add_trace(
            go.Scatter(
                x=[0.4],
                y=[e_reactants],
                mode="markers+text",
                marker=dict(size=12, color="#2563eb", symbol="circle"),
                text=[f"المتفاعلات<br>{reactants_str}"],
                textposition="bottom center",
                name="المتفاعلات (Reactants)",
            )
        )

        if has_quant_ea and activation_energy_kj is not None:
            ts_text = f"الحالة الانتقالية / المعقد النشط [‡]<br>طاقة التنشيط Ea = {activation_energy_kj:.1f} kJ"
        else:
            ts_text = "الحالة الانتقالية / المعقد النشط [‡]<br>حاجز التنشيط Ea (تخطيطي توضيحي — غير مقيس)"

        fig.add_trace(
            go.Scatter(
                x=[2.5],
                y=[e_transition],
                mode="markers+text",
                marker=dict(size=14, color="#dc2626", symbol="diamond"),
                text=[ts_text],
                textposition="top center",
                name="المعقد النشط (Transition State)",
            )
        )

        fig.add_trace(
            go.Scatter(
                x=[4.6],
                y=[e_products],
                mode="markers+text",
                marker=dict(size=12, color="#059669", symbol="square"),
                text=[f"النواتج<br>{products_str}"],
                textposition="bottom center" if eff_exothermic is True else "top center",
                name="النواتج (Products)",
            )
        )

        fig.add_shape(
            type="line",
            x0=3.8, y0=e_reactants,
            x1=4.6, y1=e_reactants,
            line=dict(color="#94a3b8", dash="dot", width=1.5),
        )

        if has_quant_dh and delta_h_kj is not None:
            enthalpy_text = f"ΔH = {delta_h_kj:.1f} kJ ({'طارد' if eff_exothermic else 'ماص'} — محسوب/موثق)"
            if has_quant_ea:
                chart_title = "مخطط مسار التفاعل ومستويات الطاقة (قيم موثقة/محسوبة)"
                yaxis_title = "طاقة الوضع النسبية (kJ)"
            else:
                chart_title = "مخطط مسار التفاعل (ΔH موثقة/محسوبة؛ حاجز التنشيط Ea توضيحي فقط)"
                yaxis_title = "طاقة الوضع النسبية (ΔH محسوبة؛ قمة Ea تخطيطية)"
            show_y_ticks = has_quant_ea
        else:
            if eff_exothermic is True:
                enthalpy_text = "ΔH < 0 (طارد كيفياً — تمثيل تخطيطي توضيحي بدون رقم)"
            elif eff_exothermic is False:
                enthalpy_text = "ΔH > 0 (ماص كيفياً — تمثيل تخطيطي توضيحي بدون رقم)"
            else:
                enthalpy_text = "ΔH غير محددة رقمياً (تمثيل تخطيطي توضيحي فقط)"
            chart_title = "مخطط توضيحي تخطيطي لمسار التفاعل (Illustrative Schematic — ليست قيماً رقمية)"
            yaxis_title = "مستوى الطاقة النسبي (تخطيطي توضيحي — غير رقمي)"
            show_y_ticks = False

        fig.add_annotation(
            x=4.2,
            y=(e_reactants + e_products) / 2.0,
            text=enthalpy_text,
            showarrow=True,
            arrowhead=2,
            arrowcolor="#0f766e",
            font=dict(color="#0f766e", size=12, family="Cairo"),
        )

        fig.update_layout(
            title=dict(
                text=chart_title,
                font=dict(family="Cairo, sans-serif", size=14, color="#0f4c5c"),
                x=0.5,
            ),
            xaxis=dict(
                title="إحداثي مسار التفاعل (Reaction Coordinate)",
                showticklabels=False,
                zeroline=False,
                gridcolor="#f1f5f9",
            ),
            yaxis=dict(
                title=yaxis_title,
                showticklabels=show_y_ticks,
                zeroline=False,
                gridcolor="#f1f5f9",
            ),
            plot_bgcolor="#ffffff",
            paper_bgcolor="#ffffff",
            showlegend=True,
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
            margin=dict(l=40, r=40, t=65, b=40),
            height=390,
        )

        return fig

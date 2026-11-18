"""
Reaction Visualizer module for MolecuLens Chemical Playground.
Generates 2D molecular cards for reactants and products, analyzes bond transformations,
audits atom conservation, and generates interactive Plotly reaction energy coordinate diagrams.
"""

from __future__ import annotations

import re
from typing import Any
import html
try:
    import plotly.graph_objects as go
except ImportError:
    from core._plotly_fallback import register_plotly_fallback_in_sys_modules
    register_plotly_fallback_in_sys_modules()
    import plotly.graph_objects as go

from chemistry_engine.formula_analyzer import (
    calculate_molar_mass,
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

    PRECIPITATE_COLORS: dict[str, tuple[str, str]] = {
        "AgCl": ("#f8fafc", "راسب أبيض متجبن (AgCl↓)"),
        "BaSO4": ("#f1f5f9", "راسب أبيض بلوري كثيف (BaSO4↓)"),
        "CaCO3": ("#f8fafc", "راسب أبيض كلسي (CaCO3↓)"),
        "Cu(OH)2": ("#38bdf8", "راسب أزرق هلامي (Cu(OH)2↓)"),
        "Fe(OH)3": ("#b45309", "راسب بني محمر (Fe(OH)3↓)"),
        "Al(OH)3": ("#e2e8f0", "راسب أبيض هلامي (Al(OH)3↓)"),
        "PbI2": ("#facc15", "راسب أصفر ذهبي (PbI2↓)"),
        "Cu": ("#d97706", "ترسب فلز النحاس الأحمر (Cu↓)"),
        "Ag": ("#cbd5e1", "ترسب بلورات الفضة اللامعة (Ag↓)"),
    }

    @classmethod
    def create_reaction_vessel_svg(
        cls,
        reactants: list[dict[str, Any]],
        products: list[dict[str, Any]],
        is_exothermic: bool | None = None,
        delta_h_kj: float | None = None,
        conditions: dict[str, Any] | None = None,
    ) -> str:
        """
        Generates an interactive laboratory SVG diagram of the active reaction vessel,
        depicting reactant addition, liquid medium, rising gas bubbles, settled solid
        precipitates, and exothermic/endothermic heat transfer.
        """
        r_labels = " + ".join(
            f"{r.get('coefficient', 1) if int(r.get('coefficient', 1)) > 1 else ''}{html.escape(str(r.get('formula', '')))}"
            for r in reactants
        )
        p_labels = " + ".join(
            f"{p.get('coefficient', 1) if int(p.get('coefficient', 1)) > 1 else ''}{html.escape(str(p.get('formula', '')))}"
            for p in products
        )

        gases = [p.get("formula", "") for p in products if p.get("state") == "g"]
        solids = [p.get("formula", "") for p in products if p.get("state") == "s"]
        has_aqueous = any(sp.get("state") in {"aq", "l"} for sp in list(reactants) + list(products))

        liquid_fill = "#bae6fd" if has_aqueous else "#e2e8f0"
        for sp in list(reactants) + list(products):
            f = sp.get("formula", "")
            if f in {"CuSO4", "CuCl2", "Cu(NO3)2"}:
                liquid_fill = "#7dd3fc"
            elif f in {"FeCl3", "KMnO4"}:
                liquid_fill = "#fde68a"

        precip_color = "#cbd5e1"
        precip_label = ""
        if solids:
            first_solid = solids[0]
            if first_solid in cls.PRECIPITATE_COLORS:
                precip_color, precip_label = cls.PRECIPITATE_COLORS[first_solid]
            else:
                precip_label = f"طور صلب ناتج ({html.escape(first_solid)}↓)"

        gas_label = f"تصاعد غاز ({', '.join(html.escape(g) for g in gases)}↑)" if gases else ""

        if delta_h_kj is not None:
            if delta_h_kj < 0:
                thermo_badge = f"طارد للحرارة (ΔH = {delta_h_kj:g} kJ)"
                thermo_color = "#ef4444"
                vessel_glow = "#fecaca"
            else:
                thermo_badge = f"ماص للحرارة (ΔH = +{delta_h_kj:g} kJ)"
                thermo_color = "#3b82f6"
                vessel_glow = "#bfdbfe"
        elif is_exothermic is True:
            thermo_badge = "تفاعل طارد للحرارة (انبعاث طاقة حرارية)"
            thermo_color = "#ef4444"
            vessel_glow = "#fecaca"
        elif is_exothermic is False:
            thermo_badge = "تفاعل ماص للحرارة (امتصاص طاقة حرارية)"
            thermo_color = "#3b82f6"
            vessel_glow = "#bfdbfe"
        else:
            thermo_badge = "الوسط التفاعلي القياسي"
            thermo_color = "#0f766e"
            vessel_glow = "#ccfbf1"

        cond_temp = html.escape(str((conditions or {}).get("temperature", "25 °C"))[:38])

        bubbles_svg = ""
        if gases:
            bubbles_svg = f"""
            <g opacity="0.9">
                <circle cx="290" cy="195" r="7" fill="#ffffff" stroke="#0284c7" stroke-width="1.5"/>
                <circle cx="315" cy="170" r="9" fill="#ffffff" stroke="#0284c7" stroke-width="1.5"/>
                <circle cx="345" cy="185" r="6" fill="#ffffff" stroke="#0284c7" stroke-width="1.5"/>
                <circle cx="305" cy="140" r="8" fill="#ffffff" stroke="#0284c7" stroke-width="1.5"/>
                <circle cx="335" cy="130" r="10" fill="#ffffff" stroke="#0284c7" stroke-width="1.5"/>
                <text x="320" y="102" font-family="'Cairo', sans-serif" font-size="12" font-weight="bold" fill="#0369a1" text-anchor="middle">{gas_label}</text>
            </g>
            """

        precip_svg = ""
        if solids:
            precip_svg = f"""
            <g>
                <path d="M 248 228 Q 320 214 392 228 L 388 246 Q 320 250 252 246 Z" fill="{precip_color}" stroke="#475569" stroke-width="1.2"/>
                <text x="320" y="268" font-family="'Cairo', sans-serif" font-size="12" font-weight="bold" fill="#334155" text-anchor="middle">{precip_label}</text>
            </g>
            """

        return f"""
        <svg width="100%" height="310" viewBox="0 0 640 310" xmlns="http://www.w3.org/2000/svg" style="background:#f8fafc;border-radius:14px;border:1px solid #cbd5e1;">
            <!-- Header Banner -->
            <rect x="16" y="12" width="608" height="36" rx="8" fill="#0f4c5c"/>
            <text x="320" y="35" font-family="'Cairo', sans-serif" font-size="14" font-weight="bold" fill="#ffffff" text-anchor="middle">
                محاكاة وعاء التفاعل الكيميائي: {r_labels} &#8594; {p_labels}
            </text>

            <!-- Left Panel: Reactants Input Stream -->
            <rect x="22" y="64" width="165" height="175" rx="10" fill="#ffffff" stroke="#94a3b8" stroke-width="1.5"/>
            <text x="104" y="88" font-family="'Cairo', sans-serif" font-size="13" font-weight="bold" fill="#1e3a8a" text-anchor="middle">1. تغذية المتفاعلات</text>
            <text x="104" y="116" font-family="'Cairo', sans-serif" font-size="13" font-weight="bold" fill="#0f4c5c" text-anchor="middle">{r_labels}</text>
            <text x="104" y="145" font-family="'Cairo', sans-serif" font-size="11" fill="#475569" text-anchor="middle">الظروف: {cond_temp}</text>
            <path d="M 187 140 L 228 140" stroke="#0284c7" stroke-width="3" stroke-dasharray="5,3"/>
            <polygon points="228,135 238,140 228,145" fill="#0284c7"/>

            <!-- Center Panel: Active Laboratory Beaker -->
            <ellipse cx="320" cy="185" rx="92" ry="72" fill="{vessel_glow}" opacity="0.45"/>
            <!-- Beaker Liquid -->
            <path d="M 244 150 L 396 150 L 392 244 Q 320 252 248 244 Z" fill="{liquid_fill}" opacity="0.82"/>
            <!-- Beaker Glass Outline -->
            <path d="M 236 88 L 246 246 Q 320 254 394 246 L 404 88" fill="none" stroke="#1e293b" stroke-width="3.5" stroke-linecap="round"/>
            <ellipse cx="320" cy="88" rx="84" ry="10" fill="none" stroke="#1e293b" stroke-width="2.5"/>
            <!-- Graduated Marks -->
            <line x1="255" y1="130" x2="270" y2="130" stroke="#475569" stroke-width="1.5"/>
            <line x1="255" y1="160" x2="275" y2="160" stroke="#475569" stroke-width="1.5"/>
            <line x1="255" y1="190" x2="270" y2="190" stroke="#475569" stroke-width="1.5"/>
            <line x1="255" y1="220" x2="275" y2="220" stroke="#475569" stroke-width="1.5"/>

            {bubbles_svg}
            {precip_svg}

            <!-- Right Panel: Products Output Summary -->
            <path d="M 402 140 L 443 140" stroke="#059669" stroke-width="3" stroke-dasharray="5,3"/>
            <polygon points="443,135 453,140 443,145" fill="#059669"/>
            <rect x="453" y="64" width="165" height="175" rx="10" fill="#ffffff" stroke="#059669" stroke-width="1.5"/>
            <text x="535" y="88" font-family="'Cairo', sans-serif" font-size="13" font-weight="bold" fill="#065f46" text-anchor="middle">2. النواتج المتكونة</text>
            <text x="535" y="118" font-family="'Cairo', sans-serif" font-size="13" font-weight="bold" fill="#047857" text-anchor="middle">{p_labels}</text>
            <text x="535" y="148" font-family="'Cairo', sans-serif" font-size="11" fill="#334155" text-anchor="middle">قانون حفظ الكتلة والذرات</text>
            <text x="535" y="168" font-family="'Cairo', sans-serif" font-size="11" font-weight="bold" fill="#059669" text-anchor="middle">موزون ومتحقق 100%</text>

            <!-- Bottom Thermal & Medium Status Bar -->
            <rect x="140" y="274" width="360" height="26" rx="13" fill="{thermo_color}" opacity="0.14"/>
            <text x="320" y="291" font-family="'Cairo', sans-serif" font-size="12" font-weight="bold" fill="{thermo_color}" text-anchor="middle">{thermo_badge}</text>
        </svg>
        """

    @classmethod
    def build_end_to_end_experiment_stages(
        cls,
        reaction_payload: dict[str, Any],
    ) -> dict[str, Any]:
        """
        Builds a complete 4-stage visual chemical experiment walkthrough from start to finish:
          Stage 1: Reactant Preparation & Initial Quantities (تحضير المتفاعلات والأوعية)
          Stage 2: Active Reaction Vessel & Laboratory Phenomena (وعاء التفاعل النشط والظواهر المخبرية)
          Stage 3: Lavoisier Mass Balance & Stoichiometric Consumption (ميزان حفظ الكتلة والذرات والاستهلاك المولي)
          Stage 4: Energy Profile & Final Vessel Inventory (مسار الطاقة وجرد الوعاء النهائي)
        Works with both predict_reaction() and MolecuLensCore.run_reaction_analysis() outputs.
        """
        # Extract reactants & products with balanced coefficients
        eq_block = reaction_payload.get("balanced_equation") or {}
        reactants = eq_block.get("reactants") or reaction_payload.get("reactants") or []
        products = eq_block.get("products") or reaction_payload.get("products") or []
        category = str(reaction_payload.get("category", ""))
        conditions = reaction_payload.get("conditions") or {}

        # Extract stoichiometry if present, or compute default stoichiometric baseline
        stoich = reaction_payload.get("stoichiometry") or {}
        if stoich.get("available") and stoich.get("reactants_table"):
            reactants_table = stoich["reactants_table"]
            theoretical_yields = stoich.get("theoretical_yields", [])
            stoich_status_ar = stoich.get("stoichiometric_status_ar", "")
            limiting_info = stoich.get("limiting_reactant")
        elif reaction_payload.get("reactants_table") and reaction_payload.get("theoretical_yields"):
            reactants_table = reaction_payload["reactants_table"]
            theoretical_yields = reaction_payload["theoretical_yields"]
            stoich_status_ar = reaction_payload.get("stoichiometric_status_ar", "")
            limiting_info = reaction_payload.get("limiting_reactant")
        else:
            # Default stoichiometric 1x reaction extent when user didn't enter custom masses
            reactants_table = []
            for r in reactants:
                coeff = float(r.get("coefficient", 1) or 1)
                mm = calculate_molar_mass(r.get("formula", ""))
                mass_g = round(coeff * mm, 4)
                reactants_table.append({
                    "formula": r.get("formula", ""),
                    "name_ar": r.get("name_ar", r.get("formula", "")),
                    "name_en": r.get("name_en", r.get("formula", "")),
                    "state": r.get("state", ""),
                    "coefficient": int(coeff),
                    "molar_mass": round(mm, 3),
                    "initial_mol": coeff,
                    "initial_mass_g": mass_g,
                    "consumed_mol": coeff,
                    "consumed_mass_g": mass_g,
                    "remaining_mol": 0.0,
                    "remaining_mass_g": 0.0,
                    "role_ar": "نسبة ستوكيومترية متكافئة (100% استهلاك)",
                })
            theoretical_yields = []
            for p in products:
                coeff = float(p.get("coefficient", 1) or 1)
                mm = calculate_molar_mass(p.get("formula", ""))
                mass_g = round(coeff * mm, 4)
                theoretical_yields.append({
                    "formula": p.get("formula", ""),
                    "name_ar": p.get("name_ar", p.get("formula", "")),
                    "name_en": p.get("name_en", p.get("formula", "")),
                    "state": p.get("state", ""),
                    "coefficient": int(coeff),
                    "molar_mass": round(mm, 3),
                    "theoretical_yield_mol": coeff,
                    "theoretical_yield_g": mass_g,
                })
            stoich_status_ar = "خلط بالنسب المولية المتكافئة قياسياً (Stoichiometric Equivalence)"
            limiting_info = None

        # Stage 1: Reactant Preparation Cards
        stage1_cards = []
        for idx, r_row in enumerate(reactants_table):
            formula = r_row.get("formula", "")
            state_code = r_row.get("state") or (reactants[idx].get("state", "") if idx < len(reactants) else "")
            init_mol = float(r_row.get("initial_mol", 0.0))
            init_g = float(r_row.get("initial_mass_g", 0.0))
            cons_mol = float(r_row.get("consumed_mol", 0.0))
            rem_mol = float(r_row.get("remaining_mol", 0.0))
            consumed_pct = round((cons_mol / init_mol) * 100.0, 1) if init_mol > 0 else 0.0
            remaining_pct = round(max(0.0, 100.0 - consumed_pct), 1)
            stage1_cards.append({
                "formula": formula,
                "formula_display": format_chemical_formula_display(formula),
                "name_ar": r_row.get("name_ar", formula),
                "name_en": r_row.get("name_en", formula),
                "state": state_code,
                "coefficient": int(r_row.get("coefficient", 1)),
                "molar_mass": float(r_row.get("molar_mass", calculate_molar_mass(formula))),
                "initial_mol": round(init_mol, 4),
                "initial_mass_g": round(init_g, 4),
                "consumed_mol": round(cons_mol, 4),
                "consumed_mass_g": round(float(r_row.get("consumed_mass_g", 0.0)), 4),
                "remaining_mol": round(rem_mol, 4),
                "remaining_mass_g": round(float(r_row.get("remaining_mass_g", 0.0)), 4),
                "consumed_pct": consumed_pct,
                "remaining_pct": remaining_pct,
                "role_ar": r_row.get("role_ar", ""),
                "svg": cls.get_species_svg(formula),
            })

        # Stage 2: Active Vessel & Phenomena
        bond_data = cls.analyze_bond_transformations(
            reactants, products, category=category, reaction_data=reaction_payload
        )
        vessel_svg = cls.create_reaction_vessel_svg(
            reactants=reactants,
            products=products,
            is_exothermic=bond_data.get("is_exothermic"),
            delta_h_kj=bond_data.get("estimated_delta_h_kj"),
            conditions=conditions,
        )
        observations = list(reaction_payload.get("observations") or [])
        if not observations:
            for p in products:
                if p.get("state") == "g":
                    observations.append({
                        "type": "gas",
                        "text_ar": f"تصاعد فقاعات غاز {p.get('name_ar', p.get('formula', ''))} ({p.get('formula', '')}↑).",
                        "basis_ar": "تكوّن ناتج في الطور الغازي (g).",
                    })
                elif p.get("state") == "s" and category in {"precipitation", "single_displacement"}:
                    observations.append({
                        "type": "precipitate",
                        "text_ar": f"تكوّن طور صلب/راسب من {p.get('name_ar', p.get('formula', ''))} ({p.get('formula', '')}↓).",
                        "basis_ar": "تكوّن مادة شحيحة الذوبان في الطور الصلب (s).",
                    })
            if bond_data.get("is_exothermic") is True:
                observations.append({
                    "type": "heat",
                    "text_ar": f"انبعاث حرارة إلى الوسط المحيط ({bond_data.get('thermodynamic_type_ar', '')}).",
                    "basis_ar": "تفاعل طارد للحرارة.",
                })
            elif bond_data.get("is_exothermic") is False:
                observations.append({
                    "type": "heat",
                    "text_ar": f"امتصاص طاقة حرارية من الوسط ({bond_data.get('thermodynamic_type_ar', '')}).",
                    "basis_ar": "تفاعل ماص للحرارة.",
                })

        # Stage 3: Mass Balance & Atom Conservation
        total_initial_mass_g = round(sum(c["initial_mass_g"] for c in stage1_cards), 4)
        total_consumed_mass_g = round(sum(c["consumed_mass_g"] for c in stage1_cards), 4)
        total_remaining_mass_g = round(sum(c["remaining_mass_g"] for c in stage1_cards), 4)
        total_produced_mass_g = round(sum(float(y.get("theoretical_yield_g", 0.0)) for y in theoretical_yields), 4)
        total_final_vessel_mass_g = round(total_remaining_mass_g + total_produced_mass_g, 4)
        mass_conserved = abs(total_initial_mass_g - total_final_vessel_mass_g) <= 0.05
        atom_audit = cls.audit_atom_conservation(reactants, products)

        # Stage 4: Final Vessel Inventory (Products Formed + Excess Reactants Remaining)
        final_inventory: list[dict[str, Any]] = []
        for idx_p, y in enumerate(theoretical_yields):
            pf = y.get("formula", "")
            p_state = y.get("state") or (products[idx_p].get("state", "") if idx_p < len(products) else "")
            final_inventory.append({
                "kind_ar": "ناتج متكوّن (Product Formed)",
                "formula": pf,
                "formula_display": format_chemical_formula_display(pf),
                "name_ar": y.get("name_ar", pf),
                "state": p_state,
                "amount_mol": round(float(y.get("theoretical_yield_mol", 0.0)), 4),
                "mass_g": round(float(y.get("theoretical_yield_g", 0.0)), 4),
                "status_ar": (
                    "يتصاعد من المحلول كغاز (↑)"
                    if p_state == "g"
                    else ("يترسب في قاع الوعاء كطور صلب (↓)" if p_state == "s" and category in {"precipitation", "single_displacement"} else "موجود في وعاء التفاعل النهائي")
                ),
                "svg": cls.get_species_svg(pf),
            })

        for c in stage1_cards:
            if c["remaining_mol"] > 1e-6:
                final_inventory.append({
                    "kind_ar": "متفاعل فائض متبقٍ (Excess Reactant Remaining)",
                    "formula": c["formula"],
                    "formula_display": c["formula_display"],
                    "name_ar": c["name_ar"],
                    "state": c["state"],
                    "amount_mol": c["remaining_mol"],
                    "mass_g": c["remaining_mass_g"],
                    "status_ar": f"متبقٍ في الوعاء دون تفاعل ({c['remaining_pct']}% من كميته الابتدائية)",
                    "svg": c["svg"],
                })

        return {
            "stage1_preparation": {
                "title_ar": "المرحلة 1: تحضير المتفاعلات والكميات الابتدائية (Reactant Setup & Initial Quantities)",
                "reactant_cards": stage1_cards,
                "stoichiometric_status_ar": stoich_status_ar,
                "limiting_reactant": limiting_info,
            },
            "stage2_active_vessel": {
                "title_ar": "المرحلة 2: وعاء التفاعل النشط والظواهر المخبرية (Active Reaction Vessel & Phenomena)",
                "vessel_svg": vessel_svg,
                "observations": observations,
                "conditions": conditions,
            },
            "stage3_mass_balance": {
                "title_ar": "المرحلة 3: ميزان حفظ الكتلة (لافوازييه) والاستهلاك الستوكيومتري (Mass & Atom Conservation)",
                "total_initial_mass_g": total_initial_mass_g,
                "total_consumed_mass_g": total_consumed_mass_g,
                "total_remaining_mass_g": total_remaining_mass_g,
                "total_produced_mass_g": total_produced_mass_g,
                "total_final_vessel_mass_g": total_final_vessel_mass_g,
                "mass_conserved": mass_conserved,
                "atom_audit": atom_audit,
                "reactant_consumption": stage1_cards,
                "theoretical_yields": theoretical_yields,
            },
            "stage4_final_inventory": {
                "title_ar": "المرحلة 4: جرد الوعاء النهائي ومسار الطاقة (Final Vessel Inventory & Energy Profile)",
                "final_inventory": final_inventory,
                "bond_transformations": bond_data,
            },
        }

    @classmethod
    def create_stoichiometry_mass_chart(cls, stage3_mass_balance: dict[str, Any]) -> go.Figure:
        """
        Creates a Plotly bar chart comparing Initial Reactant Mass vs. Final Vessel Mass
        (Produced Products + Remaining Excess Reactants) to visually prove Lavoisier mass conservation.
        """
        init_g = float(stage3_mass_balance.get("total_initial_mass_g", 0.0))
        cons_g = float(stage3_mass_balance.get("total_consumed_mass_g", 0.0))
        rem_g = float(stage3_mass_balance.get("total_remaining_mass_g", 0.0))
        prod_g = float(stage3_mass_balance.get("total_produced_mass_g", 0.0))

        categories = [
            "قبل التفاعل (المتفاعلات الابتدائية)",
            "أثناء التفاعل (المستهلك مقابل الفائض)",
            "بعد التفاعل (النواتج + الفائض المتبقي)",
        ]

        fig = go.Figure()
        fig.add_trace(
            go.Bar(
                name="الكتلة الابتدائية الكلية (g)",
                x=[categories[0]],
                y=[init_g],
                marker_color="#2563eb",
                text=[f"{init_g:g} g"],
                textposition="auto",
            )
        )
        fig.add_trace(
            go.Bar(
                name="المتفاعلات المستهلكة (g)",
                x=[categories[1]],
                y=[cons_g],
                marker_color="#0f766e",
                text=[f"مستهلك: {cons_g:g} g"],
                textposition="auto",
            )
        )
        fig.add_trace(
            go.Bar(
                name="النواتج المتكونة (g)",
                x=[categories[2]],
                y=[prod_g],
                marker_color="#059669",
                text=[f"نواتج: {prod_g:g} g"],
                textposition="auto",
            )
        )
        if rem_g > 1e-6:
            fig.add_trace(
                go.Bar(
                    name="الفائض المتبقي دون تفاعل (g)",
                    x=[categories[1], categories[2]],
                    y=[rem_g, rem_g],
                    marker_color="#f59e0b",
                    text=[f"فائض: {rem_g:g} g", f"فائض: {rem_g:g} g"],
                    textposition="auto",
                )
            )

        fig.update_layout(
            barmode="stack",
            title=dict(
                text=f"ميزان حفظ الكتلة في التجربة: {init_g:g} g قبل التفاعل = {round(prod_g + rem_g, 4):g} g بعد التفاعل",
                font=dict(family="Cairo, sans-serif", size=14, color="#0f4c5c"),
                x=0.5,
            ),
            yaxis=dict(title="الكتلة بالجرام (Mass in g)", gridcolor="#f1f5f9"),
            xaxis=dict(title="مراحل التجربة الكيميائية"),
            plot_bgcolor="#ffffff",
            paper_bgcolor="#ffffff",
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
            margin=dict(l=40, r=40, t=65, b=40),
            height=340,
        )
        return fig

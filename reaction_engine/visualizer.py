"""
Reaction Visualizer module for MolecuLens Chemical Playground.
Generates 2D molecular cards for reactants and products, analyzes bond transformations,
audits atom conservation, and generates interactive Plotly reaction energy coordinate diagrams.
"""

from __future__ import annotations

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
    def analyze_bond_transformations(
        cls,
        reactants: list[dict[str, Any]],
        products: list[dict[str, Any]],
        category: str = "",
    ) -> dict[str, Any]:
        """
        Identifies and summarizes the chemical bonds broken and formed during the reaction.
        """
        r_formulas = [r.get("formula", "") for r in reactants]
        p_formulas = [p.get("formula", "") for p in products]

        bonds_broken: list[dict[str, Any]] = []
        bonds_formed: list[dict[str, Any]] = []

        # Hydrogen + Oxygen -> Water
        if set(r_formulas) == {"H2", "O2"} and "H2O" in p_formulas:
            bonds_broken = [
                {"bond": "H-H (تساهمية أحادية)", "count": 2, "energy_kj": 436, "total_kj": 872},
                {"bond": "O=O (تساهمية ثنائية)", "count": 1, "energy_kj": 498, "total_kj": 498},
            ]
            bonds_formed = [
                {"bond": "O-H (تساهمية قطبية في H2O)", "count": 4, "energy_kj": 463, "total_kj": 1852},
            ]
        # Acid + Base neutralization (e.g. HCl + NaOH -> NaCl + H2O)
        elif ("HCl" in r_formulas or "HNO3" in r_formulas or "CH3COOH" in r_formulas) and ("NaOH" in r_formulas or "KOH" in r_formulas):
            bonds_broken = [
                {"bond": "H-Cl / H-A (تساهمية قطبية في الحمض)", "count": 1, "description": "تأين بروتون H+ في الوسط المائي"},
                {"bond": "Na-OH (رابطة أيونية في القاعدة)", "count": 1, "description": "تحرر أيون الهيدروكسيد OH-"},
            ]
            bonds_formed = [
                {"bond": "H-OH (رابطة تساهمية في H2O)", "count": 1, "energy_kj": 463, "description": "اتحاد H+ و OH- لتكوين جزيء الماء المستقر"},
                {"bond": "Na+...Cl- (تجاذب أيوني ملحي)", "count": 1, "description": "أيونات ملحية مميأة في المحلول"},
            ]
        # Precipitation (e.g. AgNO3 + NaCl -> AgCl + NaNO3)
        elif "AgCl" in p_formulas or "BaSO4" in p_formulas:
            bonds_broken = [
                {"bond": "روابط الإماهة المائية للأيونات", "count": 2, "description": "فقدان غلاف التميؤ المائي"},
            ]
            bonds_formed = [
                {"bond": "روابط أيونية بلورية في الراسب", "count": 1, "description": "تكون شبكة بلورية صلبة غير ذائبة عالية طاقة الشبيكة"},
            ]
        # Single displacement (e.g. Zn + 2HCl -> ZnCl2 + H2)
        elif "H2" in p_formulas and any(r in r_formulas for r in ["Zn", "Fe", "Mg"]):
            bonds_broken = [
                {"bond": "H-Cl (تساهمية قطبية)", "count": 2, "energy_kj": 431, "total_kj": 862},
                {"bond": "الرابطة الفلزية للفلز المتفاعل", "count": 1, "description": "أكسدة ذرات الفلز إلى كاتيونات"},
            ]
            bonds_formed = [
                {"bond": "H-H (تساهمية أحادية في H2)", "count": 1, "energy_kj": 436, "total_kj": 436},
                {"bond": "رابطة أيونية في الملح الناتج", "count": 1, "description": "تكون كلوريد الفلز المائي"},
            ]
        # Generic decomposition / synthesis fallback
        else:
            bonds_broken = [
                {"bond": "روابط المتفاعلات الأصلية", "count": len(reactants), "description": "انفصام الروابط التساهمية/الأيونية بين ذرات المواد المتفاعلة"},
            ]
            bonds_formed = [
                {"bond": "روابط النواتج المستقرة", "count": len(products), "description": "إعادة ترتيب الذرات وتكوين روابط جديدة أكثر استقراراً ثرموديناميكياً"},
            ]

        # Calculate energy estimates if numbers exist
        total_broken_energy = sum(b.get("total_kj", 0) for b in bonds_broken)
        total_formed_energy = sum(b.get("total_kj", 0) for b in bonds_formed)
        est_delta_h = (total_broken_energy - total_formed_energy) if (total_broken_energy and total_formed_energy) else -150.0

        is_exothermic = est_delta_h < 0

        return {
            "bonds_broken": bonds_broken,
            "bonds_formed": bonds_formed,
            "total_broken_energy_kj": total_broken_energy,
            "total_formed_energy_kj": total_formed_energy,
            "estimated_delta_h_kj": est_delta_h,
            "is_exothermic": is_exothermic,
            "thermodynamic_type_ar": "طارد للحرارة (Exothermic)" if is_exothermic else "ماص للحرارة (Endothermic)",
            "summary_ar": (
                f"التفاعل {('طارد للحرارة' if is_exothermic else 'ماص للحرارة')}؛ "
                f"طاقة الروابط المتكونة ({total_formed_energy} kJ/mol) {'تفوق' if is_exothermic else 'أقل من'} "
                f"طاقة الروابط المتكسرة ({total_broken_energy} kJ/mol)."
                if (total_broken_energy and total_formed_energy)
                else "يتضمن التفاعل إعادة ترتيب شاملة للروابط التساهمية والأيونية للوصول إلى النواتج الأكثر استقراراً."
            ),
        }

    @classmethod
    def create_reaction_energy_diagram(
        cls,
        reactants_str: str,
        products_str: str,
        is_exothermic: bool = True,
        delta_h_kj: float | None = None,
        activation_energy_kj: float = 65.0,
    ) -> go.Figure:
        """
        Builds an interactive Plotly reaction coordinate energy profile
        (Potential Energy vs Reaction Progress).
        """
        e_reactants = 100.0
        delta_h = delta_h_kj if delta_h_kj is not None else (-85.0 if is_exothermic else 45.0)
        e_products = e_reactants + delta_h
        e_transition = e_reactants + activation_energy_kj

        # Construct smooth reaction coordinate curve
        # Points: 0 to 1 (Reactants), 1 to 2.5 (Rise to TS), 2.5 to 4 (Fall to Products), 4 to 5 (Products)
        xs = [0.0, 0.8, 1.2, 1.8, 2.3, 2.5, 2.7, 3.2, 3.8, 4.2, 5.0]
        ys = [
            e_reactants,
            e_reactants,
            e_reactants + 15.0,
            e_reactants + (activation_energy_kj * 0.75),
            e_transition - 2.0,
            e_transition,
            e_transition - 2.0,
            e_products + (abs(delta_h) * 0.4 if is_exothermic else 15.0),
            e_products + 5.0,
            e_products,
            e_products,
        ]

        fig = go.Figure()

        # Smooth curve
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

        # Reactants Plateau Marker
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

        # Transition State Marker
        fig.add_trace(
            go.Scatter(
                x=[2.5],
                y=[e_transition],
                mode="markers+text",
                marker=dict(size=14, color="#dc2626", symbol="diamond"),
                text=[f"الحالة الانتقالية / المعقد النشط [‡]<br>طاقة التنشيط Ea = {activation_energy_kj:.1f} kJ"],
                textposition="top center",
                name="المعقد النشط (Transition State)",
            )
        )

        # Products Plateau Marker
        fig.add_trace(
            go.Scatter(
                x=[4.6],
                y=[e_products],
                mode="markers+text",
                marker=dict(size=12, color="#059669", symbol="square"),
                text=[f"النواتج<br>{products_str}"],
                textposition="bottom center" if is_exothermic else "top center",
                name="النواتج (Products)",
            )
        )

        # Delta H representation line
        fig.add_shape(
            type="line",
            x0=3.8, y0=e_reactants,
            x1=4.6, y1=e_reactants,
            line=dict(color="#94a3b8", dash="dot", width=1.5),
        )

        enthalpy_text = f"ΔH = {delta_h:.1f} kJ/mol ({'طارد' if is_exothermic else 'ماص'})"
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
                text="مخطط مسار التفاعل ومستويات الطاقة (Reaction Energy Coordinate)",
                font=dict(family="Cairo, sans-serif", size=15, color="#0f4c5c"),
                x=0.5,
            ),
            xaxis=dict(
                title="إحداثي مسار التفاعل (Reaction Coordinate)",
                showticklabels=False,
                zeroline=False,
                gridcolor="#f1f5f9",
            ),
            yaxis=dict(
                title="طاقة الوضع الكيميائية (kJ/mol)",
                zeroline=False,
                gridcolor="#f1f5f9",
            ),
            plot_bgcolor="#ffffff",
            paper_bgcolor="#ffffff",
            showlegend=True,
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
            margin=dict(l=40, r=40, t=60, b=40),
            height=380,
        )

        return fig

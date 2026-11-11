import math
from typing import Any
import plotly.graph_objects as go

KW_25C = 1.0e-14

ACID_SPECIES_DATABASE: dict[str, dict[str, Any]] = {
    "HCl": {
        "formula": "HCl",
        "name_ar": "حمض الهيدروكلوريك",
        "name_en": "Hydrochloric Acid",
        "type": "strong",
        "pKa": -6.3,
        "Ka": 1e6,
        "protons": 1,
    },
    "HNO3": {
        "formula": "HNO3",
        "name_ar": "حمض النيتريك",
        "name_en": "Nitric Acid",
        "type": "strong",
        "pKa": -1.4,
        "Ka": 25.0,
        "protons": 1,
    },
    "CH3COOH": {
        "formula": "CH3COOH",
        "name_ar": "حمض الخليك (حمض الأسيتيك)",
        "name_en": "Acetic Acid",
        "type": "weak",
        "pKa": 4.756,
        "Ka": 1.75e-5,
        "protons": 1,
    },
    "HCOOH": {
        "formula": "HCOOH",
        "name_ar": "حمض الفورميك (حمض النمليك)",
        "name_en": "Formic Acid",
        "type": "weak",
        "pKa": 3.75,
        "Ka": 1.78e-4,
        "protons": 1,
    },
    "HF": {
        "formula": "HF",
        "name_ar": "حمض الهيدروفلوريك",
        "name_en": "Hydrofluoric Acid",
        "type": "weak",
        "pKa": 3.17,
        "Ka": 6.76e-4,
        "protons": 1,
    },
    "C3H6O3": {
        "formula": "C3H6O3",
        "name_ar": "حمض اللاكتيك (حمض اللبن)",
        "name_en": "Lactic Acid",
        "type": "weak",
        "pKa": 3.86,
        "Ka": 1.38e-4,
        "protons": 1,
    },
}

BASE_SPECIES_DATABASE: dict[str, dict[str, Any]] = {
    "NaOH": {
        "formula": "NaOH",
        "name_ar": "هيدروكسيد الصوديوم",
        "name_en": "Sodium Hydroxide",
        "type": "strong",
        "pKb": 0.0,
        "Kb": 1e6,
        "hydroxides": 1,
    },
    "KOH": {
        "formula": "KOH",
        "name_ar": "هيدروكسيد البوتاسيوم",
        "name_en": "Potassium Hydroxide",
        "type": "strong",
        "pKb": 0.0,
        "Kb": 1e6,
        "hydroxides": 1,
    },
    "NH3": {
        "formula": "NH3",
        "name_ar": "الأمونيا (محلول النشادر)",
        "name_en": "Ammonia",
        "type": "weak",
        "pKb": 4.75,
        "Kb": 1.77e-5,
        "hydroxides": 1,
    },
}

INDICATORS_DATABASE: dict[str, dict[str, Any]] = {
    "phenolphthalein": {
        "name_ar": "فينولفثالين (Phenolphthalein)",
        "name_en": "Phenolphthalein",
        "pH_range": (8.2, 10.0),
        "color_acid": "عديم اللون (Colorless)",
        "color_base": "وردي قرمزي (Vivid Pink)",
        "color_hex": "#EC4899",
        "suitability_ar": "مثالي لمعايرة الأحماض الضعيفة مع القواعد القوية (نقطة التكافؤ pH > 7).",
    },
    "methyl_orange": {
        "name_ar": "ميثيل برتقالي (Methyl Orange)",
        "name_en": "Methyl Orange",
        "pH_range": (3.1, 4.4),
        "color_acid": "أحمر (Red)",
        "color_base": "أصفر برتقالي (Yellow)",
        "color_hex": "#F97316",
        "suitability_ar": "مثالي لمعايرة القواعد الضعيفة مع الأحماض القوية (نقطة التكافؤ pH < 7).",
    },
    "bromothymol_blue": {
        "name_ar": "أزرق بروموثيمول (Bromothymol Blue)",
        "name_en": "Bromothymol Blue",
        "pH_range": (6.0, 7.6),
        "color_acid": "أصفر (Yellow)",
        "color_base": "أزرق (Blue)",
        "color_hex": "#3B82F6",
        "suitability_ar": "مثالي للمعايرات المتكافئة عند نقطة التعادل (الأحماض القوية والقواعد القوية pH ≈ 7).",
    },
    "phenol_red": {
        "name_ar": "أحمر الفينول (Phenol Red)",
        "name_en": "Phenol Red",
        "pH_range": (6.8, 8.4),
        "color_acid": "أصفر (Yellow)",
        "color_base": "أحمر/أرجواني (Red)",
        "color_hex": "#EF4444",
        "suitability_ar": "ملائم للمعايرات القريبة من التعادل والوسط الفسيولوجي.",
    },
}


def solve_quadratic_positive(a: float, b: float, c: float) -> float:
    """Solves ax^2 + bx + c = 0 for the physical positive root."""
    disc = b * b - 4 * a * c
    if disc < 0:
        return 0.0
    sqrt_d = math.sqrt(disc)
    r1 = (-b + sqrt_d) / (2 * a)
    r2 = (-b - sqrt_d) / (2 * a)
    if r1 > 0:
        return r1
    if r2 > 0:
        return r2
    return 0.0


def compute_ph_from_h(h_conc: float) -> float:
    h = max(h_conc, 1e-14)
    return max(0.0, min(14.0, -math.log10(h)))


def compute_ph_from_oh(oh_conc: float) -> float:
    oh = max(oh_conc, 1e-14)
    poh = -math.log10(oh)
    return max(0.0, min(14.0, 14.0 - poh))


class TitrationSimulator:
    """
    Rigorously computes acid-base titration thermodynamics and points.
    Supports Strong/Weak Acid titrated by Strong/Weak Base and vice-versa.
    """

    @classmethod
    def simulate_titration(
        cls,
        analyte_formula: str,
        titrant_formula: str,
        analyte_volume_ml: float = 25.0,
        analyte_conc_m: float = 0.1,
        titrant_conc_m: float = 0.1,
        max_titrant_ratio: float = 2.2,
        points_count: int = 150,
    ) -> dict[str, Any]:
        """
        Simulates an acid-base titration experiment.
        Default configuration: Acid as Analyte, Base as Titrant (or vice-versa).
        """
        is_acid_analyte = analyte_formula in ACID_SPECIES_DATABASE
        if is_acid_analyte:
            acid_spec = ACID_SPECIES_DATABASE[analyte_formula]
            base_spec = BASE_SPECIES_DATABASE.get(titrant_formula, BASE_SPECIES_DATABASE["NaOH"])
            mode = "acid_analyte_titrated_by_base"
            c_analyte = max(0.001, analyte_conc_m)
            v_analyte = max(1.0, analyte_volume_ml)
            c_titrant = max(0.001, titrant_conc_m)
        else:
            base_spec = BASE_SPECIES_DATABASE.get(analyte_formula, BASE_SPECIES_DATABASE["NaOH"])
            acid_spec = ACID_SPECIES_DATABASE.get(titrant_formula, ACID_SPECIES_DATABASE["HCl"])
            mode = "base_analyte_titrated_by_acid"
            c_analyte = max(0.001, analyte_conc_m)
            v_analyte = max(1.0, analyte_volume_ml)
            c_titrant = max(0.001, titrant_conc_m)

        n_analyte_mmol = c_analyte * v_analyte
        v_eq_ml = n_analyte_mmol / c_titrant
        v_max_ml = v_eq_ml * max_titrant_ratio

        # Determine Titration Case
        acid_type = acid_spec["type"]
        base_type = base_spec["type"]
        if acid_type == "strong" and base_type == "strong":
            case_key = "strong_acid_strong_base"
            case_label_ar = "معايرة حمض قوي مع قاعدة قوية (Strong Acid - Strong Base)"
        elif acid_type == "weak" and base_type == "strong":
            case_key = "weak_acid_strong_base"
            case_label_ar = "معايرة حمض ضعيف مع قاعدة قوية (Weak Acid - Strong Base)"
        elif acid_type == "strong" and base_type == "weak":
            case_key = "strong_acid_weak_base"
            case_label_ar = "معايرة حمض قوي مع قاعدة ضعيفة (Strong Acid - Weak Base)"
        else:
            case_key = "weak_acid_weak_base"
            case_label_ar = "معايرة حمض ضعيف مع قاعدة ضعيفة (Weak Acid - Weak Base)"

        # Calculate Equivalence Point pH
        total_vol_at_eq = v_analyte + v_eq_ml
        c_salt_at_eq = n_analyte_mmol / total_vol_at_eq

        if case_key == "strong_acid_strong_base":
            ph_eq = 7.00
            ph_half_eq = None
        elif case_key == "weak_acid_strong_base":
            kb = KW_25C / acid_spec["Ka"]
            oh_eq = solve_quadratic_positive(1.0, kb, -kb * c_salt_at_eq)
            ph_eq = round(compute_ph_from_oh(oh_eq), 2)
            ph_half_eq = round(acid_spec["pKa"], 2)
        elif case_key == "strong_acid_weak_base":
            ka_conj = KW_25C / base_spec["Kb"]
            h_eq = solve_quadratic_positive(1.0, ka_conj, -ka_conj * c_salt_at_eq)
            ph_eq = round(compute_ph_from_h(h_eq), 2)
            ph_half_eq = round(14.0 - base_spec["pKb"], 2)
        else:
            ph_eq = round(7.0 + 0.5 * (acid_spec["pKa"] - base_spec["pKb"]), 2)
            ph_half_eq = round(acid_spec["pKa"], 2)

        # Generate dense volume points with refinement around V_eq
        v_points: list[float] = []
        base_step = v_max_ml / (points_count - 20)
        curr_v = 0.0
        while curr_v <= v_max_ml:
            v_points.append(round(curr_v, 4))
            curr_v += base_step

        # Add refined points near V_eq
        near_eq = [
            v_eq_ml * 0.5,
            v_eq_ml * 0.90,
            v_eq_ml * 0.95,
            v_eq_ml * 0.98,
            v_eq_ml * 0.99,
            v_eq_ml * 0.995,
            v_eq_ml * 1.0,
            v_eq_ml * 1.005,
            v_eq_ml * 1.01,
            v_eq_ml * 1.02,
            v_eq_ml * 1.05,
            v_eq_ml * 1.10,
        ]
        for ep in near_eq:
            if 0.0 <= ep <= v_max_ml:
                v_points.append(round(ep, 4))
        v_points = sorted(set(v_points))

        # Compute pH for each volume point
        curve_data: list[dict[str, Any]] = []
        for v_tit in v_points:
            v_tot = v_analyte + v_tit
            n_tit_mmol = c_titrant * v_tit

            if mode == "acid_analyte_titrated_by_base":
                # Acid analyte titrated with base titrant
                if v_tit == 0.0:
                    if acid_type == "strong":
                        ph_val = compute_ph_from_h(c_analyte)
                    else:
                        ka = acid_spec["Ka"]
                        h_val = solve_quadratic_positive(1.0, ka, -ka * c_analyte)
                        ph_val = compute_ph_from_h(h_val)
                    stage = "نقطة البداية (Initial Point)"

                elif v_tit < v_eq_ml - 1e-5:
                    # Before Equivalence Point
                    if acid_type == "strong":
                        rem_h_mmol = n_analyte_mmol - n_tit_mmol
                        h_conc = rem_h_mmol / v_tot
                        ph_val = compute_ph_from_h(h_conc)
                    else:
                        # Weak acid buffer: [HA] remaining, [A-] formed
                        ha_mmol = n_analyte_mmol - n_tit_mmol
                        a_mmol = n_tit_mmol
                        pka = acid_spec["pKa"]
                        ratio = max(1e-6, a_mmol / max(1e-6, ha_mmol))
                        ph_val = pka + math.log10(ratio)
                    stage = "منطقة المحلول المنظم / قبل التكافؤ (Buffer Region)"

                elif abs(v_tit - v_eq_ml) <= 1e-4:
                    ph_val = ph_eq
                    stage = "نقطة التكافؤ (Equivalence Point)"

                else:
                    # Past Equivalence Point
                    excess_tit_mmol = n_tit_mmol - n_analyte_mmol
                    if base_type == "strong":
                        oh_conc = excess_tit_mmol / v_tot
                        ph_val = compute_ph_from_oh(oh_conc)
                    else:
                        kb = base_spec["Kb"]
                        c_excess = excess_tit_mmol / v_tot
                        c_salt = n_analyte_mmol / v_tot
                        ratio = max(1e-6, c_excess / max(1e-6, c_salt))
                        ph_val = 14.0 - (base_spec["pKb"] - math.log10(ratio))
                    stage = "فائض المحلول القياسي / بعد التكافؤ (Post-Equivalence)"

            else:
                # Base analyte titrated with acid titrant
                if v_tit == 0.0:
                    if base_type == "strong":
                        ph_val = compute_ph_from_oh(c_analyte)
                    else:
                        kb = base_spec["Kb"]
                        oh_val = solve_quadratic_positive(1.0, kb, -kb * c_analyte)
                        ph_val = compute_ph_from_oh(oh_val)
                    stage = "نقطة البداية (Initial Point)"

                elif v_tit < v_eq_ml - 1e-5:
                    if base_type == "strong":
                        rem_oh_mmol = n_analyte_mmol - n_tit_mmol
                        oh_conc = rem_oh_mmol / v_tot
                        ph_val = compute_ph_from_oh(oh_conc)
                    else:
                        b_mmol = n_analyte_mmol - n_tit_mmol
                        bh_mmol = n_tit_mmol
                        pkb = base_spec["pKb"]
                        ratio = max(1e-6, b_mmol / max(1e-6, bh_mmol))
                        ph_val = 14.0 - (pkb - math.log10(ratio))
                    stage = "منطقة المحلول المنظم / قبل التكافؤ (Buffer Region)"

                elif abs(v_tit - v_eq_ml) <= 1e-4:
                    ph_val = ph_eq
                    stage = "نقطة التكافؤ (Equivalence Point)"

                else:
                    excess_tit_mmol = n_tit_mmol - n_analyte_mmol
                    if acid_type == "strong":
                        h_conc = excess_tit_mmol / v_tot
                        ph_val = compute_ph_from_h(h_conc)
                    else:
                        c_excess = excess_tit_mmol / v_tot
                        c_salt = n_analyte_mmol / v_tot
                        ratio = max(1e-6, c_excess / max(1e-6, c_salt))
                        ph_val = acid_spec["pKa"] - math.log10(ratio)
                    stage = "فائض المحلول القياسي / بعد التكافؤ (Post-Equivalence)"

            curve_data.append({
                "volume_added_ml": round(v_tit, 2),
                "ph": round(max(0.0, min(14.0, ph_val)), 2),
                "stage_ar": stage,
            })

        # Match Suitable Indicator
        best_indicator_key = "bromothymol_blue"
        if ph_eq > 7.5:
            best_indicator_key = "phenolphthalein"
        elif ph_eq < 6.5:
            best_indicator_key = "methyl_orange"
        indicator_info = INDICATORS_DATABASE[best_indicator_key]

        explanation_ar = (
            f"عند معايرة **{analyte_formula}** ({c_analyte} M، حجم {v_analyte} mL) مع **{titrant_formula}** ({c_titrant} M)، "
            f"يتحقق التكافؤ الكيميائي التام عند إضافة **{v_eq_ml:.2f} mL** من المحلول القياسي. "
            f"قيمة الرقم الهيدروجيني عند نقطة التكافؤ هي **pH = {ph_eq:.2f}** ({'وسط متعادل تماماً' if abs(ph_eq - 7.0) < 0.2 else ('وسط قاعدي بسبب تميه ملح الحمض الضعيف' if ph_eq > 7.0 else 'وسط حمضي بسبب تميه ملح القاعدة الضعيفة')}). "
            f"الكاشف اللوني الأنسب علمياً هو **{indicator_info['name_ar']}** لأن مدى تغير لونه ({indicator_info['pH_range'][0]} - {indicator_info['pH_range'][1]}) يتقاطع تماماً مع القفزة الحادة لمنحنى المعايرة."
        )

        return {
            "mode": mode,
            "case_key": case_key,
            "case_label_ar": case_label_ar,
            "analyte": {
                "formula": analyte_formula,
                "name_ar": acid_spec["name_ar"] if is_acid_analyte else base_spec["name_ar"],
                "concentration_m": c_analyte,
                "volume_ml": v_analyte,
                "amount_mmol": round(n_analyte_mmol, 4),
            },
            "titrant": {
                "formula": titrant_formula,
                "name_ar": base_spec["name_ar"] if is_acid_analyte else acid_spec["name_ar"],
                "concentration_m": c_titrant,
            },
            "equivalence_point": {
                "volume_ml": round(v_eq_ml, 2),
                "ph": ph_eq,
            },
            "half_equivalence_point": {
                "volume_ml": round(v_eq_ml * 0.5, 2),
                "ph": ph_half_eq,
            } if ph_half_eq is not None else None,
            "recommended_indicator": indicator_info,
            "curve_points": curve_data,
            "explanation_ar": explanation_ar,
        }

    @classmethod
    def build_plotly_titration_chart(cls, sim_result: dict[str, Any]) -> go.Figure:
        curve = sim_result["curve_points"]
        v_list = [p["volume_added_ml"] for p in curve]
        ph_list = [p["ph"] for p in curve]
        v_eq = sim_result["equivalence_point"]["volume_ml"]
        ph_eq = sim_result["equivalence_point"]["ph"]
        ind = sim_result["recommended_indicator"]
        ind_low, ind_high = ind["pH_range"]

        fig = go.Figure()

        # Indicator Shaded Band
        fig.add_hrect(
            y0=ind_low, y1=ind_high,
            fillcolor=ind["color_hex"], opacity=0.15,
            line_width=1, line_dash="dash", line_color=ind["color_hex"],
            annotation_text=f"نطاق كاشف {ind['name_ar']} (pH {ind_low} - {ind_high})",
            annotation_position="bottom right",
        )

        # Neutral line at pH 7
        fig.add_hline(y=7.0, line_dash="dot", line_color="#94a3b8", line_width=1.2)

        # Main pH Titration Curve
        fig.add_trace(
            go.Scatter(
                x=v_list,
                y=ph_list,
                mode="lines",
                name="منحنى pH (Titration Curve)",
                line=dict(color="#0f4c5c", width=3.5),
                hoverinfo="text",
                hovertext=[f"الحجم: {v} mL<br>pH: {ph}<br>{stg}" for v, ph, stg in zip(v_list, ph_list, [p["stage_ar"] for p in curve])],
            )
        )

        # Equivalence Point Marker
        fig.add_trace(
            go.Scatter(
                x=[v_eq],
                y=[ph_eq],
                mode="markers+text",
                name="نقطة التكافؤ (Equivalence Point)",
                text=[f"نقطة التكافؤ<br>V = {v_eq} mL, pH = {ph_eq}"],
                textposition="top right",
                marker=dict(size=14, color="#DC2626", symbol="diamond", line=dict(color="#ffffff", width=2)),
            )
        )

        # Half-Equivalence Point (for weak species)
        half_eq = sim_result.get("half_equivalence_point")
        if half_eq and half_eq.get("ph") is not None:
            v_half = half_eq["volume_ml"]
            ph_half = half_eq["ph"]
            fig.add_trace(
                go.Scatter(
                    x=[v_half],
                    y=[ph_half],
                    mode="markers+text",
                    name="نصف التكافؤ (pH = pKa)",
                    text=[f"نصف التكافؤ<br>pH = {ph_half}"],
                    textposition="bottom right",
                    marker=dict(size=11, color="#059669", symbol="circle", line=dict(color="#ffffff", width=2)),
                )
            )

        titrant_label = sim_result["titrant"]["name_ar"]
        titrant_form = sim_result["titrant"]["formula"]

        fig.update_layout(
            title=dict(
                text=f"منحنى المعايرة الكيميائي: {sim_result['case_label_ar']}",
                font=dict(family="Cairo, sans-serif", size=16, color="#0f4c5c"),
                x=0.5,
            ),
            xaxis=dict(
                title=f"حجم المحلول القياسي المضاف ({titrant_label} — {titrant_form}) بوحدة mL",
                gridcolor="#f1f5f9",
                zeroline=False,
            ),
            yaxis=dict(
                title="الرقم الهيدروجيني (pH)",
                range=[0, 14],
                gridcolor="#f1f5f9",
                zeroline=False,
            ),
            plot_bgcolor="#ffffff",
            paper_bgcolor="#ffffff",
            hovermode="closest",
            legend=dict(orientation="h", yanchor="bottom", y=1.02, xanchor="right", x=1),
            margin=dict(l=40, r=40, t=60, b=40),
            height=480,
        )

        return fig


class TitrationEngine:
    """Unified engine for Titration views and test suites."""

    ACID_DATABASE = ACID_SPECIES_DATABASE
    BASE_DATABASE = BASE_SPECIES_DATABASE
    INDICATORS = {
        k: {
            **v,
            "transition_range": v.get("pH_range", (7.0, 7.0)),
        }
        for k, v in INDICATORS_DATABASE.items()
    }

    @classmethod
    def simulate_titration(
        cls,
        acid_key: str,
        base_key: str,
        acid_conc_m: float,
        acid_volume_ml: float,
        base_conc_m: float,
    ) -> dict[str, Any]:
        if acid_key not in cls.ACID_DATABASE:
            raise ValueError(f"Unknown acid key: {acid_key}")
        if base_key not in cls.BASE_DATABASE:
            raise ValueError(f"Unknown base key: {base_key}")
        if acid_conc_m <= 0 or acid_volume_ml <= 0 or base_conc_m <= 0:
            raise ValueError("Concentrations and volumes must be strictly positive.")

        sim = TitrationSimulator.simulate_titration(
            analyte_formula=acid_key,
            titrant_formula=base_key,
            analyte_volume_ml=acid_volume_ml,
            analyte_conc_m=acid_conc_m,
            titrant_conc_m=base_conc_m,
        )

        eq = sim["equivalence_point"]
        half_eq = sim.get("half_equivalence_point")
        curve = sim["curve_points"]
        v_list = [p["volume_added_ml"] for p in curve]
        ph_list = [p["ph"] for p in curve]

        return {
            "success": True,
            "v_eq_ml": eq["volume_ml"],
            "ph_eq": eq["ph"],
            "initial_ph": curve[0]["ph"] if curve else 7.0,
            "half_v_eq_ml": half_eq["volume_ml"] if half_eq else None,
            "half_ph": half_eq["ph"] if half_eq else None,
            "volumes_ml": v_list,
            "ph_values": ph_list,
            "sim_result": sim,
            "recommended_indicator": sim.get("recommended_indicator"),
        }

    @classmethod
    def calculate_point_ph(
        cls,
        v_add: float,
        v_eq: float,
        v_acid: float,
        ca: float,
        cb: float,
        acid_info: dict[str, Any],
        base_info: dict[str, Any],
    ) -> dict[str, Any]:
        sim = TitrationSimulator.simulate_titration(
            analyte_formula=acid_info["formula"],
            titrant_formula=base_info["formula"],
            analyte_volume_ml=v_acid,
            analyte_conc_m=ca,
            titrant_conc_m=cb,
            max_titrant_ratio=max(2.2, (v_add / v_eq) + 0.2 if v_eq > 0 else 2.2),
        )
        curve = sim["curve_points"]
        closest = min(curve, key=lambda p: abs(p["volume_added_ml"] - v_add))
        return {
            "ph": closest["ph"],
            "region": closest["stage_ar"],
        }

    @classmethod
    def create_plotly_curve(
        cls,
        sim_data: dict[str, Any],
        current_v: float = None,
        indicator_key: str = None,
    ) -> go.Figure:
        sim_result = sim_data.get("sim_result")
        if not sim_result:
            return go.Figure()

        fig = TitrationSimulator.build_plotly_titration_chart(sim_result)

        if indicator_key and indicator_key in cls.INDICATORS:
            ind = cls.INDICATORS[indicator_key]
            ind_low, ind_high = ind["transition_range"]
            fig.add_hrect(
                y0=ind_low,
                y1=ind_high,
                fillcolor=ind.get("color_hex", "#3b82f6"),
                opacity=0.18,
                line_width=1.5,
                line_dash="dot",
                annotation_text=f"كاشف {ind['name_ar']} (pH {ind_low}-{ind_high})",
                annotation_position="top left",
            )

        if current_v is not None:
            v_list = sim_data["volumes_ml"]
            ph_list = sim_data["ph_values"]
            closest_idx = min(range(len(v_list)), key=lambda i: abs(v_list[i] - current_v))
            cur_ph = ph_list[closest_idx]
            fig.add_trace(
                go.Scatter(
                    x=[current_v],
                    y=[cur_ph],
                    mode="markers+text",
                    name="الموضع الحالي (Current)",
                    text=[f"الحالي: {current_v:.1f} mL, pH {cur_ph:.2f}"],
                    textposition="bottom right",
                    marker=dict(size=13, color="#f59e0b", symbol="circle", line=dict(color="#ffffff", width=2)),
                )
            )

        return fig

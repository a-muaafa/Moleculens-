"""
Scientific Experiment Report Generator (PDF & Printable HTML) for MolecuLens.

Generates structured, print-ready laboratory experiment reports for Reaction Lab
and Virtual Titration Lab including:
1. Materials, formulas, molar masses, and input quantities/units.
2. Experimental conditions (temperature, pressure, solvent/medium, catalyst) and compatibility.
3. Balanced chemical equation, complete ionic equation, net ionic equation, and spectator ions.
4. Stoichiometric calculations (limiting reagent, excess reagents, consumed/remaining amounts, theoretical yield).
5. Scientific explanation, laboratory observations, equilibrium status, and competing pathways.
6. Evidence tier, verification badge, and reliability classification.
7. Prediction limits, model assumptions, and safety/hazard warnings.
8. Bibliographic sources, reference IDs, and URLs.
"""

from __future__ import annotations

from datetime import datetime, timezone
import html
from io import BytesIO
from typing import Any

from chemistry_engine.formula_analyzer import format_chemical_equation_display, format_chemical_formula_display


def _esc(value: Any) -> str:
    if value is None:
        return ""
    return html.escape(str(value), quote=True)


class ExperimentReportGenerator:
    """Builds structured report payloads, printable HTML reports, and valid PDF binary files."""

    @classmethod
    def build_reaction_report_data(
        cls,
        reaction_payload: dict[str, Any],
        custom_conditions: dict[str, Any] | None = None,
    ) -> dict[str, Any]:
        ts_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
        supported = bool(reaction_payload.get("supported"))
        mode = reaction_payload.get("mode") or (
            "documented" if reaction_payload.get("documented_reaction_info") else "insufficient_evidence"
        )
        mode_label_ar = reaction_payload.get("mode_label_ar") or (
            "تفاعل موثّق في سجل التفاعلات المعتمد (Documented Reaction)"
            if mode == "documented" and supported
            else "لا تتوفر بيانات كافية / لا يحدث تفاعل"
        )

        reactants = reaction_payload.get("reactants") or []
        products = reaction_payload.get("products") or []
        reactants_table = reaction_payload.get("reactants_table") or []
        yields = reaction_payload.get("theoretical_yields") or []

        if not reactants and reactants_table:
            reactants = [
                {
                    "formula": r.get("formula", ""),
                    "name_ar": r.get("name_ar", ""),
                    "name_en": r.get("name_en", r.get("formula", "")),
                    "coefficient": r.get("coefficient", 1),
                    "state": r.get("state", ""),
                }
                for r in reactants_table
            ]
        if not products and yields:
            products = [
                {
                    "formula": y.get("formula", ""),
                    "name_ar": y.get("name_ar", ""),
                    "name_en": y.get("name_en", y.get("formula", "")),
                    "coefficient": y.get("coefficient", 1),
                    "state": y.get("state", ""),
                }
                for y in yields
            ]

        eq_text = reaction_payload.get("equation_text", "")
        eq_unicode = reaction_payload.get("equation_unicode") or (
            format_chemical_equation_display(eq_text) if eq_text else "لا توجد معادلة تفاعل (No Reaction Equation)"
        )

        conds = dict(reaction_payload.get("conditions") or {})
        if custom_conditions:
            conds.update({k: v for k, v in custom_conditions.items() if v is not None})

        ionic = reaction_payload.get("ionic_analysis") or {}
        equil = reaction_payload.get("equilibrium_info") or {}
        competing = reaction_payload.get("competing_pathways") or {}
        safety = reaction_payload.get("safety_warnings") or {}
        confidence = reaction_payload.get("confidence") or {}

        sources_list: list[dict[str, str]] = []
        for ev in reaction_payload.get("evidence") or []:
            sources_list.append({
                "source_id": str(ev.get("source_id", "")),
                "source": str(ev.get("source", "")),
                "reference_url": str(ev.get("reference_url", "")),
                "rule_or_name_ar": str(ev.get("rule_name_ar", "")),
            })
        doc_info = reaction_payload.get("documented_reaction_info") or {}
        if not sources_list and doc_info.get("source"):
            sources_list.append({
                "source_id": str(doc_info.get("source_id", "")),
                "source": str(doc_info.get("source", "")),
                "reference_url": str(doc_info.get("reference_url", "")),
                "rule_or_name_ar": str(doc_info.get("reaction_name_ar", "")),
            })

        limits_list: list[str] = list(reaction_payload.get("limitations_ar") or [])
        if reaction_payload.get("assumptions_and_limits_ar"):
            limits_list.append(str(reaction_payload["assumptions_and_limits_ar"]))
        if equil.get("yield_equilibrium_note_ar"):
            limits_list.append(str(equil["yield_equilibrium_note_ar"]))
        if not limits_list:
            limits_list.append(
                "الحسابات المعروضة تمثل النمذجة الستوكيومترية النظرية بفرض نقاء المواد؛ ولا يُعرض مردود فعلي بدون قياس تجريبي."
            )

        return {
            "report_type": "reaction_experiment",
            "title_ar": "تقرير تجربة مختبر التفاعلات الكيميائية — MolecuLens Scientific Report",
            "title_en": "MolecuLens Virtual Reaction Lab Experiment Report",
            "generated_at_utc": ts_utc,
            "supported": supported,
            "mode": mode,
            "mode_label_ar": mode_label_ar,
            "confidence_label_ar": confidence.get("label_ar") or doc_info.get("documentation_label_ar", "محسوب برمجياً"),
            "equation_text": eq_text,
            "equation_unicode": eq_unicode,
            "complete_ionic_equation": ionic.get("complete_ionic_equation") or reaction_payload.get("complete_ionic", ""),
            "net_ionic_equation": ionic.get("net_ionic_equation") or reaction_payload.get("net_ionic", ""),
            "spectator_ions": ionic.get("spectator_ions") or reaction_payload.get("spectator_ions", []),
            "reactants": reactants,
            "products": products,
            "reactants_table": reactants_table,
            "theoretical_yields": yields,
            "limiting_reagent": reaction_payload.get("limiting_reagent"),
            "unreacted_species": reaction_payload.get("unreacted_species") or [],
            "conditions": {
                "temperature": conds.get("temperature", "25 °C (الظروف القياسية)"),
                "pressure": conds.get("pressure", "1 atm"),
                "solvent": conds.get("solvent", "محلول مائي / طور قياسي"),
                "catalyst": conds.get("catalyst", "بدون حفاز إضافي"),
                "compatibility_status_ar": (
                    (reaction_payload.get("conditions_compatibility") or {}).get("status_ar")
                    or "ضمن الظروف المرجعية المعتمدة"
                ),
            },
            "atom_conservation": reaction_payload.get("atom_conservation") or {},
            "explanation_ar": reaction_payload.get("explanation_ar") or reaction_payload.get("verdict_ar", ""),
            "observations": reaction_payload.get("observations") or [],
            "equilibrium_info": equil,
            "competing_pathways": competing,
            "limitations_ar": limits_list,
            "sources": sources_list,
            "safety_warnings": safety,
        }

    @classmethod
    def build_titration_report_data(cls, sim_data: dict[str, Any]) -> dict[str, Any]:
        ts_utc = datetime.now(timezone.utc).strftime("%Y-%m-%d %H:%M:%S UTC")
        sim = sim_data.get("sim_result") or sim_data
        analyte = sim.get("analyte") or {}
        titrant = sim.get("titrant") or {}
        eq = sim.get("equivalence_point") or {"volume_ml": sim_data.get("v_eq_ml", 0), "ph": sim_data.get("ph_eq", 7)}
        endpt = sim.get("endpoint") or {"volume_ml": sim_data.get("endpoint_volume_ml", 0), "ph": sim_data.get("endpoint_ph", 7)}
        ind = sim.get("recommended_indicator") or {}

        return {
            "report_type": "titration_experiment",
            "title_ar": "تقرير تجربة مختبر المعايرة الافتراضي — MolecuLens Titration Report",
            "title_en": "MolecuLens Virtual Titration Lab Experiment Report",
            "generated_at_utc": ts_utc,
            "case_label_ar": sim.get("case_label_ar", "معايرة حمض وقاعدة"),
            "analyte": analyte,
            "titrant": titrant,
            "initial_ph": sim_data.get("initial_ph", 7.0),
            "equivalence_point": eq,
            "endpoint": endpt,
            "step_equivalence_points": sim.get("step_equivalence_points") or [],
            "buffer_analysis": sim.get("buffer_analysis") or {},
            "titration_error_analysis": sim.get("titration_error_analysis") or {},
            "recommended_indicator": ind,
            "explanation_ar": sim.get("explanation_ar", ""),
            "step_by_step_calculation": sim.get("step_by_step_calculation") or [],
            "model_limitations_ar": sim.get("model_limitations_ar") or [],
        }

    @classmethod
    def generate_printable_html_report(cls, report_data: dict[str, Any]) -> str:
        """Generates a standalone, print-ready HTML document with RTL layout and scientific tables."""
        rtype = report_data.get("report_type", "reaction_experiment")
        title_ar = _esc(report_data.get("title_ar", "تقرير تجربة MolecuLens"))
        ts = _esc(report_data.get("generated_at_utc", ""))

        if rtype == "titration_experiment":
            analyte = report_data.get("analyte") or {}
            titrant = report_data.get("titrant") or {}
            eq = report_data.get("equivalence_point") or {}
            endpt = report_data.get("endpoint") or {}
            ind = report_data.get("recommended_indicator") or {}
            err_info = report_data.get("titration_error_analysis") or {}
            lims = "".join(f"<li>{_esc(x)}</li>" for x in report_data.get("model_limitations_ar", []))
            steps = "".join(f"<li>{_esc(x)}</li>" for x in report_data.get("step_by_step_calculation", []))
            return f"""<!DOCTYPE html>
<html lang="ar" dir="rtl">
<head>
<meta charset="utf-8"/>
<title>{title_ar}</title>
<style>
  body {{ font-family: 'Cairo', 'Noto Sans Arabic', Arial, sans-serif; color: #0f172a; margin: 28px; line-height: 1.65; }}
  .header {{ border-bottom: 3px solid #0f4c5c; padding-bottom: 12px; margin-bottom: 18px; }}
  .box {{ border: 1px solid #cbd5e1; border-radius: 8px; padding: 14px; margin-bottom: 14px; background: #f8fafc; }}
  table {{ width: 100%; border-collapse: collapse; margin: 10px 0; }}
  th, td {{ border: 1px solid #cbd5e1; padding: 8px 10px; text-align: right; font-size: 0.92rem; }}
  th {{ background: #e2e8f0; color: #0f4c5c; }}
  @media print {{ body {{ margin: 12mm; }} }}
</style>
</head>
<body>
  <div class="header">
    <h1>{title_ar}</h1>
    <div>تاريخ التوليد: <span dir="ltr">{ts}</span> | التصنيف: {_esc(report_data.get('case_label_ar', ''))}</div>
  </div>
  <div class="box">
    <h3>1. المواد والتراكيز والحجوم (Materials &amp; Setup)</h3>
    <p><strong>المادة في الدورق (Analyte):</strong> {_esc(analyte.get('name_ar'))} ({_esc(analyte.get('formula'))}) — التركيز: {_esc(analyte.get('concentration_m'))} M — الحجم: {_esc(analyte.get('volume_ml'))} mL</p>
    <p><strong>المحلول القياسي في السحاحة (Titrant):</strong> {_esc(titrant.get('name_ar'))} ({_esc(titrant.get('formula'))}) — التركيز: {_esc(titrant.get('concentration_m'))} M</p>
  </div>
  <div class="box">
    <h3>2. نقطة التكافؤ مقابل نقطة النهاية اللونية (Equivalence Point vs. Indicator Endpoint)</h3>
    <p><strong>نقطة التكافؤ الستوكيومترية:</strong> V_eq = {_esc(eq.get('volume_ml'))} mL | pH_eq = {_esc(eq.get('ph'))}</p>
    <p><strong>نقطة نهاية الكاشف ({_esc(ind.get('name_ar', ''))}):</strong> V_end = {_esc(endpt.get('volume_ml'))} mL | pH_end = {_esc(endpt.get('ph'))}</p>
    <p><strong>فرق الحجم وخطأ المعايرة اللوني:</strong> {_esc(err_info.get('summary_ar', 'متطابق عملياً'))}</p>
  </div>
  <div class="box">
    <h3>3. خطوات الحساب الكيميائي والتفسير (Step-by-Step Calculations &amp; Explanation)</h3>
    <p>{_esc(report_data.get('explanation_ar', ''))}</p>
    <ul>{steps}</ul>
  </div>
  <div class="box">
    <h3>4. حدود نموذج المعايرة المبسط وإخلاء المسؤولية (Model Limitations &amp; Disclaimer)</h3>
    <ul>{lims}</ul>
    <p><em>تنبيه علمي: هذا التقرير ناتج عن محاكاة افتراضية تعليمية في منصة MolecuLens وليس بديلاً عن التجارب الواقعية أو دساتير الأدوية الرسمية.</em></p>
  </div>
</body>
</html>"""

        conds = report_data.get("conditions") or {}
        stoich_rows_html = ""
        for r in report_data.get("reactants_table") or []:
            stoich_rows_html += (
                f"<tr><td>{_esc(r.get('name_ar'))} ({_esc(r.get('formula'))})</td>"
                f"<td>{_esc(r.get('coefficient'))}</td>"
                f"<td>{_esc(r.get('molar_mass'))} g/mol</td>"
                f"<td>{_esc(r.get('role_ar'))}</td>"
                f"<td>{_esc(r.get('initial_mol'))} mol ({_esc(r.get('initial_mass_g'))} g)</td>"
                f"<td>{_esc(r.get('consumed_mol'))} mol ({_esc(r.get('consumed_mass_g'))} g)</td>"
                f"<td>{_esc(r.get('remaining_mol'))} mol ({_esc(r.get('remaining_mass_g'))} g)</td></tr>"
            )
        if not stoich_rows_html:
            for r in report_data.get("reactants") or []:
                stoich_rows_html += (
                    f"<tr><td>{_esc(r.get('name_ar'))} ({_esc(r.get('formula'))})</td>"
                    f"<td>{_esc(r.get('coefficient', 1))}</td><td>—</td><td>متفاعل</td><td>—</td><td>—</td><td>—</td></tr>"
                )

        yield_rows_html = ""
        for y in report_data.get("theoretical_yields") or []:
            yield_rows_html += (
                f"<tr><td>{_esc(y.get('name_ar'))} ({_esc(y.get('formula'))})</td>"
                f"<td>{_esc(y.get('coefficient'))}</td>"
                f"<td>{_esc(y.get('molar_mass'))} g/mol</td>"
                f"<td>{_esc(y.get('theoretical_yield_mol'))} mol</td>"
                f"<td>{_esc(y.get('theoretical_yield_g'))} g</td></tr>"
            )

        sources_html = ""
        for s in report_data.get("sources") or []:
            sources_html += (
                f"<li><strong>{_esc(s.get('source_id'))}:</strong> {_esc(s.get('source'))} "
                f"(<span dir='ltr'>{_esc(s.get('reference_url'))}</span>)</li>"
            )
        if not sources_html:
            sources_html = "<li>لم يُسجّل مرجع تفاعلي بسبب غياب الأدلة أو عدم حدوث تفاعل.</li>"

        limits_html = "".join(f"<li>{_esc(lim)}</li>" for lim in report_data.get("limitations_ar") or [])
        safety = report_data.get("safety_warnings") or {}
        safety_html = "".join(f"<li>{_esc(w)}</li>" for w in safety.get("warnings_ar") or [])
        equil = report_data.get("equilibrium_info") or {}

        return f"""<!DOCTYPE html>
<html lang="ar" dir="rtl">
<head>
<meta charset="utf-8"/>
<title>{title_ar}</title>
<style>
  body {{ font-family: 'Cairo', 'Noto Sans Arabic', Arial, sans-serif; color: #0f172a; margin: 28px; line-height: 1.65; }}
  .header {{ border-bottom: 3px solid #0f4c5c; padding-bottom: 12px; margin-bottom: 18px; }}
  .box {{ border: 1px solid #cbd5e1; border-radius: 8px; padding: 14px; margin-bottom: 14px; background: #f8fafc; }}
  .eq {{ font-size: 1.25rem; font-weight: bold; color: #0f4c5c; direction: ltr; text-align: center; padding: 10px; background: #e2e8f0; border-radius: 6px; }}
  table {{ width: 100%; border-collapse: collapse; margin: 10px 0; }}
  th, td {{ border: 1px solid #cbd5e1; padding: 8px 10px; text-align: right; font-size: 0.9rem; }}
  th {{ background: #e2e8f0; color: #0f4c5c; }}
  @media print {{ body {{ margin: 12mm; }} }}
</style>
</head>
<body>
  <div class="header">
    <h1>{title_ar}</h1>
    <div>تاريخ التقرير: <span dir="ltr">{ts}</span> | مستوى الدليل: <strong>{_esc(report_data.get('mode_label_ar'))}</strong> | الموثوقية: {_esc(report_data.get('confidence_label_ar'))}</div>
  </div>

  <div class="box">
    <h3>1. المعادلة الكيميائية الموزونة والمعادلات الأيونية (Balanced &amp; Ionic Equations)</h3>
    <div class="eq">{_esc(report_data.get('equation_unicode'))}</div>
    <p><strong>المعادلة الأيونية الكاملة (Complete Ionic):</strong> <span dir="ltr">{_esc(report_data.get('complete_ionic_equation', '—'))}</span></p>
    <p><strong>المعادلة الأيونية الصافية (Net Ionic):</strong> <span dir="ltr">{_esc(report_data.get('net_ionic_equation', '—'))}</span></p>
    <p><strong>الأيونات المتفرجة (Spectator Ions):</strong> {_esc(' ، '.join(report_data.get('spectator_ions') or []) or 'لا توجد')}</p>
    <p><strong>حالة الاتزان الكيميائي:</strong> {_esc(equil.get('equilibrium_type_ar', '—'))} — {_esc(equil.get('driving_force_ar', ''))}</p>
  </div>

  <div class="box">
    <h3>2. الظروف التجريبية المستخدمة (Experimental Conditions)</h3>
    <p><strong>درجة الحرارة:</strong> {_esc(conds.get('temperature'))} | <strong>الضغط:</strong> {_esc(conds.get('pressure'))} | <strong>المذيب/الوسط:</strong> {_esc(conds.get('solvent'))} | <strong>الحفاز:</strong> {_esc(conds.get('catalyst'))}</p>
    <p><strong>توافق الظروف:</strong> {_esc(conds.get('compatibility_status_ar'))}</p>
  </div>

  <div class="box">
    <h3>3. المواد والكميات والحسابات الستوكيومترية (Materials &amp; Stoichiometry)</h3>
    <table>
      <thead>
        <tr><th>المادة المتفاعلة</th><th>المعامل</th><th>الكتلة المولية</th><th>الدور</th><th>الابتدائية</th><th>المستهلكة</th><th>المتبقية</th></tr>
      </thead>
      <tbody>{stoich_rows_html}</tbody>
    </table>
    <h4>المردود النظري للنواتج (Theoretical Yields — بدون ادعاء مردود فعلي تجريبي)</h4>
    <table>
      <thead>
        <tr><th>الناتج</th><th>المعامل</th><th>الكتلة المولية</th><th>المردود النظري (mol)</th><th>المردود النظري (g)</th></tr>
      </thead>
      <tbody>{yield_rows_html or '<tr><td colspan="5">لا توجد نواتج محسوبة</td></tr>'}</tbody>
    </table>
  </div>

  <div class="box">
    <h3>4. التفسير العلمي والملاحظات المخبرية (Scientific Explanation &amp; Observations)</h3>
    <p>{_esc(report_data.get('explanation_ar'))}</p>
  </div>

  <div class="box">
    <h3>5. حدود التنبؤ ومستوى الدليل وتحذيرات السلامة (Prediction Limits &amp; Safety)</h3>
    <ul>{limits_html}</ul>
    <h4>تحذيرات السلامة المخبرية:</h4>
    <ul>{safety_html}</ul>
    <p><em>{_esc(safety.get('virtual_lab_disclaimer_ar', 'نتائج المختبر افتراضية تعليمية وليست بديلاً عن التجارب الواقعية.'))}</em></p>
  </div>

  <div class="box">
    <h3>6. المصادر والمراجع العلمية الموثقة (Scientific Sources &amp; References)</h3>
    <ul>{sources_html}</ul>
  </div>
</body>
</html>"""

    @classmethod
    def generate_experiment_pdf_bytes(cls, report_data: dict[str, Any]) -> bytes:
        """
        Generates a valid PDF document (bytes starting with %PDF-) containing the full
        experiment summary (materials, quantities, conditions, balanced equation, calculations,
        evidence tier, sources, and limitations).
        """
        lines = cls._build_pdf_text_lines(report_data)
        try:
            from reportlab.lib.pagesizes import A4
            from reportlab.pdfgen import canvas

            buf = BytesIO()
            c = canvas.Canvas(buf, pagesize=A4)
            width, height = A4
            c.setTitle(report_data.get("title_en", "MolecuLens Experiment Report"))
            c.setAuthor("MolecuLens Scientific Platform")

            y = height - 45
            c.setFont("Helvetica-Bold", 14)
            c.drawString(40, y, report_data.get("title_en", "MolecuLens Scientific Experiment Report"))
            y -= 18
            c.setFont("Helvetica", 9)
            c.drawString(40, y, f"Generated UTC: {report_data.get('generated_at_utc', '')}")
            y -= 14
            c.line(40, y, width - 40, y)
            y -= 18

            for section_title, section_lines in lines:
                if y < 80:
                    c.showPage()
                    y = height - 45
                c.setFont("Helvetica-Bold", 11)
                c.drawString(40, y, section_title[:95])
                y -= 15
                c.setFont("Helvetica", 9.5)
                for ln in section_lines:
                    # Wrap long lines to 95 chars
                    chunks = [ln[i : i + 92] for i in range(0, max(1, len(ln)), 92)]
                    for ch in chunks:
                        if y < 55:
                            c.showPage()
                            y = height - 45
                            c.setFont("Helvetica", 9.5)
                        c.drawString(50, y, ch)
                        y -= 13
                y -= 6

            c.save()
            return buf.getvalue()
        except Exception:
            return cls._build_minimal_pdf_bytes(lines)

    @classmethod
    def _build_pdf_text_lines(cls, report_data: dict[str, Any]) -> list[tuple[str, list[str]]]:
        rtype = report_data.get("report_type", "reaction_experiment")
        if rtype == "titration_experiment":
            analyte = report_data.get("analyte") or {}
            titrant = report_data.get("titrant") or {}
            eq = report_data.get("equivalence_point") or {}
            endpt = report_data.get("endpoint") or {}
            ind = report_data.get("recommended_indicator") or {}
            return [
                (
                    "1. Titration Setup (Materials, Concentrations & Volumes)",
                    [
                        f"Analyte in Flask: {analyte.get('formula', '')} ({analyte.get('name_ar', '')}) | Conc = {analyte.get('concentration_m', '')} M | Vol = {analyte.get('volume_ml', '')} mL",
                        f"Titrant in Burette: {titrant.get('formula', '')} ({titrant.get('name_ar', '')}) | Conc = {titrant.get('concentration_m', '')} M",
                        f"Titration Type: {report_data.get('case_label_ar', '')} | Initial pH = {report_data.get('initial_ph', '')}",
                    ],
                ),
                (
                    "2. Equivalence Point vs. Colorimetric Indicator Endpoint",
                    [
                        f"Stoichiometric Equivalence Point: V_eq = {eq.get('volume_ml', '')} mL | pH_eq = {eq.get('ph', '')}",
                        f"Indicator Endpoint ({ind.get('name_en', ind.get('name_ar', ''))}): V_end = {endpt.get('volume_ml', '')} mL | pH_end = {endpt.get('ph', '')}",
                    ],
                ),
                (
                    "3. Model Limitations & Virtual Simulation Disclaimer",
                    [
                        "Ideal aqueous equilibrium calculation (25 C, activity coefficients gamma = 1.0).",
                        "Virtual educational simulation only; not a substitute for physical laboratory titration.",
                    ],
                ),
            ]

        conds = report_data.get("conditions") or {}
        reactant_lines = []
        for r in report_data.get("reactants_table") or []:
            reactant_lines.append(
                f"{r.get('formula', '')} (coeff={r.get('coefficient', 1)}, MM={r.get('molar_mass', '')} g/mol): "
                f"Initial={r.get('initial_mol', '')} mol ({r.get('initial_mass_g', '')} g) | "
                f"Consumed={r.get('consumed_mol', '')} mol | Remaining={r.get('remaining_mol', '')} mol"
            )
        if not reactant_lines:
            for r in report_data.get("reactants") or []:
                reactant_lines.append(f"Reactant: {r.get('formula', '')} (coeff={r.get('coefficient', 1)}, state={r.get('state', '')})")

        yield_lines = []
        for y in report_data.get("theoretical_yields") or []:
            yield_lines.append(
                f"Product {y.get('formula', '')} (coeff={y.get('coefficient', 1)}, MM={y.get('molar_mass', '')} g/mol): "
                f"Theoretical Yield = {y.get('theoretical_yield_mol', '')} mol ({y.get('theoretical_yield_g', '')} g)"
            )
        if not yield_lines:
            yield_lines.append("No products formed (No Reaction or Insufficient Evidence).")

        source_lines = []
        for s in report_data.get("sources") or []:
            source_lines.append(f"[{s.get('source_id', 'REF')}] {s.get('source', '')} ({s.get('reference_url', '')})")
        if not source_lines:
            source_lines.append("No reaction reference cited (non-reacting or insufficient evidence).")

        equil = report_data.get("equilibrium_info") or {}
        return [
            (
                "1. Evidence Tier & Balanced Equations",
                [
                    f"Evidence Tier: {report_data.get('mode', '')} | Supported: {report_data.get('supported', False)}",
                    f"Balanced Equation: {report_data.get('equation_text') or report_data.get('equation_unicode', 'N/A')}",
                    f"Reversibility / Equilibrium: {'Reversible (Equilibrium)' if equil.get('is_reversible') else 'Completion-Driven'}",
                    f"Spectator Ions: {', '.join(report_data.get('spectator_ions') or []) or 'None'}",
                ],
            ),
            (
                "2. Experimental Conditions",
                [
                    f"Temperature: {conds.get('temperature', '25 C')} | Pressure: {conds.get('pressure', '1 atm')}",
                    f"Solvent / Medium: {conds.get('solvent', 'aqueous')} | Catalyst: {conds.get('catalyst', 'none')}",
                ],
            ),
            ("3. Reactants & Stoichiometric Breakdown", reactant_lines),
            ("4. Theoretical Yields (100% Limiting Conversion Upper Bound)", yield_lines),
            ("5. Scientific Sources & Citations", source_lines),
            (
                "6. Prediction Limits & Safety Disclaimer",
                [
                    "Theoretical yield assumes 100% conversion of limiting reagent; actual yield requires experimental data.",
                    "Virtual educational simulation (In-Silico); not a substitute for real laboratory experiments or SDS.",
                ],
            ),
        ]

    @staticmethod
    def _build_minimal_pdf_bytes(sections: list[tuple[str, list[str]]]) -> bytes:
        """Pure-Python PDF 1.4 fallback generator if ReportLab is absent."""
        flat_lines = ["MolecuLens Scientific Experiment Report"]
        for title, s_lines in sections:
            flat_lines.append(title)
            flat_lines.extend(s_lines)

        text_ops = ["BT", "/F1 10 Tf", "40 790 Td", "13 TL"]
        for ln in flat_lines[:52]:
            safe = ln.encode("ascii", "replace").decode("ascii").replace("\\", "\\\\").replace("(", "\\(").replace(")", "\\)")
            text_ops.append(f"({safe}) '")
        text_ops.append("ET")
        stream_data = "\n".join(text_ops).encode("latin1")

        objs = [
            b"1 0 obj\n<< /Type /Catalog /Pages 2 0 R >>\nendobj\n",
            b"2 0 obj\n<< /Type /Pages /Kids [3 0 R] /Count 1 >>\nendobj\n",
            b"3 0 obj\n<< /Type /Page /Parent 2 0 R /MediaBox [0 0 595 842] /Resources << /Font << /F1 4 0 R >> >> /Contents 5 0 R >>\nendobj\n",
            b"4 0 obj\n<< /Type /Font /Subtype /Type1 /BaseFont /Helvetica >>\nendobj\n",
            f"5 0 obj\n<< /Length {len(stream_data)} >>\nstream\n".encode("latin1") + stream_data + b"\nendstream\nendobj\n",
        ]
        out = bytearray(b"%PDF-1.4\n")
        offsets = []
        for obj in objs:
            offsets.append(len(out))
            out.extend(obj)
        xref_pos = len(out)
        out.extend(f"xref\n0 {len(objs)+1}\n0000000000 65535 f \n".encode("latin1"))
        for off in offsets:
            out.extend(f"{off:010d} 00000 n \n".encode("latin1"))
        out.extend(
            f"trailer\n<< /Size {len(objs)+1} /Root 1 0 R >>\nstartxref\n{xref_pos}\n%%EOF\n".encode("latin1")
        )
        return bytes(out)

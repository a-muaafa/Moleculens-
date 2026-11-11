"""
Unit tests for the AIScientificAssistant module in MolecuLens.
Verifies property explanations, Lipinski rule evaluations, step-by-step breakdowns,
and zero hallucination responses.
"""

from education.ai_assistant import AIScientificAssistant


def test_property_explanation_ar_and_en():
    assistant = AIScientificAssistant()

    exp_ar = assistant.explain_property("logp", value=1.5, lang="ar")
    assert "معامل التقسيم" in exp_ar["property_name"]
    assert exp_ar["optimal_range"] == "LogP <= 5.0"
    assert "سليم ومثالي" in exp_ar["context_evaluation"]

    exp_en = assistant.explain_property("mw", value=600.0, lang="en")
    assert "Molecular Weight" in exp_en["property_name"]
    assert exp_en["optimal_range"] == "MW <= 500 Da"
    assert "Exceeds" in exp_en["context_evaluation"]


def test_lipinski_evaluation_compliant():
    assistant = AIScientificAssistant()
    profile = {
        "basic_info": {"name_ar": "أسبرين", "name_en": "Aspirin"},
        "properties": {
            "molecular_weight": 180.16,
            "logp": 1.19,
            "hbd": 1,
            "hba": 3,
            "rotatable_bonds": 2,
        },
    }
    eval_res = assistant.explain_lipinski_evaluation(profile, lang="ar")
    assert eval_res["violations_count"] == 0
    assert eval_res["is_druglike"] is True
    assert "متوافق تماماً" in eval_res["summary"]


def test_step_by_step_breakdown():
    assistant = AIScientificAssistant()
    profile = {
        "basic_info": {"name_ar": "ميتفورمين", "name_en": "Metformin", "formula": "C4H11N5"},
        "properties": {
            "molecular_weight": 129.16,
            "logp": -1.43,
            "hbd": 3,
            "hba": 2,
            "tpsa": 88.0,
            "rotatable_bonds": 0,
        },
    }
    steps = assistant.generate_step_by_step_breakdown(profile, lang="ar")
    assert len(steps) >= 3
    assert any("التركيب الجزيئي" in s["title"] for s in steps)
    assert any("الذوبانية" in s["title"] for s in steps)


def test_answer_student_question():
    assistant = AIScientificAssistant()
    ans = assistant.answer_student_question("ما هي أهمية مساحة السطح القطبية TPSA؟", lang="ar")
    assert "TPSA" in ans
    assert "الحاجز الدموي الدماغي" in ans or "النفاذية" in ans


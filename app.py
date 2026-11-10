import io
import json
import os
import joblib
import numpy as np
import pandas as pd
import plotly.graph_objects as go
import streamlit as st
from rdkit import Chem, DataStructs
from rdkit.Chem import Descriptors, Draw, Lipinski, QED, rdFingerprintGenerator, rdMolDescriptors
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams
from rdkit.Chem.Scaffolds import MurckoScaffold

st.set_page_config(page_title="MolecuLens | المختبر الافتراضي للفرز الحاسوبي", layout="wide")

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
MODEL_PATH = os.path.join(BASE_DIR, "models", "dpp4_rf_model.pkl")
REF_PATH = os.path.join(BASE_DIR, "models", "training_reference.pkl")
METRICS_PATH = os.path.join(BASE_DIR, "models", "model_metrics.json")
CLEAN_REPORT_PATH = os.path.join(BASE_DIR, "data", "processed", "cleaning_report.json")
CLEAN_DATA_PATH = os.path.join(BASE_DIR, "data", "processed", "dpp4_cleaned.csv")
TEST_PRED_PATH = os.path.join(BASE_DIR, "data", "processed", "test_predictions.csv")

BENCHMARK_COMPOUNDS = {
    "سيتاجليبتين (Sitagliptin - مثبط مرجعي معتمد لإنزيم DPP4)": "N[C@@H](CC(=O)N1CCn2c(nnc2C(F)(F)F)C1)Cc1cc(F)c(F)cc1F",
    "فيلداجليبتين (Vildagliptin - مثبط مرجعي لإنزيم DPP4)": "N#C[C@@H]1CCCN1C(=O)CNC12CC3CC(CC(O)(C3)C1)C2",
    "مركب نشط من بيانات التدريب (ChEMBL Active Lead)": "N[C@@H](CC1CCCCC1)C(=O)N1CCCCC1",
    "مركب يحتوي على بنية متداخلة مخبريا (Rhodanine PAINS Example)": "O=C1CSC(=S)N1Cc1ccc(O)c(O)c1",
    "الأسبرين (Aspirin - مركب خارج نطاق التدريب لمثبطات DPP4)": "CC(=O)Oc1ccccc1C(=O)O",
}

st.markdown(
    """
<style>
@import url('https://fonts.googleapis.com/css2?family=Cairo:wght@400;500;600;700;800&family=JetBrains+Mono:wght@400;600&display=swap');

html, body, [class*="css"], .stApp {
    font-family: 'Cairo', sans-serif;
    direction: rtl;
    text-align: right;
}

.block-container {
    padding-top: 1.4rem;
    padding-bottom: 3rem;
    max-width: 1380px;
}

#MainMenu, footer {
    visibility: hidden;
}

.lims-header {
    background: linear-gradient(135deg, #0b192c 0%, #0f3460 55%, #0d5c63 100%);
    color: #f8fafc;
    padding: 1.7rem 2.2rem;
    border-radius: 14px;
    border: 1px solid rgba(255, 255, 255, 0.12);
    margin-bottom: 1rem;
    box-shadow: 0 10px 28px rgba(11, 25, 44, 0.18);
}

.lims-header-top {
    display: flex;
    justify-content: space-between;
    align-items: center;
    flex-wrap: wrap;
    gap: 0.8rem;
    margin-bottom: 0.6rem;
}

.lims-title {
    font-size: 1.9rem;
    font-weight: 800;
    margin: 0;
    color: #ffffff;
    letter-spacing: 0.3px;
}

.lims-tag {
    background: rgba(255, 255, 255, 0.14);
    border: 1px solid rgba(255, 255, 255, 0.24);
    color: #e2e8f0;
    padding: 0.25rem 0.85rem;
    border-radius: 6px;
    font-size: 0.85rem;
    font-weight: 600;
}

.lims-sub {
    margin: 0;
    font-size: 0.98rem;
    color: #cbd5e1;
    line-height: 1.65;
}

.gov-banner {
    background: #fffbeb;
    border-right: 5px solid #b45309;
    border-top: 1px solid #fde68a;
    border-left: 1px solid #fde68a;
    border-bottom: 1px solid #fde68a;
    color: #78350f;
    padding: 0.75rem 1.15rem;
    border-radius: 8px;
    margin-bottom: 1.2rem;
    font-size: 0.92rem;
    font-weight: 600;
}

.station-bar {
    display: grid;
    grid-template-columns: repeat(5, 1fr);
    gap: 0.6rem;
    margin-bottom: 1.3rem;
}

.station-step {
    background: var(--secondary-background-color);
    border: 1px solid rgba(128, 128, 128, 0.25);
    border-top: 4px solid #0d5c63;
    border-radius: 8px;
    padding: 0.65rem 0.85rem;
}

.station-step .num {
    font-size: 0.78rem;
    font-weight: 700;
    color: #0d5c63;
    text-transform: uppercase;
}

.station-step .lbl {
    font-size: 0.92rem;
    font-weight: 700;
    margin-top: 0.15rem;
}

.station-step .desc {
    font-size: 0.78rem;
    opacity: 0.8;
    margin-top: 0.1rem;
}

.panel-card {
    background: var(--secondary-background-color);
    border: 1px solid rgba(128, 128, 128, 0.22);
    border-radius: 10px;
    padding: 1.15rem 1.3rem;
    margin-bottom: 1rem;
}

.panel-card h3, .panel-card h4 {
    margin-top: 0;
    margin-bottom: 0.55rem;
    font-weight: 700;
}

.prop-card {
    background: var(--secondary-background-color);
    border: 1px solid rgba(128, 128, 128, 0.24);
    border-right: 4px solid #0f3460;
    border-radius: 8px;
    padding: 0.9rem 1rem;
    margin-bottom: 0.8rem;
    min-height: 138px;
}

.prop-title {
    font-size: 0.86rem;
    font-weight: 700;
    opacity: 0.85;
}

.prop-val {
    font-size: 1.35rem;
    font-weight: 800;
    margin: 0.25rem 0;
    direction: ltr;
    unicode-bidi: embed;
    text-align: right;
    font-family: 'JetBrains Mono', monospace;
}

.prop-explain {
    font-size: 0.82rem;
    line-height: 1.5;
    opacity: 0.9;
}

.status-pill {
    display: inline-block;
    padding: 0.22rem 0.75rem;
    border-radius: 6px;
    font-weight: 700;
    font-size: 0.86rem;
    color: #ffffff;
}

.pill-high { background: #15803d; }
.pill-mid  { background: #b45309; }
.pill-low  { background: #b91c1c; }
.pill-info { background: #0f3460; }

.ltr-box, code {
    direction: ltr;
    unicode-bidi: embed;
    font-family: 'JetBrains Mono', monospace;
    text-align: left;
    display: inline-block;
}

.stTextInput input {
    direction: ltr;
    text-align: left;
    font-family: 'JetBrains Mono', monospace;
}

.stButton > button {
    width: 100%;
    border-radius: 8px;
    font-weight: 700;
    padding: 0.55rem 0.9rem;
    background: #0f3460;
    color: #ffffff;
    border: 1px solid #1e40af;
}

.stButton > button:hover {
    background: #0d5c63;
    color: #ffffff;
    border-color: #0d5c63;
}

.conf-table {
    width: 100%;
    border-collapse: collapse;
    margin-top: 0.6rem;
    font-size: 0.9rem;
}

.conf-table th, .conf-table td {
    border: 1px solid rgba(128, 128, 128, 0.3);
    padding: 0.6rem 0.8rem;
    text-align: right;
}

.conf-table th {
    background: rgba(15, 52, 96, 0.12);
    font-weight: 700;
}
</style>
""",
    unsafe_allow_html=True,
)


@st.cache_resource
def load_artifacts():
    model = joblib.load(MODEL_PATH) if os.path.exists(MODEL_PATH) else None
    ref = joblib.load(REF_PATH) if os.path.exists(REF_PATH) else None
    metrics = None
    if os.path.exists(METRICS_PATH):
        with open(METRICS_PATH, "r", encoding="utf-8") as f:
            metrics = json.load(f)
    clean_rep = None
    if os.path.exists(CLEAN_REPORT_PATH):
        with open(CLEAN_REPORT_PATH, "r", encoding="utf-8") as f:
            clean_rep = json.load(f)
    clean_df = pd.read_csv(CLEAN_DATA_PATH) if os.path.exists(CLEAN_DATA_PATH) else None
    test_df = pd.read_csv(TEST_PRED_PATH) if os.path.exists(TEST_PRED_PATH) else None
    return model, ref, metrics, clean_rep, clean_df, test_df


@st.cache_resource
def load_pains_catalog():
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    return FilterCatalog(params)


@st.cache_resource
def get_morgan_generator():
    return rdFingerprintGenerator.GetMorganGenerator(radius=2, fpSize=1024)


def mol_to_png_bytes(mol, size=(360, 260)):
    img = Draw.MolToImage(mol, size=size)
    buf = io.BytesIO()
    img.save(buf, format="PNG")
    return buf.getvalue()


def compute_query_similarity(mol, ref_bundle):
    gen = get_morgan_generator()
    q_bv = gen.GetFingerprint(mol)
    q_arr = gen.GetFingerprintAsNumPy(mol).astype(np.float32)

    train_df = ref_bundle["train_df"]
    train_X = ref_bundle["train_fps_numpy"]

    dot_prod = np.dot(train_X, q_arr)
    q_bits = float(np.sum(q_arr))
    t_bits = np.sum(train_X, axis=1)
    denom = q_bits + t_bits - dot_prod
    sims = np.where(denom > 0, dot_prod / denom, 0.0)

    order = np.argsort(sims)[::-1]
    top_idx = order[:5]
    max_sim = float(sims[order[0]])
    mean_top5_sim = float(np.mean(sims[top_idx]))

    nearest_df = train_df.iloc[top_idx].copy()
    nearest_df["tanimoto"] = np.round(sims[top_idx], 3)
    return q_arr, sims, max_sim, mean_top5_sim, nearest_df


def evaluate_confidence(sim_value, roc_auc_value):
    if sim_value >= 0.60:
        sim_tier = 2
        sim_label = "مرتفع (>= 0.60)"
        sim_pill = "pill-high"
    elif sim_value >= 0.30:
        sim_tier = 1
        sim_label = "متوسط (0.30 - 0.59)"
        sim_pill = "pill-mid"
    else:
        sim_tier = 0
        sim_label = "منخفض (< 0.30) - خارج نطاق التدريب"
        sim_pill = "pill-low"

    if roc_auc_value >= 0.85:
        auc_tier = 2
        auc_label = "مرتفع (>= 0.85)"
        auc_pill = "pill-high"
    elif roc_auc_value >= 0.70:
        auc_tier = 1
        auc_label = "متوسط (0.70 - 0.84)"
        auc_pill = "pill-mid"
    else:
        auc_tier = 0
        auc_label = "منخفض (< 0.70)"
        auc_pill = "pill-low"

    final_tier = min(sim_tier, auc_tier)
    out_of_domain = sim_value < 0.30

    if final_tier == 2:
        final_label = "عالية (موثوق ضمن نطاق التدريب)"
        final_pill = "pill-high"
    elif final_tier == 1:
        final_label = "متوسطة (يتطلب تحققا إضافيا)"
        final_pill = "pill-mid"
    else:
        final_label = "منخفضة (خارج نطاق التدريب)" if out_of_domain else "منخفضة"
        final_pill = "pill-low"

    return {
        "sim_tier": sim_tier,
        "sim_label": sim_label,
        "sim_pill": sim_pill,
        "auc_tier": auc_tier,
        "auc_label": auc_label,
        "auc_pill": auc_pill,
        "final_tier": final_tier,
        "final_label": final_label,
        "final_pill": final_pill,
        "out_of_domain": out_of_domain,
    }


def classify_potency_educational(pic50_val):
    ic50_nm = 10.0 ** (9.0 - pic50_val)
    if pic50_val >= 7.0:
        return (
            "نشاط تثبيطي مرتفع (ضمن العتبات التعليمية)",
            "pill-high",
            "القيمة المتوقعة تقع تحت 100 nM (pIC50 >= 7.0). يمثل هذا التقدير مؤشرا حاسوبيا أوليا لارتباط قوي محتمل ولا يعد إثباتا مخبريا أو علاجيا.",
            ic50_nm,
        )
    if pic50_val >= 6.0:
        return (
            "نشاط تثبيطي متوسط (ضمن العتبات التعليمية)",
            "pill-mid",
            "القيمة المتوقعة تقع بين 100 و 1000 nM (pIC50 بين 6.0 و 7.0). يشير التقدير الأولي إلى نشاط معتدل قد يتطلب تحسينات بنائية.",
            ic50_nm,
        )
    return (
        "نشاط تثبيطي ضعيف أو غير نشط (ضمن العتبات التعليمية)",
        "pill-low",
        "القيمة المتوقعة أعلى من 1000 nM (pIC50 < 6.0). يشير التقدير الحاسوبي الأولي إلى ضعف الارتباط بالهدف البيولوجي.",
        ic50_nm,
    )


def build_pic50_gauge(pic50_val):
    clamped_val = float(np.clip(pic50_val, 3.0, 10.0))
    fig = go.Figure(
        go.Indicator(
            mode="gauge+number",
            value=clamped_val,
            number={"suffix": " pIC50", "valueformat": ".2f"},
            title={"text": "مقياس التقدير اللوغاريتمي للنشاط (pIC50 من 3.0 إلى 10.0)"},
            gauge={
                "axis": {"range": [3.0, 10.0], "tickvals": [3, 4, 5, 6, 7, 8, 9, 10]},
                "bar": {"color": "#0f3460"},
                "steps": [
                    {"range": [3.0, 6.0], "color": "#fecaca"},
                    {"range": [6.0, 7.0], "color": "#fde68a"},
                    {"range": [7.0, 10.0], "color": "#bbf7d0"},
                ],
                "threshold": {
                    "line": {"color": "#0b192c", "width": 4},
                    "thickness": 0.8,
                    "value": clamped_val,
                },
            },
        )
    )
    fig.update_layout(height=275, margin=dict(l=25, r=25, t=55, b=15))
    return fig


def build_similarity_histogram(all_sims, query_max_sim):
    fig = go.Figure()
    fig.add_trace(
        go.Histogram(
            x=all_sims,
            nbinsx=30,
            marker_color="#0d5c63",
            opacity=0.82,
            name="مركبات التدريب",
        )
    )
    fig.add_vline(
        x=0.30,
        line_dash="dash",
        line_color="#b91c1c",
        annotation_text="حد النطاق الأدنى (0.30)",
        annotation_position="top left",
    )
    fig.add_vline(
        x=0.60,
        line_dash="dash",
        line_color="#15803d",
        annotation_text="عتبة التشابه المرتفع (0.60)",
        annotation_position="top left",
    )
    fig.add_vline(
        x=query_max_sim,
        line_width=3,
        line_color="#0f3460",
        annotation_text=f"أعلى تشابه للمركب ({query_max_sim:.2f})",
        annotation_position="top right",
    )
    fig.update_layout(
        title="توزيع تشابه المركب المدخل (Tanimoto) مقارنة بمركبات التدريب فقط",
        xaxis_title="معامل تشابه Tanimoto مع مجموعة التدريب",
        yaxis_title="عدد مركبات التدريب",
        height=310,
        margin=dict(l=20, r=20, t=45, b=20),
    )
    return fig


if "experiment_history" not in st.session_state:
    st.session_state["experiment_history"] = []
if "student_guess_submitted" not in st.session_state:
    st.session_state["student_guess_submitted"] = False
if "student_potency_guess" not in st.session_state:
    st.session_state["student_potency_guess"] = "نشاط تثبيطي مرتفع (ضمن العتبات التعليمية)"
if "student_conf_guess" not in st.session_state:
    st.session_state["student_conf_guess"] = "عالية (موثوق ضمن نطاق التدريب)"

model, ref_bundle, metrics, clean_report, clean_df, test_df = load_artifacts()
pains_catalog = load_pains_catalog()

st.markdown(
    """
<div class="lims-header">
  <div class="lims-header-top">
    <h1 class="lims-title">MolecuLens | المختبر الافتراضي للفرز الحاسوبي</h1>
    <span class="lims-tag">هاكاثون شراع 2026 - مسار الصحة وجودة الحياة</span>
  </div>
  <p class="lims-sub">
    منصة تعليمية عربية تعمل كمحاكي مهام للفرز الحاسوبي الأولي للمركبات الدوائية المستهدفة لإنزيم DPP-4،
    تدمج الرحلة العلمية الموجهة، التفسير الصيدلاني العربي، وشفافية تقييم الثقة ونطاق التدريب في شاشة موحدة.
  </p>
</div>
<div class="gov-banner">
  حوكمة صارمة: مشروع تعليمي تدريبي بحثي بحت؛ لا يكتشف أدوية جديدة، ولا يصدر أحكاما نهائية على الفعالية أو الأمان، ولا يقدم استشارات علاجية أو سريرية.
</div>
""",
    unsafe_allow_html=True,
)

if model is None or ref_bundle is None or metrics is None:
    st.error("ملفات النموذج أو بصمات التدريب غير موجودة. يرجى تشغيل clean_data.py ثم train_model.py أولا.")
    st.stop()

nav_mode = st.radio(
    "وضع العمل في المختبر الافتراضي:",
    [
        "الرحلة التعليمية الموجهة (المحطات الخمس والتقرير التعليمي)",
        "صفحة التحقق ومختبر الأخطاء والاختبار النقدي",
    ],
    horizontal=True,
)

if nav_mode == "الرحلة التعليمية الموجهة (المحطات الخمس والتقرير التعليمي)":
    st.markdown(
        """
<div class="station-bar">
  <div class="station-step">
    <div class="num">المحطة 1</div>
    <div class="lbl">الهدف البيولوجي</div>
    <div class="desc">UniProt P27487 | DPP-4</div>
  </div>
  <div class="station-step">
    <div class="num">المحطة 2</div>
    <div class="lbl">الخصائص والتنبيهات</div>
    <div class="desc">RDKit &amp; PAINS Filter</div>
  </div>
  <div class="station-step">
    <div class="num">المحطة 3</div>
    <div class="lbl">الأدلة والتشابه</div>
    <div class="desc">ChEMBL &amp; Tanimoto</div>
  </div>
  <div class="station-step">
    <div class="num">المحطة 4</div>
    <div class="lbl">التوقع والتنبؤ</div>
    <div class="desc">تحدي الطالب ونموذج RF</div>
  </div>
  <div class="station-step">
    <div class="num">المحطة 5</div>
    <div class="lbl">الحكم والتقرير</div>
    <div class="desc">معادلة الثقة وملخص التعلم</div>
  </div>
</div>
""",
        unsafe_allow_html=True,
    )

    with st.container():
        st.markdown('<div class="panel-card"><h3>لوحة البداية وإدخال المركب الكيميائي</h3>', unsafe_allow_html=True)
        col_preset, col_input = st.columns([1.1, 1.4])
        with col_preset:
            preset_choice = st.selectbox(
                "اختر تجربة معيارية جاهزة أو أدخل صيغة مخصصة:",
                ["إدخال مخصص"] + list(BENCHMARK_COMPOUNDS.keys()),
                index=1,
            )
            default_smiles = (
                BENCHMARK_COMPOUNDS[preset_choice]
                if preset_choice != "إدخال مخصص"
                else "N[C@@H](CC(=O)N1CCn2c(nnc2C(F)(F)F)C1)Cc1cc(F)c(F)cc1F"
            )
        with col_input:
            smiles_input = st.text_input("الصيغة الخطية للمركب (SMILES):", value=default_smiles)

        mol = Chem.MolFromSmiles(smiles_input.strip()) if smiles_input.strip() else None
        if mol is None:
            st.error("الصيغة الكيميائية المدخلة غير صحيحة بنائيا. يرجى التحقق من تنسيق SMILES قبل بدء التجربة.")
            st.markdown("</div>", unsafe_allow_html=True)
            st.stop()

        canonical_smi = Chem.MolToSmiles(mol, canonical=True)
        scaffold_smi = MurckoScaffold.MurckoScaffoldSmiles(mol=mol, includeChirality=False)
        scaffold_mol = Chem.MolFromSmiles(scaffold_smi) if scaffold_smi else mol
        formula_str = rdMolDescriptors.CalcMolFormula(mol)

        c_prev1, c_prev2, c_prev3 = st.columns([1, 1, 1.3])
        with c_prev1:
            st.image(mol_to_png_bytes(mol, size=(360, 230)), caption="التركيب ثنائي الأبعاد للمركب المدخل")
        with c_prev2:
            st.image(mol_to_png_bytes(scaffold_mol, size=(360, 230)), caption="الهيكل الأساسي (Bemis-Murcko Scaffold)")
        with c_prev3:
            st.markdown(
                f"""
<div class="panel-card" style="margin-bottom:0;">
  <h4>بطاقة التحقق البنائي الفوري</h4>
  <p><strong>حالة التحقق:</strong> <span class="status-pill pill-high">صيغة SMILES صالحة كيميائيا</span></p>
  <p><strong>الصيغة الجزيئية:</strong> <span class="ltr-box">{formula_str}</span></p>
  <p><strong>SMILES القياسي (Canonical):</strong><br><code class="ltr-box">{canonical_smi}</code></p>
  <p><strong>هيكل Murcko الأساسي:</strong><br><code class="ltr-box">{scaffold_smi if scaffold_smi else "لا يوجد نظام حلقي"}</code></p>
</div>
""",
                unsafe_allow_html=True,
            )
        st.markdown("</div>", unsafe_allow_html=True)

    q_arr, all_train_sims, max_sim, mean_top5_sim, nearest_df = compute_query_similarity(mol, ref_bundle)
    tree_preds = np.array([est.predict([q_arr])[0] for est in model.estimators_], dtype=np.float64)
    pred_pic50 = float(np.mean(tree_preds))
    pred_std_pic50 = float(np.std(tree_preds))
    pot_label, pot_pill, pot_explain, pred_ic50_nm = classify_potency_educational(pred_pic50)
    conf_res = evaluate_confidence(max_sim, float(metrics["roc_auc"]))
    pains_match = pains_catalog.GetFirstMatch(mol)
    pains_desc = pains_match.GetDescription() if pains_match else None

    if not st.session_state["experiment_history"] or st.session_state["experiment_history"][-1]["smiles"] != canonical_smi:
        st.session_state["experiment_history"].append(
            {
                "smiles": canonical_smi,
                "formula": formula_str,
                "pred_pic50": round(pred_pic50, 2),
                "pred_ic50_nm": round(pred_ic50_nm, 2),
                "max_tanimoto": round(max_sim, 3),
                "confidence": conf_res["final_label"],
            }
        )

    view_style = st.radio(
        "طريقة عرض المحطات الخمس:",
        ["عرض الشاشة الموحدة (جميع المحطات متتابعة)", "التنقل بتبويبات المحطات الخمس"],
        horizontal=True,
    )

    def render_station_1():
        st.markdown(
            """
<div class="panel-card">
  <h3>المحطة الأولى: تعريف الهدف البيولوجي (UniProt &amp; ChEMBL Target Profile)</h3>
  <p>
    قبل فحص أي جزيء حاسوبيا، يبدأ الباحث بفهم البروتين المستهدف ووظيفته الفسيولوجية وبنية موقع الارتباط النشط.
  </p>
</div>
""",
            unsafe_allow_html=True,
        )
        c1, c2, c3, c4 = st.columns(4)
        c1.metric("الاسم العلمي للهدف", "Dipeptidyl peptidase 4")
        c2.metric("معرف UniProt المرجعي", "P27487 (DPP4_HUMAN)")
        c3.metric("معرف الهدف في ChEMBL", "CHEMBL284")
        c4.metric("التصنيف الإنزيمي", "EC 3.4.14.5")

        col_bio1, col_bio2 = st.columns(2)
        with col_bio1:
            st.markdown(
                """
<div class="panel-card">
  <h4>الآلية الفسيولوجية وعلاقتها بالسكري من النوع الثاني</h4>
  <p>
    يعمل إنزيم <strong>DPP-4</strong> كببتيداز سيريني يقوم بشطر الببتيدات الثنائية من النهاية الأمينية (N-terminus)
    للهرمونات المعوية المعروفة باسم الإنكريتينات، وعلى رأسها <code>GLP-1</code> و <code>GIP</code>.
  </p>
  <p>
    يؤدي تثبيط إنزيم DPP-4 إلى حماية هرمون <code>GLP-1</code> الداخلي من التحلل السريع، مما يعزز إفراز الإنسولين
    بطريقة معتمدة على تركيز الجلوكوز في الدم ويقلل إفراز الجلوكاجون دون التسبب بهبوط حاد في سكر الدم.
  </p>
</div>
""",
                unsafe_allow_html=True,
            )
        with col_bio2:
            st.markdown(
                """
<div class="panel-card">
  <h4>البنية الفراغية لموقع الارتباط والتحديات الانتقائية</h4>
  <p>
    يتكون الموقع النشط للإنزيم من:
  </p>
  <ul>
    <li><strong>الثالوث التحفيزي (Catalytic Triad):</strong> الأحماض الأمينية <code>Ser630</code> و <code>Asp708</code> و <code>His740</code>.</li>
    <li><strong>الجيب الكاره للماء (S1 Hydrophobic Pocket):</strong> يستوعب المجموعات الحلقية أو المفلورة في المثبطات الانتقائية.</li>
    <li><strong>منطقة الإرساء الأيونية (S2 Subsite):</strong> تحتوي على <code>Glu205</code> و <code>Glu206</code> التي تشكل جسورا ملحية مع مجموعة الأمين في المثبط.</li>
  </ul>
  <p>
    <strong>تنبيه تعليمي:</strong> الارتباط القوي بإنزيم DPP-4 وحده لا يكفي؛ إذ يجب أن يكون المركب انتقائيا ولا يثبط إنزيمي
    <code>DPP-8</code> و <code>DPP-9</code> لتجنب السمية الخلوية والمناعية.
  </p>
</div>
""",
                unsafe_allow_html=True,
            )

    def render_station_2():
        st.markdown(
            """
<div class="panel-card">
  <h3>المحطة الثانية: الخصائص الفيزيائية الكيميائية وفحص التنبيهات الهيكلية (RDKit)</h3>
  <p>
    تحسب هذه المحطة المؤشرات الجزيئية الحاسمة للامتصاص الفموي والنفاذية الخلوية، وتفحص وجود المجموعات المسببة للتداخلات المخبرية.
  </p>
</div>
""",
            unsafe_allow_html=True,
        )
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        tpsa = Descriptors.TPSA(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        rot = Lipinski.NumRotatableBonds(mol)
        qed_val = QED.qed(mol)
        fsp3 = rdMolDescriptors.CalcFractionCSP3(mol)
        rings = rdMolDescriptors.CalcNumRings(mol)

        viol = sum([mw > 500, logp > 5, hbd > 5, hba > 10])

        props = [
            (
                "الوزن الجزيئي (Molecular Weight)",
                f"{mw:.1f} Da",
                "المدى المرجعي لقاعدة ليينسكي: <= 500 Da. المركبات الأصغر حجما تعبر الأغشية الحيوية بسهولة أكبر عند الإعطاء الفموي.",
            ),
            (
                "معامل التوزع الدهني (LogP)",
                f"{logp:.2f}",
                "المدى المرجعي: <= 5.0. يقيس التوازن بين الذوبانية المائية والألفة الدهنية اللازمة لعبور جدار الأمعاء.",
            ),
            (
                "المساحة السطحية القطبية (TPSA)",
                f"{tpsa:.1f} Å²",
                "المدى المرجعي للامتصاص الفموي الجيد: < 140 Å². القيم المرتفعة جدا قد تحد من النفاذية السلبية عبر الأغشية.",
            ),
            (
                "مانحات روابط الهيدروجين (HBD)",
                f"{hbd}",
                "المدى المرجعي: <= 5. زيادة مجموعات OH و NH ترفع القطبية وتزيد طاقة نزع الماء أثناء عبور الغشاء الخلوي.",
            ),
            (
                "مستقبلات روابط الهيدروجين (HBA)",
                f"{hba}",
                "المدى المرجعي: <= 10. تعكس ذرات النيتروجين والأكسجين القادرة على تكوين تفاعلات قطبية مع الهدف والماء.",
            ),
            (
                "الروابط القابلة للدوران (Rotatable Bonds)",
                f"{rot}",
                "المدى المرجعي (Veber): <= 10. المرونة المفرطة للجزيء تزيد فقدان الإنتروبيا عند الارتباط بموقع الإنزيم.",
            ),
            (
                "مؤشر التشابه الدوائي الكمي (QED)",
                f"{qed_val:.2f}",
                "يتراوح بين 0 و 1. يدمج ثماني خصائص فيزيائية كيميائية لقياس مدى قرب الجزيء من خصائص الأدوية الفموية المعروفة.",
            ),
            (
                "نسبة ذرات الكربون المشبعة (Fraction Csp3)",
                f"{fsp3:.2f}",
                "المدى المفضل: >= 0.25. ارتفاع نسبة الكربون ثلاثي الأبعاد يحسن الذوبانية ويقلل التراكم البلوري المسطح.",
            ),
            (
                "عدد الحلقات الكيميائية (Ring Count)",
                f"{rings}",
                "توفر الأنظمة الحلقية الصلابة البنائية اللازمة لتوجيه المجموعات الوظيفية داخل جيوب S1 و S2 في إنزيم DPP-4.",
            ),
        ]

        for i in range(0, len(props), 3):
            cols = st.columns(3)
            for j in range(3):
                title, val, expl = props[i + j]
                with cols[j]:
                    st.markdown(
                        f"""
<div class="prop-card">
  <div class="prop-title">{title}</div>
  <div class="prop-val">{val}</div>
  <div class="prop-explain">{expl}</div>
</div>
""",
                        unsafe_allow_html=True,
                    )

        col_lip, col_pains = st.columns(2)
        with col_lip:
            lip_pill = "pill-high" if viol <= 1 else "pill-mid"
            lip_status = f"مخالفات قاعدة Lipinski: {viol} من 4"
            st.markdown(
                f"""
<div class="panel-card">
  <h4>تقييم قابلية الامتصاص الفموي (Lipinski Rule of Five)</h4>
  <p><span class="status-pill {lip_pill}">{lip_status}</span></p>
  <p>
    تقيم قاعدة Lipinski احتمالية الذوبان والامتصاص الفموي فقط، ولا تقيس قدرة المركب على تثبيط إنزيم DPP-4؛
    فقد يحقق جزيء خامل كافة شروط Lipinski دون أي نشاط بيولوجي على الهدف.
  </p>
</div>
""",
                unsafe_allow_html=True,
            )
        with col_pains:
            if pains_desc:
                st.markdown(
                    f"""
<div class="panel-card">
  <h4>فحص التداخلات المخبرية الهيكلي (PAINS Filter)</h4>
  <p><span class="status-pill pill-low">تم رصد تنبيه PAINS هيكلي: {pains_desc}</span></p>
  <p>
    يحتوي هذا المركب على بنية فرعية معروفة بالتفاعل غير النوعي في الفحوصات الحيوية (Pan-Assay Interference Compounds)،
    مما قد يعطي نتائج إيجابية كاذبة مخبريا نتيجة التكتل أو التفاعل الكيميائي المباشر بدلا من التثبيط الانتقائي.
  </p>
</div>
""",
                    unsafe_allow_html=True,
                )
            else:
                st.markdown(
                    """
<div class="panel-card">
  <h4>فحص التداخلات المخبرية الهيكلي (PAINS Filter)</h4>
  <p><span class="status-pill pill-info">لم ترصد تنبيهات PAINS المعروفة</span></p>
  <p>
    عدم رصد تنبيهات PAINS يعني خلو الجزيء من الأنماط الهيكلية المسجلة في الفهرس كمسببات شائعة للتداخل المخبري الكاذب،
    ولا يعني بأي حال أن المركب آمن سميا أو صالح للاستخدام الدوائي.
  </p>
</div>
""",
                    unsafe_allow_html=True,
                )

    def render_station_3():
        st.markdown(
            """
<div class="panel-card">
  <h3>المحطة الثالثة: الأدلة المخبرية المنشورة (ChEMBL) والتشابه الهيكلي مع مجموعة التدريب (Tanimoto)</h3>
  <p>
    تفصل هذه المحطة فصلا تاما بين القياسات المخبرية الحقيقية المنشورة في قاعدة بيانات ChEMBL وبين تنبؤات النموذج الحاسوبي،
    وتحسب التشابه الهيكلي حصريا مقابل بصمات مركبات التدريب لمنع أي تسرب من مجموعة الاختبار.
  </p>
</div>
""",
            unsafe_allow_html=True,
        )

        exact_matches = pd.DataFrame()
        if clean_df is not None and "canonical_smiles" in clean_df.columns:
            exact_matches = clean_df[clean_df["canonical_smiles"] == canonical_smi]

        if not exact_matches.empty:
            row_ex = exact_matches.iloc[0]
            st.markdown(
                f"""
<div class="panel-card" style="border-right: 5px solid #15803d;">
  <h4>سجل مخبري مطابق في قاعدة بيانات ChEMBL</h4>
  <p>
    المركب المدخل مسجل فعليا في البيانات المنظفة تحت المعرف
    <code class="ltr-box">{row_ex['molecule_chembl_id']}</code>
    في تجربة الربط المخبري <code class="ltr-box">{row_ex['assay_chembl_id']}</code> (النوع: <code class="ltr-box">{row_ex['assay_type']}</code>)
    بعلاقة قياس <code class="ltr-box">{row_ex['standard_relation']}</code> وقيمة مخبرية
    <code class="ltr-box">{row_ex['standard_value']:.2f} {row_ex['standard_units']}</code>
    (<code class="ltr-box">pIC50 = {row_ex['pIC50']:.2f}</code>).
  </p>
</div>
""",
                unsafe_allow_html=True,
            )
        else:
            st.markdown(
                """
<div class="panel-card">
  <h4>حالة المطابقة المباشرة في البيانات المحلية</h4>
  <p>
    هذا المركب غير مطابق تطابقا تاما لأي سجل في مجموعة التدريب، ولذلك يعتمد تقييم النطاق التطبيقي على قياس تشابه
    بصمات Morgan الجزيئية (Tanimoto Similarity) مع أقرب المركبات المقاسة مخبريا في مجموعة التدريب.
  </p>
</div>
""",
                unsafe_allow_html=True,
            )

        m1, m2, m3 = st.columns(3)
        m1.metric("أعلى تشابه مع مركبات التدريب (Max Tanimoto)", f"{max_sim:.3f}")
        m2.metric("متوسط تشابه أقرب 5 مركبات تدريب", f"{mean_top5_sim:.3f}")
        m3.metric("عدد مركبات التدريب المرجعية", f"{len(ref_bundle['train_df'])}")

        st.plotly_chart(build_similarity_histogram(all_train_sims, max_sim), use_container_width=True)

        st.markdown("<h4>أقرب 5 مركبات مرجعية في مجموعة التدريب مع بيانات التجارب والمصادر</h4>", unsafe_allow_html=True)
        display_df = pd.DataFrame(
            {
                "معرف المركب (ChEMBL ID)": nearest_df["molecule_chembl_id"].astype(str),
                "معامل التشابه (Tanimoto)": nearest_df["tanimoto"].round(3),
                "معرف التجربة (Assay ID)": nearest_df["assay_chembl_id"].astype(str),
                "نوع التجربة": nearest_df["assay_type"].astype(str),
                "نوع القياس والعلاقة": nearest_df["standard_type"].astype(str) + " (" + nearest_df["standard_relation"].astype(str) + ")",
                "القيمة المخبرية (IC50 nM)": nearest_df["standard_value"].round(2),
                "قيمة pIC50 المخبرية": nearest_df["pIC50"].round(2),
                "المصدر المرجعي (Document ID)": nearest_df["document_chembl_id"].fillna("غير محدد").astype(str),
                "المجلة العلمية": nearest_df["document_journal"].fillna("غير محدد").astype(str),
                "سنة النشر": nearest_df["document_year"].fillna(0).astype(int).astype(str).replace("0", "غير محدد"),
                "الصيغة القياسية (SMILES)": nearest_df["canonical_smiles"].astype(str),
            }
        )
        st.dataframe(display_df, use_container_width=True, hide_index=True)

        st.markdown("<h4>مقارنة بصرية مع أقرب 3 هياكل كيميائية في مجموعة التدريب</h4>", unsafe_allow_html=True)
        g_cols = st.columns(3)
        for idx_col in range(min(3, len(nearest_df))):
            r_near = nearest_df.iloc[idx_col]
            m_near = Chem.MolFromSmiles(str(r_near["canonical_smiles"]))
            with g_cols[idx_col]:
                if m_near is not None:
                    st.image(
                        mol_to_png_bytes(m_near, size=(320, 210)),
                        caption=(
                            f"{r_near['molecule_chembl_id']} | Tanimoto: {r_near['tanimoto']:.2f} | "
                            f"IC50: {r_near['standard_value']:.1f} nM (pIC50: {r_near['pIC50']:.2f})"
                        ),
                    )

    def render_station_4():
        st.markdown(
            """
<div class="panel-card">
  <h3>المحطة الرابعة: تحدي «خمّن قبل أن ترى» والتنبؤ الحاسوبي (Random Forest Regressor)</h3>
  <p>
    لتعزيز التفكير النقدي وتجنب الثقة العمياء بالذكاء الاصطناعي، سجل توقعك العلمي أولا بناء على خصائص المركب في المحطة الثانية
    ودرجة تشابهه مع الأدلة المخبرية في المحطة الثالثة، ثم قارن حكمك بتقدير النموذج الحاسوبي.
  </p>
</div>
""",
            unsafe_allow_html=True,
        )

        with st.container():
            st.markdown('<div class="panel-card"><h4>تحدي الطالب: خمّن قبل أن ترى النتيجة الحاسوبية</h4>', unsafe_allow_html=True)
            cg1, cg2, cg3 = st.columns([1.2, 1.2, 0.8])
            with cg1:
                pot_guess = st.selectbox(
                    "1. ما توقعك لمستوى النشاط التثبيطي لهذا المركب؟",
                    [
                        "نشاط تثبيطي مرتفع (ضمن العتبات التعليمية)",
                        "نشاط تثبيطي متوسط (ضمن العتبات التعليمية)",
                        "نشاط تثبيطي ضعيف أو غير نشط (ضمن العتبات التعليمية)",
                    ],
                    key="pot_guess_select",
                )
            with cg2:
                conf_guess = st.selectbox(
                    "2. ما توقعك لمستوى الثقة في تنبؤ النموذج لهذا المركب؟",
                    [
                        "عالية (موثوق ضمن نطاق التدريب)",
                        "متوسطة (يتطلب تحققا إضافيا)",
                        "منخفضة (خارج نطاق التدريب)",
                    ],
                    key="conf_guess_select",
                )
            with cg3:
                st.write("")
                st.write("")
                if st.button("تسجيل التوقع ومقارنته"):
                    st.session_state["student_guess_submitted"] = True
                    st.session_state["student_potency_guess"] = pot_guess
                    st.session_state["student_conf_guess"] = conf_guess
            st.markdown("</div>", unsafe_allow_html=True)

        col_pred1, col_pred2 = st.columns([1.15, 1.25])
        with col_pred1:
            st.plotly_chart(build_pic50_gauge(pred_pic50), use_container_width=True)
        with col_pred2:
            st.markdown(
                f"""
<div class="panel-card">
  <h4>التقدير الحاسوبي الأولي للنموذج (Random Forest على مقياس pIC50)</h4>
  <p>
    <strong>قيمة pIC50 المتوقعة:</strong>
    <span class="ltr-box" style="font-size:1.25rem;font-weight:800;">{pred_pic50:.2f} ± {pred_std_pic50:.2f}</span>
  </p>
  <p>
    <strong>القيمة المكافئة لتركيز التثبيط النصفي (IC50):</strong>
    <span class="ltr-box" style="font-size:1.2rem;font-weight:700;">{pred_ic50_nm:,.2f} nM</span>
  </p>
  <p>
    <strong>التصنيف الاسترشادي:</strong>
    <span class="status-pill {pot_pill}">{pot_label}</span>
  </p>
  <p style="margin-top:0.6rem;">{pot_explain}</p>
  <p style="font-size:0.84rem;opacity:0.85;margin-top:0.5rem;">
    <strong>ملاحظة منهجية:</strong> يعتمد المقياس على التحويل اللوغاريتمي <code class="ltr-box">pIC50 = 9 - log10(IC50_nM)</code>
    بحيث يغطي القيم الممتدة من <code class="ltr-box">0.1 nM</code> حتى <code class="ltr-box">1,000,000 nM</code> دون إخفاء القيم المرتفعة،
    مع التأكيد على أن هذه العتبات تعليمية استرشادية وليست حكما نهائيا على فعالية المركب.
  </p>
</div>
""",
                unsafe_allow_html=True,
            )

    def render_station_5():
        st.markdown(
            """
<div class="panel-card">
  <h3>المحطة الخامسة: الحكم النهائي، تطبيق معادلة الثقة المعتمدة، والتقرير التعليمي</h3>
  <p>
    تطبيقا للوثيقة المعتمدة للمشروع، لا يحدد مستوى الثقة بناء على التشابه وحده ولا بناء على أداء النموذج وحده،
    بل يحسب كحد أدنى صريح بين المعيارين مع عرض كل مؤشر بشكل منفصل.
  </p>
</div>
""",
            unsafe_allow_html=True,
        )

        roc_auc_val = float(metrics["roc_auc"])
        mae_val = float(metrics["mae_pic50"])
        rmse_val = float(metrics["rmse_pic50"])
        r2_val = float(metrics["r2_pic50"])

        st.markdown(
            f"""
<div class="panel-card">
  <h4>تفصيل حساب مستوى الثقة التشغيلي (القاعدة الحسابية الصريحة)</h4>
  <p>
    <strong>المعادلة المعتمدة:</strong>
    <code>مستوى الثقة النهائي = الحد الأدنى (معيار التشابه الهيكلي للمركب , معيار الأداء العام للنموذج بتقسيم Scaffold Split)</code>
  </p>
  <table class="conf-table">
    <thead>
      <tr>
        <th>المؤشر</th>
        <th>القيمة المحسوبة</th>
        <th>المستوى المعياري</th>
        <th>التفسير المنهجي</th>
      </tr>
    </thead>
    <tbody>
      <tr>
        <td><strong>1. معيار التشابه الهيكلي للمركب (Tanimoto)</strong></td>
        <td><span class="ltr-box">Max = {max_sim:.3f} | Top-5 Mean = {mean_top5_sim:.3f}</span></td>
        <td><span class="status-pill {conf_res['sim_pill']}">{conf_res['sim_label']}</span></td>
        <td>يقيس مدى قرب المركب المدخل من الفضاء الكيميائي لمركبات التدريب فقط.</td>
      </tr>
      <tr>
        <td><strong>2. معيار الأداء العام للنموذج (Scaffold Split ROC-AUC)</strong></td>
        <td><span class="ltr-box">ROC-AUC = {roc_auc_val:.3f} | MAE = {mae_val:.3f} | RMSE = {rmse_val:.3f} | R² = {r2_val:.3f}</span></td>
        <td><span class="status-pill {conf_res['auc_pill']}">{conf_res['auc_label']}</span></td>
        <td>قيمة عالمية ثابتة للنموذج محسوبة على هياكل Murcko غير مرئية لمنع تسرب البيانات.</td>
      </tr>
      <tr>
        <td><strong>3. مستوى الثقة النهائي (حسب القيمة الأدنى)</strong></td>
        <td><span class="ltr-box">min(Tier {conf_res['sim_tier']}, Tier {conf_res['auc_tier']}) = Tier {conf_res['final_tier']}</span></td>
        <td><span class="status-pill {conf_res['final_pill']}">{conf_res['final_label']}</span></td>
        <td>{"تنبيه حرج: المركب خارج نطاق التدريب (Tanimoto < 0.30)، التنبؤ غير موثوق." if conf_res["out_of_domain"] else "محكوم بالحد الأدنى بين المؤشرين لضمان الصدق العلمي."}</td>
      </tr>
    </tbody>
  </table>
  <p style="margin-top:0.8rem;font-size:0.86rem;opacity:0.9;">
    <strong>إقرار منهجي صريح:</strong> هذه العتبات ومستوى الثقة الناتج عنها هي معايير تشغيلية تعليمية مقترحة لتدريب الطالب على الفحص النقدي
    ونطاق التطبيق (Applicability Domain)، وليست احتمالا إحصائيا مطلقا لصحة التنبؤ.
  </p>
</div>
""",
            unsafe_allow_html=True,
        )

        st.markdown('<div class="panel-card"><h4>التقرير التعليمي الشامل للتجربة</h4>', unsafe_allow_html=True)
        student_pot = st.session_state["student_potency_guess"]
        student_conf = st.session_state["student_conf_guess"]
        pot_match = student_pot == pot_label
        conf_match = student_conf == conf_res["final_label"]

        st.markdown(
            f"""
<p><strong>1. مقارنة توقع الطالب («خمّن قبل أن ترى») بالنتيجة الحاسوبية:</strong></p>
<ul>
  <li>
    <strong>توقع النشاط:</strong> اختار الطالب (<em>{student_pot}</em>) بينما أظهر التقدير الحاسوبي (<em>{pot_label}</em>)
    بقيمة <code class="ltr-box">pIC50 = {pred_pic50:.2f}</code> (<code class="ltr-box">{pred_ic50_nm:,.2f} nM</code>) —
    <strong>{"تطابق في تقدير فئة النشاط" if pot_match else "اختلاف يوضح أهمية مراجعة الأدلة التجريبية"}</strong>.
  </li>
  <li>
    <strong>توقع مستوى الثقة:</strong> اختار الطالب (<em>{student_conf}</em>) بينما أظهرت قاعدة الحد الأدنى (<em>{conf_res['final_label']}</em>)
    بناء على <code class="ltr-box">Tanimoto = {max_sim:.3f}</code> و <code class="ltr-box">ROC-AUC = {roc_auc_val:.3f}</code> —
    <strong>{"تطابق في تقييم الموثوقية" if conf_match else "اختلاف يبرز دور قاعدة الحد الأدنى ونطاق التدريب"}</strong>.
  </li>
</ul>
<p><strong>2. الخلاصة التعليمية المستفادة من هذه التجربة:</strong></p>
<ul>
  <li>
    {"المركب يقع خارج نطاق التدريب (Tanimoto < 0.30)، ما يعني أن النموذج يستقرئ في منطقة كيميائية لم يتدرب عليها، وبالتالي لا يعتمد على الرقم المتوقع مهما بدا دقيقا." if conf_res["out_of_domain"] else "المركب يقع ضمن نطاق تشابه مقبول مع مركبات التدريب، ما يمنح التقدير الحاسوبي سياقا مرجعيا قابلا للمقارنة المخبرية."}
  </li>
  <li>
    {f"رصد تنبيه هيكلي ({pains_desc}) يستوجب الحذر من التداخلات المخبرية الكاذبة قبل إجراء أي فحص حيوي." if pains_desc else "لم ترصد تنبيهات PAINS المعروفة، مع بقاء الحاجة لاختبارات الانتقائية والسمية المخبرية."}
  </li>
</ul>
""",
            unsafe_allow_html=True,
        )

        report_md = f"""# تقرير تجربة MolecuLens التعليمية
- الصيغة القياسية (Canonical SMILES): `{canonical_smi}`
- الصيغة الجزيئية: `{formula_str}`
- هيكل Murcko الأساسي: `{scaffold_smi}`
- التقدير الحاسوبي (pIC50): `{pred_pic50:.2f} +/- {pred_std_pic50:.2f}`
- التقدير المكافئ (IC50): `{pred_ic50_nm:.2f} nM`
- فئة النشاط الاسترشادية: {pot_label}
- فحص تنبيهات PAINS: {"تم رصد تنبيه: " + pains_desc if pains_desc else "لم ترصد تنبيهات PAINS المعروفة"}
- أعلى تشابه Tanimoto مع التدريب: `{max_sim:.3f}`
- متوسط تشابه أقرب 5 مركبات تدريب: `{mean_top5_sim:.3f}`
- أداء النموذج العالمي (Scaffold Split): `ROC-AUC = {roc_auc_val:.3f}, MAE = {mae_val:.3f}, RMSE = {rmse_val:.3f}, R2 = {r2_val:.3f}`
- مستوى الثقة التشغيلي النهائي: {conf_res['final_label']}
- توقع الطالب للنشاط: {student_pot}
- توقع الطالب للثقة: {student_conf}
"""
        st.download_button(
            label="تحميل التقرير التعليمي للتجربة (ملف نصي موثق)",
            data=report_md.encode("utf-8"),
            file_name="MolecuLens_Experiment_Report.md",
            mime="text/markdown",
        )
        st.markdown("</div>", unsafe_allow_html=True)

        if st.session_state["experiment_history"]:
            st.markdown("<h4>سجل التجارب المنفذة خلال الجلسة الحالية</h4>", unsafe_allow_html=True)
            hist_df = pd.DataFrame(st.session_state["experiment_history"]).rename(
                columns={
                    "smiles": "SMILES القياسي",
                    "formula": "الصيغة الجزيئية",
                    "pred_pic50": "pIC50 المتوقع",
                    "pred_ic50_nm": "IC50 المتوقع (nM)",
                    "max_tanimoto": "أعلى تشابه (Tanimoto)",
                    "confidence": "مستوى الثقة النهائي",
                }
            )
            st.dataframe(hist_df, use_container_width=True, hide_index=True)

    if view_style == "عرض الشاشة الموحدة (جميع المحطات متتابعة)":
        render_station_1()
        render_station_2()
        render_station_3()
        render_station_4()
        render_station_5()
    else:
        t1, t2, t3, t4, t5 = st.tabs(
            [
                "المحطة 1: الهدف (UniProt)",
                "المحطة 2: الخصائص و PAINS",
                "المحطة 3: الأدلة و Tanimoto",
                "المحطة 4: التوقع والتنبؤ",
                "المحطة 5: الحكم والتقرير",
            ]
        )
        with t1:
            render_station_1()
        with t2:
            render_station_2()
        with t3:
            render_station_3()
        with t4:
            render_station_4()
        with t5:
            render_station_5()

else:
    st.markdown(
        """
<div class="panel-card">
  <h3>صفحة التحقق المنهجي، مختبر الأخطاء، والاختبار النقدي</h3>
  <p>
    تكشف هذه الصفحة بشفافية كاملة كيف تم تنظيف بيانات ChEMBL، ولماذا استخدم تقسيم الهياكل (Bemis-Murcko Scaffold Split)
    بدلا من التقسيم العشوائي، وأين يخطئ نموذج الذكاء الاصطناعي عند مواجهة المركبات البعيدة أو قفزات النشاط (Activity Cliffs).
  </p>
</div>
""",
        unsafe_allow_html=True,
    )

    if clean_report is not None:
        st.markdown('<div class="panel-card"><h4>1. سياسة التنظيف العلمي للبيانات (Data Curation Funnel)</h4>', unsafe_allow_html=True)
        c_s1, c_s2, c_s3, c_s4, c_s5 = st.columns(5)
        c_s1.metric("السجلات الخام من ChEMBL", clean_report["initial_records"])
        c_s2.metric("بعد تصفية العلاقة (= فقط)", clean_report["after_exact_relation"])
        c_s3.metric("بعد استبعاد تحذيرات الصلاحية", clean_report["after_data_validity"])
        c_s4.metric("بعد توحيد التجارب والوحدات (nM)", clean_report["after_unit_and_structure_standardization"])
        c_s5.metric("المركبات الفريدة النهائية", clean_report["final_unique_compounds"])

        st.markdown(
            f"""
<ul>
  <li><strong>معالجة العلاقات غير المحددة:</strong> تم الإبقاء حصرا على القياسات الكمية الدقيقة ذات العلاقة <code class="ltr-box">=</code> واستبعاد القياسات الحدية (<code class="ltr-box">&gt;</code> و <code class="ltr-box">&gt;&gt;</code>).</li>
  <li><strong>معالجة صلاحية البيانات (Data Validity):</strong> تم استبعاد كافة السجلات الحاملة لتحذير <code>Outside typical range</code>.</li>
  <li><strong>توحيد الوحدات وتجانس التجارب:</strong> تم الاقتصار على تجارب الربط المخبري المباشر (<code class="ltr-box">Assay Type = B</code>) وتوحيد <code class="ltr-box">{clean_report['converted_non_nm_units']}</code> سجلا من وحدة <code class="ltr-box">ug.mL-1</code> إلى <code class="ltr-box">nM</code> بالاعتماد على الوزن الجزيئي الدقيق.</li>
  <li><strong>معالجة التكرار وتضارب القياسات:</strong> تم استبعاد <code class="ltr-box">{clean_report['discordant_replicates_removed']}</code> سجلا متكررا يتجاوز فرق القياس فيها وحدة لوغاريتمية كاملة (<code class="ltr-box">&gt; 1.0 pIC50</code>)، ودمج القياسات المتوافقة بأخذ الوسيط الحسابي عبر <code class="ltr-box">{clean_report['unique_murcko_scaffolds']}</code> هيكلا أساسيا فريدا.</li>
</ul>
</div>
""",
            unsafe_allow_html=True,
        )

    st.markdown(
        '<div class="panel-card"><h4>2. تقييم النموذج: مقارنة التقسيم الهيكلي (Scaffold Split) بالتقسيم العشوائي (Random Split)</h4>',
        unsafe_allow_html=True,
    )
    rand_comp = metrics["random_split_comparison"]
    comp_df = pd.DataFrame(
        [
            {
                "منهجية التقسيم": "تقسيم الهياكل المعتمد (Bemis-Murcko Scaffold Split)",
                "عدد مركبات التدريب": metrics["train_count"],
                "عدد مركبات الاختبار": metrics["test_count"],
                "الهياكل المشتركة المتسربة": metrics["scaffold_leakage_count"],
                "خطأ MAE (pIC50)": metrics["mae_pic50"],
                "خطأ RMSE (pIC50)": metrics["rmse_pic50"],
                "معامل التحديد (R²)": metrics["r2_pic50"],
                "مساحة المنحنى (ROC-AUC)": metrics["roc_auc"],
            },
            {
                "منهجية التقسيم": "التقسيم العشوائي التقليدي (Random Split 80/20 - للمقارنة فقط)",
                "عدد مركبات التدريب": metrics["train_count"],
                "عدد مركبات الاختبار": metrics["test_count"],
                "الهياكل المشتركة المتسربة": rand_comp["shared_scaffolds_leaked"],
                "خطأ MAE (pIC50)": rand_comp["mae_pic50"],
                "خطأ RMSE (pIC50)": rand_comp["rmse_pic50"],
                "معامل التحديد (R²)": rand_comp["r2_pic50"],
                "مساحة المنحنى (ROC-AUC)": rand_comp["roc_auc"],
            },
        ]
    )
    st.dataframe(comp_df, use_container_width=True, hide_index=True)
    st.markdown(
        f"""
<p>
  <strong>لماذا لا نعرض عبارة «دقة النموذج» اعتمادا على R² وحده؟</strong>
  لأن التقسيم العشوائي يسرب <code class="ltr-box">{rand_comp['shared_scaffolds_leaked']}</code> هيكلا كيميائيا مشتركا بين التدريب والاختبار،
  مما يضخم قيمة <code class="ltr-box">R²</code> ظاهريا إلى <code class="ltr-box">{rand_comp['r2_pic50']:.3f}</code>.
  أما عند اختبار النموذج بصرامة علمية على هياكل جديدة تماما عبر <strong>Scaffold Split</strong>
  (تسرب الهياكل = <code class="ltr-box">0</code>)، نحصل على تقييم واقعي يتضمن أخطاء التنبؤ المطلقة
  (<code class="ltr-box">MAE = {metrics['mae_pic50']:.3f}</code> و <code class="ltr-box">RMSE = {metrics['rmse_pic50']:.3f}</code>)
  إلى جانب قدرة التمييز الترتيبية (<code class="ltr-box">ROC-AUC = {metrics['roc_auc']:.3f}</code>).
</p>
</div>
""",
        unsafe_allow_html=True,
    )

    if test_df is not None and not test_df.empty:
        st.markdown('<div class="panel-card"><h4>3. مختبر الأخطاء (Error Lab): تحليل أخطاء مجموعة الاختبار</h4>', unsafe_allow_html=True)
        fig_err = go.Figure()
        fig_err.add_trace(
            go.Scatter(
                x=test_df["true_pIC50"],
                y=test_df["pred_pIC50"],
                mode="markers",
                marker=dict(
                    size=9,
                    color=test_df["max_train_tanimoto"],
                    colorscale="Viridis",
                    showscale=True,
                    colorbar=dict(title="Tanimoto"),
                ),
                text=test_df["molecule_chembl_id"]
                + " | خطأ مطلق: "
                + test_df["abs_error_pIC50"].round(2).astype(str)
                + " | Tanimoto: "
                + test_df["max_train_tanimoto"].round(2).astype(str),
                hoverinfo="text",
                name="مركبات الاختبار",
            )
        )
        min_v = float(min(test_df["true_pIC50"].min(), test_df["pred_pIC50"].min()) - 0.3)
        max_v = float(max(test_df["true_pIC50"].max(), test_df["pred_pIC50"].max()) + 0.3)
        fig_err.add_trace(
            go.Scatter(
                x=[min_v, max_v],
                y=[min_v, max_v],
                mode="lines",
                line=dict(color="#b91c1c", dash="dash"),
                name="خط التطابق المثالي",
            )
        )
        fig_err.update_layout(
            title="القيم المخبرية الحقيقية مقابل المتوقعة (pIC50) لمركبات الاختبار (Scaffold Split)",
            xaxis_title="قيمة pIC50 المخبرية الحقيقية",
            yaxis_title="قيمة pIC50 المتوقعة حاسوبيا",
            height=380,
            margin=dict(l=20, r=20, t=45, b=20),
        )
        st.plotly_chart(fig_err, use_container_width=True)

        st.markdown("<h5>أبرز حالات الخطأ الحقيقي في مجموعة الاختبار (لدراسة حدود الذكاء الاصطناعي)</h5>", unsafe_allow_html=True)
        worst_df = test_df.sort_values("abs_error_pIC50", ascending=False).head(8)
        err_table = pd.DataFrame(
            {
                "معرف المركب (ChEMBL ID)": worst_df["molecule_chembl_id"].astype(str),
                "أعلى تشابه مع التدريب (Tanimoto)": worst_df["max_train_tanimoto"].round(3),
                "أقرب مركب في التدريب": worst_df["nearest_train_chembl_id"].astype(str),
                "pIC50 الحقيقي": worst_df["true_pIC50"].round(2),
                "pIC50 المتوقع": worst_df["pred_pIC50"].round(2),
                "الخطأ المطلق |ΔpIC50|": worst_df["abs_error_pIC50"].round(2),
                "IC50 الحقيقي (nM)": worst_df["true_ic50_nm"].round(1),
                "IC50 المتوقع (nM)": worst_df["pred_ic50_nm"].round(1),
                "تفسير سبب الخطأ": np.where(
                    worst_df["max_train_tanimoto"] < 0.35,
                    "ضعف التشابه مع مركبات التدريب (هيكل جديد)",
                    "قفزة نشاط هيكلية (Activity Cliff) رغم تقارب البصمة",
                ),
            }
        )
        st.dataframe(err_table, use_container_width=True, hide_index=True)
        st.markdown("</div>", unsafe_allow_html=True)

    st.markdown(
        '<div class="panel-card"><h4>4. الاختبار القياسي للفهم النقدي للثقة والتشابه (5 أسئلة تقييمية)</h4>',
        unsafe_allow_html=True,
    )
    q1 = st.radio(
        "السؤال 1: إذا حصل مركب على تشابه Tanimoto مرتفع (0.72) مع مركبات التدريب، لكن أداء النموذج العام (ROC-AUC) كان متوسطا (0.78)، فما مستوى الثقة النهائي وفق وثيقة MolecuLens؟",
        [
            "ثقة عالية لأن التشابه الهيكلي أعلى من 0.60",
            "ثقة متوسطة تطبيقا لقاعدة الحد الأدنى بين معيار التشابه ومعيار أداء النموذج",
            "ثقة مطلقة بنسبة 100%",
        ],
        key="quiz_q1",
    )
    q2 = st.radio(
        "السؤال 2: لماذا نستخدم تقسيم الهياكل (Scaffold Split) بدلا من التقسيم العشوائي عند تقييم النموذج؟",
        [
            "لتسريع وقت تدريب النموذج في بايثون",
            "لمنع تسرب المركبات المتشابهة في هيكلها الأساسي بين التدريب والاختبار وقياس قدرة النموذج على التعميم",
            "لزيادة قيمة R² وجعل دقة النموذج تبدو أعلى",
        ],
        key="quiz_q2",
    )
    q3 = st.radio(
        "السؤال 3: ماذا يعني عدم ظهور أي تنبيه في فحص PAINS للمركب المدخل؟",
        [
            "أن المركب دواء آمن تماما وخال من السمية السريرية",
            "أنه لم ترصد أنماط هيكلية معروفة بالتسبب في تداخلات مخبرية كاذبة، ولا يعني إثبات الأمان الدوائي",
            "أن المركب يثبط إنزيم DPP-4 بفعالية مؤكدة",
        ],
        key="quiz_q3",
    )
    q4 = st.radio(
        "السؤال 4: لماذا نقارن المركب المدخل ببصمات مركبات التدريب فقط عند حساب Tanimoto وليس بكامل ملف البيانات المنظفة؟",
        [
            "لتجنب المقارنة بمركبات مجموعة الاختبار وضمان قياس نطاق التطبيق الحقيقي لما تعلمه النموذج",
            "لأن مركبات الاختبار لا تحتوي على صيغ SMILES",
            "لأن حساب Tanimoto لا يعمل إلا على أقل من 500 مركب",
        ],
        key="quiz_q4",
    )
    q5 = st.radio(
        "السؤال 5: إذا كان تشابه Tanimoto للمركب المدخل يساوي 0.18 وتنبأ النموذج بنشاط 45 nM، فما الحكم العلمي الصحيح؟",
        [
            "اعتماد المركب كمثبط قوي لأن القيمة أقل من 100 nM",
            "تصنيف التنبؤ كمنخفض الثقة مع تنبيه «خارج نطاق التدريب» لعدم وجود تغطية كيميائية مشابهة في التدريب",
            "الاكتفاء بقاعدة Lipinski للحكم على فعالية المركب",
        ],
        key="quiz_q5",
    )

    if st.button("تصحيح الاختبار النقدي وعرض النتيجة"):
        correct_answers = [
            q1 == "ثقة متوسطة تطبيقا لقاعدة الحد الأدنى بين معيار التشابه ومعيار أداء النموذج",
            q2 == "لمنع تسرب المركبات المتشابهة في هيكلها الأساسي بين التدريب والاختبار وقياس قدرة النموذج على التعميم",
            q3 == "أنه لم ترصد أنماط هيكلية معروفة بالتسبب في تداخلات مخبرية كاذبة، ولا يعني إثبات الأمان الدوائي",
            q4 == "لتجنب المقارنة بمركبات مجموعة الاختبار وضمان قياس نطاق التطبيق الحقيقي لما تعلمه النموذج",
            q5 == "تصنيف التنبؤ كمنخفض الثقة مع تنبيه «خارج نطاق التدريب» لعدم وجود تغطية كيميائية مشابهة في التدريب",
        ]
        score = sum(correct_answers)
        st.markdown(
            f'<p><span class="status-pill {"pill-high" if score >= 4 else "pill-mid"}">النتيجة: {score} من 5 إجابات صحيحة</span></p>',
            unsafe_allow_html=True,
        )
    st.markdown("</div>", unsafe_allow_html=True)


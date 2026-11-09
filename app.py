import os
import numpy as np
import pandas as pd
import joblib
import plotly.graph_objects as go
import streamlit as st
from rdkit import Chem, DataStructs
from rdkit.Chem import AllChem, Draw, Descriptors, Lipinski
from rdkit.Chem.FilterCatalog import FilterCatalog, FilterCatalogParams

st.set_page_config(page_title="MolecuLens", page_icon="🧪", layout="wide")

ROOT = os.path.dirname(os.path.abspath(__file__))
MODEL_PATH = os.path.join(ROOT, "models", "dpp4_rf_model.pkl")
DATA_PATH = os.path.join(ROOT, "data", "processed", "dpp4_cleaned.csv")

# ---------------------------------------------------------------- التنسيق (CSS)
st.markdown(
    """
<style>
@import url('https://fonts.googleapis.com/css2?family=Cairo:wght@400;600;700&display=swap');

html, body, [class*="css"], .stApp {
    font-family: 'Cairo', sans-serif;
    direction: rtl;
    text-align: right;
}
.block-container { padding-top: 2rem; max-width: 1200px; }
#MainMenu, footer { visibility: hidden; }

.hero {
    background: linear-gradient(135deg, #0f766e 0%, #1d4ed8 100%);
    color: #fff; padding: 2rem 2.2rem; border-radius: 18px;
    margin-bottom: 1.5rem; box-shadow: 0 8px 24px rgba(29,78,216,.25);
}
.hero h1 { margin: 0; font-size: 2.2rem; color: #fff; }
.hero p  { margin: .5rem 0 0; opacity: .92; font-size: 1.05rem; }

.card {
    background: var(--secondary-background-color);
    border: 1px solid rgba(128,128,128,.25);
    border-radius: 14px; padding: 1.1rem 1.3rem; margin-bottom: 1rem;
}
.card h4 { margin: 0 0 .6rem; }

.badge {
    display: inline-block; padding: .25rem .8rem; border-radius: 999px;
    font-weight: 700; font-size: .95rem; color: #fff;
}
.b-high { background: #16a34a; } .b-mid { background: #f59e0b; }
.b-low  { background: #dc2626; } .b-ok  { background: #16a34a; }
.b-warn { background: #f59e0b; }

[data-testid="stMetric"] {
    background: var(--secondary-background-color);
    border: 1px solid rgba(128,128,128,.25);
    border-radius: 14px; padding: .8rem 1rem;
}
.stButton > button {
    width: 100%; border-radius: 12px; font-weight: 700; padding: .6rem 0;
    background: linear-gradient(135deg, #0f766e, #1d4ed8); color: #fff; border: 0;
}
.stButton > button:hover { filter: brightness(1.1); color: #fff; }
.stTextInput input { direction: ltr; text-align: left; font-family: monospace; }
.smiles-ltr, code { direction: ltr; unicode-bidi: embed; }
</style>
""",
    unsafe_allow_html=True,
)

st.markdown(
    """
<div class="hero">
  <h1>🧪 MolecuLens</h1>
  <p>منصة الفحص الحاسوبي الافتراضي لتقدير نشاط مثبطات DPP-4،
     مع فحص PAINS وخصائص الدواء وقياس التشابه مع المركبات المعروفة.</p>
</div>
""",
    unsafe_allow_html=True,
)


# ---------------------------------------------------------------- الموارد
@st.cache_resource
def load_model():
    return joblib.load(MODEL_PATH) if os.path.exists(MODEL_PATH) else None


@st.cache_resource
def load_pains():
    params = FilterCatalogParams()
    params.AddCatalog(FilterCatalogParams.FilterCatalogs.PAINS)
    return FilterCatalog(params)


@st.cache_resource
def load_reference():
    """بصمات مركبات التدريب لحساب تشابه Tanimoto."""
    if not os.path.exists(DATA_PATH):
        return None
    df = pd.read_csv(DATA_PATH).dropna(subset=["canonical_smiles"])
    df = df.drop_duplicates("canonical_smiles").reset_index(drop=True)
    fps, keep = [], []
    for i, smi in enumerate(df["canonical_smiles"]):
        m = Chem.MolFromSmiles(str(smi))
        if m is not None:
            fps.append(AllChem.GetMorganFingerprintAsBitVect(m, 2, nBits=1024))
            keep.append(i)
    return df.iloc[keep].reset_index(drop=True), fps


def potency_info(v):
    if v < 100:
        return "نشاط عالٍ", "b-high", "مرشح واعد؛ يستحق المتابعة المخبرية."
    if v < 1000:
        return "نشاط متوسط", "b-mid", "يفضّل إجراء تحسينات هيكلية."
    return "نشاط ضعيف", "b-low", "قد لا يكون مثبطاً فعالاً لهذا الهدف."


def gauge(v):
    fig = go.Figure(
        go.Indicator(
            mode="gauge+number",
            value=v,
            number={"suffix": " nM"},
            title={"text": "IC50 المتوقع"},
            gauge={
                "axis": {"range": [0, 2000]},
                "bar": {"color": "#1d4ed8"},
                "steps": [
                    {"range": [0, 100], "color": "#bbf7d0"},
                    {"range": [100, 1000], "color": "#fde68a"},
                    {"range": [1000, 2000], "color": "#fecaca"},
                ],
            },
        )
    )
    fig.update_layout(height=260, margin=dict(l=20, r=20, t=50, b=10))
    return fig


# ---------------------------------------------------------------- الإدخال
model = load_model()
pains = load_pains()

left, right = st.columns([2, 1])
with left:
    smiles = st.text_input(
        "أدخل الصيغة الكيميائية (SMILES):", "CC(=O)OC1=CC=CC=C1C(=O)O"
    )
with right:
    st.write("")
    st.write("")
    run = st.button("🔬 بدء الفحص والتنبؤ")

if run:
    if model is None:
        st.error("الرجاء تشغيل سكربت التدريب أولاً لتوليد ملف النموذج.")
        st.stop()
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        st.error("الصيغة المدخلة غير صحيحة، يرجى التحقق من تنسيق SMILES.")
        st.stop()

    fp = AllChem.GetMorganFingerprintAsBitVect(mol, 2, nBits=1024)
    arr = np.zeros((1024,))
    DataStructs.ConvertToNumpyArray(fp, arr)
    pred = float(model.predict([arr])[0])
    label, cls, advice = potency_info(pred)

    tab1, tab2, tab3 = st.tabs(["📊 النتيجة", "💊 الخصائص الدوائية", "🔗 المركبات المشابهة"])

    # ---- النتيجة
    with tab1:
        c1, c2 = st.columns([1, 1.3])
        with c1:
            st.image(Draw.MolToImage(mol, size=(400, 400)), caption="الهيكل الكيميائي")
        with c2:
            st.plotly_chart(gauge(pred), use_container_width=True)
            st.markdown(
                f'<div class="card"><h4>التصنيف</h4>'
                f'<span class="badge {cls}">{label}</span><p style="margin-top:.7rem">{advice}</p></div>',
                unsafe_allow_html=True,
            )
            hit = pains.GetFirstMatch(mol)
            if hit:
                st.markdown(
                    f'<div class="card"><h4>فحص PAINS</h4><span class="badge b-warn">⚠️ تحذير</span>'
                    f'<p style="margin-top:.7rem">يحتوي على مجموعة مزعجة قد تسبب تداخلات مخبرية: '
                    f'<code>{hit.GetDescription()}</code></p></div>',
                    unsafe_allow_html=True,
                )
            else:
                st.markdown(
                    '<div class="card"><h4>فحص PAINS</h4><span class="badge b-ok">✅ آمن</span>'
                    '<p style="margin-top:.7rem">خالٍ من تنبيهات PAINS المعروفة.</p></div>',
                    unsafe_allow_html=True,
                )

    # ---- الخصائص (Lipinski)
    with tab2:
        mw = Descriptors.MolWt(mol)
        logp = Descriptors.MolLogP(mol)
        hbd = Lipinski.NumHDonors(mol)
        hba = Lipinski.NumHAcceptors(mol)
        tpsa = Descriptors.TPSA(mol)
        rot = Lipinski.NumRotatableBonds(mol)
        viol = sum([mw > 500, logp > 5, hbd > 5, hba > 10])

        m = st.columns(3)
        m[0].metric("الوزن الجزيئي", f"{mw:.1f}")
        m[1].metric("LogP", f"{logp:.2f}")
        m[2].metric("TPSA", f"{tpsa:.1f}")
        m = st.columns(3)
        m[0].metric("مانحات H", hbd)
        m[1].metric("مستقبلات H", hba)
        m[2].metric("روابط دوّارة", rot)
        if viol == 0:
            st.success("✅ يحقق قاعدة Lipinski بالكامل.")
        else:
            st.warning(f"⚠️ عدد مخالفات قاعدة Lipinski: {viol}")

    # ---- التشابه
    with tab3:
        ref = load_reference()
        if ref is None:
            st.info("ملف البيانات المعالجة غير موجود لحساب التشابه.")
        else:
            df, fps = ref
            sims = np.array(DataStructs.BulkTanimotoSimilarity(fp, fps))
            top = np.argsort(sims)[::-1][:5]
            best = sims[top[0]]
            if best >= 0.6:
                st.success(f"أعلى تشابه Tanimoto = {best:.2f} ← مستوى الثقة: **عالية** (موثوق ضمن النطاق).")
            elif best > 0.3:
                st.warning(f"أعلى تشابه Tanimoto = {best:.2f} ← مستوى الثقة: **متوسطة**.")
            else:
                st.error(f"أعلى تشابه Tanimoto = {best:.2f} ← مستوى الثقة: **منخفضة**. تنبيه: خارج نطاق التدريب.")
            out = pd.DataFrame(
                {
                    "SMILES": df.loc[top, "canonical_smiles"].values,
                    "Tanimoto": sims[top].round(2),
                    "IC50 الفعلي (nM)": df.loc[top, "standard_value"].values,
                }
            )
            st.dataframe(out, use_container_width=True, hide_index=True)
else:
    st.info("أدخل SMILES واضغط «بدء الفحص والتنبؤ» لعرض النتائج.")

st.markdown(
    '<div class="card" style="margin-top:1.5rem"><b>⚠️ حكمة صارمة:</b> '
    "مشروع تعليمي تدريبي بحثي؛ لا يكتشف أدوية جديدة ولا يقدم استشارات علاجية.</div>",
    unsafe_allow_html=True,
)

import json
import os
import joblib
import numpy as np
import pandas as pd
import streamlit as st
from chemistry.descriptors import PROPERTY_EDUCATIONAL_GUIDE, calculate_descriptors
from chemistry.molecular_structure import get_murcko_scaffold, render_mol_png, validate_smiles
from config.settings import (
    APP_NAME,
    APP_SLOGAN,
    CLEAN_DATA_PATH,
    CLEAN_REPORT_PATH,
    DEMO_COMPOUNDS_PATH,
    METRICS_PATH,
    MODEL_PATH,
    TARGET_INFO,
    TEAM_NAME,
    TEST_PRED_PATH,
    TRAINING_REF_PATH,
)
from model.predict import predict_dpp4_activity
from reports.generate_report import generate_educational_report
from search.autocomplete import (
    build_instant_search_catalog,
    format_suggestion_label,
    get_autocomplete_suggestions,
)
from search.compound_search import (
    SEARCH_ONLINE_STATE_KEY,
    SEARCH_QUERY_STATE_KEY,
    SEARCH_TYPE_STATE_KEY,
    search_compound,
    update_search_session_state,
    verify_compound_identity,
)
from search.local_dictionary import get_all_local_compounds
from visualization.charts import (
    build_nearest_compounds_chart,
    build_pic50_gauge,
    build_properties_chart,
    build_published_activity_distribution_chart,
    build_similarity_chart,
    build_split_comparison_chart,
    build_test_scatter_chart,
)

st.set_page_config(page_title=f"{APP_NAME} | {APP_SLOGAN}", layout="wide")

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
    padding-top: 1.2rem;
    padding-bottom: 3rem;
    max-width: 1400px;
}

#MainMenu, footer {
    visibility: hidden;
}

.hero-banner {
    background: linear-gradient(135deg, #0b192c 0%, #0f3460 55%, #0d5c63 100%);
    color: #f8fafc;
    padding: 1.7rem 2.2rem;
    border-radius: 14px;
    border: 1px solid rgba(255, 255, 255, 0.14);
    margin-bottom: 1rem;
    box-shadow: 0 10px 28px rgba(11, 25, 44, 0.16);
}

.hero-top {
    display: flex;
    justify-content: space-between;
    align-items: center;
    flex-wrap: wrap;
    gap: 0.8rem;
    margin-bottom: 0.55rem;
}

.hero-brand {
    font-size: 2rem;
    font-weight: 800;
    margin: 0;
    color: #ffffff;
}

.hero-slogan {
    direction: ltr;
    unicode-bidi: embed;
    font-family: 'JetBrains Mono', monospace;
    background: rgba(255, 255, 255, 0.14);
    border: 1px solid rgba(255, 255, 255, 0.25);
    color: #e2e8f0;
    padding: 0.28rem 0.9rem;
    border-radius: 6px;
    font-size: 0.86rem;
    font-weight: 600;
}

.hero-desc {
    margin: 0;
    font-size: 0.98rem;
    color: #cbd5e1;
    line-height: 1.7;
}

.gov-notice {
    background: #fffbeb;
    border-right: 5px solid #b45309;
    border-top: 1px solid #fde68a;
    border-left: 1px solid #fde68a;
    border-bottom: 1px solid #fde68a;
    color: #78350f;
    padding: 0.75rem 1.15rem;
    border-radius: 8px;
    margin-bottom: 1.15rem;
    font-size: 0.91rem;
    font-weight: 600;
}

.stage-grid {
    display: grid;
    grid-template-columns: repeat(auto-fit, minmax(210px, 1fr));
    gap: 0.65rem;
    margin-bottom: 1.2rem;
}

.stage-box {
    background: var(--secondary-background-color);
    border: 1px solid rgba(128, 128, 128, 0.25);
    border-top: 4px solid #0d5c63;
    border-radius: 8px;
    padding: 0.7rem 0.9rem;
}

.stage-box .st-num {
    font-size: 0.78rem;
    font-weight: 700;
    color: #0d5c63;
}

.stage-box .st-title {
    font-size: 0.94rem;
    font-weight: 700;
    margin-top: 0.15rem;
}

.stage-box .st-sub {
    font-size: 0.78rem;
    opacity: 0.82;
    margin-top: 0.1rem;
}

.panel-card {
    background: var(--secondary-background-color);
    border: 1px solid rgba(128, 128, 128, 0.22);
    border-radius: 10px;
    padding: 1.15rem 1.35rem;
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
    padding: 0.95rem 1.05rem;
    margin-bottom: 0.85rem;
}

.prop-val {
    font-size: 1.3rem;
    font-weight: 800;
    margin: 0.25rem 0 0.45rem 0;
    direction: ltr;
    unicode-bidi: embed;
    text-align: right;
    font-family: 'JetBrains Mono', monospace;
}

.status-pill {
    display: inline-block;
    padding: 0.22rem 0.78rem;
    border-radius: 6px;
    font-weight: 700;
    font-size: 0.86rem;
    color: #ffffff;
}

.pill-high { background: #15803d; }
.pill-mid  { background: #b45309; }
.pill-low  { background: #b91c1c; }
.pill-info { background: #0f3460; }

.ltr-inline, code {
    direction: ltr;
    unicode-bidi: embed;
    font-family: 'JetBrains Mono', monospace;
    text-align: left;
    display: inline-block;
}

.sep-col-chem {
    border-top: 4px solid #0f3460;
    background: var(--secondary-background-color);
    padding: 1rem;
    border-radius: 8px;
    border-right: 1px solid rgba(128,128,128,0.2);
    border-left: 1px solid rgba(128,128,128,0.2);
    border-bottom: 1px solid rgba(128,128,128,0.2);
}

.sep-col-exp {
    border-top: 4px solid #15803d;
    background: var(--secondary-background-color);
    padding: 1rem;
    border-radius: 8px;
    border-right: 1px solid rgba(128,128,128,0.2);
    border-left: 1px solid rgba(128,128,128,0.2);
    border-bottom: 1px solid rgba(128,128,128,0.2);
}

.sep-col-ai {
    border-top: 4px solid #b45309;
    background: var(--secondary-background-color);
    padding: 1rem;
    border-radius: 8px;
    border-right: 1px solid rgba(128,128,128,0.2);
    border-left: 1px solid rgba(128,128,128,0.2);
    border-bottom: 1px solid rgba(128,128,128,0.2);
}

.stButton > button {
    width: 100%;
    border-radius: 8px;
    font-weight: 700;
    padding: 0.5rem 0.85rem;
    background: #0f3460;
    color: #ffffff;
    border: 1px solid #1e40af;
}

.stButton > button:hover {
    background: #0d5c63;
    color: #ffffff;
    border-color: #0d5c63;
}

.suggestion-card {
    background: var(--secondary-background-color);
    border: 1px solid rgba(15, 52, 96, 0.25);
    border-right: 4px solid #0d5c63;
    border-radius: 10px;
    padding: 0.75rem 0.95rem;
    margin-bottom: 0.65rem;
    direction: rtl;
    text-align: right;
}

.match-badge {
    display: inline-block;
    padding: 0.18rem 0.55rem;
    border-radius: 999px;
    font-size: 0.78rem;
    font-weight: 700;
    background: #0f3460;
    color: #ffffff;
    margin-left: 0.4rem;
}
</style>
""",
    unsafe_allow_html=True,
)


@st.cache_resource
def load_all_resources():
    model = joblib.load(MODEL_PATH) if os.path.exists(MODEL_PATH) else None
    ref_bundle = joblib.load(TRAINING_REF_PATH) if os.path.exists(TRAINING_REF_PATH) else None
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
    demo_compounds = []
    if os.path.exists(DEMO_COMPOUNDS_PATH):
        with open(DEMO_COMPOUNDS_PATH, "r", encoding="utf-8") as f:
            demo_compounds = json.load(f)
    return model, ref_bundle, metrics, clean_rep, clean_df, test_df, demo_compounds


model, ref_bundle, metrics, clean_report, clean_df, test_df, demo_compounds = load_all_resources()

update_search_session_state(st.session_state, default_query="سيتاجليبتين")
if "confirmed_compound" not in st.session_state:
    initial_res = search_compound("سيتاجليبتين", allow_network=False)
    st.session_state["confirmed_compound"] = initial_res.get("selected")
if "confirmed_compound_key" not in st.session_state:
    comp0 = st.session_state.get("confirmed_compound")
    st.session_state["confirmed_compound_key"] = comp0.get("key") if isinstance(comp0, dict) else None
if "search_verification_ok" not in st.session_state:
    st.session_state["search_verification_ok"] = True
if "selected_override_compound" not in st.session_state:
    st.session_state["selected_override_compound"] = None
if "student_guesses_by_smiles" not in st.session_state:
    st.session_state["student_guesses_by_smiles"] = {}


def _set_query_callback(new_query: str, override_comp: dict | None = None) -> None:
    st.session_state[SEARCH_QUERY_STATE_KEY] = new_query
    st.session_state["selected_override_compound"] = override_comp


def _on_text_query_change() -> None:
    st.session_state["selected_override_compound"] = None


def _on_instant_catalog_change() -> None:
    chosen_lbl = st.session_state.get("instant_catalog_selector", "")
    if not chosen_lbl or chosen_lbl.startswith("--"):
        return
    for item in build_instant_search_catalog():
        if item["instant_label"] == chosen_lbl:
            st.session_state[SEARCH_QUERY_STATE_KEY] = item["name_ar"]
            st.session_state["selected_override_compound"] = item
            break

st.markdown(
    f"""
<div class="hero-banner">
  <div class="hero-top">
    <div>
      <h1 class="hero-brand">{APP_NAME} | المختبر الافتراضي لاستكشاف المركبات الدوائية</h1>
      <div style="font-size:0.88rem;color:#94a3b8;font-weight:600;">الفريق التطويري: {TEAM_NAME}</div>
    </div>
    <span class="hero-slogan">{APP_SLOGAN}</span>
  </div>
  <p class="hero-desc">
    مختبر افتراضي تعليمي باللغة العربية يساعد طلاب الصيدلة والعلوم الصحية والكيمياء على البحث الذكي عن المركبات،
    وفهم خصائصها الفيزيائية الكيميائية، والاطلاع على الأدلة العلمية المنشورة عنها في ChEMBL، والتعرف على كيفية استخدام
    الذكاء الاصطناعي في التنبؤ بالنشاط الحيوي تجاه الهدف البروتيني البشري DPP4 مع تقييم نطاق تطبيق النموذج.
  </p>
</div>
<div class="gov-notice">
  حدود المشروع وحوكمته العلمية: MolecuLens مختبر تعليمي للفرز الحاسوبي الأولي، وليس أداة لاكتشاف أدوية جديدة أو إثبات فعاليتها العلاجية أو تقديم نصائح طبية.
</div>
""",
    unsafe_allow_html=True,
)

if model is None or ref_bundle is None or metrics is None:
    st.error("ملفات النموذج أو المرجع غير موجودة. يرجى تشغيل clean_data.py ثم train_model.py أولا.")
    st.stop()

main_section = st.radio(
    "أقسام منصة MolecuLens:",
    [
        "1. البحث الذكي والمختبر الافتراضي (المراحل الخمس)",
        "2. مكتبة المركبات التجريبية لـ DPP4 (Demo Library)",
        "3. التحقق من النموذج ومقارنة Random Split و Scaffold Split وتنظيف البيانات",
    ],
    horizontal=True,
)

if main_section == "1. البحث الذكي والمختبر الافتراضي (المراحل الخمس)":
    st.markdown(
        """
<div class="stage-grid">
  <div class="stage-box">
    <div class="st-num">نقطة الدخول</div>
    <div class="st-title">البحث الذكي والتحقق</div>
    <div class="st-sub">عربي / إنجليزي / صيغة / اختصار / PubChem</div>
  </div>
  <div class="stage-box">
    <div class="st-num">المرحلة الأولى</div>
    <div class="st-title">الهدف البروتيني</div>
    <div class="st-sub">Human DPP4 (UniProt P27487)</div>
  </div>
  <div class="stage-box">
    <div class="st-num">المرحلة الثانية</div>
    <div class="st-title">الخصائص الكيميائية</div>
    <div class="st-sub">RDKit, Lipinski &amp; PAINS</div>
  </div>
  <div class="stage-box">
    <div class="st-num">المرحلة الثالثة</div>
    <div class="st-title">الأدلة المنشورة</div>
    <div class="st-sub">ChEMBL Assays &amp; Documents</div>
  </div>
  <div class="stage-box">
    <div class="st-num">المرحلة الرابعة</div>
    <div class="st-title">تحدي التوقع وتنبؤ AI</div>
    <div class="st-sub">Guess Before You See &amp; AD</div>
  </div>
  <div class="stage-box">
    <div class="st-num">المرحلة الخامسة</div>
    <div class="st-title">التقييم والتقرير التعليمي</div>
    <div class="st-sub">فصل الخصائص والأدلة والتنبؤ</div>
  </div>
</div>
""",
        unsafe_allow_html=True,
    )

    @st.fragment
    def render_smart_search_fragment():
        st.markdown('<div class="panel-card"><h3>مربع البحث الذكي عن المركبات (Smart Compound Search &amp; Autocomplete)</h3>', unsafe_allow_html=True)
        st.markdown(
            "<p>ابحث تلقائيا باستخدام الاسم العربي الشائع، الاسم الإنجليزي، الاسم العلمي (IUPAC)، الصيغة الكيميائية (مثل CO2 أو H2O أو O2 أو H2O2 أو C6H12O6)، الاختصار العلمي (مثل ATP أو NADH)، الاسم التجاري (مثل Januvia أو Galvus)، أو صيغة SMILES.</p>",
            unsafe_allow_html=True,
        )

        quick_examples = [
            ("سيتاجليبتين", "سيتاجليبتين"),
            ("Vildagliptin", "Vildagliptin"),
            ("CO2", "CO2"),
            ("الماء", "الماء"),
            ("H2O", "H2O"),
            ("O2", "O2"),
            ("H2O2", "H2O2"),
            ("ATP", "ATP"),
            ("NADH", "NADH"),
            ("C6H12O6 (تعدد مصاوغات)", "C6H12O6"),
        ]
        ex_cols = st.columns(len(quick_examples))
        for idx_ex, (btn_lbl, q_val) in enumerate(quick_examples):
            with ex_cols[idx_ex]:
                st.button(
                    btn_lbl,
                    key=f"q_ex_{idx_ex}",
                    on_click=_set_query_callback,
                    args=(q_val, None),
                )

        c_search1, c_search2, c_search3 = st.columns([1.8, 1.0, 0.9])
        with c_search1:
            query_text = st.text_input(
                "أدخل اسم المركب أو صيغته الكيميائية أو اختصاره:",
                key=SEARCH_QUERY_STATE_KEY,
                on_change=_on_text_query_change,
            )
        with c_search2:
            type_options = [
                ("auto", "تلقائي شامل (Auto)"),
                ("name_ar", "الاسم العربي (Arabic Name)"),
                ("name_en", "الاسم الإنجليزي (English Name)"),
                ("iupac", "الاسم العلمي (IUPAC Name)"),
                ("formula", "الصيغة الكيميائية (Chemical Formula)"),
                ("abbreviation", "الاختصار العلمي (Abbreviation)"),
                ("synonym", "المرادفات (Synonyms)"),
                ("brand", "الاسم التجاري (Brand Name)"),
                ("smiles", "صيغة SMILES"),
            ]
            type_codes = [code for code, _ in type_options]
            type_labels_map = dict(type_options)
            selected_search_type = st.selectbox(
                "تصفية نوع المدخل (دون مسح النص المكتوب):",
                options=type_codes,
                format_func=lambda code: type_labels_map.get(code, code),
                key=SEARCH_TYPE_STATE_KEY,
            )
        with c_search3:
            allow_online = st.checkbox(
                "تفعيل الاتصال بـ PubChem عند عدم توفر المركب محليا",
                key=SEARCH_ONLINE_STATE_KEY,
            )

        catalog_items = build_instant_search_catalog()
        catalog_labels = ["-- اكتب هنا للفلترة الحية الفورية أثناء الكتابة (عربي / إنجليزي / صيغة / اختصار / اسم تجاري) --"] + [
            item["instant_label"] for item in catalog_items
        ]
        st.selectbox(
            "البحث التفاعلي الفوري أثناء الكتابة (Live Autocomplete Palette - لا يتطلب الضغط على Enter):",
            options=catalog_labels,
            index=0,
            key="instant_catalog_selector",
            on_change=_on_instant_catalog_change,
        )

        suggestions = get_autocomplete_suggestions(query_text, limit=6, search_type=selected_search_type)
        if suggestions:
            st.markdown("<strong>قائمة الاقتراحات التفاعلية المطابقة للمدخل الحالي (اضغط على أي مركب لاختياره فورا):</strong>", unsafe_allow_html=True)
            sugg_cols = st.columns(min(len(suggestions), 3))
            for idx_s, sugg in enumerate(suggestions[:6]):
                with sugg_cols[idx_s % min(len(suggestions), 3)]:
                    s_mol, _, _ = validate_smiles(sugg.get("canonical_smiles", ""))
                    s_cid = sugg.get("pubchem_cid") or "N/A"
                    s_chembl = sugg.get("chembl_id") or "N/A"
                    st.markdown(
                        f"""
<div class="suggestion-card">
  <div><span class="match-badge">{sugg.get('match_category', 'Match')}</span> <strong>{sugg['name_ar']}</strong> | <span class="ltr-inline">{sugg['name_en']}</span></div>
  <div style="font-size:0.86rem;margin-top:0.25rem;">
    <strong>الصيغة:</strong> <code class="ltr-inline">{sugg.get('molecular_formula', '')}</code> |
    <strong>CID:</strong> <code class="ltr-inline">{s_cid}</code> |
    <strong>ChEMBL:</strong> <code class="ltr-inline">{s_chembl}</code>
  </div>
  <div style="font-size:0.82rem;color:#475569;margin-top:0.2rem;"><strong>سبب الظهور:</strong> {sugg.get('match_reason_ar', '')}</div>
</div>
""",
                        unsafe_allow_html=True,
                    )
                    if s_mol is not None and idx_s < 3:
                        with st.expander(f"معاينة البنية: {sugg['name_ar']} ({sugg.get('molecular_formula', '')})", expanded=False):
                            st.image(render_mol_png(s_mol, size=(220, 150)), caption=f"{sugg['name_ar']} | {sugg['name_en']}")
                    st.button(
                        f"اختيار: {sugg['name_ar']} ({sugg.get('molecular_formula', '')})",
                        key=f"pick_sugg_{idx_s}_{sugg.get('key', idx_s)}",
                        on_click=_set_query_callback,
                        args=(sugg["name_ar"], sugg),
                    )

        override_comp = st.session_state.get("selected_override_compound")
        if override_comp is not None and query_text.strip() in {override_comp.get("name_ar"), override_comp.get("name_en"), override_comp.get("molecular_formula")}:
            search_res = {
                "found": True,
                "match_type": "exact",
                "is_ambiguous": False,
                "ambiguity_message": None,
                "candidates": [override_comp],
                "selected": override_comp,
                "did_you_mean": [],
                "data_source": "local_dictionary",
                "connection_status": override_comp.get("source", "القاموس المحلي الموثق"),
                "message": f"تم اختيار المركب تفاعليا ({override_comp.get('match_reason_ar', 'مطابقة مباشرة')}).",
            }
        else:
            search_res = search_compound(query_text, search_type=selected_search_type, allow_network=allow_online)

        st.markdown(
            f"<p><strong>مصدر البحث وحالة الاتصال:</strong> <span class='status-pill pill-info'>{search_res['connection_status']}</span> — {search_res['message']}</p>",
            unsafe_allow_html=True,
        )

        if not search_res["found"]:
            st.warning(search_res["message"])
            if search_res["did_you_mean"]:
                st.markdown("<strong>مركبات مقترحة يمكنك اختيارها مباشرة (Did You Mean):</strong>", unsafe_allow_html=True)
                dym_cols = st.columns(min(len(search_res["did_you_mean"]), 5))
                for idx_d, dym_item in enumerate(search_res["did_you_mean"][:5]):
                    clean_dym = dym_item.split("(")[0].strip()
                    with dym_cols[idx_d]:
                        st.button(
                            dym_item,
                            key=f"dym_btn_{idx_d}",
                            on_click=_set_query_callback,
                            args=(clean_dym, None),
                        )
            st.markdown("</div>", unsafe_allow_html=True)
            prev_ok = st.session_state.get("search_verification_ok", True)
            st.session_state["search_verification_ok"] = False
            st.session_state["confirmed_compound"] = None
            st.session_state["confirmed_compound_key"] = None
            if prev_ok:
                st.rerun()
            return

        if search_res["is_ambiguous"]:
            st.warning(search_res["ambiguity_message"])
            cand_list = search_res["candidates"]
            cand_options = [
                (
                    f"{idx_c + 1}. {c['name_ar']} | {c['name_en']} | الصيغة: {c.get('molecular_formula', '')} "
                    f"| CID: {c.get('pubchem_cid') or 'غير متوفر'} | ChEMBL: {c.get('chembl_id') or 'غير متوفر'} "
                    f"— [{c.get('match_reason_ar', 'تطابق الصيغة')}]"
                )
                for idx_c, c in enumerate(cand_list)
            ]
            chosen_cand_label = st.radio(
                "نظرا لتعدد المركبات أو المصاوغات المشتركة، يرجى اختيار المركب المقصود لتأكيد هويته (لا يتم اختيار مركب عشوائي):",
                cand_options,
                key=f"ambig_choice_{query_text}_{selected_search_type}",
            )
            chosen_idx = cand_options.index(chosen_cand_label)
            active_compound = cand_list[chosen_idx]
        else:
            active_compound = search_res["selected"]

        verification = verify_compound_identity(active_compound)
        if not verification["is_verified"]:
            st.error(verification["error"])
            st.markdown("</div>", unsafe_allow_html=True)
            prev_ok = st.session_state.get("search_verification_ok", True)
            st.session_state["search_verification_ok"] = False
            st.session_state["confirmed_compound"] = None
            st.session_state["confirmed_compound_key"] = None
            if prev_ok:
                st.rerun()
            return

        verified_comp = verification["compound"]
        v_steps = verification["verification_steps"]
        mol_obj_frag, canonical_smi_frag, _ = validate_smiles(verified_comp["canonical_smiles"])
        scaf_mol_frag, _ = get_murcko_scaffold(mol_obj_frag)
        desc_data_frag = calculate_descriptors(mol_obj_frag)

        col_id1, col_id2, col_id3 = st.columns([1.0, 1.0, 1.5])
        with col_id1:
            st.image(render_mol_png(mol_obj_frag, size=(340, 230)), caption=f"البنية ثنائية الأبعاد: {verified_comp['name_ar']} ({verified_comp['name_en']})")
        with col_id2:
            st.image(render_mol_png(scaf_mol_frag, size=(340, 230)), caption="الهيكل الأساسي (Bemis-Murcko Scaffold)")
        with col_id3:
            syns_ar_str = "، ".join(verified_comp.get("synonyms_ar", [])[:4]) or "غير محدد"
            syns_en_str = ", ".join(verified_comp.get("synonyms_en", [])[:4]) or "N/A"
            brands_str = "، ".join(verified_comp.get("brand_names", [])) or "غير متوفر"
            match_reason_str = verified_comp.get("match_reason_ar", "تطابق مباشر")
            match_cat_str = verified_comp.get("match_category", "Exact Match")
            st.markdown(
                f"""
<div class="panel-card" style="margin-bottom:0;">
  <h4>بطاقة التحقق من هوية المركب (Compound Identity Verification)</h4>
  <p>
    <span class="status-pill pill-high">تم التحقق من الهوية والبنية الجزيئية بنجاح</span>
    <span class="match-badge">{match_cat_str}</span>
  </p>
  <p><strong>سبب ظهور النتيجة:</strong> {match_reason_str}</p>
  <p><strong>الاسم العربي:</strong> {verified_comp['name_ar']} | <strong>الاسم الإنجليزي:</strong> <span class="ltr-inline">{verified_comp['name_en']}</span></p>
  <p><strong>الصيغة الجزيئية المستخرجة (RDKit):</strong> <code class="ltr-inline">{desc_data_frag['molecular_formula']}</code> (مطابقة للصيغة المسجلة: <code class="ltr-inline">{verified_comp.get('molecular_formula', '')}</code>)</p>
  <p><strong>الاسم العلمي (IUPAC Name):</strong> <code class="ltr-inline">{verified_comp.get('iupac_name', 'N/A')}</code></p>
  <p><strong>المرادفات (Synonyms):</strong> {syns_ar_str} | <span class="ltr-inline">{syns_en_str}</span></p>
  <p><strong>الاسم التجاري الموثق:</strong> {brands_str}</p>
  <p><strong>المعرفات المرجعية:</strong> PubChem CID: <code class="ltr-inline">{verified_comp.get('pubchem_cid') or 'غير متوفر'}</code> | ChEMBL ID: <code class="ltr-inline">{verified_comp.get('chembl_id') or 'غير متوفر'}</code></p>
  <p><strong>صيغة SMILES الموحدة:</strong> <code class="ltr-inline">{canonical_smi_frag}</code></p>
  <p><strong>فحوصات التحقق الأربع:</strong> صلاحية SMILES: <code class="ltr-inline">{v_steps['smiles_valid']}</code> | استخراج الصيغة: <code class="ltr-inline">{v_steps['formula_extracted']}</code> | تطابق الصيغة: <code class="ltr-inline">{v_steps['formula_verified']}</code> | المعرف المرجعي: <code class="ltr-inline">{v_steps['has_reference_id']}</code></p>
  <p><strong>نبذة تعريفية:</strong> {verified_comp.get('info_ar', '')}</p>
</div>
""",
                unsafe_allow_html=True,
            )
        st.markdown("</div>", unsafe_allow_html=True)

        new_comp_key = f"{verified_comp.get('key', '')}_{canonical_smi_frag}"
        prev_comp_key = st.session_state.get("confirmed_compound_key")
        prev_ok = st.session_state.get("search_verification_ok", False)
        st.session_state["confirmed_compound"] = verified_comp
        st.session_state["confirmed_compound_key"] = new_comp_key
        st.session_state["search_verification_ok"] = True
        if prev_comp_key != new_comp_key or not prev_ok:
            st.rerun()

    render_smart_search_fragment()

    if not st.session_state.get("search_verification_ok") or not st.session_state.get("confirmed_compound"):
        st.stop()

    verified_comp = st.session_state["confirmed_compound"]
    mol_obj, canonical_smi, _ = validate_smiles(verified_comp["canonical_smiles"])
    scaf_mol, scaf_smi = get_murcko_scaffold(mol_obj)
    desc_data = calculate_descriptors(mol_obj)

    pred_data = predict_dpp4_activity(mol_obj, model, ref_bundle)

    published_records = []
    if clean_df is not None and "canonical_smiles" in clean_df.columns:
        matched_rows = clean_df[clean_df["canonical_smiles"] == canonical_smi]
        if not matched_rows.empty:
            published_records.extend(matched_rows.to_dict(orient="records"))
    if not published_records and demo_compounds:
        for d_item in demo_compounds:
            _, d_canon, _ = validate_smiles(d_item["canonical_smiles"])
            if d_canon == canonical_smi or (verified_comp.get("chembl_id") and d_item.get("chembl_id") == verified_comp.get("chembl_id")):
                published_records.append(d_item)

    display_mode = st.radio(
        "طريقة تصفح مراحل المختبر الخمس:",
        ["عرض الرحلة الكاملة المتتابعة (المراحل الخمس في شاشة موحدة)", "التنقل عبر تبويبات المراحل الخمس"],
        horizontal=True,
    )

    def render_stage_1():
        st.markdown(
            f"""
<div class="panel-card">
  <h3>المرحلة الأولى: اختيار الهدف البروتيني (Target Protein Selection)</h3>
  <p>يركز المختبر في إصداره الحالي على الهدف البروتيني البشري المعتمد <strong>{TARGET_INFO['name_en']}</strong>.</p>
</div>
""",
            unsafe_allow_html=True,
        )
        m1, m2, m3, m4 = st.columns(4)
        m1.metric("اسم البروتين المستهدف", "DPP4 (Dipeptidyl Peptidase-4)")
        m2.metric("المعرف المرجعي (UniProt)", TARGET_INFO["uniprot_id"])
        m3.metric("معرف الهدف (ChEMBL)", TARGET_INFO["target_id"])
        m4.metric("الكائن الحي (Organism)", TARGET_INFO["organism"])

        c_t1, c_t2 = st.columns(2)
        with c_t1:
            st.markdown(
                f"""
<div class="panel-card">
  <h4>الوظيفة الحيوية والعلاقة بأبحاث السكري</h4>
  <p><strong>وصف الوظيفة:</strong> {TARGET_INFO['function_ar']}</p>
  <p><strong>علاقته بأبحاث السكري:</strong> {TARGET_INFO['diabetes_relation_ar']}</p>
  <p>
    <strong>الروابط المرجعية المعتمدة:</strong>
    <a href="{TARGET_INFO['uniprot_url']}" target="_blank">صفحة البروتين في UniProt (P27487)</a> |
    <a href="{TARGET_INFO['chembl_url']}" target="_blank">صفحة الهدف في ChEMBL (CHEMBL284)</a>
  </p>
</div>
""",
                unsafe_allow_html=True,
            )
        with c_t2:
            st.markdown(
                f"""
<div class="panel-card">
  <h4>التحقق من هوية الهدف في بيانات التدريب وحدود التفسير</h4>
  <ul>
    <li><strong>تأكيد الكائن الحي والهدف:</strong> تم التحقق برمجيا في مرحلة تنظيف البيانات من أن جميع سجلات التدريب والاختبار البالغ عددها <code class="ltr-inline">{clean_report['final_unique_compounds']}</code> مركبا فريدا تتبع حصريا الهدف البشري <code class="ltr-inline">CHEMBL284</code> والكائن الحي <code class="ltr-inline">Homo sapiens</code>.</li>
    <li><strong>بنية موقع الارتباط:</strong> يضم الموقع النشط الثالوث التحفيزي (<code class="ltr-inline">Ser630, Asp708, His740</code>) والجيب الكاره للماء <code class="ltr-inline">S1</code> ومنطقة الإرساء الأيونية <code class="ltr-inline">Glu205/Glu206</code>.</li>
    <li><strong>تنبيه علمي:</strong> وجود نشاط تثبيطي تجاه DPP4 حاسوبيا أو مخبريا لا يعني إثبات فعالية علاجية لدى مرضى السكري؛ إذ يتطلب تطوير الدواء إثبات الانتقائية عدم تثبيط إنزيمي <code class="ltr-inline">DPP8</code> و <code class="ltr-inline">DPP9</code> ودراسات الحرائك والسمية السريرية.</li>
  </ul>
</div>
""",
                unsafe_allow_html=True,
            )

    def render_stage_2():
        st.markdown(
            """
<div class="panel-card">
  <h3>المرحلة الثانية: تحليل الخصائص الكيميائية (Chemical Properties Analysis - RDKit)</h3>
  <p>
    استخراج الخصائص الفيزيائية الكيميائية للمركب باستخدام مكتبة RDKit مع شرح تعليمي مفصل لكل خاصية (معناها، طريقة حسابها، أهميتها الدوائية، وحدود تفسيرها).
  </p>
</div>
""",
            unsafe_allow_html=True,
        )

        prop_keys = ["molecular_weight", "logp", "hbd", "hba", "rotatable_bonds", "tpsa"]
        for i in range(0, len(prop_keys), 2):
            cols = st.columns(2)
            for j in range(2):
                pk = prop_keys[i + j]
                guide = PROPERTY_EDUCATIONAL_GUIDE[pk]
                val = desc_data[pk]
                with cols[j]:
                    st.markdown(
                        f"""
<div class="prop-card">
  <div style="font-weight:700;font-size:0.95rem;color:#0f3460;">{guide['title_ar']}</div>
  <div class="prop-val">{val} {guide['unit']} <span style="font-size:0.8rem;font-weight:500;opacity:0.75;">(المدى الاسترشادي: {guide['reference_range']})</span></div>
  <p style="margin:0.2rem 0;font-size:0.84rem;"><strong>ما معناها؟</strong> {guide['meaning']}</p>
  <p style="margin:0.2rem 0;font-size:0.84rem;"><strong>كيف يتم حسابها؟</strong> {guide['calculation']}</p>
  <p style="margin:0.2rem 0;font-size:0.84rem;"><strong>لماذا تهم في دراسة المركبات؟</strong> {guide['importance']}</p>
  <p style="margin:0.2rem 0;font-size:0.84rem;color:#78350f;"><strong>ما حدود تفسيرها؟</strong> {guide['limitations']}</p>
</div>
""",
                        unsafe_allow_html=True,
                    )

        fig_props, exp_props = build_properties_chart(desc_data)
        st.plotly_chart(fig_props, use_container_width=True)
        st.caption(exp_props)

        c_lip, c_pains = st.columns(2)
        with c_lip:
            lip = desc_data["lipinski"]
            lip_pill = "pill-high" if lip["passes_rule_of_five"] else "pill-mid"
            st.markdown(
                f"""
<div class="panel-card">
  <h4>قاعدة ليبينسكي الخماسية (Lipinski's Rule of Five)</h4>
  <p><span class="status-pill {lip_pill}">{lip['summary_ar']}</span></p>
  <ul>
    <li>الوزن الجزيئي &lt;= 500 Da: <strong>{"مطابق" if lip['checks']['mw_ok'] else "غير مطابق"}</strong> (<code class="ltr-inline">{desc_data['molecular_weight']} Da</code>)</li>
    <li>معامل التوزع LogP &lt;= 5.0: <strong>{"مطابق" if lip['checks']['logp_ok'] else "غير مطابق"}</strong> (<code class="ltr-inline">{desc_data['logp']}</code>)</li>
    <li>مانحات الهيدروجين HBD &lt;= 5: <strong>{"مطابق" if lip['checks']['hbd_ok'] else "غير مطابق"}</strong> (<code class="ltr-inline">{desc_data['hbd']}</code>)</li>
    <li>مستقبلات الهيدروجين HBA &lt;= 10: <strong>{"مطابق" if lip['checks']['hba_ok'] else "غير مطابق"}</strong> (<code class="ltr-inline">{desc_data['hba']}</code>)</li>
  </ul>
  <p style="font-size:0.84rem;color:#78350f;"><strong>حدود التفسير:</strong> {lip['disclaimer_ar']}</p>
</div>
""",
                unsafe_allow_html=True,
            )
        with c_pains:
            pains = desc_data["pains"]
            pains_pill = "pill-low" if pains["has_pains_alert"] else "pill-info"
            st.markdown(
                f"""
<div class="panel-card">
  <h4>فحص الأنماط المتداخلة مخبريا (PAINS Filter)</h4>
  <p><span class="status-pill {pains_pill}">{pains['status_ar']}</span></p>
  <p style="font-size:0.86rem;">
    يكشف فحص PAINS (Pan-Assay Interference Compounds) عن المجموعات الكيميائية التي تميل للتفاعل غير النوعي
    أو التكتل أو التلوين الكاذب في الفحوصات الحيوية.
  </p>
  <p style="font-size:0.84rem;color:#78350f;"><strong>حدود التفسير:</strong> {pains['disclaimer_ar']}</p>
</div>
""",
                unsafe_allow_html=True,
            )

    def render_stage_3():
        st.markdown(
            """
<div class="panel-card">
  <h3>المرحلة الثالثة: استكشاف الأدلة العلمية المنشورة (Published Experimental Evidence - ChEMBL)</h3>
  <p>
    تعرض هذه المرحلة النتائج التجريبية المنشورة في الأدبيات العلمية وقاعدة بيانات ChEMBL بمعزل تام عن تنبؤات الذكاء الاصطناعي،
    مع توثيق معرف التجربة (Assay ID) ومعرف المرجع العلمي (Document ID) وعدد القياسات المتوافقة.
  </p>
</div>
""",
            unsafe_allow_html=True,
        )

        st.markdown(
            """
<div class="panel-card">
  <h4>توضيح تعليمي: الفرق العلمي بين قياس IC50 وقياس Ki</h4>
  <ul>
    <li><strong>تركيز التثبيط النصفي (IC50):</strong> هو تركيز المثبط اللازم لخفض نشاط الإنزيم بنسبة 50% في ظروف تجربة محددة؛ وتعتمد قيمته على تركيز الإنزيم وتركيز الركيزة (Substrate) في التجربة.</li>
    <li><strong>ثابت التثبيط (Ki):</strong> يعبر عن ألفة الارتباط الجوهرية (Intrinsic Binding Affinity) بين المثبط والإنزيم بمعزل عن تركيز الركيزة (يرتبط بـ IC50 عبر معادلة Cheng-Prusoff).</li>
    <li><strong>القاعدة العلمية في MolecuLens:</strong> يتم التعامل مع كل نوع قياس بصورة منفصلة تماما ولا يتم دمج قيم IC50 مع Ki في تدريب النموذج أو حساب الوسيط.</li>
  </ul>
</div>
""",
            unsafe_allow_html=True,
        )

        highlight_val = None
        highlight_lbl = verified_comp["name_ar"]
        if published_records:
            st.markdown(
                f"<div class='panel-card' style='border-right:5px solid #15803d;'><h4>الأدلة التجريبية المنشورة للمركب المختار ({verified_comp['name_ar']} | {verified_comp['name_en']})</h4>",
                unsafe_allow_html=True,
            )
            ev_rows = []
            for rec in published_records:
                med_val = float(rec.get("median_activity", rec.get("standard_value", 0.0)))
                pic50_v = float(rec.get("pIC50", 9.0 - np.log10(med_val) if med_val > 0 else 0.0))
                highlight_val = pic50_v
                ev_rows.append(
                    {
                        "اسم المركب": f"{verified_comp['name_ar']} ({verified_comp['name_en']})",
                        "معرف المركب (ChEMBL ID)": str(rec.get("compound_id") or rec.get("molecule_chembl_id") or verified_comp.get("chembl_id") or "N/A"),
                        "اسم الهدف البروتيني": f"{TARGET_INFO['name_en']} ({rec.get('target_id', 'CHEMBL284')})",
                        "الكائن الحي": str(rec.get("target_organism", "Homo sapiens")),
                        "نوع التجربة (Assay Type)": str(rec.get("assay_type", "B")),
                        "نوع القياس (Standard Type)": str(rec.get("standard_type", "IC50")),
                        "وسيط القيمة (Median Activity)": round(med_val, 2),
                        "وحدة القياس": str(rec.get("standard_units", "nM")),
                        "قيمة pIC50": round(pic50_v, 2),
                        "عدد القياسات المجمعة (Measurement Count)": int(rec.get("measurement_count", 1)),
                        "معرف التجربة (Assay ID)": str(rec.get("assay_id") or rec.get("assay_chembl_id") or "N/A"),
                        "معرف المرجع (Document ID)": str(rec.get("document_id") or rec.get("document_chembl_id") or "N/A"),
                        "جميع مراجع القياسات المجمعة (Source IDs)": str(rec.get("source_ids") or rec.get("reference_citation") or f"{rec.get('assay_id')} ({rec.get('document_id')})"),
                    }
                )
            st.dataframe(pd.DataFrame(ev_rows), use_container_width=True, hide_index=True)
            st.caption(
                "ملاحظة منهجية: عند وجود أكثر من قياس متوافق للمركب نفسه، تعرض القيمة الوسيطة (Median Activity) مع ذكر عدد القياسات (Measurement Count) وكافة معرفات التجارب والمراجع الأصلية، ولا تعرض كأنها تجربة منفردة."
            )
            st.markdown("</div>", unsafe_allow_html=True)
        else:
            st.markdown(
                """
<div class="panel-card">
  <h4>حالة الأدلة التجريبية المباشرة للمركب المختار</h4>
  <p>
    لا يوجد سجل تجريبي منشور مطابق تماما لهذا المركب في مجموعة بيانات ChEMBL المنظفة لهدف DPP4 البشري.
    يمكنك الاطلاع أدناه على توزيع القياسات المنشورة في ChEMBL وأقرب المركبات المقاسة مخبريا في المرحلة الرابعة.
  </p>
</div>
""",
                unsafe_allow_html=True,
            )

        if clean_df is not None and not clean_df.empty:
            fig_pub, exp_pub = build_published_activity_distribution_chart(clean_df, highlight_val, highlight_lbl)
            st.plotly_chart(fig_pub, use_container_width=True)
            st.caption(exp_pub)

    def render_stage_4():
        st.markdown(
            """
<div class="panel-card">
  <h3>المرحلة الرابعة: تحدي «خمّن قبل أن ترى» والتنبؤ بالنشاط الحيوي باستخدام الذكاء الاصطناعي</h3>
  <p>
    تبدأ هذه المرحلة بتجربة تفاعلية تطلب منك صياغة توقعك العلمي بناء على خصائص المركب وأدلته المنشورة قبل مشاهدة نتائج نموذج Random Forest،
    يليها عرض تصنيف نطاق تطبيق النموذج (Applicability Domain) وأقرب مركبات التدريب.
  </p>
</div>
""",
            unsafe_allow_html=True,
        )

        current_guess = st.session_state["student_guesses_by_smiles"].get(canonical_smi, {"submitted": False})

        st.markdown('<div class="panel-card"><h4>تحدي الطالب التفاعلي: خمّن قبل أن ترى (Guess Before You See)</h4>', unsafe_allow_html=True)
        g_col1, g_col2, g_col3 = st.columns([1.2, 1.2, 0.8])
        act_placeholder = "-- اختر تخمينك لمستوى النشاط الحيوي تجاه DPP4 --"
        dom_placeholder = "-- اختر تخمينك لتصنيف نطاق تطبيق النموذج --"
        act_options = [
            act_placeholder,
            "نشاط تثبيطي مرتفع (ضمن العتبات التعليمية)",
            "نشاط تثبيطي متوسط (ضمن العتبات التعليمية)",
            "نشاط تثبيطي ضعيف أو غير نشط (ضمن العتبات التعليمية)",
        ]
        dom_options = [
            dom_placeholder,
            "داخل النطاق – In-Domain",
            "حدّي – Borderline",
            "خارج النطاق – Out-of-Domain",
        ]

        with g_col1:
            sel_act = st.selectbox(
                "1. بناء على المراحل السابقة، ما تخمينك لمستوى نشاط المركب؟",
                act_options,
                index=0,
                key=f"guess_act_{canonical_smi}",
            )
        with g_col2:
            sel_dom = st.selectbox(
                "2. ما تخمينك لمدى انتماء المركب إلى نطاق بيانات تدريب النموذج؟",
                dom_options,
                index=0,
                key=f"guess_dom_{canonical_smi}",
            )
        with g_col3:
            st.write("")
            st.write("")
            if st.button("إرسال التخمين وإظهار النتيجة", key=f"btn_guess_{canonical_smi}"):
                if sel_act == act_placeholder or sel_dom == dom_placeholder:
                    st.warning("يرجى اختيار تخمينك لكلا السؤالين أولا قبل الإرسال (لا توجد قيم افتراضية مسبقة).")
                else:
                    ad_actual = pred_data["applicability_domain"]["status_ar"]
                    act_actual = pred_data["category_ar"]
                    act_ok = sel_act == act_actual
                    dom_ok = sel_dom == ad_actual
                    comp_exp = (
                        f"{'أصبت في تقدير فئة النشاط؛ ' if act_ok else 'اختلف تخمينك للنشاط عن تقدير النموذج؛ '}"
                        f"{'وأصبت في تحديد نطاق التطبيق.' if dom_ok else 'واختلف تقديرك لنطاق التطبيق (حيث يعتمد النطاق على تشابه Tanimoto مع مركبات التدريب وفحص البنية الكيميائية).'}"
                    )
                    st.session_state["student_guesses_by_smiles"][canonical_smi] = {
                        "submitted": True,
                        "activity_guess": sel_act,
                        "domain_guess": sel_dom,
                        "activity_matched": act_ok,
                        "domain_matched": dom_ok,
                        "comparison_explanation": comp_exp,
                    }
                    current_guess = st.session_state["student_guesses_by_smiles"][canonical_smi]

        if current_guess.get("submitted"):
            st.success(
                f"تم تسجيل تخمينك لهذا المركب ({verified_comp['name_ar']}): "
                f"النشاط المختار: [{current_guess['activity_guess']}] | النطاق المختار: [{current_guess['domain_guess']}]. "
                f"النتيجة التعليمية: {current_guess['comparison_explanation']}"
            )
        else:
            st.info("أدخل تخمينك أعلاه واضغط «إرسال التخمين وإظهار النتيجة» لتفعيل مقارنة إجابتك بنتيجة النموذج.")
        st.markdown("</div>", unsafe_allow_html=True)

        ad_info = pred_data["applicability_domain"]
        st.markdown(
            f"""
<div class="panel-card">
  <h4>تصنيف نطاق تطبيق النموذج (Prediction Applicability Domain)</h4>
  <p>
    <strong>حالة النطاق للمركب الحالي:</strong>
    <span class="status-pill {ad_info['pill_class']}">{ad_info['status_ar']}</span>
    &nbsp;|&nbsp;
    <strong>أقصى تشابه Tanimoto مع مركبات التدريب:</strong> <code class="ltr-inline">{pred_data['max_tanimoto']:.3f}</code>
    &nbsp;|&nbsp;
    <strong>متوسط أقرب 5 مركبات تدريب:</strong> <code class="ltr-inline">{pred_data['mean_top5_tanimoto']:.3f}</code>
  </p>
  <p>{ad_info['summary_ar']}</p>
  <ul>
    {"".join(f"<li>{r}</li>" for r in ad_info['reasons_ar'])}
  </ul>
  <p style="font-size:0.84rem;opacity:0.9;"><strong>المبرر العلمي للحدود:</strong> {ad_info['thresholds_justification_ar']}</p>
  <p style="font-size:0.84rem;color:#78350f;"><strong>إقرار منهجي:</strong> {ad_info['scientific_disclaimer_ar']}</p>
</div>
""",
            unsafe_allow_html=True,
        )

        show_ood_raw = False
        if not ad_info["can_show_reliable_prediction"]:
            st.error(
                "تنبيه علمي: هذا المركب يقع خارج نطاق تطبيق النموذج (Out-of-Domain). "
                "لذلك لا يقدم النظام تنبؤا عدديا معتمدا له، منعا للتضليل العلمي."
            )
            show_ood_raw = st.checkbox(
                "إظهار القيمة الحاسوبية الخام لأغراض المقارنة والتعلم فقط (مع الإقرار بعدم موثوقيتها)",
                value=False,
                key=f"ood_toggle_{canonical_smi}",
            )

        if ad_info["can_show_reliable_prediction"] or show_ood_raw:
            c_p1, c_p2 = st.columns([1.1, 1.3])
            with c_p1:
                fig_g, exp_g = build_pic50_gauge(pred_data["pred_pic50"])
                st.plotly_chart(fig_g, use_container_width=True)
                st.caption(exp_g)
            with c_p2:
                warn_banner = (
                    "<p><span class='status-pill pill-low'>تحذير: قيمة خام خارج نطاق التدريب (للتوضيح التعليمي فقط)</span></p>"
                    if not ad_info["can_show_reliable_prediction"]
                    else (
                        "<p><span class='status-pill pill-mid'>تنبيه: المركب في المنطقة الحدية (Borderline) - يفسر التنبؤ بحذر</span></p>"
                        if ad_info["status_code"] == "Borderline"
                        else ""
                    )
                )
                st.markdown(
                    f"""
<div class="panel-card">
  <h4>مخرجات نموذج الذكاء الاصطناعي (Random Forest Regressor)</h4>
  {warn_banner}
  <p><strong>نوع القياس المتنبأ به:</strong> {pred_data['measurement_type']} (خاص بـ IC50 فقط ولا يخلط مع Ki)</p>
  <p><strong>وحدة القياس:</strong> {pred_data['units']}</p>
  <p><strong>قيمة التنبؤ اللوغاريتمي (Predicted pIC50):</strong> <code class="ltr-inline" style="font-size:1.2rem;font-weight:800;">{pred_data['pred_pic50']:.2f}</code></p>
  <p><strong>القيمة المكافئة لتركيز التثبيط النصفي (Predicted IC50):</strong> <code class="ltr-inline" style="font-size:1.15rem;font-weight:700;">{pred_data['pred_ic50_nm']:,.2f} nM</code></p>
  <p><strong>تفسير النتيجة:</strong> <span class="status-pill {pred_data['pill_class']}">{pred_data['category_ar']}</span></p>
  <p>{pred_data['explanation_ar']}</p>
</div>
""",
                    unsafe_allow_html=True,
                )

        st.markdown('<div class="panel-card"><h4>مقاييس التحقق العالمية للنموذج على مجموعة الاختبار (Bemis-Murcko Scaffold Split)</h4>', unsafe_allow_html=True)
        mc1, mc2, mc3, mc4 = st.columns(4)
        mc1.metric("متوسط الخطأ المطلق (MAE)", f"{metrics['mae_pic50']:.3f} pIC50")
        mc2.metric("جذر متوسط مربع الخطأ (RMSE)", f"{metrics['rmse_pic50']:.3f} pIC50")
        mc3.metric("معامل التحديد العددي (R²)", f"{metrics['r2_pic50']:.3f}")
        mc4.metric("أداء التمييز التصنيفي (ROC-AUC)", f"{metrics['roc_auc']:.3f}")
        st.caption(metrics.get("roc_auc_note_ar", ""))
        st.markdown("</div>", unsafe_allow_html=True)

        fig_sim, exp_sim = build_similarity_chart(pred_data["all_train_similarities"], pred_data["max_tanimoto"])
        st.plotly_chart(fig_sim, use_container_width=True)
        st.caption(exp_sim)

        nearest_df = pred_data["nearest_train_compounds"]
        fig_near, exp_near = build_nearest_compounds_chart(nearest_df)
        st.plotly_chart(fig_near, use_container_width=True)
        st.caption(exp_near)

        st.markdown("<h4>جدول أقرب 5 مركبات في بيانات التدريب فقط (مع المعرفات والمصادر)</h4>", unsafe_allow_html=True)
        near_table = pd.DataFrame(
            {
                "معرف المركب (ChEMBL ID)": nearest_df["molecule_chembl_id"].astype(str),
                "معامل التشابه (Tanimoto)": nearest_df["tanimoto"].round(3),
                "نوع القياس": nearest_df["standard_type"].astype(str),
                "وسيط النشاط (IC50 nM)": nearest_df["standard_value"].round(2),
                "قيمة pIC50": nearest_df["pIC50"].round(2),
                "عدد القياسات": nearest_df.get("measurement_count", pd.Series(1, index=nearest_df.index)).astype(int),
                "معرف التجربة (Assay ID)": nearest_df["assay_chembl_id"].astype(str),
                "معرف المرجع (Document ID)": nearest_df["document_chembl_id"].astype(str),
                "مراجع القياسات (Source IDs)": nearest_df.get("source_ids", nearest_df["assay_chembl_id"]).astype(str),
            }
        )
        st.dataframe(near_table, use_container_width=True, hide_index=True)

    def render_stage_5():
        st.markdown(
            """
<div class="panel-card">
  <h3>المرحلة الخامسة: التقييم النهائي والتقرير التعليمي (Final Evaluation &amp; Educational Report)</h3>
  <p>
    يفصل هذا الملخص النهائي بصريا وعلميا بين ثلاث ركائز مستقلة: الخصائص الكيميائية المحسوبة، الأدلة التجريبية المنشورة، وتوقعات نموذج الذكاء الاصطناعي.
  </p>
</div>
""",
            unsafe_allow_html=True,
        )

        ad_info = pred_data["applicability_domain"]
        c_sep1, c_sep2, c_sep3 = st.columns(3)
        with c_sep1:
            st.markdown(
                f"""
<div class="sep-col-chem">
  <h4>1. الخصائص الكيميائية المحسوبة (RDKit)</h4>
  <p><em>حسابات بنيوية حتمية مشتقة من الصيغة الجزيئية:</em></p>
  <ul>
    <li><strong>الصيغة الجزيئية:</strong> <code class="ltr-inline">{desc_data['molecular_formula']}</code></li>
    <li><strong>الوزن الجزيئي:</strong> <code class="ltr-inline">{desc_data['molecular_weight']} Da</code></li>
    <li><strong>معامل LogP:</strong> <code class="ltr-inline">{desc_data['logp']}</code></li>
    <li><strong>مانحات / مستقبلات H:</strong> <code class="ltr-inline">HBD={desc_data['hbd']} | HBA={desc_data['hba']}</code></li>
    <li><strong>الروابط الدوارة و TPSA:</strong> <code class="ltr-inline">Rot={desc_data['rotatable_bonds']} | TPSA={desc_data['tpsa']} Å²</code></li>
    <li><strong>قاعدة Lipinski:</strong> {desc_data['lipinski']['summary_ar']}</li>
    <li><strong>فحص PAINS:</strong> {desc_data['pains']['status_ar']}</li>
  </ul>
</div>
""",
                unsafe_allow_html=True,
            )
        with c_sep2:
            if published_records:
                rec0 = published_records[0]
                med_v = float(rec0.get("median_activity", rec0.get("standard_value", 0.0)))
                pic_v = float(rec0.get("pIC50", 9.0 - np.log10(med_v) if med_v > 0 else 0.0))
                exp_html = f"""
  <p><em>قياسات مخبرية حقيقية موثقة في قاعدة بيانات ChEMBL:</em></p>
  <ul>
    <li><strong>معرف المركب:</strong> <code class="ltr-inline">{rec0.get('compound_id') or rec0.get('molecule_chembl_id')}</code></li>
    <li><strong>الهدف البشري:</strong> <code class="ltr-inline">{rec0.get('target_id', 'CHEMBL284')} ({rec0.get('target_organism', 'Homo sapiens')})</code></li>
    <li><strong>نوع القياس:</strong> <code class="ltr-inline">{rec0.get('standard_type', 'IC50')}</code></li>
    <li><strong>الوسيط المخبري:</strong> <code class="ltr-inline">{med_v:.2f} nM (pIC50 = {pic_v:.2f})</code></li>
    <li><strong>عدد القياسات:</strong> <code class="ltr-inline">{rec0.get('measurement_count', 1)}</code></li>
    <li><strong>المراجع (Assay / Doc):</strong> <code class="ltr-inline">{rec0.get('source_ids') or (str(rec0.get('assay_id')) + ' / ' + str(rec0.get('document_id')))}</code></li>
  </ul>
"""
            else:
                exp_html = """
  <p><em>قياسات مخبرية حقيقية موثقة في قاعدة بيانات ChEMBL:</em></p>
  <p>لا تتوفر تجربة مخبرية منشورة مطابقة تماما لهذا المركب في مجموعة بيانات DPP4 المحلية المنظفة.</p>
"""
            st.markdown(
                f"""
<div class="sep-col-exp">
  <h4>2. الأدلة المنشورة من التجارب (ChEMBL)</h4>
  {exp_html}
</div>
""",
                unsafe_allow_html=True,
            )
        with c_sep3:
            pred_val_html = (
                f"<li><strong>التنبؤ اللوغاريتمي (pIC50):</strong> <code class='ltr-inline'>{pred_data['pred_pic50']:.2f}</code></li>"
                f"<li><strong>التركيز المكافئ (IC50):</strong> <code class='ltr-inline'>{pred_data['pred_ic50_nm']:,.2f} nM</code></li>"
                f"<li><strong>تفسير النشاط:</strong> {pred_data['category_ar']}</li>"
                if ad_info["can_show_reliable_prediction"]
                else "<li><strong>التنبؤ العددي:</strong> محجوب كقيمة معتمدة لوقوع المركب خارج نطاق تطبيق النموذج (Out-of-Domain).</li>"
            )
            st.markdown(
                f"""
<div class="sep-col-ai">
  <h4>3. توقعات الذكاء الاصطناعي ونطاق التطبيق</h4>
  <p><em>تقديرات حاسوبية لنموذج Random Forest (ليست قياسا مخبريا):</em></p>
  <ul>
    <li><strong>نطاق التطبيق (AD):</strong> <span class="status-pill {ad_info['pill_class']}">{ad_info['status_ar']}</span></li>
    <li><strong>أقصى تشابه Tanimoto:</strong> <code class="ltr-inline">{pred_data['max_tanimoto']:.3f}</code></li>
    {pred_val_html}
    <li><strong>أداء النموذج (Scaffold Split):</strong> <code class="ltr-inline">MAE={metrics['mae_pic50']:.2f}, RMSE={metrics['rmse_pic50']:.2f}, R²={metrics['r2_pic50']:.2f}</code></li>
  </ul>
</div>
""",
                unsafe_allow_html=True,
            )

        student_guess_state = st.session_state["student_guesses_by_smiles"].get(canonical_smi)
        report_markdown = generate_educational_report(
            compound_info=verified_comp,
            descriptors=desc_data,
            published_records=published_records,
            prediction=pred_data,
            metrics=metrics,
            student_guess=student_guess_state,
        )

        with st.expander("معاينة النص الكامل للتقرير التعليمي الموثق", expanded=False):
            st.markdown(report_markdown)

        st.download_button(
            label="تحميل التقرير التعليمي الكامل للمركب (Markdown Report)",
            data=report_markdown.encode("utf-8"),
            file_name=f"MolecuLens_Report_{verified_comp.get('id', 'compound')}.md",
            mime="text/markdown",
        )

    if display_mode == "عرض الرحلة الكاملة المتتابعة (المراحل الخمس في شاشة موحدة)":
        render_stage_1()
        render_stage_2()
        render_stage_3()
        render_stage_4()
        render_stage_5()
    else:
        tab1, tab2, tab3, tab4, tab5 = st.tabs(
            [
                "المرحلة 1: الهدف البروتيني",
                "المرحلة 2: الخصائص الكيميائية",
                "المرحلة 3: الأدلة المنشورة",
                "المرحلة 4: التوقع وتنبؤ AI",
                "المرحلة 5: التقييم والتقرير",
            ]
        )
        with tab1:
            render_stage_1()
        with tab2:
            render_stage_2()
        with tab3:
            render_stage_3()
        with tab4:
            render_stage_4()
        with tab5:
            render_stage_5()

elif main_section == "2. مكتبة المركبات التجريبية لـ DPP4 (Demo Library)":
    st.markdown(
        """
<div class="panel-card">
  <h3>مكتبة المركبات التجريبية المرجعية المثبطة لإنزيم DPP4 البشري (DPP4 Demo Compounds)</h3>
  <p>
    تضم هذه المكتبة خمسة مثبطات قياسية معروفة لإنزيم DPP4 البشري، تم التحقق من سجلاتها المرجعية في ChEMBL و PubChem
    وارتباط قياساتها بالهدف البشري <code class="ltr-inline">CHEMBL284</code> ونوع التجربة والوحدة والمرجع العلمي.
  </p>
</div>
""",
        unsafe_allow_html=True,
    )

    for idx_d, demo in enumerate(demo_compounds):
        d_mol, d_canon, _ = validate_smiles(demo["canonical_smiles"])
        with st.container():
            st.markdown('<div class="panel-card">', unsafe_allow_html=True)
            cd1, cd2 = st.columns([1.0, 2.2])
            with cd1:
                if d_mol is not None:
                    st.image(render_mol_png(d_mol, size=(340, 220)), caption=f"{demo['name_ar']} ({demo['name_en']})")
            with cd2:
                st.markdown(
                    f"""
<h4>{demo['name_ar']} | <span class="ltr-inline">{demo['name_en']}</span></h4>
<p><strong>الأسماء التجارية والرموز:</strong> {", ".join(demo.get("brand_names", []))}</p>
<p><strong>الصيغة الكيميائية:</strong> <code class="ltr-inline">{demo['molecular_formula']}</code> | <strong>PubChem CID:</strong> <code class="ltr-inline">{demo['pubchem_cid']}</code> | <strong>ChEMBL ID:</strong> <code class="ltr-inline">{demo['chembl_id']}</code></p>
<p><strong>الهدف البروتيني:</strong> <code class="ltr-inline">{demo['target_id']}</code> ({demo['target_name']} - <code class="ltr-inline">{demo['target_organism']}</code>)</p>
<p><strong>بيانات النشاط الحيوي الموثقة:</strong> نوع التجربة: <code class="ltr-inline">{demo['assay_type']}</code> | نوع القياس: <code class="ltr-inline">{demo['standard_type']}</code> | الوسيط: <code class="ltr-inline">{demo['median_activity']} {demo['standard_units']}</code> (<code class="ltr-inline">pIC50 = {demo['pIC50']}</code>) | عدد القياسات: <code class="ltr-inline">{demo['measurement_count']}</code></p>
<p><strong>مصدر القياس المرجعي:</strong> Assay ID: <code class="ltr-inline">{demo['assay_id']}</code> | Document ID: <code class="ltr-inline">{demo['document_id']}</code> (<span class="ltr-inline">{demo['reference_citation']}</span>)</p>
<p><strong>ملاحظة علمية:</strong> {demo['notes_ar']}</p>
""",
                    unsafe_allow_html=True,
                )
                if st.button(f"نسخ اسم {demo['name_ar']} إلى مربع البحث لتحليله في المختبر", key=f"load_demo_{idx_d}"):
                    st.session_state[SEARCH_QUERY_STATE_KEY] = demo["name_en"]
                    st.session_state["selected_override_compound"] = None
                    st.success(f"تم تحديد {demo['name_ar']} ({demo['name_en']}). انتقل إلى تبويب «1. البحث الذكي والمختبر الافتراضي» لاستعراض المراحل الخمس.")
            st.markdown("</div>", unsafe_allow_html=True)

else:
    st.markdown(
        """
<div class="panel-card">
  <h3>الدرس التعليمي: مقارنة Random Split و Scaffold Split وتوثيق تنظيف بيانات ChEMBL</h3>
  <p>
    يوضح هذا القسم كيف تم تنظيف بيانات ChEMBL بصرامة علمية، ولماذا نعتمد تقسيم الهياكل الأساسية (Bemis-Murcko Scaffold Split)
    بدلا من التقسيم العشوائي عند اختبار قدرة نموذج الذكاء الاصطناعي على التعميم على مركبات جديدة.
  </p>
</div>
""",
        unsafe_allow_html=True,
    )

    st.markdown('<div class="panel-card"><h4>أولا: مقارنة منهجية التقسيم العشوائي (Random Split) والتقسيم الهيكلي (Scaffold Split)</h4>', unsafe_allow_html=True)
    rand_m = metrics["random_split_comparison"]
    comp_table = pd.DataFrame(
        [
            {
                "المقياس": "MAE (متوسط الخطأ المطلق في pIC50 - الأقل أفضل)",
                "Random Split (تقسيم عشوائي 80/20)": rand_m["mae_pic50"],
                "Scaffold Split (تقسيم الهياكل المعتمد)": metrics["mae_pic50"],
            },
            {
                "المقياس": "RMSE (جذر متوسط مربع الخطأ في pIC50 - الأقل أفضل)",
                "Random Split (تقسيم عشوائي 80/20)": rand_m["rmse_pic50"],
                "Scaffold Split (تقسيم الهياكل المعتمد)": metrics["rmse_pic50"],
            },
            {
                "المقياس": "R² (معامل التحديد للتنبؤ العددي - الأعلى أفضل)",
                "Random Split (تقسيم عشوائي 80/20)": rand_m["r2_pic50"],
                "Scaffold Split (تقسيم الهياكل المعتمد)": metrics["r2_pic50"],
            },
            {
                "المقياس": "ROC-AUC (للتمييز التصنيفي عند عتبة pIC50 >= 6.0)",
                "Random Split (تقسيم عشوائي 80/20)": rand_m["roc_auc"],
                "Scaffold Split (تقسيم الهياكل المعتمد)": metrics["roc_auc"],
            },
            {
                "المقياس": "عدد الهياكل الأساسية المشتركة المتسربة بين التدريب والاختبار",
                "Random Split (تقسيم عشوائي 80/20)": rand_m["shared_scaffolds_leaked"],
                "Scaffold Split (تقسيم الهياكل المعتمد)": metrics["scaffold_leakage_count"],
            },
        ]
    )
    st.dataframe(comp_table, use_container_width=True, hide_index=True)

    fig_split, exp_split = build_split_comparison_chart(metrics)
    st.plotly_chart(fig_split, use_container_width=True)
    st.caption(exp_split)

    st.markdown(
        f"""
<h5>الشرح التعليمي للمقارنة (بذرة عشوائية ثابتة Seed = {metrics['random_seed']} ونفس إعدادات النموذج):</h5>
<ul>
  <li><strong>لماذا يبدو Random Split أسهل؟</strong> لأن التقسيم العشوائي يوزع المركبات المتشابهة التي تشترك في نفس النواة الكيميائية (حيث تسرب <code class="ltr-inline">{rand_m['shared_scaffolds_leaked']}</code> هيكلا أساسيا مشتركا) بين التدريب والاختبار، مما يجعل الاختبار أقرب إلى استيفاء داخلي ويعطي تقديرا متفائلا (<code class="ltr-inline">R² = {rand_m['r2_pic50']:.3f}, MAE = {rand_m['mae_pic50']:.3f}</code>).</li>
  <li><strong>لماذا يعتبر Scaffold Split اختبارا أصعب وأكثر صدقا علميا؟</strong> لأنه يعزل كل هيكل Bemis-Murcko بالكامل بحيث لا يظهر أي هيكل أساسي مشترك بين التدريب والاختبار (<code class="ltr-inline">تسرب الهياكل = {metrics['scaffold_leakage_count']}</code>)، مما يحاكي الواقع الفعلي عندما يختبر الباحث عائلة كيميائية جديدة تماما.</li>
  <li><strong>ماذا تعني هذه النتائج للطالب؟</strong> تعلمنا ألا ننخدع بارتفاع <code class="ltr-inline">R²</code> في التقسيم العشوائي، وأن نعتمد مقاييس الخطأ الواقعية على الهياكل الجديدة مع فحص نطاق التطبيق (Applicability Domain) لكل مركب.</li>
</ul>
</div>
""",
        unsafe_allow_html=True,
    )

    if test_df is not None and not test_df.empty:
        st.markdown('<div class="panel-card"><h4>ثانيا: تحليل تنبؤات مجموعة الاختبار (Scaffold Split Test Set)</h4>', unsafe_allow_html=True)
        fig_sc, exp_sc = build_test_scatter_chart(test_df)
        st.plotly_chart(fig_sc, use_container_width=True)
        st.caption(exp_sc)
        st.markdown("</div>", unsafe_allow_html=True)

    if clean_report is not None:
        st.markdown('<div class="panel-card"><h4>ثالثا: تقرير تنظيف بيانات ChEMBL المفصل (Advanced ChEMBL Data Cleaning Report)</h4>', unsafe_allow_html=True)
        cr1, cr2, cr3, cr4, cr5 = st.columns(5)
        cr1.metric("السجلات قبل التنظيف", clean_report["initial_records"])
        cr2.metric("السجلات الصالحة بعد التنظيف", clean_report["records_after_cleaning"])
        cr3.metric("المركبات الفريدة النهائية", clean_report["unique_compounds_count"])
        cr4.metric("عدد التجارب الفريدة (Assays)", clean_report["unique_assays_count"])
        cr5.metric("عدد المراجع العلمية (Documents)", clean_report["unique_documents_count"])

        excl = clean_report.get("exclusion_reasons", {})
        excl_df = pd.DataFrame(
            [
                {"سبب الاستبعاد العلمي": "قيم مفقودة أو غير موجبة", "عدد السجلات المستبعدة": excl.get("missing_or_nonpositive_value", 0)},
                {"سبب الاستبعاد العلمي": "هدف بروتيني مختلف أو كائن حي غير بشري", "عدد السجلات المستبعدة": excl.get("wrong_target_or_non_human_organism", 0)},
                {"سبب الاستبعاد العلمي": "نوع قياس غير مطابق لـ IC50", "عدد السجلات المستبعدة": excl.get("non_ic50_measurement_type", 0)},
                {"سبب الاستبعاد العلمي": "علاقات قياس مقيدة أو حدية (> أو < أو >>)", "عدد السجلات المستبعدة": excl.get("censored_or_inequality_relation", 0)},
                {"سبب الاستبعاد العلمي": "ملاحظات عدم صلاحية البيانات (Outside typical range)", "عدد السجلات المستبعدة": excl.get("data_validity_warning_flagged", 0)},
                {"سبب الاستبعاد العلمي": "اختلاف نوع التجربة (غير تجارب الربط المباشر B)", "عدد السجلات المستبعدة": excl.get("non_binding_assay_type", 0)},
                {"سبب الاستبعاد العلمي": "تكرارات متضاربة للمركب نفسه (فرق يتجاوز 1.0 pIC50)", "عدد السجلات المستبعدة": excl.get("discordant_replicate_measurements", 0)},
            ]
        )
        st.dataframe(excl_df, use_container_width=True, hide_index=True)
        st.caption(
            f"تم توحيد {clean_report['converted_non_nm_units']} سجلا من وحدة ug.mL-1 إلى nM باستخدام الوزن الجزيئي، "
            f"وتجميع القياسات المتوافقة بحساب الوسيط (Median Activity) مع حفظ عدد القياسات وجميع معرفات التجارب والمراجع في حقل source_ids."
        )
        st.markdown("</div>", unsafe_allow_html=True)

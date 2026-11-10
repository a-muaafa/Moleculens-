from typing import Any
import pandas as pd
import streamlit as st
from config.settings import (
    ENTITY_BIOMOLECULE_PROTEIN,
    ENTITY_CHEMICAL_COMPOUND,
    ENTITY_CHEMICAL_ELEMENT,
    ENTITY_DRUG_SUBSTANCE,
)
from core.orchestrator import MolecuLensCore
from frontend.structure_viewer import render_structure_viewer
from search.normalizer import detect_query_input_type

EXPLORE_QUERY_KEY = "ml_explore_query"


def _set_explore_query(q_val: str) -> None:
    st.session_state[EXPLORE_QUERY_KEY] = q_val


def _badge_class_for_type(entity_type: str) -> str:
    if entity_type == ENTITY_CHEMICAL_ELEMENT:
        return "ml-badge ml-badge-element"
    if entity_type == ENTITY_DRUG_SUBSTANCE:
        return "ml-badge ml-badge-drug"
    if entity_type == ENTITY_BIOMOLECULE_PROTEIN:
        return "ml-badge ml-badge-protein"
    return "ml-badge ml-badge-compound"


def _render_suggestions(query_text: str) -> None:
    suggestions = MolecuLensCore.fetch_suggestions(query_text, limit=6)
    if not suggestions:
        return
    st.markdown("<strong>اقتراحات فورية (اختر لتعبئة خانة البحث):</strong>", unsafe_allow_html=True)
    for idx, s in enumerate(suggestions):
        c_text, c_btn = st.columns([5, 1])
        with c_text:
            st.markdown(
                f"""
<div class="ml-suggestion">
  <span class="ml-sug-name">{s['highlighted_name_ar']}</span>
  <span class="ml-sug-sep">|</span>
  <span class="ml-sug-name ltr-inline">{s['highlighted_name_en']}</span>
  <span class="ml-sug-formula ltr-inline">{s['highlighted_formula']}</span>
  <span class="ml-sug-class">{s['subclass_label_ar']}</span>
</div>
""",
                unsafe_allow_html=True,
            )
        with c_btn:
            st.button("اختيار", key=f"exp_sug_{idx}_{s['key']}", on_click=_set_explore_query, args=(s["symbol_or_formula"] if s["entity_type"] != ENTITY_BIOMOLECULE_PROTEIN else s["name_en"],))


def render_explore_section(on_open_reaction_lab) -> None:
    if EXPLORE_QUERY_KEY not in st.session_state:
        st.session_state[EXPLORE_QUERY_KEY] = "CH3COOH"

    st.markdown(
        """
<div class="ml-card">
  <h3>Explore &amp; Understand — استكشاف المادة وفهم تركيبها وخصائصها</h3>
  <p>
    اكتب في خانة البحث الواحدة أي مدخل: الاسم العربي، الاسم الإنجليزي، الاسم العلمي (IUPAC)، الصيغة الكيميائية (H2O أو H₂O أو CH3COOH)،
    الاختصار، الاسم التجاري، أو المرادفات. يتعرف النظام على نوع المدخل تلقائيا.
  </p>
</div>
""",
        unsafe_allow_html=True,
    )

    query_text = st.text_input(
        "خانة البحث العلمي الموحدة:",
        key=EXPLORE_QUERY_KEY,
        placeholder="مثال: حمض الخليك، Acetic acid، CH3COOH، Aspirin، Fe",
    )
    detected = detect_query_input_type(query_text)
    if detected["type"] != "empty":
        st.caption(f"نوع المدخل المكتشف تلقائيا: {detected['label_ar']}")

    _render_suggestions(query_text)

    res = MolecuLensCore.explore_entity(query_text, allow_network=True)

    if not res["found"]:
        st.warning(res["message"])
        if res["did_you_mean"]:
            st.markdown("<strong>مقترحات يمكنك اختيارها مباشرة:</strong>", unsafe_allow_html=True)
            d_cols = st.columns(min(len(res["did_you_mean"]), 5))
            for idx_d, d_txt in enumerate(res["did_you_mean"][:5]):
                clean_q = d_txt.split("(")[-1].replace(")", "").strip() if "(" in d_txt else d_txt
                with d_cols[idx_d]:
                    st.button(d_txt, key=f"exp_dym_{idx_d}", on_click=_set_explore_query, args=(clean_q,))
        return

    if res["is_ambiguous"]:
        st.warning(res["ambiguity_message"])
        cand_list = res["candidates"]
        cand_labels = [
            f"{idx_c + 1}. {c['name_ar']} | {c['name_en']} ({c['symbol_or_formula']}) — [{c['entity_type_label_ar']}]"
            for idx_c, c in enumerate(cand_list)
        ]
        chosen_lbl = st.radio(
            "توجد أكثر من مادة تشترك في هذا المدخل؛ اختر المادة المقصودة لعرض ملفها العلمي:",
            options=cand_labels,
            key=f"exp_ambig_{query_text}",
        )
        chosen_cand = cand_list[cand_labels.index(chosen_lbl)]
        res = MolecuLensCore.explore_entity(query_text, selected_candidate=chosen_cand)

    entity: dict[str, Any] = res["entity"]
    profile: dict[str, Any] = res["profile"]
    entity_type = entity["entity_type"]
    badge_cls = _badge_class_for_type(entity_type)

    st.markdown(
        f"""
<div class="ml-card">
  <div>
    <span class="{badge_cls}">{entity['entity_type_label_ar']} | {entity['entity_type_label_en']}</span>
    <span class="ml-badge" style="background:#475569;">{entity.get('subclass_label_ar', '')}</span>
    <span class="ml-badge" style="background:#334155;">{entity.get('match_category', 'Exact Match')}</span>
  </div>
  <h3 style="margin-top:0.65rem;margin-bottom:0.25rem;">
    {entity['name_ar']} — <span class="ltr-inline">{entity['name_en']} ({entity['symbol_or_formula']})</span>
  </h3>
  <p style="margin-bottom:0;color:#64748b;font-size:0.92rem;">
    <strong>سبب المطابقة:</strong> {entity.get('match_reason_ar', '')} |
    <strong>المصادر العلمية الموحدة:</strong> <span class="ltr-inline">{", ".join(profile.get('sources', []))}</span>
  </p>
</div>
""",
        unsafe_allow_html=True,
    )

    if entity_type == ENTITY_CHEMICAL_ELEMENT:
        _render_element_profile_ui(profile)
    elif entity_type == ENTITY_DRUG_SUBSTANCE:
        _render_drug_profile_ui(profile)
    elif entity_type == ENTITY_BIOMOLECULE_PROTEIN:
        _render_protein_profile_ui(profile)
    else:
        _render_compound_profile_ui(profile)

    if entity_type != ENTITY_CHEMICAL_ELEMENT:
        render_structure_viewer(entity)

    _render_explanation_and_applications_ui(profile)

    st.markdown('<div class="ml-card">', unsafe_allow_html=True)
    st.markdown(
        f"<h4>الانتقال إلى مختبر التفاعلات (Explore Reactions / Open Reaction Lab)</h4>"
        f"<p>هل ترغب في اختبار تفاعلات <strong>{entity['name_ar']} ({entity['symbol_or_formula']})</strong> وحساب النواتج والكميات المولية؟</p>",
        unsafe_allow_html=True,
    )
    if st.button(
        f"Explore Reactions / Open Reaction Lab — نقل ({entity['symbol_or_formula']}) إلى مختبر التفاعلات",
        key="btn_pass_to_reaction_lab",
    ):
        on_open_reaction_lab(profile["reaction_lab_entity"])
    st.markdown("</div>", unsafe_allow_html=True)


def _render_element_profile_ui(profile: dict[str, Any]) -> None:
    basic = profile["basic_info"]
    atomic = profile["atomic_structure"]
    props = profile["properties"]

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("الاسم والرمز الكيميائي", f"{basic['name_ar']} ({basic['symbol']})")
    m2.metric("العدد الذري (Atomic Number)", str(basic["atomic_number"]))
    m3.metric("الكتلة الذرية (Atomic Mass)", str(basic["atomic_mass"]))
    m4.metric("الحالة الفيزيائية", props["physical_state"].split(" ")[0])

    c1, c2 = st.columns(2)
    with c1:
        st.markdown(
            f"""
<div class="ml-card">
  <h4>1. المعلومات الأساسية والتركيب الذري والتوزيع الإلكتروني</h4>
  <p><strong>الاسم الإنجليزي:</strong> <span class="ltr-inline">{basic['name_en']}</span> | <strong>الاسم العربي:</strong> {basic['name_ar']}</p>
  <p><strong>الرمز الكيميائي:</strong> <code class="ltr-inline">{basic['symbol']}</code> | <strong>العدد الذري (Z):</strong> <code class="ltr-inline">{basic['atomic_number']}</code></p>
  <p><strong>عدد البروتونات والإلكترونات:</strong> <code class="ltr-inline">{atomic['protons_count']} بروتون / {atomic['electrons_count']} إلكترون</code></p>
  <p><strong>التوزيع الإلكتروني (Electron Configuration):</strong> <code class="ltr-inline">{atomic['electron_configuration']}</code></p>
  <p><strong>توزيع الأغلفة الإلكترونية:</strong> <code class="ltr-inline">{atomic.get('electron_shells') or 'غير محدد'}</code></p>
  <p><strong>الموقع في الجدول الدوري:</strong> {atomic.get('group', '')} — {atomic.get('period', '')} ({atomic.get('block', '')})</p>
  <p><strong>حالات الأكسدة الشائعة:</strong> <span class="ltr-inline">{atomic.get('oxidation_states') or 'غير محدد'}</span></p>
  <p><strong>البنية البلورية الذرية:</strong> {atomic.get('crystal_structure') or 'شبكة بلورية عنصرية'}</p>
</div>
""",
            unsafe_allow_html=True,
        )
    with c2:
        st.markdown(
            f"""
<div class="ml-card">
  <h4>2. الخصائص الفيزيائية والكيميائية للعنصر</h4>
  <p><strong>الحالة الفيزيائية:</strong> {props['physical_state']}</p>
  <p><strong>الكثافة:</strong> <span class="ltr-inline">{props['density']}</span></p>
  <p><strong>درجة الانصهار:</strong> <span class="ltr-inline">{props['melting_point']}</span></p>
  <p><strong>درجة الغليان:</strong> <span class="ltr-inline">{props['boiling_point']}</span></p>
  <p><strong>السلوك والتفاعلية الكيميائية:</strong> {props['chemical_behavior']}</p>
</div>
""",
            unsafe_allow_html=True,
        )


def _render_compound_profile_ui(profile: dict[str, Any]) -> None:
    basic = profile["basic_info"]
    comp = profile["composition"]
    struct = profile["structure"]
    props = profile["properties"]

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("الاسم", f"{basic['name_ar']} ({basic['name_en']})")
    m2.metric("الصيغة الكيميائية", basic["formula"])
    m3.metric("الكتلة المولية", basic["molar_mass"])
    m4.metric("إجمالي عدد الذرات في الجزيء", str(comp.get("total_atoms", 0)))

    c1, c2 = st.columns(2)
    with c1:
        st.markdown(
            f"""
<div class="ml-card">
  <h4>1. المعلومات الأساسية والتركيب والمكونات</h4>
  <p><strong>الاسم الإنجليزي:</strong> <span class="ltr-inline">{basic['name_en']}</span> | <strong>الاسم العربي:</strong> {basic['name_ar']}</p>
  <p><strong>الاسم النظامي (IUPAC):</strong> <code class="ltr-inline">{basic.get('iupac_name') or 'N/A'}</code></p>
  <p><strong>الصيغة الجزيئية:</strong> <code class="ltr-inline">{basic['formula']}</code> | <strong>الكتلة المولية:</strong> <code class="ltr-inline">{basic['molar_mass']}</code></p>
  <p><strong>ملخص التركيب الذري:</strong> {comp.get('summary_ar', '')}</p>
</div>
""",
            unsafe_allow_html=True,
        )
        if comp.get("components"):
            comp_df = pd.DataFrame(
                [
                    {
                        "العنصر": f"{item['name_ar']} ({item['symbol']})",
                        "عدد الذرات": item["count"],
                        "الوصف": item["description_ar"],
                        "الكتلة الذرية (u)": item["atomic_mass"],
                        "النسبة الكتلية (%)": f"{item['mass_percentage']}%",
                    }
                    for item in comp["components"]
                ]
            )
            st.dataframe(comp_df, use_container_width=True, hide_index=True)

    with c2:
        st.markdown(
            f"""
<div class="ml-card">
  <h4>2. البنية الجزيئية والروابط الكيميائية</h4>
  <p><strong>شكل الجزيء (Molecular Geometry):</strong> {struct['molecular_geometry']}</p>
  <p><strong>التهجين المداري:</strong> {struct.get('hybridization', 'غير محدد')}</p>
  <p><strong>الروابط الكيميائية:</strong> {struct['chemical_bonds']}</p>
  <p><strong>القطبية وعزم ثنائي القطب:</strong> {struct.get('dipole_moment', 'غير محدد')}</p>
  <p><strong>التمثيل ثلاثي الأبعاد:</strong> {struct['representation_3d_note']}</p>
</div>
""",
            unsafe_allow_html=True,
        )
        if struct.get("pubchem_3d_url"):
            st.markdown(
                f"<p><a href='{struct['pubchem_3d_url']}' target='_blank'>فتح عارض البنية ثلاثية الأبعاد التفاعلي في PubChem 3D</a></p>",
                unsafe_allow_html=True,
            )

    st.markdown(
        f"""
<div class="ml-card">
  <h4>3. الخصائص الفيزيائية والكيميائية</h4>
  <p><strong>الحالة الفيزيائية:</strong> {props['physical_state']}</p>
  <p><strong>الكثافة:</strong> <span class="ltr-inline">{props['density']}</span> | <strong>درجة الانصهار:</strong> <span class="ltr-inline">{props['melting_point']}</span> | <strong>درجة الغليان:</strong> <span class="ltr-inline">{props['boiling_point']}</span></p>
  <p><strong>الخصائص والتفاعلية الكيميائية:</strong> {props['chemical_properties_summary']} {props.get('reactivity', '')}</p>
</div>
""",
        unsafe_allow_html=True,
    )


def _render_drug_profile_ui(profile: dict[str, Any]) -> None:
    _render_compound_profile_ui(profile)
    pharm = profile.get("pharmacology", {})
    brands = "، ".join(profile.get("brand_names", [])) or "غير متوفر"
    desc = pharm.get("calculated_drug_descriptors") or {}

    st.markdown(
        f"""
<div class="ml-card">
  <h4>4. الملف الدوائي والنشاط الحيوي الموثق (Pharmacology &amp; Bioactivity)</h4>
  <p><strong>الأسماء التجارية المسجلة:</strong> {brands} | <strong>حالة الاعتماد:</strong> {pharm.get('approval_status_ar', '')}</p>
  <p><strong>الفئة العلاجية:</strong> {pharm.get('therapeutic_class_ar', '')}</p>
  <p><strong>الاستطباب الدوائي:</strong> {pharm.get('indication_ar', '')}</p>
  <p><strong>الهدف البروتيني الجزيئي:</strong> {pharm.get('target_protein_ar', '')}</p>
  <p><strong>آلية العمل الدوائية (Mechanism of Action):</strong> {pharm.get('mechanism_of_action_ar', '')}</p>
  <p><strong>النشاط الحيوي التجريبي المنشور (ChEMBL):</strong> <span class="ltr-inline">{pharm.get('documented_bioactivity_ar', '')}</span></p>
</div>
""",
        unsafe_allow_html=True,
    )
    if desc:
        d1, d2, d3, d4, d5 = st.columns(5)
        d1.metric("الوزن الجزيئي (MW)", f"{desc['molecular_weight']} Da")
        d2.metric("معامل الدهنية (LogP)", str(desc["logp"]))
        d3.metric("مانحات / مستقبلات H", f"HBD={desc['hbd']} | HBA={desc['hba']}")
        d4.metric("السطح القطبي (TPSA)", f"{desc['tpsa']} Å²")
        d5.metric("قاعدة Lipinski", "متوافق" if desc["lipinski_compliant"] else "مخالف")


def _render_protein_profile_ui(profile: dict[str, Any]) -> None:
    basic = profile["basic_info"]
    bio = profile["biology"]
    struct = profile["structure"]

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("اسم البروتين", basic["name_ar"])
    m2.metric("معرف UniProt", str(bio.get("uniprot_id", "N/A")))
    m3.metric("الكائن الحي", str(bio.get("organism", "Homo sapiens")))
    m4.metric("الكتلة الجزيئية", str(basic.get("molar_mass", "N/A")))

    c1, c2 = st.columns(2)
    with c1:
        st.markdown(
            f"""
<div class="ml-card">
  <h4>1. الهوية البيولوجية والوظيفة الحيوية (UniProt)</h4>
  <p><strong>الاسم الإنجليزي:</strong> <span class="ltr-inline">{basic['name_en']}</span> | <strong>الاسم العربي:</strong> {basic['name_ar']}</p>
  <p><strong>رمز الجين (Gene):</strong> <code class="ltr-inline">{bio.get('gene_name', '')}</code> | <strong>معرف UniProt:</strong> <code class="ltr-inline">{bio.get('uniprot_id', '')}</code></p>
  <p><strong>طول التسلسل الأميني:</strong> {bio.get('sequence_length', '')}</p>
  <p><strong>الموقع الخلوي/الحيوي:</strong> {bio.get('subcellular_location', '')}</p>
  <p><strong>الوظيفة الحيوية والفسيولوجية:</strong> {bio.get('biological_function_ar', '')}</p>
</div>
""",
            unsafe_allow_html=True,
        )
    with c2:
        pdb_items = ", ".join(bio.get("pdb_structures", [])) or "غير متوفر"
        st.markdown(
            f"""
<div class="ml-card">
  <h4>2. البنية الفراغية ثلاثية الأبعاد (RCSB PDB)</h4>
  <p><strong>البنية الجزيئية:</strong> {struct.get('molecular_geometry', '')}</p>
  <p><strong>الروابط والجسور المثبتة للبنية:</strong> {struct.get('chemical_bonds', '')}</p>
  <p><strong>التراكيب البلورية المسجلة في PDB:</strong> <code class="ltr-inline">{pdb_items}</code></p>
  <p><strong>ملاحظة تركيبية:</strong> {struct.get('representation_3d_note', '')}</p>
</div>
""",
            unsafe_allow_html=True,
        )
        for link in struct.get("pdb_links", []):
            st.markdown(
                f"<p><a href='{link['viewer_3d_url']}' target='_blank'>فتح التمثيل ثلاثي الأبعاد للبنية ({link['pdb_id']}) في RCSB PDB 3D Viewer</a></p>",
                unsafe_allow_html=True,
            )


def _render_explanation_and_applications_ui(profile: dict[str, Any]) -> None:
    apps = profile.get("applications", {})
    exp = profile.get("explanation", {})

    c_app, c_exp = st.columns(2)
    with c_app:
        st.markdown(
            f"""
<div class="ml-card">
  <h4>التطبيقات والاستخدامات العلمية والعملية</h4>
  <p><strong>1. الاستخدامات اليومية:</strong> {apps.get('daily', '')}</p>
  <p><strong>2. الاستخدامات الصناعية:</strong> {apps.get('industrial', '')}</p>
  <p><strong>3. الاستخدامات العلمية والبحثية:</strong> {apps.get('scientific', '')}</p>
</div>
""",
            unsafe_allow_html=True,
        )
    with c_exp:
        st.markdown('<div class="ml-card"><h4>التفسير العلمي (Scientific Explanation)</h4>', unsafe_allow_html=True)
        level_choice = st.radio(
            "اختر مستوى الشرح العلمي:",
            ["شرح مبسط للطالب", "شرح جامعي متقدم"],
            horizontal=True,
            key=f"exp_lvl_{profile['reaction_lab_entity']['key']}",
        )
        if level_choice == "شرح مبسط للطالب":
            st.info(exp.get("student_level_ar", ""))
        else:
            st.info(exp.get("advanced_level_ar", ""))
        st.caption(exp.get("methodology_distinction_ar", ""))
        st.markdown("</div>", unsafe_allow_html=True)


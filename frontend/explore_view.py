from typing import Any
import pandas as pd
import streamlit as st
from chemistry_engine.formula_analyzer import format_chemical_formula_display
from config.settings import (
    ENTITY_BIOMOLECULE_PROTEIN,
    ENTITY_CHEMICAL_COMPOUND,
    ENTITY_CHEMICAL_ELEMENT,
    ENTITY_DRUG_SUBSTANCE,
)
from core.orchestrator import MolecuLensCore
from frontend.structure_viewer import render_structure_viewer
from frontend.ai_assistant_view import render_ai_assistant_component
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
    st.markdown("<strong>اقتراحات فورية (Live Suggestions — اختر لتعبئة خانة البحث):</strong>", unsafe_allow_html=True)
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
            st.button(
                "اختيار",
                key=f"exp_sug_{idx}_{s['key']}",
                on_click=_set_explore_query,
                args=(s["symbol_or_formula"] if s["entity_type"] != ENTITY_BIOMOLECULE_PROTEIN else s["name_en"],),
            )


def render_explore_section(on_open_reaction_lab, on_open_prediction_lab=None, lang: str = "ar") -> None:
    if EXPLORE_QUERY_KEY not in st.session_state:
        st.session_state[EXPLORE_QUERY_KEY] = "CH3COOH"

    st.markdown(
        """
<div class="ml-card">
  <h3>Explore &amp; Understand — استكشاف المادة وفهم تركيبها وخصائصها</h3>
  <p>
    اكتب في خانة البحث الموحدة أي مدخل: الاسم العربي، الاسم الإنجليزي، الاسم النظامي (IUPAC Name)، الصيغة الكيميائية (H₂O أو CH₃COOH)،
    الرمز الكيميائي، التمثيل البنيوي (SMILES)، أو الاسم التجاري؛ وسيتعرف المحرك العلمي على نوع المدخل تلقائياً.
  </p>
</div>
""",
        unsafe_allow_html=True,
    )

    query_text = st.text_input(
        "خانة البحث العلمي الموحدة (Unified Scientific Search):",
        key=EXPLORE_QUERY_KEY,
        placeholder="مثال: حمض الخليك، Acetic Acid، CH3COOH، Aspirin، Fe",
    )
    detected = detect_query_input_type(query_text)
    if detected["type"] != "empty":
        st.caption(f"نوع المدخل المكتشف تلقائياً (Detected Input Type): {detected['label_ar']}")

    _render_suggestions(query_text)

    res = MolecuLensCore.explore_entity(query_text, allow_network=True)

    if not res["found"]:
        st.warning(res["message"])
        if res["did_you_mean"]:
            st.markdown("<strong>مقترحات بديلة يمكنك اختيارها مباشرة (Suggested Alternatives):</strong>", unsafe_allow_html=True)
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
            f"{idx_c + 1}. {c['name_ar']} | {c['name_en']} ({format_chemical_formula_display(c['symbol_or_formula'])}) — [{c['entity_type_label_ar']}]"
            for idx_c, c in enumerate(cand_list)
        ]
        chosen_lbl = st.radio(
            "توجد أكثر من مادة تشترك في هذا المدخل؛ اختر المادة المقصودة لعرض ملفها العلمي الكامل:",
            options=cand_labels,
            key=f"exp_ambig_{query_text}",
        )
        chosen_cand = cand_list[cand_labels.index(chosen_lbl)]
        res = MolecuLensCore.explore_entity(query_text, selected_candidate=chosen_cand)

    entity: dict[str, Any] = res["entity"]
    profile: dict[str, Any] = res["profile"]
    entity_type = entity["entity_type"]
    badge_cls = _badge_class_for_type(entity_type)
    disp_formula = format_chemical_formula_display(entity["symbol_or_formula"])
    id_ver = profile.get("identity_verification") or entity.get("identity_verification") or {}

    st.markdown(
        f"""
<div class="ml-card">
  <div>
    <span class="{badge_cls}">{entity['entity_type_label_ar']} | {entity['entity_type_label_en']}</span>
    <span class="ml-badge" style="background:#475569;">{entity.get('subclass_label_ar', '')}</span>
    <span class="ml-badge" style="background:#334155;">{entity.get('match_category', 'Exact Match')}</span>
  </div>
  <h3 style="margin-top:0.65rem;margin-bottom:0.35rem;">
    {entity['name_ar']} — <span class="ltr-inline">{entity['name_en']}</span> <span class="ml-formula-pill">{disp_formula}</span>
  </h3>
  <p style="margin-bottom:0.2rem;color:#64748b;font-size:0.92rem;">
    <strong>سبب المطابقة (Match Reason):</strong> {entity.get('match_reason_ar', '')} |
    <strong>المصادر العلمية الموحدة (Verified Sources):</strong> <span class="ltr-inline">{", ".join(profile.get('sources', []))}</span>
  </p>
  <p style="margin-bottom:0;color:{'#0f766e' if id_ver.get('is_verified') else '#b45309'};font-size:0.88rem;font-weight:700;">
    {id_ver.get('status_ar', '')}
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

    render_ai_assistant_component(compound_data=profile, lang=lang)

    st.markdown('<div class="ml-card">', unsafe_allow_html=True)
    st.markdown(
        f"<h4>الانتقال إلى مختبر التفاعلات (Explore Reactions / Open Reaction Lab)</h4>"
        f"<p>هل ترغب في اختبار تفاعلات <strong>{entity['name_ar']} ({disp_formula})</strong> وحساب المتفاعل المحدد (Limiting Reagent) والناتج النظري (Theoretical Yield)؟</p>",
        unsafe_allow_html=True,
    )
    nav_cols = st.columns(2)
    with nav_cols[0]:
        if st.button(
            f"Explore Reactions / Open Reaction Lab — نقل ({disp_formula}) إلى مختبر التفاعلات",
            key="btn_pass_to_reaction_lab",
        ):
            on_open_reaction_lab(profile["reaction_lab_entity"])
    if on_open_prediction_lab is not None:
        with nav_cols[1]:
            if st.button(
                f"Reaction Prediction — نقل ({disp_formula}) إلى مختبر التوقع",
                key="btn_pass_to_prediction_lab",
            ):
                on_open_prediction_lab(profile["reaction_lab_entity"])
    st.markdown("</div>", unsafe_allow_html=True)


def _render_element_profile_ui(profile: dict[str, Any]) -> None:
    basic = profile["basic_info"]
    atomic = profile["atomic_structure"]
    props = profile["properties"]

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("الاسم والرمز (Element Symbol)", f"{basic['name_ar']} ({basic['symbol']})")
    m2.metric("العدد الذري (Atomic Number)", str(basic["atomic_number"]))
    m3.metric("الكتلة الذرية (Atomic Mass)", str(basic["atomic_mass"]))
    m4.metric("الحالة الفيزيائية (Physical State)", props["physical_state"].split(" ")[0])

    c1, c2 = st.columns(2)
    with c1:
        st.markdown(
            f"""
<div class="ml-card">
  <h4>1. المعلومات الأساسية والتركيب الذري (Atomic Structure &amp; Configuration)</h4>
  <p><strong>الاسم العلمي (Scientific Name):</strong> {basic['name_ar']} | <span class="ltr-inline">{basic['name_en']}</span></p>
  <p><strong>الرمز الكيميائي (Symbol):</strong> <code class="ltr-inline">{basic['symbol']}</code> | <strong>العدد الذري (Atomic Number Z):</strong> <code class="ltr-inline">{basic['atomic_number']}</code></p>
  <p><strong>عدد البروتونات والإلكترونات:</strong> <code class="ltr-inline">{atomic['protons_count']} Protons / {atomic['electrons_count']} Electrons</code></p>
  <p><strong>التوزيع الإلكتروني (Electron Configuration):</strong> <code class="ltr-inline">{atomic['electron_configuration']}</code></p>
  <p><strong>توزيع الأغلفة الإلكترونية (Electron Shells):</strong> <code class="ltr-inline">{atomic.get('electron_shells') or 'غير محدد'}</code></p>
  <p><strong>الموقع في الجدول الدوري:</strong> {atomic.get('group', '')} — {atomic.get('period', '')} ({atomic.get('block', '')})</p>
  <p><strong>حالات الأكسدة الشائعة (Oxidation States):</strong> <span class="ltr-inline">{atomic.get('oxidation_states') or 'غير محدد'}</span></p>
  <p><strong>البنية البلورية الذرية (Crystal Structure):</strong> {atomic.get('crystal_structure') or 'شبكة بلورية عنصرية'}</p>
</div>
""",
            unsafe_allow_html=True,
        )
    with c2:
        st.markdown(
            f"""
<div class="ml-card">
  <h4>2. الخصائص الفيزيائية والكيميائية (Physical &amp; Chemical Properties)</h4>
  <p><strong>الحالة الفيزيائية (Physical State):</strong> {props['physical_state']}</p>
  <p><strong>الكثافة (Density):</strong> <span class="ltr-inline">{props['density']}</span></p>
  <p><strong>درجة الانصهار (Melting Point):</strong> <span class="ltr-inline">{props['melting_point']}</span></p>
  <p><strong>درجة الغليان (Boiling Point):</strong> <span class="ltr-inline">{props['boiling_point']}</span></p>
  <p><strong>السلوك والتفاعلية الكيميائية (Chemical Reactivity):</strong> {props['chemical_behavior']}</p>
</div>
""",
            unsafe_allow_html=True,
        )


def _render_compound_profile_ui(profile: dict[str, Any]) -> None:
    basic = profile["basic_info"]
    comp = profile["composition"]
    struct = profile["structure"]
    props = profile["properties"]
    disp_formula = basic.get("formatted_formula") or format_chemical_formula_display(basic["formula"])

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("الاسم العلمي (Compound Name)", f"{basic['name_ar']} ({basic['name_en']})")
    m2.metric("الصيغة الجزيئية (Molecular Formula)", disp_formula)
    m3.metric("الكتلة المولية (Molar Mass)", basic["molar_mass"])
    m4.metric("إجمالي الذرات (Total Atoms)", str(comp.get("total_atoms", 0)))

    c1, c2 = st.columns(2)
    with c1:
        st.markdown(
            f"""
<div class="ml-card">
  <h4>1. المعلومات الأساسية والتركيب الذري (Composition &amp; Formula)</h4>
  <p><strong>الاسم العلمي:</strong> {basic['name_ar']} | <span class="ltr-inline">{basic['name_en']}</span></p>
  <p><strong>الاسم النظامي (IUPAC Name):</strong> <code class="ltr-inline">{basic.get('iupac_name') or 'N/A'}</code></p>
  <p><strong>الصيغة الجزيئية (Molecular Formula):</strong> <code class="ltr-inline">{disp_formula}</code> | <strong>الكتلة المولية (Molar Mass):</strong> <code class="ltr-inline">{basic['molar_mass']}</code></p>
  <p><strong>ملخص التركيب الذري:</strong> {comp.get('summary_ar', '')}</p>
</div>
""",
            unsafe_allow_html=True,
        )
        if comp.get("components"):
            comp_df = pd.DataFrame(
                [
                    {
                        "العنصر (Element)": f"{item['name_ar']} ({item['symbol']})",
                        "عدد الذرات (Count)": int(item["count"]),
                        "الوصف الذري": item["description_ar"],
                        "الكتلة الذرية (Atomic Mass u)": float(item["atomic_mass"]),
                        "النسبة الكتلية (Mass %)": f"{item['mass_percentage']:.2f}%",
                    }
                    for item in comp["components"]
                ]
            )
            st.dataframe(
                comp_df,
                use_container_width=True,
                hide_index=True,
                column_config={
                    "عدد الذرات (Count)": st.column_config.NumberColumn(format="%d"),
                    "الكتلة الذرية (Atomic Mass u)": st.column_config.NumberColumn(format="%.4f u"),
                },
            )

    with c2:
        st.markdown(
            f"""
<div class="ml-card">
  <h4>2. البنية الجزيئية والروابط الكيميائية (Molecular Structure &amp; Bonding)</h4>
  <p><strong>الشكل الهندسي للجزيء (Molecular Geometry):</strong> {struct['molecular_geometry']}</p>
  <p><strong>التهجين المداري (Orbital Hybridization):</strong> {struct.get('hybridization', 'غير محدد')}</p>
  <p><strong>الروابط الكيميائية (Chemical Bonds):</strong> {struct['chemical_bonds']}</p>
  <p><strong>القطبية وعزم ثنائي القطب (Dipole Moment):</strong> {struct.get('dipole_moment', 'غير محدد')}</p>
  <p><strong>التمثيل الفراغي ثلاثي الأبعاد (3D Representation):</strong> {struct['representation_3d_note']}</p>
</div>
""",
            unsafe_allow_html=True,
        )
        if struct.get("pubchem_3d_url"):
            st.markdown(
                f"<p><a href='{struct['pubchem_3d_url']}' target='_blank'>فتح عارض البنية ثلاثية الأبعاد التفاعلي في PubChem 3D Conformer</a></p>",
                unsafe_allow_html=True,
            )

    st.markdown(
        f"""
<div class="ml-card">
  <h4>3. الخصائص الفيزيائية والكيميائية (Physicochemical Properties)</h4>
  <p><strong>الحالة الفيزيائية (Physical State):</strong> {props['physical_state']}</p>
  <p><strong>الكثافة (Density):</strong> <span class="ltr-inline">{props['density']}</span> | <strong>درجة الانصهار (Melting Point):</strong> <span class="ltr-inline">{props['melting_point']}</span> | <strong>درجة الغليان (Boiling Point):</strong> <span class="ltr-inline">{props['boiling_point']}</span></p>
  <p><strong>السلوك والتفاعلية الكيميائية (Chemical Reactivity):</strong> {props['chemical_properties_summary']} {props.get('reactivity', '')}</p>
</div>
""",
        unsafe_allow_html=True,
    )

    safety = profile.get("lab_safety")
    if safety and profile.get("profile_kind") != ENTITY_DRUG_SUBSTANCE:
        st.markdown(
            f"""
<div class="ml-safety-card">
  <h4 style="margin-top:0;color:#b45309;">إرشادات السلامة المخبرية لما قبل المختبر (Pre-Lab Safety &amp; Handling)</h4>
  <p><strong>المخاطر الكيميائية (Hazard Profile):</strong> {safety.get('hazard_ar', '')}</p>
  <p><strong>أدوات الوقاية الشخصية (Required PPE):</strong> {safety.get('ppe_ar', '')}</p>
  <p style="margin-bottom:0;"><strong>احتياطات التداول المخبري (Safe Lab Handling):</strong> {safety.get('handling_ar', '')}</p>
</div>
""",
            unsafe_allow_html=True,
        )


def _render_drug_profile_ui(profile: dict[str, Any]) -> None:
    _render_compound_profile_ui(profile)
    pharm = profile.get("pharmacology", {})
    brands = "، ".join(profile.get("brand_names", [])) or "غير متوفر"
    desc = pharm.get("calculated_drug_descriptors") or {}
    prelab = pharm.get("prelab_safety_and_qc") or {}

    st.markdown(
        f"""
<div class="ml-card">
  <h4>4. الملف الدوائي والنشاط الحيوي الموثق (Pharmacology &amp; Bioactivity — ChEMBL)</h4>
  <p><strong>الأسماء التجارية (Brand Names):</strong> {brands} | <strong>حالة الاعتماد (Approval Status):</strong> {pharm.get('approval_status_ar', '')}</p>
  <p><strong>الفئة العلاجية (Therapeutic Class):</strong> {pharm.get('therapeutic_class_ar', '')}</p>
  <p><strong>الاستطباب الدوائي (Clinical Indication):</strong> {pharm.get('indication_ar', '')}</p>
  <p><strong>الهدف البروتيني الجزيئي (Molecular Target):</strong> {pharm.get('target_protein_ar', '')}</p>
  <p><strong>آلية العمل الدوائية (Mechanism of Action):</strong> {pharm.get('mechanism_of_action_ar', '')}</p>
  <p><strong>النشاط الحيوي التجريبي المنشور (Documented Bioactivity):</strong> <span class="ltr-inline">{pharm.get('documented_bioactivity_ar', '')}</span></p>
</div>
""",
        unsafe_allow_html=True,
    )
    if desc:
        d1, d2, d3, d4, d5 = st.columns(5)
        d1.metric("الوزن الجزيئي (Molecular Weight)", f"{desc['molecular_weight']} Da")
        d2.metric("معامل التوزع الدهني (LogP)", str(desc["logp"]))
        d3.metric("الروابط الهيدروجينية (HBD / HBA)", f"HBD={desc['hbd']} | HBA={desc['hba']}")
        d4.metric("مساحة السطح القطبي (TPSA)", f"{desc['tpsa']} Å²")
        d5.metric("قاعدة ليبينسكي (Lipinski Rule)", "متوافق (Compliant)" if desc["lipinski_compliant"] else "مخالف (Non-compliant)")
    else:
        st.info(
            pharm.get(
                "descriptors_status_ar",
                "الخصائص الجزيئية الدوائية المحسوبة (MW, LogP, HBD, HBA, TPSA) غير متاحة حالياً لعدم إمكانية حسابها فعلياً، ويلتزم النظام بعدم عرض قيم افتراضية.",
            )
        )

    if prelab:
        st.markdown(
            f"""
<div class="ml-safety-card">
  <h4 style="margin-top:0;color:#b45309;">5. بطاقة تهيئة طالب الصيدلة قبل المختبر (Pharmacy Pre-Lab Safety &amp; QC Protocol)</h4>
  <p><strong>تقييم المخاطر الدوائية (Hazard Assessment):</strong> {prelab.get('hazard_ar', '')}</p>
  <p><strong>معدات الوقاية الشخصية (Required PPE):</strong> {prelab.get('ppe_ar', '')}</p>
  <p><strong>اختبار الهوية والنقاوة الدستوري (Pharmacopoeial QC Test):</strong> {prelab.get('quality_test_ar', '')}</p>
  <p style="margin-bottom:0;"><strong>التداول والحفظ في المختبر (Handling &amp; Storage):</strong> {prelab.get('handling_ar', '')}</p>
</div>
""",
            unsafe_allow_html=True,
        )


def _render_protein_profile_ui(profile: dict[str, Any]) -> None:
    basic = profile["basic_info"]
    bio = profile["biology"]
    struct = profile["structure"]

    m1, m2, m3, m4 = st.columns(4)
    m1.metric("اسم البروتين (Protein Name)", basic["name_ar"])
    m2.metric("معرف UniProt", str(bio.get("uniprot_id", "N/A")))
    m3.metric("الكائن الحي (Organism)", str(bio.get("organism", "Homo sapiens")))
    m4.metric("الكتلة الجزيئية (Molecular Mass)", str(basic.get("molar_mass", "N/A")))

    c1, c2 = st.columns(2)
    with c1:
        st.markdown(
            f"""
<div class="ml-card">
  <h4>1. الهوية البيولوجية والوظيفة الحيوية (Biological Identity — UniProt)</h4>
  <p><strong>الاسم العلمي:</strong> {basic['name_ar']} | <span class="ltr-inline">{basic['name_en']}</span></p>
  <p><strong>رمز الجين (Gene Symbol):</strong> <code class="ltr-inline">{bio.get('gene_name', '')}</code> | <strong>معرف UniProtKB:</strong> <code class="ltr-inline">{bio.get('uniprot_id', '')}</code></p>
  <p><strong>طول التسلسل الأميني (Sequence Length):</strong> {bio.get('sequence_length', '')}</p>
  <p><strong>الموقع الخلوي (Subcellular Location):</strong> {bio.get('subcellular_location', '')}</p>
  <p><strong>الوظيفة الحيوية والفسيولوجية (Biological Function):</strong> {bio.get('biological_function_ar', '')}</p>
</div>
""",
            unsafe_allow_html=True,
        )
    with c2:
        pdb_items = ", ".join(bio.get("pdb_structures", [])) or "غير متوفر"
        st.markdown(
            f"""
<div class="ml-card">
  <h4>2. البنية الفراغية ثلاثية الأبعاد (3D Crystal Structure — RCSB PDB)</h4>
  <p><strong>البنية الجزيئية (Quaternary/Tertiary Geometry):</strong> {struct.get('molecular_geometry', '')}</p>
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
  <h4>التطبيقات والاستخدامات العلمية والعملية (Scientific &amp; Practical Applications)</h4>
  <p><strong>1. الاستخدامات اليومية والطبية (Daily &amp; Medical):</strong> {apps.get('daily', '')}</p>
  <p><strong>2. الاستخدامات الصناعية والدوائية (Industrial &amp; Pharmaceutical):</strong> {apps.get('industrial', '')}</p>
  <p><strong>3. الاستخدامات العلمية والمخبرية (Laboratory &amp; Research):</strong> {apps.get('scientific', '')}</p>
</div>
""",
            unsafe_allow_html=True,
        )
    with c_exp:
        st.markdown('<div class="ml-card"><h4>التفسير العلمي المتدرج (Scientific Explanation)</h4>', unsafe_allow_html=True)
        level_choice = st.radio(
            "اختر مستوى الشرح العلمي:",
            ["شرح مبسط للطالب (Student Level)", "شرح جامعي متقدم (Advanced University Level)"],
            horizontal=True,
            key=f"exp_lvl_{profile['reaction_lab_entity']['key']}",
        )
        if "مبسط" in level_choice:
            st.info(exp.get("student_level_ar", ""))
        else:
            st.info(exp.get("advanced_level_ar", ""))
        st.caption(exp.get("methodology_distinction_ar", ""))
        st.markdown("</div>", unsafe_allow_html=True)

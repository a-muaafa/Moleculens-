import json
from typing import Any

import plotly.graph_objects as go
import streamlit as st
import streamlit.components.v1 as components

from chemistry_engine.formula_analyzer import format_chemical_formula_display
from chemistry_engine.structure_analyzer import (
    analyze_molecular_structure,
    generate_3d_conformer_data,
    render_molecule_2d_svg,
    structure_matches_formula,
    verify_substance_structure_identity,
)

ZOOM_LEVELS = [0.6, 0.8, 1.0, 1.25, 1.5, 2.0, 2.5]
BASE_2D_SIZE = (440, 300)
BOND_COLOR = "#64748b"


@st.cache_data(show_spinner=False)
def _cached_svg(smiles: str, zoom: float, show_indices: bool) -> str | None:
    return render_molecule_2d_svg(smiles, int(BASE_2D_SIZE[0] * zoom), int(BASE_2D_SIZE[1] * zoom), show_indices)


@st.cache_data(show_spinner=False)
def _cached_analysis(smiles: str) -> dict[str, Any] | None:
    return analyze_molecular_structure(smiles)


@st.cache_data(show_spinner="جارٍ توليد البنية ثلاثية الأبعاد عبر RDKit...")
def _cached_conformer(smiles: str) -> dict[str, Any] | None:
    return generate_3d_conformer_data(smiles)


def build_plotly_3d_figure(conf: dict[str, Any]) -> go.Figure:
    atoms = conf["atoms"]
    fig = go.Figure()
    for b in conf["bonds"]:
        a1, a2 = atoms[b["begin"]], atoms[b["end"]]
        fig.add_trace(
            go.Scatter3d(
                x=[a1["x"], a2["x"]], y=[a1["y"], a2["y"]], z=[a1["z"], a2["z"]],
                mode="lines", line={"color": BOND_COLOR, "width": 3 + 3 * (b["order"] - 1)},
                hoverinfo="skip", showlegend=False,
            )
        )
    fig.add_trace(
        go.Scatter3d(
            x=[a["x"] for a in atoms], y=[a["y"] for a in atoms], z=[a["z"] for a in atoms],
            mode="markers+text" if len(atoms) <= 40 else "markers",
            text=[a["symbol"] for a in atoms], textposition="top center",
            marker={"size": [a["radius"] for a in atoms], "color": [a["color"] for a in atoms], "line": {"color": "#0f172a", "width": 1}},
            hovertext=[f"{a['symbol']} #{a['index']}" for a in atoms], hoverinfo="text", showlegend=False,
        )
    )
    axis = {"visible": False}
    fig.update_layout(scene={"xaxis": axis, "yaxis": axis, "zaxis": axis, "aspectmode": "data"}, margin={"l": 0, "r": 0, "t": 0, "b": 0}, height=460)
    return fig


def _3dmol_html(mol_block: str | None, pdb_id: str | None) -> str:
    if pdb_id:
        loader = f"$3Dmol.download('pdb:{pdb_id}', viewer, {{}}, function(){{viewer.setStyle({{}},{{cartoon:{{color:'spectrum'}}}}); viewer.zoomTo(); viewer.render();}});"
    else:
        loader = (
            f"viewer.addModel({json.dumps(mol_block)}, 'mol'); viewer.setStyle({{}}, {{stick:{{radius:0.15}}, sphere:{{scale:0.25}}}}); viewer.zoomTo(); viewer.render();"
        )
    return f"""
<div id="mol3d" style="width:100%;height:440px;position:relative;border:1px solid #cbd5e1;border-radius:10px;"></div>
<script src="https://3dmol.org/build/3Dmol-min.js"></script>
<script>
  var viewer = $3Dmol.createViewer(document.getElementById('mol3d'), {{backgroundColor:'white'}});
  {loader}
  viewer.spin(false);
</script>
"""


def render_structure_viewer(entity: dict[str, Any], key_prefix: str = "sv") -> None:
    smiles = (entity.get("canonical_smiles") or "").strip()
    formula = entity.get("symbol_or_formula", "")
    disp_formula = format_chemical_formula_display(formula)
    pdb_ids = entity.get("identifiers", {}).get("pdb_ids") or entity.get("pdb_ids") or []
    ekey = entity.get("key", "x")
    zoom_key = f"{key_prefix}_zoom_{ekey}"

    st.markdown('<div class="ml-card"><h4>عارض التراكيب الكيميائية ثنائي وثلاثي الأبعاد (2D &amp; 3D Molecular Structure Viewer)</h4>', unsafe_allow_html=True)

    if not smiles and not pdb_ids:
        st.info("لا تتوفر بنية جزيئية مرجعية لهذه المادة في البيانات الحالية، لذلك لا يعرض النظام أي تركيب بديل أو عشوائي.")
        st.markdown("</div>", unsafe_allow_html=True)
        return

    if smiles:
        id_check = entity.get("identity_verification") or verify_substance_structure_identity(entity)
        verified = structure_matches_formula(smiles, formula) or entity.get("entity_type") == "chemical_element"
        badge = (
            id_check.get("status_ar")
            if id_check.get("status_ar")
            else ("تحقق RDKit: البنية تطابق الصيغة المسجلة" if verified else "تنبيه: الصيغة المحسوبة من البنية تختلف عن الصيغة المسجلة")
        )
        st.caption(f"{badge} | المادة المرتبطة: {entity.get('name_ar')} ({disp_formula}) | Canonical SMILES: {smiles}")

    tab_labels = []
    if smiles:
        tab_labels.append("2D")
    if smiles or pdb_ids:
        tab_labels.append("3D")
    tabs = st.tabs(tab_labels)
    tab_map = dict(zip(tab_labels, tabs))

    if smiles and "2D" in tab_map:
        with tab_map["2D"]:
            if zoom_key not in st.session_state:
                st.session_state[zoom_key] = ZOOM_LEVELS.index(1.0)

            def _zoom(delta: int) -> None:
                st.session_state[zoom_key] = max(0, min(len(ZOOM_LEVELS) - 1, st.session_state[zoom_key] + delta))

            c1, c2, c3, c4 = st.columns([1, 1, 1.4, 2])
            c1.button("تكبير (Zoom In +)", key=f"{zoom_key}_in", on_click=_zoom, args=(1,))
            c2.button("تصغير (Zoom Out -)", key=f"{zoom_key}_out", on_click=_zoom, args=(-1,))
            c3.markdown(f"مستوى التكبير (Zoom): **{int(ZOOM_LEVELS[st.session_state[zoom_key]] * 100)}%**")
            show_idx = c4.checkbox("إظهار أرقام الذرات (Show Atom Indices)", key=f"{key_prefix}_idx_{ekey}")
            svg = _cached_svg(smiles, ZOOM_LEVELS[st.session_state[zoom_key]], show_idx)
            if svg:
                st.markdown(
                    f'<div style="overflow:auto;direction:ltr;background:#fff;border:1px solid #e2e8f0;border-radius:10px;padding:6px;text-align:center;">{svg}</div>',
                    unsafe_allow_html=True,
                )
            analysis = _cached_analysis(smiles)
            if analysis:
                b = analysis["bonds"]
                m1, m2, m3, m4 = st.columns(4)
                m1.metric("روابط أحادية (Single)", b["single"])
                m2.metric("روابط ثنائية (Double)", b["double"])
                m3.metric("روابط ثلاثية (Triple)", b["triple"])
                m4.metric("روابط عطرية (Aromatic)", b["aromatic"])
                st.markdown(
                    f"**الحلقات (Rings):** {analysis['rings']['summary_ar']}  \n"
                    f"**المجموعات الوظيفية (Functional Groups):** {analysis['functional_groups_summary_ar']}  \n"
                    f"**الكيمياء الفراغية (Stereochemistry):** {analysis['stereochemistry']['summary_ar']}"
                )

    if "3D" in tab_map:
        with tab_map["3D"]:
            conf = _cached_conformer(smiles) if smiles else None
            if conf:
                st.caption("نموذج كرات وعصي ثلاثي الأبعاد (3D Ball-and-Stick) قابل للتدوير والتكبير، مولد عبر خوارزمية ETKDGv3 وتحسين الطاقة MMFF94/UFF.")
                st.plotly_chart(build_plotly_3d_figure(conf), use_container_width=True, key=f"{key_prefix}_plotly_{ekey}")
                with st.expander("عارض 3Dmol.js التفاعلي (Interactive WebGL 3D Viewer)"):
                    components.html(_3dmol_html(conf["mol_block"], None), height=470)
            elif pdb_ids:
                st.caption(f"عرض البنية البلورية المسجلة في RCSB PDB ({pdb_ids[0]}).")
                components.html(_3dmol_html(None, pdb_ids[0]), height=470)
            else:
                st.info("تعذر توليد بنية ثلاثية الأبعاد موثوقة لهذه المادة.")
    st.markdown("</div>", unsafe_allow_html=True)

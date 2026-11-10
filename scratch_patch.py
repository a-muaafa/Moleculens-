p="frontend/reaction_lab_view.py"
t=open(p,encoding="utf-8").read()
block='''    doc = rxn_result["documented_reaction_info"]
    calc = rxn_result["calculated_stoichiometry_result"]
    theo = rxn_result["theoretical_prediction_notes"]
    d1, d2, d3 = st.columns(3)
    d1.info(f"**{doc['label_ar']}**  \\nالمصدر: {doc['source']}  \\nالنوع: {doc['reaction_type_ar']}")
    d2.success(f"**{calc['label_ar']}**  \\nالمعادلة موزونة وحفظ الذرات: {'متحقق' if calc['atom_conservation_verified'] else 'غير متحقق'}")
    d3.warning(f"**{theo['label_ar']}**  \\n{theo['notes_ar']}")
    if rxn_result.get("alternative_documented_reactions"):
        alt_names = " ، ".join(a["reaction_name_ar"] for a in rxn_result["alternative_documented_reactions"])
        st.caption(f"توجد تفاعلات موثقة أخرى لنفس المواد (غيّر نوع التفاعل لاختيارها): {alt_names}")

    primary_yield = rxn_result["theoretical_yields"][0]'''
t=t.replace("    primary_yield = rxn_result[\"theoretical_yields\"][0]",block,1)
open(p,"w",encoding="utf-8").write(t)

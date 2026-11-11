import ast
import os
import re
import tempfile

try:
    import pytest
except ImportError:
    class _DummyPytest:
        @staticmethod
        def approx(val, rel=1e-6, abs=1e-6):
            return val

        class mark:
            @staticmethod
            def parametrize(argnames, argvalues):
                def decorator(func):
                    func._parametrize = (argnames, argvalues)
                    return func
                return decorator

    pytest = _DummyPytest()  # type: ignore

from chemistry_engine.formula_analyzer import parse_chemical_formula
from chemistry_engine.structure_analyzer import (
    HAS_RDKIT,
    Chem,
    analyze_molecular_structure,
    generate_3d_conformer_data,
    render_molecule_2d_svg,
    structure_matches_formula,
    validate_and_parse_smiles,
)
from config.settings import CLOSED_DROPDOWN_OPTIONS
from data_integrations.local_db import SmartLocalDatabase
from reaction_engine.equation_balancer import verify_atom_conservation
from reaction_engine.reaction_registry import (
    VERIFIED_REACTIONS_REGISTRY,
    find_verified_reactions,
    search_reactions_by_substance,
)
from reaction_engine.reaction_service import analyze_chemical_reaction
from search.local_catalog import LOCAL_SCIENTIFIC_CATALOG, get_catalog_entity_by_key
from search.normalizer import detect_query_input_type, highlight_matched_text
from search.search_engine import get_live_suggestions, get_reaction_substance_options, search_scientific_entity

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))


def _sel(q):
    res = search_scientific_entity(q, allow_network=False)
    assert res["found"], q
    assert not res["is_ambiguous"], (q, [c["key"] for c in res["candidates"]])
    return res["selected"]


@pytest.mark.parametrize("query", ["CH3COOH", "ch3cooh", "حمض الخليك", "حمض الأسيتيك", "Acetic acid", "acetic acid", "Ethanoic acid", "AcOH", "الخل"])
def test_acetic_acid_recognized_everywhere(query):
    assert _sel(query)["key"] == "acetic_acid"


def test_acetic_acid_formula_not_confused_with_isomer():
    assert _sel("CH3COOH")["key"] == "acetic_acid"
    assert _sel("HCOOCH3")["key"] == "methyl_formate"
    res = search_scientific_entity("C2H4O2", allow_network=False)
    assert res["is_ambiguous"]
    assert {"acetic_acid", "methyl_formate"} <= {c["key"] for c in res["candidates"]}


@pytest.mark.parametrize("a,b", [("H2O", "H₂O"), ("CO2", "CO₂"), ("H2SO4", "H₂SO₄")])
def test_subscript_variants_resolve_identically(a, b):
    assert _sel(a)["key"] == _sel(b)["key"]


@pytest.mark.parametrize(
    "query,key",
    [
        ("NaOH", "sodium_hydroxide_naoh"), ("Ca(OH)2", "calcium_hydroxide"), ("هيدروكسيد الكالسيوم", "calcium_hydroxide"),
        ("NaHCO3", "sodium_bicarbonate"), ("بيكربونات الصوديوم", "sodium_bicarbonate"), ("Baking soda", "sodium_bicarbonate"),
        ("C6H6", "benzene"), ("Phenol", "phenol"), ("C6H5OH", "phenol"), ("Ethanol", "ethanol"), ("الأسيتون", "acetone"),
        ("CH3COOCH2CH3", "ethyl_acetate"), ("Glycerol", "glycerol"), ("Metformin", "metformin"), ("Paracetamol", "paracetamol"),
        ("Acetaminophen", "paracetamol"), ("ATP", "atp"), ("Urea", "urea"), ("Hemoglobin", "human_hemoglobin"), ("NH4+", "ammonium_ion"),
        ("SO4^2-", "sulfate_ion"), ("Cl-", "chloride_ion"), ("Au", "element_au"), ("الذهب", "element_au"), ("Oganesson", "element_og"),
    ],
)
def test_family_coverage_search(query, key):
    assert _sel(query)["key"] == key


def test_case_sensitive_formula_priority():
    assert _sel("CO")["key"] == "carbon_monoxide"
    assert _sel("Co")["key"] == "element_co"
    res = search_scientific_entity("co", allow_network=False)
    assert res["is_ambiguous"]


def test_isopropanol_and_propanol_are_separated():
    assert _sel("(CH3)2CHOH")["key"] == "isopropanol"
    assert _sel("CH3CH2CH2OH")["key"] == "propanol"
    res = search_scientific_entity("C3H8O", allow_network=False)
    assert res["is_ambiguous"]
    assert {"isopropanol", "propanol"} <= {c["key"] for c in res["candidates"]}


def test_smiles_search_resolves_to_exact_substance():
    assert _sel("CC(=O)O")["key"] == "acetic_acid"


def test_input_type_auto_detection():
    assert detect_query_input_type("حمض الخليك")["type"] == "arabic_name"
    assert detect_query_input_type("CH3COOH")["type"] == "condensed_formula"
    assert detect_query_input_type("H2O")["type"] == "chemical_formula"
    assert detect_query_input_type("Fe")["type"] == "element_symbol"
    assert detect_query_input_type("CC(=O)O")["type"] == "smiles"
    assert detect_query_input_type("aspirin")["type"] == "english_name"


def test_yellow_highlighting_markup():
    out = highlight_matched_text("Acetic Acid", "acet")
    assert out == '<mark class="ml-highlight">Acet</mark>ic Acid'
    assert '<mark class="ml-highlight">خلي</mark>' in highlight_matched_text("حمض الخليك", "خلي")
    assert "<mark" not in highlight_matched_text("Water", "")
    assert "&lt;" in highlight_matched_text("<b>", "x")


def test_suggestions_have_name_formula_class_and_highlight():
    suggs = get_live_suggestions("acet")
    assert suggs and suggs[0]["key"] == "acetic_acid"
    top = suggs[0]
    assert "ml-highlight" in top["highlighted_name_en"]
    assert top["symbol_or_formula"] == "CH3COOH"
    assert top["subclass_label_ar"] == "حمض كربوكسيلي عضوي"
    assert [s["match_score"] for s in suggs] == sorted((s["match_score"] for s in suggs), reverse=True)


def test_suggestions_by_formula_prefix_and_typo():
    assert any(s["key"] == "acetic_acid" for s in get_live_suggestions("CH3C"))
    assert any(s["key"] == "sitagliptin" for s in get_live_suggestions("Sitaglptin"))


def test_results_ranked_by_match_quality():
    res = search_scientific_entity("ace", allow_network=False)
    scores = [c["match_score"] for c in res["candidates"]]
    assert scores == sorted(scores, reverse=True)


def test_catalog_coverage_and_integrity():
    keys = [e["key"] for e in LOCAL_SCIENTIFIC_CATALOG]
    assert len(keys) == len(set(keys))
    assert len(LOCAL_SCIENTIFIC_CATALOG) == 241
    elements = [e for e in LOCAL_SCIENTIFIC_CATALOG if e["entity_type"] == "chemical_element"]
    assert len(elements) == 118
    subclasses = {e.get("subclass") for e in LOCAL_SCIENTIFIC_CATALOG}
    for needed in ["element", "acid", "base", "salt", "ion", "inorganic_compound", "organic_compound", "aromatic_compound", "alcohol", "aldehyde", "ketone", "ester", "drug", "biomolecule"]:
        assert needed in subclasses


def test_every_structure_is_valid_and_matches_its_substance():
    for e in LOCAL_SCIENTIFIC_CATALOG:
        smiles = e.get("canonical_smiles")
        if not smiles:
            continue
        if HAS_RDKIT and Chem is not None:
            assert Chem.MolFromSmiles(smiles) is not None, e["key"]
        else:
            mol, _, _ = validate_and_parse_smiles(smiles)
            assert mol is not None, e["key"]
        assert structure_matches_formula(smiles, e["symbol_or_formula"]), (e["key"], e["symbol_or_formula"])


def test_structure_analysis_bonds_rings_groups():
    acetic = analyze_molecular_structure(get_catalog_entity_by_key("acetic_acid")["canonical_smiles"])
    assert acetic["bonds"]["double"] == 1
    assert any(g["key"] == "carboxylic_acid" for g in acetic["functional_groups"])
    benzene = analyze_molecular_structure("c1ccccc1")
    assert benzene["rings"]["aromatic_rings"] == 1 and benzene["bonds"]["aromatic"] == 6
    assert analyze_molecular_structure("C#C")["bonds"]["triple"] == 1
    assert analyze_molecular_structure("C[C@H](N)C(=O)O")["stereochemistry"]["chiral_centers_count"] == 1


def test_2d_svg_and_3d_conformer_generation():
    svg = render_molecule_2d_svg("CC(=O)O", 400, 300)
    assert svg and "<svg" in svg
    bigger = render_molecule_2d_svg("CC(=O)O", 800, 600)
    assert "width='800px'" in bigger or 'width="800px"' in bigger or "800" in bigger
    conf = generate_3d_conformer_data("CC(=O)O")
    assert conf["atom_count"] == 8
    assert conf["bond_count"] == 7
    assert all({"x", "y", "z"} <= set(a) for a in conf["atoms"])


def test_local_database_dedup_update_and_timestamps():
    with tempfile.TemporaryDirectory() as tmp:
        db = SmartLocalDatabase(os.path.join(tmp, "t.db"))
        entry = dict(get_catalog_entity_by_key("acetic_acid"))
        db.upsert_substance(entry)
        first = db.lookup_by_term("حمض الخليك")[0]
        entry_dup = dict(entry, key="another_key_same_cid")
        db.upsert_substance(entry_dup)
        assert db.get_database_statistics()["substances"] == 1
        changed = dict(entry, molar_mass=60.052)
        db.upsert_substance(changed)
        after = db.lookup_by_formula("CH3COOH")[0]
        assert after["molar_mass"] == 60.052
        assert after["fetched_at"] == first["fetched_at"] and after["updated_at"] >= first["updated_at"]
        assert db.lookup_by_term("Ethanoic acid")[0]["key"] == "acetic_acid"
        assert db.get_database_statistics()["synonyms"] > 5


def test_local_database_external_records_never_overwrite_curated():
    with tempfile.TemporaryDirectory() as tmp:
        db = SmartLocalDatabase(os.path.join(tmp, "t.db"))
        curated = dict(get_catalog_entity_by_key("acetic_acid"))
        db.upsert_substance(curated, source="Curated Scientific Catalog")
        db.upsert_substance(dict(curated, name_en="Hacked"), source="External: PubChem")
        assert db.lookup_by_formula("CH3COOH")[0]["name_en"] == "Acetic Acid"


def test_local_database_reactions_search_by_formula():
    with tempfile.TemporaryDirectory() as tmp:
        db = SmartLocalDatabase(os.path.join(tmp, "t.db"))
        db.seed(LOCAL_SCIENTIFIC_CATALOG[:5], VERIFIED_REACTIONS_REGISTRY)
        stats = db.get_database_statistics()
        assert stats["reactions"] == len(VERIFIED_REACTIONS_REGISTRY)
        assert any("CH3COOH" in r["equation_text"] or any(x["formula"] == "CH3COOH" for x in r["reactants"]) for r in db.search_reactions_by_formula("CH3COOH"))


def _walk_widgets(path):
    tree = ast.parse(open(path, encoding="utf-8").read())
    return [n.func.attr for n in ast.walk(tree) if isinstance(n, ast.Call) and isinstance(n.func, ast.Attribute)]


def test_exactly_one_free_text_input_in_whole_frontend():
    total = 0
    for fname in os.listdir(os.path.join(ROOT, "frontend")):
        if fname.endswith(".py"):
            calls = _walk_widgets(os.path.join(ROOT, "frontend", fname))
            total += calls.count("text_input") + calls.count("text_area")
    assert total == 1
    assert "text_input" in _walk_widgets(os.path.join(ROOT, "frontend", "explore_view.py"))
    assert "text_input" not in _walk_widgets(os.path.join(ROOT, "frontend", "reaction_lab_view.py"))


def test_explore_has_no_search_type_dropdown_and_dropdowns_are_closed():
    src = open(os.path.join(ROOT, "frontend", "explore_view.py"), encoding="utf-8").read()
    assert "selectbox" not in src and "نوع البحث" not in src
    styles = open(os.path.join(ROOT, "frontend", "styles.py"), encoding="utf-8").read()
    assert "pointer-events: none" in styles and "readonly" in styles
    lab = open(os.path.join(ROOT, "frontend", "reaction_lab_view.py"), encoding="utf-8").read()
    assert "accept_new_options" not in lab
    for key in ["formula_types", "substance_categories", "quantity_types", "amount_units", "reaction_types", "temperature_conditions", "pressure_conditions", "solvent_conditions", "catalyst_conditions"]:
        assert CLOSED_DROPDOWN_OPTIONS[key]


def test_reaction_substance_options_filtering():
    assert all(o["subclass"] == "acid" or o["subclass"] == "carboxylic_acid" for o in get_reaction_substance_options("acid"))
    assert len(get_reaction_substance_options("element")) == 118
    assert any(o["formula"] == "CH3COOH" for o in get_reaction_substance_options("acid"))


def test_registry_size_unique_and_atom_conservation_for_every_reaction():
    assert len(VERIFIED_REACTIONS_REGISTRY) == 87
    ids = [r["reaction_id"] for r in VERIFIED_REACTIONS_REGISTRY]
    assert len(ids) == len(set(ids))
    categories = {r["category"] for r in VERIFIED_REACTIONS_REGISTRY}
    assert {"neutralization", "combustion", "precipitation", "single_displacement", "redox_synthesis", "decomposition", "organic_reaction", "biochemical_reaction"} <= categories
    for r in VERIFIED_REACTIONS_REGISTRY:
        assert r.get("source"), r["reaction_id"]
        for s in r["reactants"] + r["products"]:
            assert parse_chemical_formula(s["formula"]), (r["reaction_id"], s["formula"])
        assert verify_atom_conservation(r["reactants"], r["products"])["is_balanced"], r["reaction_id"]


def _run(formulas_amounts, reaction_type="all", domain="all_domains"):
    inputs = []
    for f, amt in formulas_amounts:
        ent = {"symbol_or_formula": f, "name_ar": f}
        inputs.append({"substance": ent, "formula": f, "amount": amt, "unit": "mol"})
    return analyze_chemical_reaction(inputs, domain=domain, custom_conditions={"reaction_type": reaction_type})


def test_acetic_acid_neutralization_and_limiting_reagent():
    res = _run([("CH3COOH", 2.0), ("NaOH", 1.0)])
    assert res["supported"]
    assert res["equation_text"] == "CH3COOH + NaOH -> CH3COONa + H2O"
    assert res["limiting_reagent"]["formula"] == "NaOH"
    yields = {y["formula"]: y["theoretical_yield_mol"] for y in res["theoretical_yields"]}
    assert abs(yields["CH3COONa"] - 1.0) < 1e-6
    assert res["excess_reagents"][0]["formula"] == "CH3COOH"


def test_acetic_acid_bicarbonate_produces_co2():
    res = _run([("CH3COOH", 1.0), ("NaHCO3", 1.0)])
    assert res["supported"] and "CO2" in res["equation_text"]


def test_stoichiometry_mass_units_and_nonunit_coefficients():
    res = _run([("H2SO4", 1.0), ("NaOH", 1.0)])
    assert res["supported"] and res["limiting_reagent"]["formula"] == "NaOH"
    y = {p["formula"]: p["theoretical_yield_mol"] for p in res["theoretical_yields"]}
    assert abs(y["Na2SO4"] - 0.5) < 1e-6


@pytest.mark.parametrize(
    "reactants,fragment",
    [
        ([("AgNO3", 1.0), ("NaCl", 1.0)], "AgCl"), ([("Zn", 1.0), ("CuSO4", 1.0)], "ZnSO4"), ([("CaCO3", 1.0)], "CaO + CO2"),
        ([("C2H5OH", 1.0), ("O2", 3.0)], "2CO2"), ([("CH3COOH", 1.0), ("C2H5OH", 1.0)], "CH3COOCH2CH3"),
        ([("Na", 2.0), ("H2O", 2.0)], "H2"), ([("C6H12O6", 1.0)], "C2H5OH"),
    ],
)
def test_diverse_documented_reactions(reactants, fragment):
    res = _run(reactants)
    assert res["supported"], reactants
    assert fragment in res["equation_text"]
    assert res["calculated_stoichiometry_result"]["atom_conservation_verified"]


def test_refuses_undocumented_inert_and_isomer_confusion():
    assert not _run([("Fe", 1.0), ("He", 1.0)])["supported"]
    assert not _run([("NaCl", 1.0), ("KCl", 1.0)])["supported"]
    assert not _run([("Cu", 1.0), ("H2O", 1.0)])["supported"]
    assert not _run([("HCOOCH3", 1.0), ("NaOH", 1.0)])["supported"]


def test_three_way_distinction_fields_and_alternatives():
    res = _run([("CO2", 1.0), ("H2O", 1.0)])
    assert res["documented_reaction_info"]["label_ar"] == "التفاعل الموثق"
    assert res["calculated_stoichiometry_result"]["label_ar"] == "النتيجة المحسوبة"
    assert res["theoretical_prediction_notes"]["label_ar"] == "التوقع النظري"
    assert res["alternative_documented_reactions"]
    bio = _run([("CO2", 6.0), ("H2O", 6.0)], reaction_type="biochemical_reaction")
    assert bio["supported"] and "C6H12O6" in bio["equation_text"]
    mismatch = _run([("HCl", 1.0), ("NaOH", 1.0)], reaction_type="combustion")
    assert not mismatch["supported"] and mismatch["error_type"] == "reaction_type_mismatch"


def test_search_reactions_by_substance_formula():
    found = search_reactions_by_substance("CH3COOH")
    assert len(found) >= 4
    assert all(any(s["formula"] == "CH3COOH" for s in r["reactants"] + r["products"]) for r in found)
    assert find_verified_reactions(["HCl", "NaOH"])
    assert not search_reactions_by_substance("HCOOCH3")
    # Domain filtering on substance reaction search
    org_only = search_reactions_by_substance("CH3COOH", domain="organic_chemistry")
    assert len(org_only) == 2
    assert all(r["domain"] == "organic_chemistry" for r in org_only)


def test_no_emoji_in_project_sources():
    emoji = re.compile("[\U0001F300-\U0001FAFF\u2600-\u27BF]")
    for folder in ["frontend", "search", "reaction_engine", "data_integrations", "core", "chemistry_engine"]:
        for fname in os.listdir(os.path.join(ROOT, folder)):
            if fname.endswith(".py"):
                text = open(os.path.join(ROOT, folder, fname), encoding="utf-8-sig").read()
                assert not emoji.search(text), (folder, fname)

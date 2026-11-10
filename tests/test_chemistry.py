import pytest
from chemistry.descriptors import calculate_descriptors
from chemistry.molecular_structure import get_molecular_formula, render_mol_png, validate_smiles


def test_validate_smiles_valid_and_invalid():
    mol, canon, err = validate_smiles("CC(=O)Oc1ccccc1C(=O)O")
    assert mol is not None
    assert canon is not None
    assert err is None

    mol_bad, canon_bad, err_bad = validate_smiles("INVALID_SMILES_STRING_###")
    assert mol_bad is None
    assert canon_bad is None
    assert err_bad is not None

    mol_empty, _, err_empty = validate_smiles("   ")
    assert mol_empty is None
    assert err_empty is not None


def test_calculate_descriptors_and_formula():
    desc = calculate_descriptors("CC(=O)Oc1ccccc1C(=O)O")
    assert desc["molecular_formula"] == "C9H8O4"
    assert 179.0 <= desc["molecular_weight"] <= 181.0
    assert desc["hbd"] == 1
    assert desc["hba"] == 3
    assert desc["lipinski"]["passes_rule_of_five"] is True
    assert desc["pains"]["has_pains_alert"] is False
    assert desc["pains"]["status_ar"] == "لم ترصد تنبيهات PAINS المعروفة"


def test_pains_detection():
    desc_pains = calculate_descriptors("O=C1CSC(=S)N1Cc1ccc(O)c(O)c1")
    assert desc_pains["pains"]["has_pains_alert"] is True
    assert desc_pains["pains"]["alert_description"] is not None


def test_render_2d_structure_png():
    mol, _, _ = validate_smiles("N#C[C@@H]1CCCN1C(=O)CNC12CC3CC(CC(O)(C3)C1)C2")
    png_bytes = render_mol_png(mol)
    assert isinstance(png_bytes, bytes)
    assert png_bytes[:8] == b"\x89PNG\r\n\x1a\n"
    assert get_molecular_formula(mol) == "C17H25N3O2"

    with pytest.raises(ValueError):
        render_mol_png(None)

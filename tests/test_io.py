import pytest

from autocycle.io_spec import from_reaction_smiles, load_yaml
from autocycle.spec import SpecError, drop_side

SPEC = "examples/formose_gain.yaml"


def test_missing_key(tmp_path):
    p = tmp_path / "bad.yaml"
    p.write_text("title: x\nnodes: [C]\n")
    with pytest.raises(SpecError, match="missing 'steps'"):
        load_yaml(p)


def test_not_a_mapping(tmp_path):
    p = tmp_path / "bad.yaml"
    p.write_text("- a\n- b\n")
    with pytest.raises(SpecError, match="mapping"):
        load_yaml(p)


def test_from_reaction_smiles_assigns_side_species():
    c = from_reaction_smiles(
        ["OCC=O.C=O>>OCC(O)C=O", "OCC(O)C=O>>OCC=O.C=O"],
        ["OCC=O", "OCC(O)C=O"],
    )
    assert [s.smiles for s in c.steps[0].consumes] == ["C=O"]
    assert [s.smiles for s in c.steps[1].produces] == ["C=O"]


def test_from_reaction_smiles_reports_the_missing_link():
    with pytest.raises(SpecError, match="no reaction converts"):
        from_reaction_smiles(["OCC=O>>OCC(O)C=O"], ["OCC=O", "OCC(O)C=O"])


def test_reaction_needs_double_arrow():
    with pytest.raises(SpecError, match="needs '>>'"):
        from_reaction_smiles(["OCC=O -> C"], ["OCC=O", "C"])


def test_drop_side_removes_water_everywhere():
    c = drop_side(load_yaml(SPEC), ["O"])
    for group in [c.steps] + [s.steps for s in c.subs]:
        for st in group:
            assert all(x.smiles != "O" for x in st.consumes + st.produces)


NODES = "nodes: ['OCC=O','OCC(O)C=O']\nsteps:\n"


def test_a_mistyped_key_is_rejected(tmp_path):
    """`flux` and `count` decide the stoichiometry, so a typo must not pass silently."""
    p = tmp_path / "bad.yaml"
    p.write_text(NODES + "  - {id: a, flx: 2}\n  - {id: b}\n")
    with pytest.raises(SpecError, match="unknown key"):
        load_yaml(p)


def test_a_step_without_an_id_is_named_not_a_traceback(tmp_path):
    p = tmp_path / "bad.yaml"
    p.write_text(NODES + "  - {rule: x}\n  - {id: b}\n")
    with pytest.raises(SpecError, match="needs an 'id'"):
        load_yaml(p)

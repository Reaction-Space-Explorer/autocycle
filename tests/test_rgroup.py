import pytest

from autocycle.depict import mol_svg
from autocycle.io_spec import load_yaml
from autocycle.render import render
from autocycle.spec import Mol, Side, SpecError, parse_smiles

ARM = "examples/canonical/malyl_coa_arm.yaml"


def test_pseudo_atom_is_parsed_and_named():
    smi, label = parse_smiles("CC(=O)S[CoA]")
    assert label == "CoA"
    assert smi == "*SC(C)=O"


def test_plain_smiles_reports_no_rgroup():
    assert parse_smiles("OCC=O")[1] is None


def test_invalid_smiles_still_raises():
    with pytest.raises(SpecError, match="bad SMILES"):
        parse_smiles("not a molecule")


def test_a_bracket_that_is_a_real_atom_is_untouched():
    """The fallback must only fire when the SMILES genuinely fails to parse."""
    smi, label = parse_smiles("C[N+](C)(C)C")
    assert label is None
    assert "N+" in smi


def test_mol_and_side_carry_the_label():
    assert Mol("CC(=O)S[CoA]").rgroup == "CoA"
    assert Side("CC(=O)S[CoA]").rgroup == "CoA"


def test_an_explicit_rgroup_wins():
    assert Mol("CC(=O)S[CoA]", rgroup="R").rgroup == "R"


def test_the_label_reaches_the_depiction():
    with_label = mol_svg("*SC(C)=O", rgroup="CoA")
    without = mol_svg("*SC(C)=O")
    assert with_label != without


def test_the_paper_example_labels_every_thioester():
    c = load_yaml(ARM)
    named = {m.label: m.rgroup for m in c.nodes}
    assert named["acetyl-CoA"] == "CoA"
    assert named["malyl-CoA"] == "CoA"
    assert named["oxaloacetate"] is None
    assert render(c, style="annotated").startswith("<svg")


def test_one_arm_alone_is_not_a_net_gain():
    """The published cycle is n = 2 over 11 steps; this arm alone is n = 1."""
    from autocycle.verify import SIMPLE, verify

    assert verify(load_yaml(ARM)).status == SIMPLE

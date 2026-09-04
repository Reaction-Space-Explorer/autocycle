import pytest

from autocycle.spec import Cycle, Mol, Side, SpecError, Step


def _steps(n):
    return [Step(f"r{i}") for i in range(n)]


def _nodes(n):
    return [Mol("C" * (i + 1)) for i in range(n)]


def test_count_mismatch():
    with pytest.raises(SpecError, match="one step per node"):
        Cycle(nodes=_nodes(3), steps=_steps(2))


def test_duplicate_reaction_ids():
    with pytest.raises(SpecError, match="duplicate reaction ids"):
        Cycle(nodes=_nodes(2), steps=[Step("r1"), Step("r1")])


def test_bad_smiles():
    with pytest.raises(SpecError, match="bad SMILES"):
        Mol("this is not smiles")


def test_smiles_canonicalised():
    assert Mol("OCC=O").smiles == Mol("O=CCO").smiles


def test_total_dg_is_none_if_any_step_lacks_one():
    c = Cycle(nodes=_nodes(2), steps=[Step("a", dg=-1.0), Step("b")])
    assert c.total_dg is None


def test_total_dg_sums_when_complete():
    c = Cycle(nodes=_nodes(2), steps=[Step("a", dg=-1.5), Step("b", dg=2.0)])
    assert c.total_dg == pytest.approx(0.5)


def test_gain_steps():
    c = Cycle(nodes=_nodes(3), steps=[Step("a"), Step("b", gain=True), Step("c")])
    assert c.gain_steps == [1]


def test_reversible_flag():
    assert not Step("a").reversible
    assert Step("a", rev_mag=0.2).reversible


def test_side_count_must_be_positive():
    with pytest.raises(SpecError, match="count must be"):
        Side("O", count=0)


def test_seed_out_of_range():
    with pytest.raises(SpecError, match="seed"):
        Cycle(nodes=_nodes(2), steps=_steps(2), seed=5)

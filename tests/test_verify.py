import pytest

from autocycle.io_spec import load_yaml
from autocycle.spec import Cycle, Mol, Side, Step
from autocycle.verify import AUTOCATALYTIC, CANDIDATE, INCOMPLETE, SIMPLE, verify


def _cycle(seed=0, consumes=(), produces=()):
    nodes = [Mol("OCC=O"), Mol("OCC(O)C=O"), Mol("OCC(O)C(O)C=O")]
    steps = [Step(f"r{i}") for i in range(3)]
    for i, group in enumerate(consumes):
        steps[i].consumes = [Side(s) for s in group]
    for i, group in enumerate(produces):
        steps[i].produces = [Side(*g) if isinstance(g, tuple) else Side(g) for g in group]
    return Cycle(nodes=nodes, steps=steps, seed=seed)


def test_extra_yield_of_the_seed_is_autocatalytic():
    c = _cycle(consumes=[["C=O"], [], []], produces=[[], [], ["OCC=O"]])
    v = verify(c)
    assert v.conditions["extra_yield"] == "yes"
    assert v.status == AUTOCATALYTIC
    assert v.seed_yield == 2.0


def test_no_extra_copy_is_a_candidate_not_a_simple_cycle():
    """Absent coefficients the yield is unknown, so it must not be called simple."""
    c = _cycle(consumes=[["C=O"], [], []])
    v = verify(c)
    assert v.conditions["extra_yield"] == "unknown"
    assert v.status == CANDIDATE
    assert v.status != SIMPLE


def test_a_cycle_with_no_feeder_is_incomplete():
    v = verify(_cycle())
    assert v.conditions["feeder"] == "no"
    assert v.status == INCOMPLETE


def test_no_seed_is_incomplete():
    c = _cycle(consumes=[["C=O"], [], []])
    object.__setattr__(c, "seed", None)
    v = verify(c)
    assert v.conditions["seed_identified"] == "no"
    assert v.status == INCOMPLETE
    assert v.seed_yield is None


def test_stoichiometric_count_is_honoured():
    c = _cycle(consumes=[["C=O"], [], []], produces=[[], [], [("OCC=O", 2)]])
    assert verify(c).seed_yield == 3.0


def test_ring_species_do_not_count_as_feeders():
    """A species already on the ring is not an external feeder."""
    c = _cycle(consumes=[["OCC(O)C=O"], [], []])
    assert verify(c).conditions["feeder"] == "no"


def test_outlet_is_reported_separately():
    c = _cycle(consumes=[["C=O"], [], []], produces=[[], ["O=C=O"], []])
    v = verify(c)
    assert v.conditions["outlet"] == "yes"
    assert verify(_cycle(consumes=[["C=O"], [], []])).conditions["outlet"] == "no"


def test_declared_gain_without_support_is_flagged():
    c = _cycle(consumes=[["C=O"], [], []])
    c.steps[2].gain = True
    v = verify(c)
    assert v.declared_gain
    assert v.disagrees_with_declaration


def test_declared_gain_with_support_is_not_flagged():
    c = _cycle(consumes=[["C=O"], [], []], produces=[[], [], ["OCC=O"]])
    c.steps[2].gain = True
    assert not verify(c).disagrees_with_declaration


def test_hand_written_example_is_autocatalytic():
    v = verify(load_yaml("examples/formose_gain.yaml"))
    assert v.status == AUTOCATALYTIC
    assert "extra_yield=yes" in v.summary()


def test_summary_lists_every_condition():
    s = verify(load_yaml("examples/formose_gain.yaml")).summary()
    for cond in ("seed_identified", "feeder", "outlet", "seed_regenerated", "extra_yield"):
        assert cond in s


def test_a_shunt_carries_the_topological_criterion():
    from autocycle.spec import Shunt
    from autocycle.verify import TOPOLOGICAL

    c = _cycle(consumes=[["C=O"], [], []])
    c.shunt = Shunt(from_node=1, steps=[Step("s1")])
    v = verify(c)
    assert v.conditions["shunt"] == "yes"
    assert v.conditions["extra_yield"] == "unknown"
    assert v.status == TOPOLOGICAL


def test_a_stated_coefficient_outranks_a_shunt():
    from autocycle.spec import Shunt

    c = _cycle(consumes=[["C=O"], [], []], produces=[[], [], ["OCC=O"]])
    c.shunt = Shunt(from_node=1, steps=[Step("s1")])
    assert verify(c).status == AUTOCATALYTIC


def test_no_shunt_is_reported_as_such():
    assert verify(_cycle(consumes=[["C=O"], [], []])).conditions["shunt"] == "no"


def test_shunt_out_of_range_rejected():
    from autocycle.spec import Cycle, Shunt, SpecError

    with pytest.raises(SpecError, match="shunt from_node"):
        Cycle(
            nodes=[Mol("OCC=O"), Mol("OCC(O)C=O")],
            steps=[Step("a"), Step("b")],
            shunt=Shunt(from_node=9, steps=[Step("s")]),
        )

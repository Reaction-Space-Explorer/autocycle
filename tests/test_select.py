import pytest

from autocycle.io_spec import load_yaml
from autocycle.select import (
    consumers,
    feeders,
    min_feeder_weight,
    rank_key,
    role_violations,
)
from autocycle.spec import Cycle, Mol, Side, Step, canonical


def _cycle(nodes, consumes=()):
    steps = [Step(f"r{i}") for i in range(len(nodes))]
    for i, group in enumerate(consumes):
        steps[i].consumes = [Side(s) for s in group]
    return Cycle(nodes=[Mol(n) for n in nodes], steps=steps, seed=0)


def test_feeders_and_consumers_are_distinct_sets():
    c = load_yaml("examples/formose_gain.yaml")
    assert feeders(c) == {canonical("C=O")}
    # the main ring only; the sub-cycle's water is not counted here
    assert consumers(c) == {canonical("O=C=O"), canonical("OCC=O")}


def test_repeated_feeder_counts_once():
    c = _cycle(["OCC=O", "OCC(O)C=O"], consumes=[["C=O"], ["C=O"]])
    assert len(feeders(c)) == 1


def test_min_feeder_weight_picks_the_lightest():
    c = _cycle(["OCC=O", "OCC(O)C=O"], consumes=[["C=O", "OCC=O"], []])
    assert min_feeder_weight(c) == pytest.approx(30.026, abs=0.01)  # formaldehyde


def test_no_feeders_ranks_last_on_weight():
    assert min_feeder_weight(_cycle(["OCC=O", "OCC(O)C=O"])) == float("inf")


def test_rank_prefers_fewer_feeders_then_lighter():
    one = _cycle(["OCC=O", "OCC(O)C=O"], consumes=[["OCC=O"], []])
    two = _cycle(["OCC=O", "OCC(O)C=O"], consumes=[["C=O"], ["OCC=O"]])
    light = _cycle(["OCC=O", "OCC(O)C=O"], consumes=[["C=O"], []])
    assert sorted([two, one, light], key=rank_key) == [light, one, two]


def test_restricted_molecule_on_the_ring_is_flagged():
    c = _cycle(["C=O", "OCC=O", "OCC(O)C=O"])
    assert role_violations(c) == [canonical("C=O")]


def test_restricted_molecule_as_a_feeder_is_fine():
    c = _cycle(["OCC=O", "OCC(O)C=O"], consumes=[["C=O"], []])
    assert role_violations(c) == []


def test_methanol_on_the_ring_is_flagged():
    assert role_violations(_cycle(["CO", "OCC=O"])) == [canonical("CO")]


def test_restricted_set_is_configurable():
    c = _cycle(["OCC=O", "OCC(O)C=O"])
    assert role_violations(c, restricted=("OCC=O",)) == [canonical("OCC=O")]

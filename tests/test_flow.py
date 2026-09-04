"""Coefficient-based verdicts read from MOD flow summaries."""

import pytest

from autocycle.io_flow import AUTOCATALYTIC, FlowSolution, read_flow_summary

SAMPLE = "examples/flow_sample/acetyl_coa_flow.txt"


@pytest.fixture(scope="module")
def sols():
    return read_flow_summary(SAMPLE)


def test_reads_every_solution(sols):
    assert len(sols) == 2
    assert all(s.objective == 11 for s in sols)


def test_target_is_the_molecule_fed_in_and_returned_in_excess(sols):
    for s in sols:
        assert s.target == "ACETYL COA"
        assert s.balance["ACETYL COA"][:2] == (1, 2)
        assert s.gain == 1


def test_cofactors_are_not_mistaken_for_the_target(sols):
    """Several cofactors are produced net, but are not fed in."""
    net = sols[0].net_produced
    assert "ACETYL COA" in net
    assert any(m in net for m in ("Pi", "AMP", "NAD+"))
    assert sols[0].target == "ACETYL COA"


def test_coefficients_settle_the_criterion(sols):
    assert all(s.verdict() == AUTOCATALYTIC for s in sols)


def test_a_balanced_solution_is_not_autocatalytic():
    s = FlowSolution(objective=3, balance={"A": (1, 1, 0), "B": (0, 1, 1)})
    assert s.target is None
    assert s.gain is None
    assert s.verdict() != AUTOCATALYTIC


def test_two_candidate_targets_are_not_resolved():
    s = FlowSolution(objective=3, balance={"A": (1, 2, 1), "B": (1, 3, 2)})
    assert s.target is None


def test_empty_text_yields_nothing(tmp_path):
    p = tmp_path / "x.txt"
    p.write_text("no tables here\n")
    assert read_flow_summary(p) == []

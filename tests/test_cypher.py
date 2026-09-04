import contextlib

import pytest

from autocycle.io_cypher import from_cypher_row, read_rows
from autocycle.render import render
from autocycle.spec import SpecError

CSV = "examples/cypher_sample.csv"


@pytest.fixture(scope="module")
def rows():
    return read_rows(CSV)


def _cycles(rows):
    out = []
    for r in rows:
        with contextlib.suppress(SpecError):
            out.append(from_cypher_row(r))
    return out


def test_sample_splits_into_simple_and_pinched(rows):
    good = _cycles(rows)
    assert len(good) == 4
    assert len(rows) - len(good) == 2


def test_pinched_rings_are_rejected_by_name(rows):
    bad = [r for r in rows if _is_pinched(r)]
    assert bad
    with pytest.raises(SpecError, match="pinched ring path"):
        from_cypher_row(bad[0])


def _is_pinched(row):
    try:
        from_cypher_row(row)
    except SpecError:
        return True
    return False


def test_seed_is_the_begin_molecule(rows):
    import ast

    for r in rows:
        try:
            c = from_cypher_row(r)
        except SpecError:
            continue
        begin = ast.literal_eval(r["beginMol"])["smiles_str"]
        from autocycle.spec import canonical

        assert c.nodes[c.seed].smiles == canonical(begin)


def test_gain_step_arrives_at_the_seed(rows):
    for c in _cycles(rows):
        assert c.gain_steps == [(c.seed - 1) % len(c.nodes)]
        assert c.steps[c.gain_steps[0]].mag > c.steps[c.seed].mag


def test_side_species_are_attached(rows):
    cycles = _cycles(rows)
    assert any(
        s.consumes or s.produces for c in cycles for s in c.steps
    ), "feeder/consumer/export should land on some step"


def test_no_thermodynamics_in_this_format(rows):
    for c in _cycles(rows):
        assert all(s.dg is None for s in c.steps)
        assert c.total_dg is None  # None, not a partial sum


def test_every_sample_renders(rows):
    for c in _cycles(rows):
        assert render(c, style="paper").startswith("<svg")


def test_missing_column_is_named():
    with pytest.raises(SpecError, match="missing column 'ringMols'"):
        from_cypher_row({"beginMol": "{}"})

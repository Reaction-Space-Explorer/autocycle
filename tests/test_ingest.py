import pytest

from autocycle import ingest
from autocycle.spec import SpecError

CSV = "examples/edges.csv"


def test_missing_required_column(tmp_path):
    p = tmp_path / "e.csv"
    p.write_text("from,to\nC,CC\n")
    with pytest.raises(SpecError, match="missing column"):
        ingest.read_edges(p)


def test_empty_file(tmp_path):
    p = tmp_path / "e.csv"
    p.write_text("source,target\n")
    with pytest.raises(SpecError, match="no rows"):
        ingest.read_edges(p)


def test_unknown_dg_is_none_not_false():
    g = ingest.spontaneous(ingest.read_edges(CSV))
    for _, _, d in g.edges(data=True):
        assert d["spontaneous"] is (None if d["dg"] is None else d["dg"] < 0)


def test_cycle_length_filter():
    g = ingest.read_edges(CSV)
    assert ingest.find_cycles(g, 3, 12) == [c for c in ingest.find_cycles(g, 2, 12) if len(c) >= 3]
    assert ingest.find_cycles(g, 99, 100) == []


def test_to_cycle_rejects_a_non_cycle():
    g = ingest.read_edges(CSV)
    with pytest.raises(SpecError, match="not a cycle"):
        ingest.to_cycle(g, ["OCC=O", "O=C=O"])


def test_to_cycle_marks_uphill_steps_filtered_not_dropped():
    g = ingest.spontaneous(ingest.read_edges(CSV))
    ring = ingest.find_cycles(g, 3, 12)[0]
    c = ingest.to_cycle(g, ring)
    assert len(c.steps) == len(ring)  # nothing dropped
    assert any(s.filtered for s in c.steps)  # the +5.3 step is flagged


def test_summary_csv_has_a_row_per_cycle(tmp_path):
    g = ingest.spontaneous(ingest.read_edges(CSV))
    cycles = ingest.find_cycles(g, 3, 12)
    out = tmp_path / "s.csv"
    ingest.summary_csv(g, cycles, out)
    lines = out.read_text().strip().splitlines()
    assert len(lines) == len(cycles) + 1
    assert lines[0].startswith("index,length,total_dg")


def test_gain_at_sets_the_gain_step():
    g = ingest.read_edges(CSV)
    ring = ingest.find_cycles(g, 3, 12)[0]
    c = ingest.to_cycle(g, ring, gain_at=2)
    assert c.gain_steps == [2]
    assert c.steps[2].mag > c.steps[0].mag

from autocycle import layout as L
from autocycle.check import boxes, collisions
from autocycle.io_cypher import from_cypher_row, read_rows
from autocycle.io_spec import load_yaml
from autocycle.style import PAPER


def test_example_figure_has_no_overlapping_depictions():
    assert collisions(load_yaml("examples/formose_gain.yaml")) == []


def test_real_cycles_have_no_overlapping_depictions():
    import contextlib

    from autocycle.spec import SpecError

    for row in read_rows("examples/cypher_sample.csv"):
        with contextlib.suppress(SpecError):
            c = from_cypher_row(row)
            assert collisions(c) == [], c.nodes[0].smiles


def test_fused_molecule_is_counted_once():
    c = load_yaml("examples/formose_gain.yaml")
    n_side = sum(len(s.consumes) + len(s.produces) for s in c.steps)
    n_side += sum(len(s.consumes) + len(s.produces) for sub in c.subs for s in sub.steps)
    expected = len(c.nodes) + sum(len(s.nodes) - 1 for s in c.subs) + n_side
    assert len(boxes(c)) == expected


def test_collisions_are_found_when_the_ring_is_too_tight(monkeypatch):
    """The checker must actually be able to fail."""
    c = load_yaml("examples/formose_gain.yaml")
    monkeypatch.setattr(L, "MOL_FRAC", 1.2)  # depictions far wider than the gap
    assert collisions(c, PAPER)


def test_verify_command_reports_status_and_balance(capsys):
    from autocycle.cli import main

    assert main(["verify", "examples/canonical/formose_core.yaml"]) == 0
    out = capsys.readouterr().out
    assert "autocatalytic" in out
    assert "atoms      balanced" in out
    assert "one extreme current" in out


def test_verify_command_fails_on_an_unreturned_intermediate(tmp_path, capsys):
    from autocycle.cli import main

    spec = tmp_path / "bad.yaml"
    spec.write_text(
        "title: broken\nseed: 0\n"
        "nodes:\n  - {smiles: 'OCC=O'}\n  - {smiles: 'OCC(O)C=O'}\n"
        "steps:\n  - {id: r0}\n"
        "  - {id: r1, consumes: [{smiles: 'OCC(O)C=O'}]}\n"
    )
    assert main(["verify", str(spec)]) == 1
    assert "imbalanced" in capsys.readouterr().out

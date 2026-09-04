"""Solution 0 of Fig. 3 of Abel et al., transcribed from the authors' MOD flow output."""

import pytest

from autocycle.io_spec import load_yaml
from autocycle.render import render
from autocycle.spec import SpecError
from autocycle.verify import AUTOCATALYTIC, verify

SPEC = "examples/canonical/acetyl_coa_sol0.yaml"


@pytest.fixture(scope="module")
def cycle():
    return load_yaml(SPEC)


def test_step_count_matches_the_published_objective_value(cycle):
    assert len(cycle.steps) + len(cycle.shunt.steps) == 11


def test_seed_yield_matches_the_published_in_out(cycle):
    """The flow output gives ACETYL COA In 1, Out 2."""
    v = verify(cycle)
    assert v.seed_yield == 2.0
    assert v.status == AUTOCATALYTIC


def test_the_extra_copy_comes_from_the_shunt(cycle):
    """Without counting shunt products the yield would look like 1."""
    ring_only = sum(
        sp.count for st in cycle.steps for sp in st.produces
        if sp.smiles == cycle.nodes[cycle.seed].smiles
    )
    shunt_only = sum(
        sp.count for st in cycle.shunt.steps for sp in st.produces
        if sp.smiles == cycle.nodes[cycle.seed].smiles
    )
    assert ring_only == 0
    assert shunt_only == 1


def test_shunt_leaves_malyl_coa(cycle):
    assert cycle.nodes[cycle.shunt.from_node].label == "malyl-CoA"
    assert [m.label for m in cycle.shunt.nodes] == [
        "p0,12", "succinyl-CoA", "p0,5", "p0,6", "citryl-CoA"
    ]


def test_the_doubled_flux_segment_is_marked(cycle):
    """Oxaloacetate to malyl-CoA carries flux 2 in the source."""
    by_id = {s.rid: s for s in cycle.steps}
    assert by_id["r2"].mag == 2.0
    assert by_id["r23"].mag == 2.0
    assert by_id["r16"].mag == 1.0


def test_every_thioester_is_labelled_coa(cycle):
    thio = [m for m in cycle.nodes + cycle.shunt.nodes if "CoA" in (m.label or "")]
    assert len(thio) == 5
    assert all(m.rgroup == "CoA" for m in thio)


def test_glyoxylate_is_internal_not_a_feeder(cycle):
    """It is produced in the ring and consumed in the shunt."""
    gly = "O=CC(=O)O"
    assert any(sp.smiles == gly for st in cycle.steps for sp in st.produces)
    assert any(sp.smiles == gly for st in cycle.shunt.steps for sp in st.consumes)


def test_the_shunt_intermediates_are_drawn(cycle):
    svg = render(cycle, style="annotated")
    assert svg.startswith("<svg")
    for mol in cycle.shunt.nodes:
        assert f">{mol.label}<" in svg
    for step in cycle.shunt.steps:
        assert f">r:{step.rid}<" in svg
    # 11 depictions: 5 ring molecules, 5 shunt intermediates, plus side species
    assert svg.count("<svg x='") >= 10


def test_a_shunt_without_from_node_is_named(tmp_path):
    p = tmp_path / "s.yaml"
    p.write_text(
        "nodes: [OCC=O, OCC(O)C=O]\nsteps:\n  - {id: a}\n  - {id: b}\nshunt: {steps: []}\n"
    )
    with pytest.raises(SpecError, match="shunt needs 'from_node'"):
        load_yaml(p)

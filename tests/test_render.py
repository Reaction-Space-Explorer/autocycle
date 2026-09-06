import re
import xml.etree.ElementTree as ET

import pytest

from autocycle import layout as L
from autocycle.io_spec import load_yaml
from autocycle.render import render
from autocycle.style import PAPER, RICH, get

SPEC = "examples/formose_gain.yaml"


@pytest.fixture(scope="module")
def cycle():
    return load_yaml(SPEC)


def test_output_is_well_formed_xml(cycle):
    ET.fromstring(render(cycle))


def test_one_flow_arrow_per_step(cycle):
    svg = render(cycle, style="paper")
    n_steps = len(cycle.steps) + sum(len(s.steps) for s in cycle.subs)
    n_rev = sum(1 for s in cycle.steps if s.reversible)
    n_gain = len(cycle.gain_steps)
    # each step draws one arrow; reversible steps add a second, gain steps a halo band
    assert svg.count("<path d='M") >= n_steps + n_rev + n_gain


def test_every_molecule_is_depicted_once(cycle):
    svg = render(cycle, style="paper")
    n_side = sum(len(s.consumes) + len(s.produces) for s in cycle.steps)
    n_side += sum(len(s.consumes) + len(s.produces) for sub in cycle.subs for s in sub.steps)
    # ring molecules + side species; the fused molecule is shared, so drawn once
    expected = len(cycle.nodes) + sum(len(s.nodes) - 1 for s in cycle.subs) + n_side
    assert len(re.findall(r"<svg x='", svg)) == expected


def test_nothing_is_drawn_outside_the_viewbox(cycle):
    """Regression: sub-rings and side species used to be cropped."""
    svg = render(cycle, style="paper")
    m = re.search(r"translate\((-?[\d.]+),(-?[\d.]+)\)", svg)
    tx, ty = float(m.group(1)), float(m.group(2))
    w, h = (float(v) for v in re.search(r"viewBox='0 0 ([\d.]+) ([\d.]+)'", svg).groups())
    scale = float(re.search(r"scale\(([\d.]+)\)", svg).group(1))
    # nested depictions: check the whole box, not just its top-left corner
    boxes = [
        (float(x), float(y), float(bw), float(bh))
        for x, y, bw, bh in re.findall(
            r"<svg x='(-?[\d.]+)' y='(-?[\d.]+)' width='([\d.]+)' height='([\d.]+)'", svg
        )
    ]
    assert boxes
    for x, y, bw, bh in boxes:
        assert (x + tx) * scale >= 0 and (x + bw + tx) * scale <= w + 1e-6
        assert (y + ty) * scale >= 0 and (y + bh + ty) * scale <= h + 1e-6


def test_paper_style_has_no_node_circles_and_grey_arrows(cycle):
    svg = render(cycle, style="paper")
    assert "<circle" not in svg
    assert PAPER.ring_grey in svg


def test_rich_style_adds_circles_and_a_legend(cycle):
    svg = render(cycle, style="rich")
    assert "<circle" in svg
    assert "cycle dG" in svg
    assert RICH.ring_colour == "dg"


def test_centre_label_names_the_seed(cycle):
    centred = r"<text x='-?0\.0+' y='-?0\.0+'[^>]*>glycolaldehyde</text>"
    # paper puts the autocatalyst's name in the middle of the ring
    assert re.search(centred, render(cycle, style="paper"))
    # rich labels the node instead, so the name is present but not at the centre
    rich = render(cycle, style="rich")
    assert "glycolaldehyde" in rich
    assert not re.search(centred, rich)


def test_legend_override_beats_the_style(cycle):
    assert "cycle dG" in render(cycle, style="paper", legend=True)
    assert "cycle dG" not in render(cycle, style="rich", legend=False)


def test_unknown_style_rejected():
    with pytest.raises(ValueError, match="unknown style"):
        get("neon")


def test_figure_grows_with_cycle_length():
    small = L.lay_out(4).radius
    big = L.lay_out(12).radius
    assert big > small


def test_no_legend_drops_the_rule_list_too():
    cycle = load_yaml("examples/canonical/formose_core.yaml")
    assert "Retro-aldol cleavage" in render(cycle, style="annotated")
    assert "Retro-aldol cleavage" not in render(cycle, style="annotated", legend=False)

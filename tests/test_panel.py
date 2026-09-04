import re
import xml.etree.ElementTree as ET

import pytest

from autocycle.panel import Series, bar_panel, compose, svg


def _series():
    return [
        Series("Total", {4: 300, 5: 1800}),
        Series("Two feeders", {4: 40, 5: 200}),
        Series("One feeder", {4: 12, 5: 90}),
    ]


def test_panel_is_well_formed_svg():
    ET.fromstring(svg(bar_panel(_series()), 520, 300))


def test_one_bar_per_series_and_length():
    out = bar_panel(_series())
    # 3 series x 2 lengths, plus the legend swatches
    assert len(re.findall(r"<rect ", out)) == 6 + 3 + 1


def test_zero_counts_draw_no_bar():
    out = bar_panel([Series("Total", {4: 0, 5: 10})])
    assert len(re.findall(r"<rect x='", out)) == 1 + 1 + 1  # one bar, legend box, swatch


def test_log_axis_is_monotonic_in_value():
    out = bar_panel([Series("Total", {4: 10, 5: 1000})])
    widths = [float(w) for w in re.findall(r"<rect x='[\d.]+' y='[\d.]+' width='([\d.]+)'", out)]
    bars = sorted(widths)[-2:]
    assert bars[1] > bars[0]


def test_axis_labels_are_present():
    out = bar_panel(_series(), x_label="Frequency", y_label="Cycle length")
    assert "Frequency" in out
    assert "Cycle length" in out
    for s in _series():
        assert s.label in out


def test_tag_is_drawn_when_given():
    assert ">A<" in bar_panel(_series(), tag="A")
    assert ">A<" not in bar_panel(_series())


def test_empty_input_rejected():
    with pytest.raises(ValueError, match="no series"):
        bar_panel([])
    with pytest.raises(ValueError, match="no data"):
        bar_panel([Series("Total", {})])


def test_compose_grids_and_tags_cells():
    cells = [svg("<g/>", 100, 100) for _ in range(3)]
    out = compose(cells, 100, 100, cols=2)
    ET.fromstring(out)
    assert len(re.findall(r"<svg x='", out)) == 3
    for tag in ("A", "B", "C"):
        assert f">{tag}<" in out


def test_compose_size_follows_the_grid():
    out = compose([svg("<g/>", 100, 100)] * 4, 100, 100, cols=2, gap=10)
    w, h = (float(v) for v in re.search(r"viewBox='0 0 (\d+) (\d+)'", out).groups())
    assert w == 2 * 100 + 3 * 10
    assert h == 2 * 100 + 3 * 10


def test_compose_strips_the_nested_white_background():
    """A nested cell must not paint over the panel to its left."""
    out = compose([svg("<g/>", 100, 100)] * 2, 100, 100)
    assert out.count("width='100%' height='100%' fill='white'") == 1


def test_compose_rejects_nothing():
    with pytest.raises(ValueError, match="nothing to compose"):
        compose([], 100, 100)

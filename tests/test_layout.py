import math

import pytest

from autocycle import layout as L


def test_chord_between_adjacent_vertices_equals_pitch():
    for n_mol in (2, 3, 5, 9, 20):
        ring = L.lay_out(n_mol)
        a, b = ring.verts[0], ring.verts[1]
        assert math.hypot(a.x - b.x, a.y - b.y) == pytest.approx(L.PITCH, rel=1e-9)


def test_vertices_alternate_and_count_is_even():
    ring = L.lay_out(6)
    assert ring.n == 12
    assert ring.n % 2 == 0
    assert [v.kind for v in ring.verts[:4]] == ["mol", "rxn", "mol", "rxn"]
    assert len(ring.of("mol")) == len(ring.of("rxn")) == 6


def test_all_vertices_lie_on_the_ring():
    ring = L.lay_out(7)
    for v in ring.verts:
        assert math.hypot(v.x - ring.cx, v.y - ring.cy) == pytest.approx(ring.radius)


def test_first_molecule_is_at_top():
    assert L.lay_out(4).verts[0].angle == pytest.approx(90.0)


def test_step_arc_spans_two_vertex_steps_less_gaps():
    ring = L.lay_out(5)
    a0, a1 = L.step_arc(ring, 0, cw=True, gap=5.0)
    span = 2 * 360.0 / ring.n
    assert a0 == pytest.approx(90.0 - 5.0)
    assert a1 == pytest.approx(90.0 - span + 5.0)


def test_step_arc_direction_flips_with_cw():
    ring = L.lay_out(5)
    cw = L.step_arc(ring, 0, cw=True, gap=5.0)
    ccw = L.step_arc(ring, 0, cw=False, gap=5.0)
    assert cw[1] < cw[0]
    assert ccw[1] > ccw[0]


def test_node_gap_clears_the_circle():
    ring = L.lay_out(4)
    gap = L.node_gap(ring)
    # a point at `gap` degrees from a node must be at least a node radius away
    x, y = L.polar(ring.radius, ring.verts[0].angle - gap)
    d = math.hypot(x - ring.verts[0].x, y - ring.verts[0].y)
    assert d >= L.node_radius()


def test_subring_shares_the_vertex_and_stays_outside_parent():
    parent = L.lay_out(4)
    sub = L.lay_out_sub(parent, at_step=1, n_mol=3)
    shared = parent.verts[(2 * 2) % parent.n]
    # the shared vertex lies on the sub-ring
    assert math.hypot(shared.x - sub.cx, shared.y - sub.cy) == pytest.approx(sub.radius)
    # and the sub-ring never dips inside the parent circle
    for v in sub.verts:
        assert math.hypot(v.x - parent.cx, v.y - parent.cy) >= parent.radius - 1e-6


def test_subring_starts_at_the_shared_vertex():
    parent = L.lay_out(5)
    sub = L.lay_out_sub(parent, at_step=0, n_mol=4)
    shared = parent.verts[2]
    assert (sub.verts[0].x, sub.verts[0].y) == pytest.approx((shared.x, shared.y))


def test_side_anchor_respects_keep_out_disc():
    ring = L.lay_out(3, centre=(0.0, -3.0))
    v = ring.verts[1]
    avoid = (0.0, 0.0, 2.6)
    x, y = L.side_anchor(ring, v, "in", avoid=avoid)
    assert math.hypot(x, y) >= 2.6


def test_side_anchor_rejects_bad_side():
    with pytest.raises(ValueError, match="'in' or 'out'"):
        L.side_anchor(L.lay_out(3), L.lay_out(3).verts[1], "sideways")


def test_ring_radius_needs_two_molecules():
    with pytest.raises(ValueError, match=">= 2"):
        L.ring_radius(1)


def test_a_shunt_clears_the_ring_side_species():
    """The shunt is an outer arc, so it must not run through the ring's own feeders."""
    from autocycle.io_spec import load_yaml

    c = load_yaml("examples/canonical/acetyl_coa_sol0.yaml")
    ring = L.lay_out(len(c.nodes))
    half = L.mol_half(ring)
    reach = L.side_reach(ring, c.steps, half)
    r, _, _ = L.shunt_arc(ring, c.shunt.from_node, c.seed, len(c.shunt.nodes), half, reach)
    assert reach > 0
    assert r - half > reach

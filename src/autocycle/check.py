"""Geometric checks, used by the tests and the bench commands."""

from __future__ import annotations

from autocycle import layout as L
from autocycle import tree as T
from autocycle.spec import Cycle, Pathway
from autocycle.style import PAPER, Style


def boxes(cycle, st: Style = PAPER) -> list[tuple[str, float, float, float]]:
    """(label, x, y, half) per depiction."""
    if isinstance(cycle, Pathway):
        return _route_boxes(cycle, st)
    ring = L.lay_out(len(cycle.nodes))
    subs = [(L.lay_out_sub(ring, s.at_step, len(s.nodes)), s) for s in cycle.subs]
    keep_out = (ring.cx, ring.cy, ring.radius + L.node_radius() * 2.4)

    out = []
    groups = [(ring, cycle.steps, cycle.nodes, None)]
    groups += [(sr, s.steps, s.nodes, keep_out) for sr, s in subs]
    for r, steps, nodes, avoid in groups:
        half = L.mol_half(r) * st.mol_scale
        skip = 0 if avoid is not None else -1  # a fused molecule is drawn once
        for v in r.of("mol"):
            if v.index == skip:
                continue
            out.append((nodes[v.index].smiles, v.x, v.y, half))
        for _, _, sp, (x, y) in L.side_points(r, steps, avoid=avoid):
            out.append((sp.smiles, x, y, half))
    return out


def collisions(cycle: Cycle, st: Style = PAPER, tol: float = 0.25):
    """Depiction pairs overlapping by more than `tol` of the smaller box."""
    bs = boxes(cycle, st)
    hits = []
    for i in range(len(bs)):
        for j in range(i + 1, len(bs)):
            a, b = bs[i], bs[j]
            ox = min(a[1] + a[3], b[1] + b[3]) - max(a[1] - a[3], b[1] - b[3])
            oy = min(a[2] + a[3], b[2] + b[3]) - max(a[2] - a[3], b[2] - b[3])
            if ox <= 0 or oy <= 0:
                continue
            frac = (ox * oy) / (2 * min(a[3], b[3])) ** 2
            if frac > tol:
                hits.append((a[0], b[0], round(frac, 3)))
    return sorted(hits, key=lambda h: -h[2])


def _route_boxes(pw: Pathway, st: Style) -> list[tuple[str, float, float, float]]:
    lay = T.lay_out_pathway(pw)
    half = T.MOL_HALF * st.mol_scale
    side = T.MOL_HALF * 0.62 * st.mol_scale
    out = [(n.mol.smiles, lay.mol(n).x, lay.mol(n).y, half) for n in pw.nodes]
    for node in pw.nodes:
        if not node.step:
            continue
        r = lay.rxn(node)
        for s, group in (("in", node.step.consumes), ("out", node.step.produces)):
            for k, sp in enumerate(group):
                x, y = T.side_anchor(r, s, k)
                out.append((sp.smiles, x, y, side))
    return out


def collisions_route(pw: Pathway, st: Style = PAPER, tol: float = 0.25):
    return collisions(pw, st, tol)


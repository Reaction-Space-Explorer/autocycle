"""Assemble a cycle into one SVG."""

from __future__ import annotations

from autocycle import draw as _draw
from autocycle import layout as L
from autocycle import tree as T
from autocycle.depict import fit_bond, formula
from autocycle.draw import (
    GAIN,
    _p,
    draw_nodes,
    draw_ring,
    draw_route,
    draw_route_nodes,
    draw_route_sides,
    draw_shunt,
    draw_sides,
    text,
)
from autocycle.encode import dg_colour, dg_span
from autocycle.spec import Cycle, Pathway
from autocycle.style import Style, get

PAD = 0.55
TARGET_PX = 1100.0
ROUTE_PX_PER_COLUMN = 250.0  # a deep route gets a wider canvas, not smaller structures


def render(obj, mode: str = "linear", style: str | Style = "paper", legend=None) -> str:
    """Render a Cycle or a Pathway."""
    if isinstance(obj, Pathway):
        return _render_pathway(obj, mode, style, legend)
    return _render_cycle(obj, mode, style, legend)


def _render_pathway(pw: Pathway, mode: str, style: str | Style, legend) -> str:
    st = get(style) if isinstance(style, str) else style
    show_legend = st.legend if legend is None else legend
    _draw.BOND[0] = fit_bond(_route_smiles(pw))

    lay = T.lay_out_pathway(pw)
    span = dg_span([s.dg for s in pw.steps])

    body = draw_route(pw, lay, span, st, mode)
    body += draw_route_sides(pw, lay, st)
    body += draw_route_nodes(pw, lay, st)

    half = T.MOL_HALF * 0.62 * st.mol_scale
    pts = []
    for node in pw.nodes:
        if not node.step:
            continue
        r = lay.rxn(node)
        for side, group in (("in", node.step.consumes), ("out", node.step.produces)):
            for k in range(len(group)):
                x, y = T.side_anchor(r, side, k)
                pts += [(x + half, y + half), (x - half, y - half)]

    x0, y0, x1, y1 = T.bounds(lay, pts, PAD + 0.35)
    if show_legend:
        y0 -= 0.72
    target = max(TARGET_PX, ROUTE_PX_PER_COLUMN * (pw.root.depth + 1))
    return _wrap(body, x0, y0, x1, y1, pw, span, mode, show_legend, target)


def _render_cycle(cycle: Cycle, mode: str, style: str | Style, legend) -> str:
    st = get(style) if isinstance(style, str) else style
    show_legend = st.legend if legend is None else legend
    _draw.BOND[0] = fit_bond(_cycle_smiles(cycle))

    ring = L.lay_out(len(cycle.nodes))
    subs = [(L.lay_out_sub(ring, s.at_step, len(s.nodes)), s) for s in cycle.subs]
    span = dg_span([s.dg for s in list(cycle.steps) + [x for sc in cycle.subs for x in sc.steps]])

    keep_out = (ring.cx, ring.cy, ring.radius + L.node_radius() * 2.4)
    out0 = L.side_out(ring, st.side_out)
    anchors = [(ring, cycle.steps, L.side_points(ring, cycle.steps, out0=out0))]
    anchors += [
        (sr, s.steps, L.side_points(sr, s.steps, avoid=keep_out, out0=L.side_out(sr, st.side_out)))
        for sr, s in subs
    ]

    body = draw_ring(ring, cycle.steps, span, st, mode)
    body += draw_shunt(ring, cycle, st)
    for sr, sub in subs:
        body += draw_ring(sr, sub.steps, span, st, mode, cw=False)
    for r, steps, anc in anchors:
        body += draw_sides(r, steps, st, anc)
    body += draw_nodes(ring, cycle.nodes, cycle.steps, cycle.seed, st)
    for sr, sub in subs:
        body += draw_nodes(sr, sub.nodes, sub.steps, None, st, cw=False, skip_mols=(0,))

    if st.centre_label and cycle.seed is not None:
        seed_mol = cycle.nodes[cycle.seed]
        # a raw SMILES is not a label
        name = seed_mol.label or formula(seed_mol.smiles)
        cx, cy = _p(ring.cx, ring.cy)
        body.append(text(cx, cy, name, st.centre_size, "middle", "#111"))

    # pad by the depiction each anchor carries
    pad = []
    for r, _, anc in anchors:
        h = L.mol_half(r) * st.mol_scale
        pad += [(x + dx * h, y + dy * h) for _, _, _, (x, y) in anc for dx, dy in ((1, 1), (-1, -1))]
    if cycle.shunt is not None and cycle.seed is not None:
        h = L.mol_half(ring) * st.mol_scale
        rr, _, _ = L.shunt_arc(ring, cycle.shunt.from_node, cycle.seed,
                               len(cycle.shunt.nodes), h)
        rr += h * 1.2
        pad += [
            (ring.cx + L.polar(rr, ang)[0], ring.cy + L.polar(rr, ang)[1])
            for ang in range(0, 360, 10)
        ]
    x0, y0, x1, y1 = L.bounds([ring] + [sr for sr, _ in subs], pad, PAD)
    rules = _rule_lines(cycle) if st.rule_legend else []
    if rules:
        y0 -= 0.30 + 0.155 * len(rules)
        body += [
            text(x0 + 0.3, -(y0 + 0.30 + 0.155 * (len(rules) - 1 - i)), line,
                 st.label_size, "start", "#222")
            for i, line in enumerate(rules)
        ]
    if show_legend:
        y0 -= 0.72
    return _wrap(body, x0, y0, x1, y1, cycle, span, mode, show_legend)


def _cycle_smiles(cycle) -> list[str]:
    out = [m.smiles for m in cycle.nodes if getattr(m, "structure", True)]
    groups = [cycle.steps] + [s.steps for s in cycle.subs]
    if cycle.shunt:
        groups.append(cycle.shunt.steps)
        out += [m.smiles for m in cycle.shunt.nodes if getattr(m, "structure", True)]
    for sub in cycle.subs:
        out += [m.smiles for m in sub.nodes if getattr(m, "structure", True)]
    for g in groups:
        for st_ in g:
            out += [sp.smiles for sp in st_.consumes + st_.produces
                    if getattr(sp, "structure", True)]
    return out


def _route_smiles(pw) -> list[str]:
    out = [n.mol.smiles for n in pw.nodes if getattr(n.mol, "structure", True)]
    for st_ in pw.steps:
        out += [sp.smiles for sp in st_.consumes + st_.produces
                if getattr(sp, "structure", True)]
    return out


def _rule_lines(cycle) -> list[str]:
    seen, out = set(), []
    every = cycle.steps + [x for sub in cycle.subs for x in sub.steps]
    every += cycle.shunt.steps if cycle.shunt else []
    for st_ in every:
        if st_.rule and st_.rid not in seen:
            seen.add(st_.rid)
            out.append(f"{st_.rid}  :  {st_.rule}")
    return out


def _wrap(body, x0, y0, x1, y1, obj, span, mode, show_legend, target_px=TARGET_PX) -> str:
    w, h = x1 - x0, y1 - y0
    scale = target_px / w
    head = [f"<g transform='scale({scale:.4f}) translate({-x0:.4f},{y1:.4f})'>"]
    if show_legend:
        head += _legend(obj, x0 + 0.15, y0 + 0.2, w, span, mode)
    return "\n".join(
        [
            f"<svg xmlns='http://www.w3.org/2000/svg' width='{w * scale:.0f}' "
            f"height='{h * scale:.0f}' viewBox='0 0 {w * scale:.2f} {h * scale:.2f}'>",
            "<rect width='100%' height='100%' fill='white'/>",
            *head,
            *body,
            "</g></svg>",
        ]
    )


def _legend(cycle, x: float, y: float, width: float, span: float, mode: str) -> list[str]:
    out = []
    if cycle.title:
        out.append(text(x, -(y + 0.72), cycle.title, 0.11, "start", "#111", "bold"))

    bx, by, bw, bh = x, y + 0.30, min(2.4, width * 0.28), 0.11
    stops = "".join(
        f"<stop offset='{i / 10:.2f}' stop-color='{dg_colour(-span + 2 * span * i / 10, span)}'/>"
        for i in range(11)
    )
    out.append(f"<defs><linearGradient id='dgbar'>{stops}</linearGradient></defs>")
    out.append(
        f"<rect x='{bx:.4f}' y='{-(by + bh):.4f}' width='{bw:.4f}' height='{bh:.4f}' "
        f"fill='url(#dgbar)' stroke='#999' stroke-width='0.006'/>"
    )
    out.append(text(bx, -(by - 0.09), f"{-span:.0f}", 0.062, "start", "#444"))
    out.append(text(bx + bw / 2, -(by - 0.09), "dG (kJ/mol)   0", 0.062, "middle", "#444"))
    out.append(text(bx + bw, -(by - 0.09), f"+{span:.0f}", 0.062, "end", "#444"))

    rules = [f"{s.rid}  {s.rule}" for s in cycle.steps if s.rule]
    per = 5
    for i, line in enumerate(rules[: per * 3]):
        out.append(
            text(bx + bw + 0.5 + (i // per) * 1.9, -(y + 0.44 - (i % per) * 0.088),
                 line, 0.058, "start", "#333")
        )

    tot = cycle.total_dg
    what = "route" if isinstance(cycle, Pathway) else "cycle"
    note = f"{what} dG: not computed" if tot is None else f"{what} dG = {tot:+.1f} kJ/mol"
    out.append(text(x, -(y + 0.06), note, 0.068, "start", "#222"))
    if cycle.gain_steps:
        ids = ", ".join(cycle.steps[i].rid for i in cycle.gain_steps)
        out.append(text(x + 2.9, -(y + 0.06), f"gain step: {ids}", 0.068, "start", GAIN))
    out.append(text(x + 6.0, -(y + 0.06), f"width scale: {mode}", 0.068, "start", "#666"))
    if isinstance(cycle, Pathway) and cycle.dead_ends:
        out.append(
            text(x + 8.4, -(y + 0.06),
                 f"{len(cycle.dead_ends)} leaf/leaves not traced to a seed",
                 0.068, "start", "#a06000")
        )
    return out

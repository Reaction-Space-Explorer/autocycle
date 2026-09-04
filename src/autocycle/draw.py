"""SVG emission."""

from __future__ import annotations

import math

from autocycle import layout as L
from autocycle import tree as T
from autocycle.depict import embed, formula
from autocycle.encode import GREY, dg_colour, widths
from autocycle.style import Style

FONT = "Helvetica, Arial, sans-serif"
NODE_EDGE = "#8a8a8a"
WATER = "O"
RXN_FILL = "#8a6d3b"
FLOW = "#4a4a4a"
GAIN = "#1a7f37"


def _p(x: float, y: float) -> tuple[float, float]:
    return x, -y  # layout is y-up, SVG is y-down


def _arc(cx: float, cy: float, r: float, a0: float, a1: float, res: float = 2.0):
    n = max(2, int(abs(a1 - a0) / res) + 1)
    out = []
    for i in range(n):
        a = a0 + (a1 - a0) * i / (n - 1)
        dx, dy = L.polar(r, a)
        out.append(_p(cx + dx, cy + dy))
    return out


def _path(pts) -> str:
    return "M" + " L".join(f"{x:.4f},{y:.4f}" for x, y in pts) + " Z"


def annular_arrow(cx, cy, r, a0, a1, w, head_w=2.4, head_deg=None) -> str:
    """Annular arrow a0 -> a1 on radius r."""
    total = a1 - a0
    if head_deg is None:
        head_deg = min(abs(total) * 0.4, 10.0)
    base = a1 - math.copysign(head_deg, total)
    tipx, tipy = L.polar(r, a1)
    box, boy = L.polar(r + w * head_w / 2, base)
    bix, biy = L.polar(r - w * head_w / 2, base)
    pts = (
        _arc(cx, cy, r + w / 2, a0, base)
        + [_p(cx + box, cy + boy), _p(cx + tipx, cy + tipy), _p(cx + bix, cy + biy)]
        + _arc(cx, cy, r - w / 2, base, a0)
    )
    return _path(pts)


def annular_band(cx, cy, r, a0, a1, w) -> str:
    """Annular segment, no head."""
    pts = _arc(cx, cy, r + w / 2, a0, a1) + _arc(cx, cy, r - w / 2, a1, a0)
    return _path(pts)


def glyph(v: L.Vertex, size: float, kind: str, flow_cw: bool = True) -> str:
    """A square on the ring, or a triangle along the flow."""
    heading = v.angle - 90 if flow_cw else v.angle + 90
    if kind == "square":
        pts = []
        for k in (45, 135, 225, 315):
            dx, dy = L.polar(size * 1.414 / 2 * 1.42, heading + k)
            pts.append(_p(v.x + dx, v.y + dy))
        return _path(pts)
    if kind == "triangle":
        pts = []
        for k in (0, 120, 240):
            dx, dy = L.polar(size, heading + k)
            pts.append(_p(v.x + dx, v.y + dy))
        return _path(pts)
    raise ValueError(f"unknown glyph {kind!r}")


def bezier_arrow(x0, y0, x1, y1, bow=0.25, trim=0.0, head=0.09) -> tuple[str, str]:
    """Curved arrow for side flux. Returns (curve, head); `trim` shortens the target end."""
    ax, ay = _p(x0, y0)
    bx, by = _p(x1, y1)
    vx, vy = bx - ax, by - ay
    d = math.hypot(vx, vy) or 1.0
    if trim:
        bx, by = bx - vx / d * trim, by - vy / d * trim
        vx, vy = bx - ax, by - ay
        d = math.hypot(vx, vy) or 1.0
    mx, my = (ax + bx) / 2, (ay + by) / 2
    cx, cy = mx - vy / d * bow * d * 0.5, my + vx / d * bow * d * 0.5
    # head along the incoming tangent
    tx, ty = bx - cx, by - cy
    t = math.hypot(tx, ty) or 1.0
    tx, ty = tx / t, ty / t
    px, py = -ty, tx
    h = [
        (bx, by),
        (bx - tx * head + px * head * 0.42, by - ty * head + py * head * 0.42),
        (bx - tx * head - px * head * 0.42, by - ty * head - py * head * 0.42),
    ]
    curve = f"M{ax:.4f},{ay:.4f} Q{cx:.4f},{cy:.4f} {bx:.4f},{by:.4f}"
    return curve, _path(h)



def text(x, y, s, size=0.075, anchor="middle", fill="#222", weight="normal") -> str:
    s = s.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    return (
        f"<text x='{x:.4f}' y='{y:.4f}' font-family='{FONT}' font-size='{size}' "
        f"text-anchor='{anchor}' fill='{fill}' font-weight='{weight}'>{s}</text>"
    )


def draw_ring(ring: L.Ring, steps, span: float, st_: Style, mode: str, cw: bool = True) -> list[str]:
    out: list[str] = []
    if st_.uniform_width is None:
        ws = widths([s.mag for s in steps], mode=mode)
        revs = widths([s.rev_mag or 0.0 for s in steps], mode=mode)
    else:
        ws = [st_.uniform_width] * len(steps)
        revs = [st_.uniform_width * 0.6] * len(steps)

    for i, step in enumerate(steps):
        a0, a1 = L.step_arc(ring, i, cw)
        if st_.ring_colour == "dg":
            col = GREY if step.filtered else dg_colour(step.dg, span)
            edge, ew = "#555", "0.008"
        else:
            col = GREY if step.filtered else st_.ring_grey
            edge, ew = "none", "0"
        op = 0.35 if step.filtered else 1.0
        if step.gain:
            out.append(
                f"<path d='{annular_band(ring.cx, ring.cy, ring.radius, a0, a1, ws[i] * 2.6)}' "
                f"fill='{GAIN}' opacity='0.15'/>"
            )
        if step.reversible:
            r_out = ring.radius + ws[i] * 0.6
            r_in = ring.radius - max(revs[i], 0.05) * 0.6
            out.append(
                f"<path d='{annular_arrow(ring.cx, ring.cy, r_out, a0, a1, ws[i])}' "
                f"fill='{col}' opacity='{op}' stroke='{edge}' stroke-width='{ew}'/>"
            )
            rcol = col if st_.ring_colour == "dg" else "#9a9a9a"
            out.append(
                f"<path d='{annular_arrow(ring.cx, ring.cy, r_in, a1, a0, max(revs[i], 0.05))}' "
                f"fill='{rcol}' opacity='{op}' stroke='{edge}' stroke-width='{ew}'/>"
            )
        else:
            out.append(
                f"<path d='{annular_arrow(ring.cx, ring.cy, ring.radius, a0, a1, ws[i])}' "
                f"fill='{col}' opacity='{op}' stroke='{edge}' stroke-width='{ew}'/>"
            )
    return out


def draw_sides(ring: L.Ring, steps, st_: Style, anchors=None) -> list[str]:
    out: list[str] = []
    half = L.mol_half(ring) * st_.mol_scale
    anchors = L.side_points(ring, steps) if anchors is None else anchors
    for i, side, sp, (ax, ay) in anchors:
        if st_.water_as_text and sp.smiles == WATER:
            continue
        v = ring.verts[(2 * i + 1) % ring.n]
        if side == "in":
            curve, head = bezier_arrow(ax, ay, v.x, v.y, bow=0.16, trim=st_.rxn_size * 2.4)
        else:
            curve, head = bezier_arrow(v.x, v.y, ax, ay, bow=0.16, trim=half * 1.05)
        out.append(f"<path d='{curve}' fill='none' stroke='{st_.side_grey}' stroke-width='0.026'/>")
        out.append(f"<path d='{head}' fill='{st_.side_grey}'/>")
        px, py = _p(ax, ay)
        if st_.node_circle:
            out.append(
                f"<circle cx='{px:.4f}' cy='{py:.4f}' r='{half * 0.95:.4f}' fill='{st_.side_fill}' "
                f"fill-opacity='{st_.side_alpha}' stroke='none'/>"
            )
        out.append(embed(sp.smiles, px - half, py - half, half * 2, st_.backend))
        if sp.count > 1:
            out.append(text(px + half * 1.05, py - half * 0.9, f"{sp.count}x", 0.07, "start"))
    return out


def draw_nodes(
    ring: L.Ring, nodes, steps, seed: int | None, st_: Style, cw: bool = True, skip_mols=()
) -> list[str]:
    out: list[str] = []
    half = L.mol_half(ring) * st_.mol_scale
    for v in ring.verts:
        px, py = _p(v.x, v.y)
        if v.kind == "mol":
            if v.index in skip_mols:
                continue
            m = nodes[v.index]
            is_seed = seed is not None and v.index == seed
            if st_.node_circle and (is_seed or st_.circle_ring_nodes):
                edge = (
                    f"stroke='{GAIN}' stroke-width='0.022'"
                    if is_seed and st_.seed_ring
                    else "stroke='none'"
                )
                out.append(
                    f"<circle cx='{px:.4f}' cy='{py:.4f}' r='{half * 0.95:.4f}' fill='{st_.node_fill}' "
                    f"fill-opacity='{st_.node_alpha}' {edge}/>"
                )
            out.append(embed(m.smiles, px - half, py - half, half * 2, st_.backend))
            if m.label and not st_.centre_label:
                out.append(text(px, py + half + 0.1, m.label, 0.072))
        else:
            step = steps[v.index]
            fill = GREY if step.filtered else st_.rxn_fill
            out.append(
                f"<path d='{glyph(v, st_.rxn_size, st_.rxn_glyph, cw)}' fill='{fill}' "
                f"stroke='#7a5c10' stroke-width='0.006'/>"
            )
            lines = _step_lines(step, st_)
            reach = st_.rxn_size + 0.26 + 0.10 * len(lines)
            dx, dy = L.polar(reach, v.angle)
            sign = -1.0 if st_.label_inward else 1.0
            lx, ly = _p(v.x + sign * dx, v.y + sign * dy)
            for j, line in enumerate(lines):
                out.append(
                    text(lx, ly + j * st_.label_size * 1.2, line, st_.label_size, "middle",
                         GAIN if line == "gain" else "#222")
                )
    return out


def _water_balance(step) -> str | None:
    """Water consumed or produced, annotated rather than drawn."""
    n_in = sum(sp.count for sp in step.consumes if sp.smiles == WATER)
    n_out = sum(sp.count for sp in step.produces if sp.smiles == WATER)
    net = n_in - n_out
    return None if net == 0 else f"{net:+d} H2O"


def _step_lines(step, st_: Style) -> list[str]:
    if st_.step_label == "none":
        return ["gain"] if step.gain else []
    lines = []
    if st_.step_label in ("id_dg", "id_dg_units") or step.dg is None:
        lines.append(f"r:{step.rid}" if st_.step_label == "id_dg_units" else step.rid)
    if step.dg is not None:
        if st_.step_label == "dg":
            lines.append(f"{step.dg:.1f}")
        elif st_.step_label == "id_dg_units":
            lines.append(f"{step.dg:.2f} kJ/mol")
        else:
            lines.append(f"{step.dg:+.1f} kJ/mol")
    if st_.water_as_text and (w := _water_balance(step)):
        lines.append(w)
    if step.gain:
        lines.append("gain")
    return lines


def straight_arrow(x0, y0, x1, y1, w, head_w=2.5, head_len=None) -> str:
    """Filled tapered arrow, the route counterpart of the ring arc."""
    ax, ay = _p(x0, y0)
    bx, by = _p(x1, y1)
    vx, vy = bx - ax, by - ay
    d = math.hypot(vx, vy)
    if d == 0:
        return ""
    ux, uy = vx / d, vy / d
    px, py = -uy, ux
    hl = min(w * 2.6, d * 0.45) if head_len is None else head_len
    kx, ky = bx - ux * hl, by - uy * hl
    pts = [
        (ax + px * w / 2, ay + py * w / 2),
        (kx + px * w / 2, ky + py * w / 2),
        (kx + px * w * head_w / 2, ky + py * w * head_w / 2),
        (bx, by),
        (kx - px * w * head_w / 2, ky - py * w * head_w / 2),
        (kx - px * w / 2, ky - py * w / 2),
        (ax - px * w / 2, ay - py * w / 2),
    ]
    return _path(pts)


def _trim(a, b, ta, tb, square=False):
    """Shorten the segment a->b by ta at the start and tb at the end.

    With `square`, the trim reaches the edge of a box of half-width ta/tb rather than
    its inscribed circle, so a diagonal arrow does not start inside the box.
    """
    vx, vy = b[0] - a[0], b[1] - a[1]
    d = math.hypot(vx, vy) or 1.0
    ux, uy = vx / d, vy / d
    if square:
        scale = max(abs(ux), abs(uy)) or 1.0
        ta, tb = ta / scale, tb / scale
    return (a[0] + ux * ta, a[1] + uy * ta), (b[0] - ux * tb, b[1] - uy * tb)


def draw_route(pw, lay: T.TreeLayout, span: float, st_: Style, mode: str) -> list[str]:
    out: list[str] = []
    steps = pw.steps
    if st_.uniform_width is None:
        ws = dict(zip([id(s) for s in steps], widths([s.mag for s in steps], mode=mode), strict=True))
    else:
        ws = {id(s): st_.uniform_width for s in steps}

    for node in pw.nodes:
        if not node.step:
            continue
        r = lay.rxn(node)
        prod = lay.mol(node)
        w = ws[id(node.step)]
        if st_.ring_colour == "dg":
            col = GREY if node.step.filtered else dg_colour(node.step.dg, span)
        else:
            col = GREY if node.step.filtered else st_.ring_grey
        op = 0.35 if node.step.filtered else 1.0

        boxed = st_.node_circle
        for pre in node.precursors:
            m = lay.mol(pre)
            a, b = _trim((m.x, m.y), (r.x, r.y), m.half * 0.95, r.half * 2.4, square=boxed)
            out.append(f"<path d='{straight_arrow(*a, *b, w)}' fill='{col}' opacity='{op}'/>")
        a, b = _trim((r.x, r.y), (prod.x, prod.y), r.half * 2.4, prod.half * 0.95, square=boxed)
        out.append(f"<path d='{straight_arrow(*a, *b, w)}' fill='{col}' opacity='{op}'/>")
    return out


def draw_route_sides(pw, lay: T.TreeLayout, st_: Style) -> list[str]:
    out: list[str] = []
    half = T.MOL_HALF * 0.62 * st_.mol_scale
    for node in pw.nodes:
        if not node.step:
            continue
        r = lay.rxn(node)
        for side, group in (("in", node.step.consumes), ("out", node.step.produces)):
            for k, sp in enumerate(group):
                if st_.water_as_text and sp.smiles == WATER:
                    continue
                ax, ay = T.side_anchor(r, side, k)
                if side == "in":
                    curve, head = bezier_arrow(ax, ay, r.x, r.y, bow=0.1, trim=r.half * 2.6)
                else:
                    curve, head = bezier_arrow(r.x, r.y, ax, ay, bow=0.1, trim=half * 1.1)
                out.append(
                    f"<path d='{curve}' fill='none' stroke='{st_.side_grey}' stroke-width='0.026'/>"
                )
                out.append(f"<path d='{head}' fill='{st_.side_grey}'/>")
                px, py = _p(ax, ay)
                out.append(embed(sp.smiles, px - half, py - half, half * 2, st_.backend))
                if sp.count > 1:
                    out.append(text(px + half * 1.05, py - half * 0.9, f"{sp.count}x", 0.07, "start"))
    return out


def draw_route_nodes(pw, lay: T.TreeLayout, st_: Style) -> list[str]:
    out: list[str] = []
    for node in pw.nodes:
        m = lay.mol(node)
        px, py = _p(m.x, m.y)
        half = m.half * st_.mol_scale
        is_target = node is pw.root
        if st_.node_circle:
            out.append(
                f"<rect x='{px - half:.4f}' y='{py - half:.4f}' width='{2 * half:.4f}' "
                f"height='{2 * half:.4f}' rx='0.12' fill='{st_.node_fill}' "
                f"fill-opacity='{st_.node_alpha}' "
                f"stroke='{GAIN if is_target else NODE_EDGE}' "
                f"stroke-width='{0.022 if is_target else 0.009}'/>"
            )
        out.append(embed(node.mol.smiles, px - half, py - half, half * 2, st_.backend))

        tag = None
        if is_target:
            tag = (node.mol.label or formula(node.mol.smiles), "#111", "bold")
        elif not node.precursors:
            tag = {
                "seed": ("seed", GAIN, "normal"),
                "untraced": ("not traced", "#a06000", "normal"),
                "unknown": None,  # do not invent a label
            }.get(node.terminal)
        if tag:
            off = (half + st_.label_size * 1.3) if st_.node_circle else (half * 0.80 + st_.label_size * 1.3)
            out.append(text(px, py + off, tag[0], st_.label_size, "middle", tag[1], tag[2]))

        if node.step:
            r = lay.rxn(node)
            out.append(
                f"<path d='{glyph(T_vertex(r), st_.rxn_size, st_.rxn_glyph)}' "
                f"fill='{GREY if node.step.filtered else st_.rxn_fill}' "
                f"stroke='#7a5c10' stroke-width='0.006'/>"
            )
            lines = _step_lines(node.step, st_)
            # centre on the visible arrow, or the label lands inside the product box
            seg_a = r.x + r.half * 2.4
            seg_b = m.x - half * 0.95
            lx, ly0 = _p((seg_a + seg_b) / 2, (r.y + m.y) / 2)
            ly = ly0 - st_.label_size * 1.6
            for j, line in enumerate(lines):
                out.append(
                    text(lx, ly - (len(lines) - 1 - j) * st_.label_size * 1.2, line,
                         st_.label_size, "middle", "#222")
                )
    return out


def T_vertex(p: T.Placed) -> L.Vertex:
    """Adapt a placed square to the glyph helper."""
    return L.Vertex("rxn", 0, p.x, p.y, 90.0, p.half)


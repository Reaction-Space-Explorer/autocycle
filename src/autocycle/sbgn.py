"""SBGN Process Description export (SBGN-ML).

Molecules become `simple chemical` glyphs, reactions `process` glyphs, and the ring
becomes consumption and production arcs. Coordinates come from the ring layout, so the
map opens in an SBGN editor already arranged.

Lossy by design: SBGN-PD has no glyph for a shunt or for stoichiometric gain, so the two
things this package exists to show are recorded only as notes.
"""

from __future__ import annotations

from xml.etree.ElementTree import Element, SubElement, indent, tostring

from autocycle import layout as L
from autocycle.spec import Cycle

NS = "http://sbgn.org/libsbgn/0.2"
SCALE = 90.0
BOX = 70.0
PROC = 24.0


def _glyph(parent, gid: str, cls: str, x: float, y: float, w: float, h: float, label=None):
    g = SubElement(parent, "glyph", {"id": gid, "class": cls})
    if label is not None:
        SubElement(g, "label", {"text": str(label)})
    SubElement(g, "bbox", {"x": f"{x:.1f}", "y": f"{y:.1f}", "w": f"{w:.1f}", "h": f"{h:.1f}"})
    return g


def _arc(parent, aid: str, cls: str, source: str, target: str, a, b):
    arc = SubElement(parent, "arc", {"id": aid, "class": cls, "source": source, "target": target})
    SubElement(arc, "start", {"x": f"{a[0]:.1f}", "y": f"{a[1]:.1f}"})
    SubElement(arc, "end", {"x": f"{b[0]:.1f}", "y": f"{b[1]:.1f}"})
    return arc


def to_sbgn(cycle: Cycle) -> str:
    ring = L.lay_out(len(cycle.nodes))
    root = Element("sbgn", {"xmlns": NS})
    m = SubElement(root, "map", {"language": "process description"})

    notes = SubElement(m, "notes")
    body = SubElement(notes, "body")
    body.text = (
        "Exported by autocycle. SBGN-PD has no glyph for a shunt or for stoichiometric "
        "gain; " + (f"shunt from ring position {cycle.shunt.from_node}; " if cycle.shunt else "")
        + (f"autocatalyst at ring position {cycle.seed}." if cycle.seed is not None else "")
    )

    def px(x: float, y: float) -> tuple[float, float]:
        return x * SCALE, -y * SCALE

    mids: dict[int, str] = {}
    for v in ring.of("mol"):
        mol = cycle.nodes[v.index]
        gid = f"m{v.index}"
        mids[v.index] = gid
        cx, cy = px(v.x, v.y)
        _glyph(m, gid, "simple chemical", cx - BOX / 2, cy - BOX / 2, BOX, BOX,
               mol.label or mol.smiles)

    for v in ring.of("rxn"):
        step = cycle.steps[v.index]
        pid = f"p{v.index}"
        cx, cy = px(v.x, v.y)
        _glyph(m, pid, "process", cx - PROC / 2, cy - PROC / 2, PROC, PROC)
        n = len(cycle.nodes)
        src, tgt = mids[v.index % n], mids[(v.index + 1) % n]
        _arc(m, f"a{v.index}c", "consumption", src, pid, px(*_at(ring, v.index, 0)), (cx, cy))
        _arc(m, f"a{v.index}p", "production", pid, tgt, (cx, cy), px(*_at(ring, v.index, 1)))

        for k, sp in enumerate(step.consumes):
            sid = f"s{v.index}i{k}"
            sx, sy = px(*L.side_anchor(ring, v, "in", 1.5 + 1.9 * k))
            _glyph(m, sid, "simple chemical", sx - BOX / 2, sy - BOX / 2, BOX, BOX,
                   sp.label if getattr(sp, "label", None) else sp.smiles)
            _arc(m, f"{sid}a", "consumption", sid, pid, (sx, sy), (cx, cy))
        for k, sp in enumerate(step.produces):
            sid = f"s{v.index}o{k}"
            sx, sy = px(*L.side_anchor(ring, v, "out", 1.5 + 1.9 * k))
            _glyph(m, sid, "simple chemical", sx - BOX / 2, sy - BOX / 2, BOX, BOX,
                   sp.label if getattr(sp, "label", None) else sp.smiles)
            _arc(m, f"{sid}a", "production", pid, sid, (cx, cy), (sx, sy))

    indent(root)
    return '<?xml version="1.0" encoding="UTF-8"?>\n' + tostring(root, encoding="unicode")


def _at(ring: L.Ring, step: int, which: int) -> tuple[float, float]:
    v = ring.verts[(2 * step + 2 * which) % ring.n]
    return v.x, v.y

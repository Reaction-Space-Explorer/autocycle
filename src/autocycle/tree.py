"""Route geometry. Columns are route depth, target on the right."""

from __future__ import annotations

import math
from dataclasses import dataclass

from autocycle.spec import PathNode, Pathway

# Column and row pitch are a molecule plus a fixed gap, so a style that draws larger
# molecules gets a proportionally larger figure rather than a more crowded one.
X_GAP = 1.4
Y_GAP = 0.55
MOL_HALF = 0.95
X_PITCH = 2 * MOL_HALF + X_GAP
Y_PITCH = 2 * MOL_HALF + Y_GAP


def pitches(mol_scale: float = 1.0) -> tuple[float, float]:
    return 2 * MOL_HALF * mol_scale + X_GAP, 2 * MOL_HALF * mol_scale + Y_GAP


@dataclass(frozen=True)
class Placed:
    kind: str  # "mol" or "rxn"
    x: float
    y: float
    half: float


@dataclass
class TreeLayout:
    mols: dict[int, Placed]
    rxns: dict[int, Placed]
    order: list[PathNode]

    def mol(self, node: PathNode) -> Placed:
        return self.mols[id(node)]

    def rxn(self, node: PathNode) -> Placed:
        return self.rxns[id(node)]


#: Side species hang below their reaction while the next row's step label sits above
#: its own, and the two meet in the gap between rows. A route that draws side species
#: gets a taller row so they clear; one that does not keeps the compact pitch.
SIDE_ROW_EXTRA = 0.95


def has_side_species(pw: Pathway) -> bool:
    return any(n.step and (n.step.consumes or n.step.produces) for n in pw.nodes)


def lay_out_pathway(pw: Pathway, x_pitch: float = X_PITCH, y_pitch: float = Y_PITCH) -> TreeLayout:
    if has_side_species(pw):
        y_pitch += SIDE_ROW_EXTRA
    depths: dict[int, int] = {}

    def depth(n: PathNode, d: int) -> None:
        depths[id(n)] = d
        for p in n.precursors:
            depth(p, d + 1)

    depth(pw.root, 0)

    ys: dict[int, float] = {}
    slot = [0.0]
    order: list[PathNode] = []

    def assign(n: PathNode) -> float:
        order.append(n)
        if not n.precursors:
            y = slot[0] * y_pitch
            slot[0] += 1
        else:
            kids = [assign(p) for p in n.precursors]
            y = sum(kids) / len(kids)
        ys[id(n)] = y
        return y

    assign(pw.root)

    mols = {
        id(n): Placed("mol", -depths[id(n)] * x_pitch, ys[id(n)], MOL_HALF) for n in pw.nodes
    }
    # between the precursor column and the product
    rxns = {
        id(n): Placed("rxn", mols[id(n)].x - x_pitch / 2, ys[id(n)], MOL_HALF * 0.12)
        for n in pw.nodes
        if n.step
    }
    return TreeLayout(mols=mols, rxns=rxns, order=order)


def bounds(layout: TreeLayout, points, pad: float):
    xs, ys = [], []
    for p in list(layout.mols.values()) + list(layout.rxns.values()):
        xs += [p.x - p.half, p.x + p.half]
        ys += [p.y - p.half, p.y + p.half]
    for x, y in points:
        xs.append(x)
        ys.append(y)
    return min(xs) - pad, min(ys) - pad, max(xs) + pad, max(ys) + pad


# How far the anchor may swing off vertical to clear an arrow, in radians.
SWING = (0.0, 0.38, -0.38, 0.72, -0.72)


def incident_angles(lay: TreeLayout, node: PathNode) -> tuple[float, ...]:
    """Directions of the arrows meeting a node's reaction square."""
    r = lay.rxn(node)
    around = [lay.mol(p) for p in node.precursors] + [lay.mol(node)]
    return tuple(math.atan2(m.y - r.y, m.x - r.x) for m in around)


def side_anchor(
    rxn: Placed, side: str, k: int = 0, gap: float = 1.35, avoid: tuple[float, ...] = ()
) -> tuple[float, float]:
    """Consumed above the reaction, produced below, swung aside to clear `avoid`."""
    if side not in ("in", "out"):
        raise ValueError(f"side must be 'in' or 'out', got {side!r}")
    sign = 1.0 if side == "in" else -1.0
    base = sign * math.pi / 2
    dist = gap + k * 2.1
    theta = base
    if avoid:
        def clearance(off):
            a = base + off
            gaps = [abs(math.remainder(a - b, 2 * math.pi)) for b in avoid]
            return min(gaps) - 0.1 * abs(off)
        theta = base + max(SWING, key=clearance)
    return rxn.x + dist * math.cos(theta), rxn.y + dist * math.sin(theta)

"""Route geometry. Columns are route depth, target on the right."""

from __future__ import annotations

from dataclasses import dataclass

from autocycle.spec import PathNode, Pathway

X_PITCH = 3.5
Y_PITCH = 3.0
MOL_HALF = 0.95


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


def lay_out_pathway(pw: Pathway, x_pitch: float = X_PITCH, y_pitch: float = Y_PITCH) -> TreeLayout:
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


def side_anchor(rxn: Placed, side: str, k: int = 0, gap: float = 1.35) -> tuple[float, float]:
    """Consumed above the reaction, produced below."""
    if side not in ("in", "out"):
        raise ValueError(f"side must be 'in' or 'out', got {side!r}")
    sign = 1.0 if side == "in" else -1.0
    return rxn.x, rxn.y + sign * (gap + k * 2.1)

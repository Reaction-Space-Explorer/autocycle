"""A cycle unrolled into a linear diagram, from the seed round to the seed."""

from __future__ import annotations

from autocycle.spec import SEED, Cycle, PathNode, Pathway, SpecError


def linearise(cycle: Cycle) -> Pathway:
    """The ring read once through, starting and ending at the seed.

    The seed appears twice, as the molecule the turn consumes and as the one it
    returns, and a step producing an extra copy shows it as a side product, so the
    same specification reads as a sequence without changing what it states.
    """
    if cycle.shunt is not None or cycle.subs:
        raise SpecError(
            "a cycle with a shunt or a subcycle is not a single chain; "
            "decompose it into extreme currents first"
        )
    n = len(cycle.nodes)
    first = cycle.seed or 0
    node = PathNode(mol=cycle.nodes[first], terminal=SEED)
    for k in range(n):
        node = PathNode(
            mol=cycle.nodes[(first + k + 1) % n],
            step=cycle.steps[(first + k) % n],
            precursors=[node],
        )
    return Pathway(root=node, title=cycle.title, closed=True)

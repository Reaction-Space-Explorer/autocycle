"""Pick a maximally distinct subset, farthest-first over token sets."""

from __future__ import annotations

from autocycle.spec import Pathway


def tokens(obj) -> frozenset[str]:
    if isinstance(obj, Pathway):
        return _route_tokens(obj)
    cycle = obj
    t = {f"len:{len(cycle.nodes)}"}
    t |= {f"mol:{m.smiles}" for m in cycle.nodes}
    t |= {f"rule:{s.rule}" for s in cycle.steps if s.rule}
    n_side = sum(len(s.consumes) + len(s.produces) for s in cycle.steps)
    t.add(f"side:{n_side}")
    t.add(f"seed:{cycle.nodes[cycle.seed].smiles}" if cycle.seed is not None else "seed:none")
    if cycle.subs:
        t.add(f"subs:{len(cycle.subs)}")
    return frozenset(t)


def distance(a: frozenset[str], b: frozenset[str]) -> float:
    union = a | b
    return 1.0 - len(a & b) / len(union) if union else 0.0


def farthest_first(items: list[tuple[frozenset[str], object]], k: int) -> list[object]:
    """Greedy farthest-point subset."""
    if k <= 0 or not items:
        return []
    if k >= len(items):
        return [payload for _, payload in items]
    start = max(range(len(items)), key=lambda i: len(items[i][0]))
    chosen = [start]
    best = [distance(items[start][0], t) for t, _ in items]
    while len(chosen) < k:
        nxt = max(range(len(items)), key=lambda i: (best[i], i) if i not in chosen else (-1, i))
        chosen.append(nxt)
        for i, (t, _) in enumerate(items):
            best[i] = min(best[i], distance(items[nxt][0], t))
    return [items[i][1] for i in chosen]


def _route_tokens(pw: Pathway) -> frozenset[str]:
    t = {
        f"depth:{pw.root.depth}",
        f"leaves:{len(pw.leaves)}",
        f"steps:{len(pw.steps)}",
        f"complete:{pw.complete}",
        f"target:{pw.root.mol.smiles}",
    }
    t |= {f"mol:{n.mol.smiles}" for n in pw.nodes}
    t |= {f"rule:{s.rule}" for s in pw.steps if s.rule}
    return frozenset(t)

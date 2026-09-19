"""Reflexively autocatalytic food-generated sets, after Hordijk and Steel.

A RAF is a set of reactions in which every reaction is catalysed by a species the
set itself produces from the food, and every species it needs is reachable from
the food using only reactions in the set. The maximal RAF is the union of all of
them and is reached by removing whatever fails, repeatedly, until nothing does.

A CAF is the constructive version: a catalyst has to be present already, not
produced later by the set. Every CAF is a RAF; the converse fails, and the
difference is what Golnik's theorem turns on.
"""

from __future__ import annotations

from autocycle.io_crs import System, catalysed


def closure(system: System, names: set[str]) -> set[str]:
    """Species reachable from the food using `names`, ignoring catalysis."""
    by_name = {r.name: r for r in system.reactions}
    have, changed = set(system.food), True
    while changed:
        changed = False
        for n in names:
            r = by_name[n]
            if all(s in have for s in r.reactants) and not set(r.products) <= have:
                have |= set(r.products)
                changed = True
    return have


def max_raf(system: System) -> set[str]:
    """The largest RAF, or the empty set when the system contains none."""
    by_name = {r.name: r for r in system.reactions}
    names = set(by_name)
    while names:
        have = closure(system, names)
        keep = {n for n in names
                if all(s in have for s in by_name[n].reactants)
                and catalysed(by_name[n].catalyst_groups, have)}
        if keep == names:
            return names
        names = keep
    return names


def max_caf(system: System, names: set[str] | None = None) -> set[str]:
    """The largest CAF: a catalyst must be present before its reaction is used."""
    by_name = {r.name: r for r in system.reactions}
    names = set(by_name) if names is None else set(names)
    have, keep, changed = set(system.food), set(), True
    while changed:
        changed = False
        for n in names - keep:
            r = by_name[n]
            if all(s in have for s in r.reactants) and catalysed(r.catalyst_groups, have):
                keep.add(n)
                have |= set(r.products)
                changed = True
    return keep


def stoichiometry(system: System, names, food=None) -> dict[str, dict[str, int]]:
    """Net stoichiometry over non-food species. Catalysts cancel and drop out."""
    food = system.food if food is None else food
    by_name = {r.name: r for r in system.reactions}
    out = {}
    for n in names:
        r = by_name[n]
        col: dict[str, int] = {}
        for s in r.reactants:
            col[s] = col.get(s, 0) - 1
        for s in r.products:
            col[s] = col.get(s, 0) + 1
        col = {s: c for s, c in col.items() if c and s not in food}
        if col:
            out[n] = col
    return out

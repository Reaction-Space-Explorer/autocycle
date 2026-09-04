"""CatReNet `.crs` catalytic reaction systems.

Grammar, from the shipped examples:

    # comment
    r1: a1+b1 [c3] -> c1
    Food: a1,a2,b1

Species are abstract names, not structures, so figures from a `.crs` carry no chemistry;
molecules are drawn as their names.
"""

from __future__ import annotations

import re
from dataclasses import dataclass, field
from pathlib import Path

from autocycle.spec import Mol, Side, SpecError, Step

_RXN = re.compile(
    r"^\s*(?P<name>[^:]+?)\s*:\s*(?P<lhs>[^\[\]]*?)\s*"
    r"(?:\[(?P<cats>[^\]]*)\])?\s*(?P<arrow><->|->|<-)\s*(?P<rhs>.*?)\s*$"
)
_FOOD = re.compile(r"^\s*Food\s*:\s*(?P<items>.*?)\s*$", re.IGNORECASE)


@dataclass
class Reaction:
    name: str
    reactants: list[str]
    products: list[str]
    catalysts: list[str] = field(default_factory=list)
    reversible: bool = False


@dataclass
class System:
    reactions: list[Reaction]
    food: set[str]

    @property
    def species(self) -> set[str]:
        out: set[str] = set(self.food)
        for r in self.reactions:
            out |= set(r.reactants) | set(r.products) | set(r.catalysts)
        return out


def _split(text: str, sep: str) -> list[str]:
    return [x.strip() for x in text.split(sep) if x.strip()]


def read_crs(path: str | Path) -> System:
    reactions: list[Reaction] = []
    food: set[str] = set()
    for raw in Path(path).read_text().splitlines():
        line = raw.split("#", 1)[0].strip()
        if not line:
            continue
        if m := _FOOD.match(line):
            food |= set(_split(m.group("items"), ","))
            continue
        m = _RXN.match(line)
        if not m:
            raise SpecError(f"{Path(path).name}: cannot parse line {raw!r}")
        lhs, rhs = _split(m.group("lhs"), "+"), _split(m.group("rhs"), "+")
        if m.group("arrow") == "<-":
            lhs, rhs = rhs, lhs
        reactions.append(
            Reaction(
                name=m.group("name").strip(),
                reactants=lhs,
                products=rhs,
                catalysts=_split(m.group("cats") or "", ","),
                reversible=m.group("arrow") == "<->",
            )
        )
    if not reactions:
        raise SpecError(f"{Path(path).name}: no reactions found")
    return System(reactions=reactions, food=food)


def to_graph(system: System, mode: str = "catalysis"):
    """Species graph.

    `catalysis` draws catalyst -> product, which is where a catalytic reaction system's
    self-sustaining loops live: in `r: a+b [c] -> d`, c enables the production of d.
    `flow` draws reactant -> product instead, the mass-flow view.
    """
    import networkx as nx

    if mode not in ("catalysis", "flow"):
        raise ValueError(f"unknown mode {mode!r}")
    g = nx.MultiDiGraph()
    for r in system.reactions:
        sources = r.catalysts if mode == "catalysis" else r.reactants
        for a in sources:
            for b in r.products:
                g.add_edge(
                    a, b, reaction=r.name, rule=None, dg=None,
                    consumes=list(r.reactants),
                    produces=[x for x in r.products if x != b],
                    catalysts=r.catalysts,
                )
    return g


def find_cycles(system: System, mode: str = "catalysis", min_len: int = 2, max_len: int = 12):
    import networkx as nx

    g = nx.DiGraph(to_graph(system, mode))
    out = [c for c in nx.simple_cycles(g) if min_len <= len(c) <= max_len]
    return sorted(out, key=lambda c: (len(c), c))


def to_cycle(system: System, ring: list[str], title: str | None = None,
             mode: str = "catalysis"):
    """A Cycle over abstract species, with food molecules as feeders."""
    from autocycle.spec import Cycle

    g = to_graph(system, mode)
    steps = []
    for i, src in enumerate(ring):
        tgt = ring[(i + 1) % len(ring)]
        data = g.get_edge_data(src, tgt)
        if not data:
            raise SpecError(f"no reaction {src} -> {tgt}; ring is not a cycle in this system")
        d = next(iter(data.values()))
        steps.append(
            Step(
                rid=d["reaction"],
                consumes=[
                    Side(x, structure=False) for x in d["consumes"] if x in system.food
                ],
                produces=[Side(x, structure=False) for x in d["produces"]],
            )
        )
    return Cycle(
        nodes=[Mol(s, label=s, structure=False) for s in ring],
        steps=steps,
        title=title,
        seed=0,
    )

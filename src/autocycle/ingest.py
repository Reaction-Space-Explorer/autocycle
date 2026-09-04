"""Edge list -> networkx -> cycles."""

from __future__ import annotations

import csv
from pathlib import Path

import networkx as nx

from autocycle.spec import Cycle, Mol, Side, SpecError, Step, canonical

FIELDS = ("source", "target", "reaction", "dg", "rule", "consumes", "produces")


def read_edges(path: str | Path) -> nx.MultiDiGraph:
    """Read a source/target CSV."""
    g = nx.MultiDiGraph()
    with open(path, newline="") as fh:
        rows = list(csv.DictReader(fh))
    if not rows:
        raise SpecError(f"{path}: no rows")
    missing = {"source", "target"} - set(rows[0])
    if missing:
        raise SpecError(f"{path}: missing column(s) {sorted(missing)}")
    for i, row in enumerate(rows):
        src, tgt = canonical(row["source"]), canonical(row["target"])
        dg = row.get("dg", "").strip()
        g.add_edge(
            src,
            tgt,
            reaction=row.get("reaction") or f"e{i}",
            dg=float(dg) if dg else None,
            rule=(row.get("rule") or "").strip() or None,
            consumes=_species(row.get("consumes")),
            produces=_species(row.get("produces")),
        )
    return g


def _species(field: str | None) -> list[str]:
    if not field or not field.strip():
        return []
    return [canonical(s) for s in field.replace(".", " ").split() if s]


def spontaneous(g: nx.MultiDiGraph, cut: float = 0.0) -> nx.MultiDiGraph:
    """Tag dg < cut. Unknown dg stays None, not False."""
    for _, _, d in g.edges(data=True):
        d["spontaneous"] = None if d["dg"] is None else d["dg"] < cut
    return g


def find_cycles(g: nx.MultiDiGraph, min_len: int = 3, max_len: int = 12) -> list[list[str]]:
    out = [c for c in nx.simple_cycles(nx.DiGraph(g)) if min_len <= len(c) <= max_len]
    return sorted(out, key=lambda c: (len(c), c))


def to_cycle(
    g: nx.MultiDiGraph,
    ring: list[str],
    title: str | None = None,
    gain_at: int | None = None,
) -> Cycle:
    """Build a Cycle from a ring of molecule SMILES, taking step data off the edges."""
    steps = []
    for i, src in enumerate(ring):
        tgt = ring[(i + 1) % len(ring)]
        data = g.get_edge_data(src, tgt)
        if not data:
            raise SpecError(f"no edge {src} -> {tgt}; ring is not a cycle in this graph")
        d = min(data.values(), key=lambda e: (e["dg"] is None, e["dg"] or 0.0))
        steps.append(
            Step(
                rid=d["reaction"],
                rule=d["rule"],
                dg=d["dg"],
                mag=2.0 if gain_at == i else 1.0,
                gain=gain_at == i,
                consumes=[Side(s) for s in d["consumes"]],
                produces=[Side(s) for s in d["produces"]],
                filtered=d.get("spontaneous") is False,
            )
        )
    return Cycle(nodes=[Mol(s) for s in ring], steps=steps, title=title, seed=0)


def summary_csv(g: nx.MultiDiGraph, cycles: list[list[str]], path: str | Path) -> None:
    """One row per cycle. Unknown dg is listed, not dropped."""
    with open(path, "w", newline="") as fh:
        w = csv.writer(fh)
        w.writerow(["index", "length", "total_dg", "all_spontaneous", "reactions", "molecules"])
        for i, ring in enumerate(cycles):
            c = to_cycle(g, ring)
            dgs = [s.dg for s in c.steps]
            if any(d is None for d in dgs):
                total, spont = "", ""
            else:
                total, spont = f"{sum(dgs):.2f}", all(d < 0 for d in dgs)
            w.writerow(
                [i, len(ring), total, spont,
                 ";".join(s.rid for s in c.steps), ";".join(m.smiles for m in c.nodes)]
            )

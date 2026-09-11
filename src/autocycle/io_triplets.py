"""A stoichiometric matrix in triplet form -> networkx -> cycles.

One row per non-zero entry: reaction, species, coefficient, and optionally the
rule that produced it. Negative is consumed, positive is produced, and a species
appearing twice in one reaction is summed, so multiplicity may be written either
as a coefficient or as repeated rows. This is the sparse form stoichiometry is
usually stored in, whatever wrote it.
"""

from __future__ import annotations

import collections
import csv
from pathlib import Path

import networkx as nx

from autocycle.spec import SpecError, canonical

FIELDS = ("reaction", "species", "coefficient", "rule")


def _rows(path: str | Path):
    text = Path(path).read_text().splitlines()
    if not text:
        raise SpecError(f"{path}: no rows")
    sep = "\t" if text[0].count("\t") >= text[0].count(",") else ","
    rdr = list(csv.reader(text, delimiter=sep))
    head = [h.strip().lower() for h in rdr[0]]
    if set(FIELDS[:3]) <= set(head):
        idx = [head.index(f) for f in FIELDS[:3]]
        rule = head.index("rule") if "rule" in head else None
        body = rdr[1:]
    else:                                  # no header: reaction, species, coefficient, rule
        idx, rule, body = [0, 1, 2], (3 if len(rdr[0]) > 3 else None), rdr
    for n, row in enumerate(body, 1):
        if len(row) <= max(idx):
            continue
        try:
            coef = int(float(row[idx[2]]))
        except ValueError:
            raise SpecError(
                f"{Path(path).name}: row {n} has a non-numeric coefficient {row[idx[2]]!r}"
            ) from None
        yield row[idx[0]], row[idx[1]], coef, (row[rule].strip() if rule is not None else None)


def read_triplets(path: str | Path) -> nx.MultiDiGraph:
    """Build the species graph, carrying each reaction's other species as sides."""
    nu: dict[str, dict[str, int]] = collections.defaultdict(lambda: collections.defaultdict(int))
    rules: dict[str, str | None] = {}
    for rid, smiles, coef, rule in _rows(path):
        nu[rid][canonical(smiles)] += coef
        rules.setdefault(rid, rule)
    if not nu:
        raise SpecError(f"{path}: no reactions")

    g = nx.MultiDiGraph()
    for rid, d in nu.items():
        left = [s for s, c in d.items() if c < 0]
        right = [s for s, c in d.items() if c > 0]
        if not left or not right:
            continue                        # nothing is transformed; not an edge
        for src in left:
            for tgt in right:
                if src == tgt:
                    continue
                g.add_edge(
                    src, tgt, reaction=rid, rule=rules[rid], dg=None,
                    consumes=[s for s in left for _ in range(-d[s]) if s != src]
                             + [src] * (-d[src] - 1),
                    produces=[s for s in right for _ in range(d[s]) if s != tgt]
                             + [tgt] * (d[tgt] - 1),
                )
    if not g:
        raise SpecError(f"{path}: no reaction transforms one species into another")
    return g


def gain_step(g: nx.MultiDiGraph, ring: list[str]) -> int | None:
    """The ring step that returns an extra copy of the molecule it leads to."""
    for i, src in enumerate(ring):
        tgt = ring[(i + 1) % len(ring)]
        for d in (g.get_edge_data(src, tgt) or {}).values():
            if tgt in d["produces"]:
                return i
    return None


def orient(g: nx.MultiDiGraph, ring: list[str]) -> tuple[list[str], int | None]:
    """Rotate the ring so the species produced in excess comes first.

    A cycle has no inherent starting point and a cycle search returns whichever
    rotation it reached first, so the seed has to be chosen from the coefficients
    rather than from position. Returns the rotated ring and the index of its gain
    step, which is then always the last one.
    """
    i = gain_step(g, ring)
    if i is None:
        return ring, None
    j = (i + 1) % len(ring)
    return ring[j:] + ring[:j], len(ring) - 1

"""Cypher ring-query results -> Cycle.

`ringMols` is the ring in path order, `ringRels` the joining reactions. A path that
pinches (a molecule visited twice) is not a simple cycle and is rejected.
"""

from __future__ import annotations

import ast
import csv
from pathlib import Path

from autocycle.spec import Cycle, Mol, Shunt, Side, SpecError, Step, canonical


def read_rows(path: str | Path) -> list[dict]:
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


def _lit(row: dict, key: str):
    if key not in row:
        raise SpecError(f"missing column {key!r}")
    try:
        return ast.literal_eval(row[key])
    except (ValueError, SyntaxError) as exc:
        raise SpecError(f"column {key!r} is not a Python literal: {exc}") from exc


def _generations(row: dict) -> dict[str, int]:
    """SMILES -> generation_formed, from every node list the row carries."""
    out: dict[str, int] = {}
    keys = ("ringPathNodes", "attachedPathNodes", "branchedBeginMolPathNodes", "autocatPathNodes")
    for key in keys:
        if key not in row:
            continue
        try:
            nodes = _lit(row, key)
        except SpecError:
            continue
        for n in nodes if isinstance(nodes, list) else []:
            smi = _smiles(n)
            gen = n.get("generation_formed") if isinstance(n, dict) else None
            if smi and gen is not None:
                out.setdefault(canonical(smi), int(gen))
    return out


def _smiles(node) -> str | None:
    if isinstance(node, dict):
        return node.get("smiles_str")
    return node if isinstance(node, str) else None


def from_cypher_row(row: dict, title: str | None = None) -> Cycle:
    mols = _lit(row, "ringMols")
    rels = _lit(row, "ringRels")
    if len(set(mols)) != len(mols) or len(mols) != len(rels):
        raise SpecError(
            f"not a simple cycle: {len(mols)} ring molecules "
            f"({len(set(mols))} distinct) and {len(rels)} relationships - pinched ring path"
        )

    n = len(mols)
    gens = _generations(row)
    steps = [
        Step(
            rid=str(r.get("rxn_id") or f"rel{i}"),
            rule=r.get("rule"),
            note=None if r.get("generation_formed") is None else f"gen {r['generation_formed']}",
        )
        for i, r in enumerate(rels)
    ]

    begin = _smiles(_lit(row, "beginMol"))
    if begin not in mols:
        raise SpecError(f"beginMol {begin!r} is not on the ring")
    seed = mols.index(begin)
    # No gain flag: these rows carry no stoichiometric coefficients, so an extra yield
    # of the seed cannot be established. autocycle.verify reports it as unknown.

    inter = _smiles(_lit(row, "intermediateMol")) if "intermediateMol" in row else None
    if inter in mols:
        j = mols.index(inter)
        feeder = _smiles(_lit(row, "feederMol")) if "feederMol" in row else None
        consumer = _smiles(_lit(row, "consumerMol")) if "consumerMol" in row else None
        if feeder:
            steps[(j - 1) % n].consumes.append(
                Side(feeder, generation=gens.get(canonical(feeder)))
            )
        if consumer:
            steps[j].produces.append(Side(consumer, generation=gens.get(canonical(consumer))))

    shunt = None
    cat = _smiles(_lit(row, "catMolInRing")) if (row.get("catMolInRing") or "").strip() else None
    if cat in mols and (row.get("autocatPathRels") or "").strip():
        rels_s = _lit(row, "autocatPathRels")
        nodes_s = _lit(row, "autocatPathNodes") if (row.get("autocatPathNodes") or "").strip() else []
        inner = [_smiles(x) for x in nodes_s][1:-1] if len(nodes_s) > 2 else []
        shunt = Shunt(
            from_node=mols.index(cat),
            steps=[
                Step(rid=str(r.get("rxn_id") or f"shunt{i}"), rule=r.get("rule"))
                for i, r in enumerate(rels_s)
            ],
            nodes=[Mol(m, generation=gens.get(canonical(m))) for m in inner if m],
        )

    exported = _smiles(_lit(row, "beginMolConsumer")) if "beginMolConsumer" in row else None
    if exported:
        steps[seed].produces.append(Side(exported, generation=gens.get(canonical(exported))))

    return Cycle(
        nodes=[Mol(m, generation=gens.get(canonical(m))) for m in mols],
        steps=steps,
        title=title,
        seed=seed,
        shunt=shunt,
    )

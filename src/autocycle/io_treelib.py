"""Read treelib-rendered route files.

Two traps in the format:
- `Pathway Length` is tree depth, not the reaction count.
- `Reaction IDs` cannot be mapped to nodes by position: treelib sorts children when
  printing, an id list is in insertion order. Reactions are matched by content instead.
"""

from __future__ import annotations

import csv
import re
from pathlib import Path

from autocycle.spec import (
    SEED,
    UNKNOWN,
    UNTRACED,
    Mol,
    PathNode,
    Pathway,
    SpecError,
    Step,
    canonical,
)

RECORD_SEP = "-" * 21
_BRANCH = re.compile(r"^((?:(?:│|\|)\s{3}|\s{4})*)(?:├──|└──)\s(.+?)\s*$")
_HEAD = re.compile(r"^(Network Smiles|Generation|INCHIKEY|Analogue Smiles):\s*(.*)$")
_BLOCK = re.compile(r"^\s*(?:Pathway\s+\d+|Shortest Pathway)\s*$")
_TRAIL = re.compile(r"^\s*(Pathway Length|Energy Change|Reaction IDs):\s*(.*)$")


def read_seeds(products_tsv: str | Path) -> set[str]:
    """Generation-0 molecules."""
    out = set()
    with open(products_tsv, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            if row.get("Generation") == "G0":
                out.add(canonical(row["Smiles"]))
    return out


def read_rels(paths) -> dict[str, dict]:
    """rxn_id -> {rule, reactants, products}."""
    rels: dict[str, dict] = {}
    for path in paths:
        with open(path, newline="") as fh:
            for row in csv.DictReader(fh, delimiter="\t"):
                rid = row.get("Index")
                if not rid:
                    continue
                r = rels.setdefault(rid, {"rule": row.get("Rule"), "reactants": set(), "products": set()})
                side = "products" if str(row.get("Formed/Produced", "")).strip() == "1" else "reactants"
                r[side].add(canonical(row["Reagent"]))
    return rels


def parse_file(
    path: str | Path,
    seeds: set[str] | None = None,
    rels: dict[str, dict] | None = None,
) -> list[Pathway]:
    text = Path(path).read_text()
    out: list[Pathway] = []
    for record in text.split(RECORD_SEP):
        if record.strip():
            out += _parse_record(record, Path(path).stem, seeds, rels)
    return out


def _parse_record(record: str, stem: str, seeds, rels) -> list[Pathway]:
    head: dict[str, str] = {}
    blocks: list[tuple[list[str], dict[str, str]]] = []
    lines, trail = [], {}
    in_block = False

    for raw in record.splitlines():
        if not raw.strip():
            continue
        if (m := _HEAD.match(raw)) and not in_block:
            head[m.group(1)] = m.group(2)
            continue
        if _BLOCK.match(raw):
            if in_block:
                blocks.append((lines, trail))
            lines, trail, in_block = [], {}, True
            continue
        if m := _TRAIL.match(raw):
            trail[m.group(1)] = m.group(2)
            continue
        if in_block:
            lines.append(raw)
    if in_block:
        blocks.append((lines, trail))

    return [_build(ls, tr, head, stem, seeds, rels) for ls, tr in blocks if ls]


def _build(lines: list[str], trail: dict, head: dict, stem: str, seeds, rels) -> Pathway:
    root_smi = lines[0].strip()
    kids: dict[int, list] = {}
    nodes: dict[int, tuple[str, int]] = {0: (root_smi, -1)}
    stack = [0]
    nxt = 1
    for raw in lines[1:]:
        m = _BRANCH.match(raw)
        if not m:
            raise SpecError(f"{stem}: cannot parse tree line {raw!r}")
        depth = len(m.group(1)) // 4 + 1
        if depth > len(stack):
            raise SpecError(f"{stem}: tree line jumps from depth {len(stack) - 1} to {depth}")
        del stack[depth:]
        parent = stack[-1]
        nodes[nxt] = (m.group(2), parent)
        kids.setdefault(parent, []).append(nxt)
        stack.append(nxt)
        nxt += 1

    ids = []
    if raw_ids := trail.get("Reaction IDs"):
        ids = [x.strip().strip("'\"") for x in raw_ids.strip("[]").split(",") if x.strip()]
    internal = [i for i in nodes if kids.get(i)]
    if ids and len(ids) != len(internal):
        raise SpecError(
            f"{stem}: {len(ids)} reaction ids but {len(internal)} internal nodes"
        )

    def resolve(i: int) -> Step:
        mol = canonical(nodes[i][0])
        pres = {canonical(nodes[c][0]) for c in kids[i]}
        if rels:
            hits = sorted(
                {
                    rid
                    for rid in ids
                    if rid in rels
                    and mol in rels[rid]["products"]
                    and rels[rid]["reactants"] <= pres
                }
            )
            if len(hits) == 1:
                return Step(rid=hits[0], rule=rels[hits[0]]["rule"])
            if hits:
                return Step(rid=hits[0], rule=rels[hits[0]]["rule"], note="ambiguous")
        return Step(rid=f"?{i}", note="unresolved")

    def make(i: int) -> PathNode:
        smi = nodes[i][0]
        if kids.get(i):
            return PathNode(mol=Mol(smi), step=resolve(i), precursors=[make(c) for c in kids[i]])
        term = UNKNOWN if seeds is None else (SEED if canonical(smi) in seeds else UNTRACED)
        return PathNode(mol=Mol(smi), terminal=term, generation=head.get("Generation"))

    dg = trail.get("Energy Change")
    return Pathway(
        root=make(0),
        title=head.get("Network Smiles"),
        reported_dg=float(dg) if dg else None,
        meta={
            "inchikey": head.get("INCHIKEY"),
            "generation": head.get("Generation"),
            "depth": trail.get("Pathway Length"),
            "source": stem,
        },
    )

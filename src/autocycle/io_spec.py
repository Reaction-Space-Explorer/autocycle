"""YAML and reaction-SMILES front ends."""

from __future__ import annotations

from pathlib import Path

import yaml

from autocycle.spec import (
    SEED,
    UNTRACED,
    Cycle,
    Mol,
    PathNode,
    Pathway,
    Shunt,
    Side,
    SpecError,
    Step,
    Sub,
    canonical,
)


# a mistyped key would otherwise be dropped in silence, and `flux` or `count` decide
# the reported stoichiometry
def _keys(d, allowed: set[str], what: str) -> dict:
    if not isinstance(d, dict):
        raise SpecError(f"{what}: expected a mapping, got {type(d).__name__}")
    unknown = set(d) - allowed
    if unknown:
        raise SpecError(f"{what}: unknown key(s) {sorted(unknown)}; allowed {sorted(allowed)}")
    return d


def _index(v, what: str) -> int:
    try:
        return int(v)
    except (TypeError, ValueError):
        raise SpecError(f"{what} must be a node index, got {v!r}") from None


def _sides(items) -> list[Side]:
    out = []
    for it in items or []:
        if isinstance(it, str):
            out.append(Side(it))
            continue
        _keys(it, {"smiles", "count", "generation", "structure"}, "side species")
        out.append(
            Side(
                it["smiles"],
                int(it.get("count", 1)),
                generation=it.get("generation"),
                structure=bool(it.get("structure", True)),
            )
        )
    return out


def _step(d: dict) -> Step:
    _keys(d, {"id", "rule", "dg", "mag", "flux", "rev_mag", "consumes", "produces",
              "gain", "filtered", "note"}, f"step {d.get('id', '?')}")
    if "id" not in d:
        raise SpecError(f"step needs an 'id': {d}")
    return Step(
        rid=str(d["id"]),
        rule=d.get("rule"),
        dg=None if d.get("dg") is None else float(d["dg"]),
        mag=float(d.get("mag", 1.0)),
        flux=float(d.get("flux", 1.0)),
        rev_mag=None if d.get("rev_mag") is None else float(d["rev_mag"]),
        consumes=_sides(d.get("consumes")),
        produces=_sides(d.get("produces")),
        gain=bool(d.get("gain", False)),
        filtered=bool(d.get("filtered", False)),
        note=d.get("note"),
    )


def _mol(d) -> Mol:
    if isinstance(d, str):
        return Mol(d)
    _keys(d, {"smiles", "label", "generation", "structure"}, "molecule")
    if "smiles" not in d:
        raise SpecError(f"molecule needs 'smiles': {d}")
    return Mol(
        d["smiles"],
        d.get("label"),
        generation=d.get("generation"),
        structure=bool(d.get("structure", True)),
    )


def load_yaml(path: str | Path) -> Cycle:
    raw = yaml.safe_load(Path(path).read_text())
    if not isinstance(raw, dict):
        raise SpecError(f"{path}: expected a mapping at the top level")
    for key in ("nodes", "steps"):
        if key not in raw:
            raise SpecError(f"{path}: missing '{key}'")
        if not isinstance(raw[key], list):
            raise SpecError(f"{path}: '{key}' must be a list")
    _keys(raw, {"title", "seed", "stoichiometry_complete", "nodes", "steps",
                "subcycles", "shunt"}, str(path))
    subs = [
        Sub(
            at_step=int(_keys(s, {"at_step", "nodes", "steps", "label"}, "subcycle")["at_step"]),
            nodes=[_mol(m) for m in s["nodes"]],
            steps=[_step(x) for x in s["steps"]],
            label=s.get("label"),
        )
        for s in raw.get("subcycles", [])
    ]
    shunt = None
    if "shunt" in raw:
        sh = _keys(raw["shunt"], {"from_node", "steps", "nodes"}, "shunt")
        if "from_node" not in sh:
            raise SpecError(f"{path}: shunt needs 'from_node'")
        shunt = Shunt(
            from_node=_index(sh["from_node"], "shunt from_node"),
            steps=[_step(x) for x in sh.get("steps", [])],
            nodes=[_mol(m) for m in sh.get("nodes", [])],
        )

    return Cycle(
        nodes=[_mol(m) for m in raw["nodes"]],
        steps=[_step(s) for s in raw["steps"]],
        subs=subs,
        shunt=shunt,
        title=raw.get("title"),
        seed=None if raw.get("seed") is None else _index(raw["seed"], "seed"),
        stoichiometry_complete=bool(raw.get("stoichiometry_complete", False)),
    )


def from_reaction_smiles(reactions: list[str], order: list[str], title: str | None = None) -> Cycle:
    """Match each consecutive node pair to a reaction."""
    parsed = []
    for r in reactions:
        if ">>" not in r:
            raise SpecError(f"not a reaction SMILES (needs '>>'): {r!r}")
        lhs, rhs = r.split(">>", 1)
        parsed.append(
            (
                [canonical(x) for x in lhs.split(".") if x],
                [canonical(x) for x in rhs.split(".") if x],
                r,
            )
        )
    ring = [canonical(o) for o in order]
    steps = []
    for i, src in enumerate(ring):
        tgt = ring[(i + 1) % len(ring)]
        hit = next(((L, R, r) for L, R, r in parsed if src in L and tgt in R), None)
        if hit is None:
            raise SpecError(f"no reaction converts {src} -> {tgt}")
        lhs, rhs, _ = hit
        steps.append(
            Step(
                rid=f"s{i + 1}",
                consumes=[Side(s) for s in lhs if s != src],
                produces=[Side(s) for s in rhs if s != tgt],
            )
        )
    return Cycle(nodes=[Mol(s) for s in ring], steps=steps, title=title, seed=0)


def _node(d) -> PathNode:
    if isinstance(d, str):
        return PathNode(mol=Mol(d), terminal=SEED)
    if "mol" not in d and "smiles" not in d:
        raise SpecError(f"route node needs 'mol' or 'smiles': {d}")
    _keys(d, {"mol", "smiles", "from", "reaction", "terminal", "seed", "generation"},
          "route node")
    mol = _mol(d.get("mol", d.get("smiles")))
    pre = [_node(x) for x in d.get("from", [])]
    step = None
    if "reaction" in d:
        step = _step(d["reaction"])
    elif pre:
        raise SpecError(f"{mol.smiles}: has 'from' precursors but no 'reaction'")
    if "terminal" in d:
        term = d["terminal"]
    elif "seed" in d:
        term = SEED if d["seed"] else UNTRACED
    else:
        term = SEED if not pre else None
    return PathNode(
        mol=mol,
        step=step,
        precursors=pre,
        terminal=term,
        generation=d.get("generation"),
    )


def load_pathway_yaml(path: str | Path) -> Pathway:
    raw = yaml.safe_load(Path(path).read_text())
    if not isinstance(raw, dict) or "target" not in raw:
        raise SpecError(f"{path}: a route needs a 'target' mapping at the top level")
    _keys(raw, {"target", "title"}, str(path))
    return Pathway(root=_node(raw["target"]), title=raw.get("title"))


def from_route_smiles(reactions: list[str], target: str, title: str | None = None) -> Pathway:
    """Walk back from `target`. Co-products are kept, not dropped."""
    parsed = []
    for r in reactions:
        if ">>" not in r:
            raise SpecError(f"not a reaction SMILES (needs '>>'): {r!r}")
        lhs, rhs = r.split(">>", 1)
        parsed.append(
            ([canonical(x) for x in lhs.split(".") if x], [canonical(x) for x in rhs.split(".") if x])
        )

    def build(smi: str, seen: tuple[str, ...], idx: list[int]) -> PathNode:
        for lhs, rhs in parsed:
            if smi in rhs and not any(x in seen for x in lhs):
                idx[0] += 1
                step = Step(
                    rid=f"s{idx[0]}",
                    produces=[Side(x) for x in rhs if x != smi],
                )
                return PathNode(
                    mol=Mol(smi),
                    step=step,
                    precursors=[build(x, (*seen, smi), idx) for x in lhs],
                )
        return PathNode(mol=Mol(smi), terminal=UNTRACED)

    return Pathway(root=build(canonical(target), (), [0]), title=title)


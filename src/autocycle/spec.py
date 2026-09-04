from __future__ import annotations

from dataclasses import dataclass, field

from rdkit import Chem, RDLogger

RDLogger.DisableLog("rdApp.*")


class SpecError(ValueError):
    pass


SEED = "seed"          # stopped deliberately: a network seed
UNTRACED = "untraced"  # ran out: nothing produces this
UNKNOWN = "unknown"    # the source did not say


def canonical(smi: str) -> str:
    m = Chem.MolFromSmiles(smi)
    if m is None:
        raise SpecError(f"bad SMILES: {smi!r}")
    return Chem.MolToSmiles(m)


@dataclass
class Mol:
    smiles: str
    label: str | None = None
    generation: int | None = None
    structure: bool = True   # False for abstract species, drawn as a name

    def __post_init__(self):
        if self.structure:
            self.smiles = canonical(self.smiles)


@dataclass
class Side:
    smiles: str
    count: int = 1
    generation: int | None = None
    structure: bool = True

    def __post_init__(self):
        if self.structure:
            self.smiles = canonical(self.smiles)
        if self.count < 1:
            raise SpecError(f"count must be >= 1, got {self.count}")


@dataclass
class Step:
    rid: str
    rule: str | None = None
    dg: float | None = None
    mag: float = 1.0
    rev_mag: float | None = None
    consumes: list[Side] = field(default_factory=list)
    produces: list[Side] = field(default_factory=list)
    gain: bool = False
    filtered: bool = False
    note: str | None = None

    def __post_init__(self):
        if self.mag < 0 or (self.rev_mag is not None and self.rev_mag < 0):
            raise SpecError(f"{self.rid}: magnitudes must be >= 0")

    @property
    def reversible(self) -> bool:
        return self.rev_mag is not None


@dataclass
class Shunt:
    """A bridging path from a ring molecule back to the seed.

    The published search calls this the shunt, and treats its presence as the minimal
    topological evidence that the seed is produced in quantity greater than one.
    """

    from_node: int
    steps: list[Step] = field(default_factory=list)
    nodes: list[Mol] = field(default_factory=list)


@dataclass
class Sub:
    at_step: int
    nodes: list[Mol]
    steps: list[Step]
    label: str | None = None

    def __post_init__(self):
        _check(self.nodes, self.steps)


@dataclass
class Cycle:
    nodes: list[Mol]
    steps: list[Step]
    subs: list[Sub] = field(default_factory=list)
    shunt: Shunt | None = None
    title: str | None = None
    seed: int | None = None
    # set only when the source states full stoichiometry, so that the absence of an
    # extra copy of the seed means n = 1 rather than "not recorded"
    stoichiometry_complete: bool = False

    def __post_init__(self):
        _check(self.nodes, self.steps)
        for s in self.subs:
            if not 0 <= s.at_step < len(self.steps):
                raise SpecError(f"sub at_step {s.at_step} out of range")
        if self.seed is not None and not 0 <= self.seed < len(self.nodes):
            raise SpecError(f"seed {self.seed} out of range")
        if self.shunt is not None and not 0 <= self.shunt.from_node < len(self.nodes):
            raise SpecError(f"shunt from_node {self.shunt.from_node} out of range")

    @property
    def gain_steps(self) -> list[int]:
        return [i for i, s in enumerate(self.steps) if s.gain]

    @property
    def total_dg(self) -> float | None:
        # None rather than a partial sum
        if any(s.dg is None for s in self.steps):
            return None
        return sum(s.dg for s in self.steps)


def _check(nodes: list[Mol], steps: list[Step]) -> None:
    if len(nodes) < 2:
        raise SpecError(f"need >= 2 nodes, got {len(nodes)}")
    if len(nodes) != len(steps):
        raise SpecError(f"{len(nodes)} nodes but {len(steps)} steps; need one step per node")
    ids = [s.rid for s in steps]
    if len(set(ids)) != len(ids):
        raise SpecError(f"duplicate reaction ids: {sorted({i for i in ids if ids.count(i) > 1})}")


def drop_side(obj, smiles: list[str]):
    """Remove named side species from every step."""
    drop = {canonical(s) for s in smiles}
    steps = list(obj.steps)
    for sub in getattr(obj, "subs", []):
        steps += sub.steps
    for st in steps:
        st.consumes = [x for x in st.consumes if x.smiles not in drop]
        st.produces = [x for x in st.produces if x.smiles not in drop]
    return obj


@dataclass
class PathNode:
    """A molecule in a route, with the reaction producing it from `precursors`."""

    mol: Mol
    step: Step | None = None
    precursors: list[PathNode] = field(default_factory=list)
    terminal: str | None = None
    generation: str | None = None

    def __post_init__(self):
        if self.step is None and self.precursors:
            raise SpecError(f"{self.mol.smiles}: has precursors but no reaction producing it")
        if self.step is not None and not self.precursors:
            raise SpecError(f"{self.mol.smiles}: reaction {self.step.rid} has no precursors")
        if self.precursors:
            self.terminal = None
        elif self.terminal is None:
            # unknown, never assumed to be a seed
            self.terminal = UNKNOWN
        if self.terminal not in (None, SEED, UNTRACED, UNKNOWN):
            raise SpecError(
                f"{self.mol.smiles}: terminal must be one of "
                f"{SEED!r}, {UNTRACED!r}, {UNKNOWN!r}; got {self.terminal!r}"
            )

    @property
    def seed(self) -> bool:
        return self.terminal == SEED

    def walk(self):
        yield self
        for p in self.precursors:
            yield from p.walk()

    @property
    def depth(self) -> int:
        return 1 + max((p.depth for p in self.precursors), default=-1)


@dataclass
class Pathway:
    """A route: target at the root, precursors branching back."""

    root: PathNode
    title: str | None = None
    reported_dg: float | None = None   # stated by the source, not recomputed here
    meta: dict = field(default_factory=dict)

    def __post_init__(self):
        # one reaction may legitimately appear at several places in a route
        _check_acyclic(self.root, ())

    @property
    def nodes(self) -> list[PathNode]:
        return list(self.root.walk())

    @property
    def steps(self) -> list[Step]:
        return [n.step for n in self.nodes if n.step]

    @property
    def leaves(self) -> list[PathNode]:
        return [n for n in self.nodes if not n.precursors]

    @property
    def total_dg(self) -> float | None:
        dgs = [s.dg for s in self.steps]
        if not dgs or any(d is None for d in dgs):
            return None
        return sum(dgs)

    @property
    def gain_steps(self) -> list[int]:
        return [i for i, s in enumerate(self.steps) if s.gain]

    @property
    def seeds(self) -> list[PathNode]:
        return [n for n in self.leaves if n.terminal == SEED]

    @property
    def dead_ends(self) -> list[PathNode]:
        """Leaves known not to be seeds."""
        return [n for n in self.leaves if n.terminal == UNTRACED]

    @property
    def unknown_leaves(self) -> list[PathNode]:
        """Leaves whose status the source never stated."""
        return [n for n in self.leaves if n.terminal == UNKNOWN]

    @property
    def complete(self) -> bool:
        """Every leaf is a known seed."""
        return bool(self.leaves) and all(n.terminal == SEED for n in self.leaves)


def _check_acyclic(node: PathNode, seen: tuple[str, ...]) -> None:
    """A molecule may not reappear on its own root-to-leaf path."""
    if node.mol.smiles in seen:
        chain = " <- ".join([*seen, node.mol.smiles])
        raise SpecError(f"route is not acyclic: {chain}")
    for p in node.precursors:
        _check_acyclic(p, (*seen, node.mol.smiles))


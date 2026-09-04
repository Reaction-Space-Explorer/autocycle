"""Check the conditions that make a cycle autocatalytic rather than simple.

Orgel's distinction, as stated by Pereto (Chem. Soc. Rev. 2012): a *simple* cycle
regenerates one reactant stoichiometrically, so the molecule produced replaces the one
consumed; an *autocatalytic* cycle "exhibits an additional yield of the feeder", n > 1.

The distinction is therefore a stoichiometric coefficient. Where a source records no
coefficients the condition is reported `unknown`, never assumed either way.
"""

from __future__ import annotations

from dataclasses import dataclass

from autocycle.spec import Cycle

YES, NO, UNKNOWN = "yes", "no", "unknown"

AUTOCATALYTIC = "autocatalytic"   # an extra copy is stated stoichiometrically
TOPOLOGICAL = "topological"       # a shunt gives the published topological criterion
SIMPLE = "simple"
CANDIDATE = "candidate"
INCOMPLETE = "incomplete"


@dataclass(frozen=True)
class Verdict:
    conditions: dict[str, str]
    seed_yield: float | None
    declared_gain: bool

    @property
    def status(self) -> str:
        c = self.conditions
        structural = (c["seed_identified"], c["feeder"], c["seed_regenerated"])
        if any(x != YES for x in structural):
            return INCOMPLETE
        if c["extra_yield"] == YES:
            return AUTOCATALYTIC
        if c["shunt"] == YES:
            return TOPOLOGICAL
        if c["extra_yield"] == NO:
            return SIMPLE
        return CANDIDATE

    @property
    def disagrees_with_declaration(self) -> bool:
        """A `gain` flag that the conditions do not support."""
        return self.declared_gain and self.conditions["extra_yield"] != YES

    def summary(self) -> str:
        bits = ", ".join(f"{k}={v}" for k, v in self.conditions.items())
        return f"{self.status} ({bits})"


def verify(cycle: Cycle) -> Verdict:
    ring = {m.smiles for m in cycle.nodes}
    # the shunt is part of the cycle, and is usually where the extra copy is made
    steps = list(cycle.steps) + (list(cycle.shunt.steps) if cycle.shunt else [])
    ring |= {m.smiles for m in (cycle.shunt.nodes if cycle.shunt else [])}
    consumed = [sp for st in steps for sp in st.consumes]
    produced = [sp for st in steps for sp in st.produces]

    seed = cycle.nodes[cycle.seed].smiles if cycle.seed is not None else None
    off_ring_in = [sp for sp in consumed if sp.smiles not in ring]
    off_ring_out = [sp for sp in produced if sp.smiles not in ring]

    # the ring closure returns one seed; a further copy must be produced explicitly
    extra = sum(sp.count for sp in produced if seed is not None and sp.smiles == seed)
    seed_yield = 1.0 + extra if seed is not None else None

    conditions = {
        "seed_identified": YES if seed is not None else NO,
        "feeder": YES if off_ring_in else NO,
        "outlet": YES if off_ring_out else NO,
        "seed_regenerated": YES if seed is not None else NO,
        # without coefficients the yield cannot be pinned to 1
        "extra_yield": (
            YES if extra > 0 else (NO if cycle.stoichiometry_complete else UNKNOWN)
        ),
        # a bridging path back to the seed: topological, not mass-balanced
        "shunt": YES if cycle.shunt is not None else NO,
    }
    return Verdict(
        conditions=conditions,
        seed_yield=seed_yield,
        declared_gain=bool(cycle.gain_steps),
    )

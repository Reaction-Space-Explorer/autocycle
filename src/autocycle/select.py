"""Rank and vet cycles for presentation."""

from __future__ import annotations

from rdkit import Chem
from rdkit.Chem.Descriptors import MolWt

from autocycle.spec import Cycle, canonical

# molecules that are chemically implausible as ring intermediates without a catalyst,
# and belong on a cycle only as a feeder or a product
ROLE_RESTRICTED = ("C=O", "CO")


def feeders(cycle: Cycle) -> set[str]:
    """Distinct consumed species."""
    return {sp.smiles for st in cycle.steps for sp in st.consumes}


def consumers(cycle: Cycle) -> set[str]:
    return {sp.smiles for st in cycle.steps for sp in st.produces}


def min_feeder_weight(cycle: Cycle) -> float:
    ws = [MolWt(Chem.MolFromSmiles(s)) for s in feeders(cycle)]
    return min(ws) if ws else float("inf")


def rank_key(cycle: Cycle) -> tuple:
    """Fewest distinct feeders first, then lowest feeder weight, then shortest ring."""
    return (len(feeders(cycle)), min_feeder_weight(cycle), len(cycle.nodes))


def role_violations(cycle: Cycle, restricted=ROLE_RESTRICTED) -> list[str]:
    """Restricted molecules sitting on the ring rather than feeding it."""
    bad = {canonical(s) for s in restricted}
    return [m.smiles for m in cycle.nodes if m.smiles in bad]

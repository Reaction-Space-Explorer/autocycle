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


def min_feeder_generation(cycle: Cycle) -> float:
    """Earliest generation any feeder first appears in, or inf if unrecorded."""
    gens = [sp.generation for st in cycle.steps for sp in st.consumes if sp.generation is not None]
    return min(gens) if gens else float("inf")


def rank_key(cycle: Cycle) -> tuple:
    """Fewest feeders, then lightest feeder, then earliest feeder generation."""
    return (
        len(feeders(cycle)),
        min_feeder_weight(cycle),
        min_feeder_generation(cycle),
        len(cycle.nodes),
    )


def role_violations(cycle: Cycle, restricted=ROLE_RESTRICTED) -> list[str]:
    """Restricted molecules sitting on the ring rather than feeding it."""
    bad = {canonical(s) for s in restricted}
    return [m.smiles for m in cycle.nodes if m.smiles in bad]


def centrality(graph, measure: str = "degree") -> dict[str, float]:
    """Node centrality over the whole network. Degree is cheap; betweenness is not."""
    import networkx as nx

    g = nx.DiGraph(graph)
    if measure == "degree":
        return nx.degree_centrality(g)
    if measure == "betweenness":
        return nx.betweenness_centrality(g)
    raise ValueError(f"unknown measure {measure!r}")


def cycle_centrality(cycle: Cycle, cent: dict[str, float]) -> float | None:
    """Lowest centrality among the ring molecules (Zubarev's cycle centrality)."""
    vals = [cent[m.smiles] for m in cycle.nodes if m.smiles in cent]
    return min(vals) if len(vals) == len(cycle.nodes) else None

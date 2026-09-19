"""Classify a core by what its steps do, not by what its rules are called.

Rule names alias: the same transformation appears as Aldol Condensation and as
Knoevenagel C, so counting motifs by rule name inflates them. A step is described
here by the formula change it makes to the ring plus the formulas it takes in and
lets out, which is substrate-agnostic and alias-proof: every aldol addition of
formaldehyde reads the same whatever sugar it happens on, and whatever the rule
that generated it was called.
"""

from __future__ import annotations

import collections

from rdkit import Chem
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

_F: dict[str, str] = {}


def formula(smiles: str) -> str:
    if smiles not in _F:
        m = Chem.MolFromSmiles(smiles)
        _F[smiles] = CalcMolFormula(m) if m else smiles
    return _F[smiles]


def _counts(f: str) -> collections.Counter:
    return collections.Counter(
        {e: int(n or 1) for e, n in __import__("re").findall(r"([A-Z][a-z]?)(\d*)", f) if e}
    )


def delta(before: str, after: str) -> str:
    """The formula change from one ring species to the next, as a signed string."""
    d = _counts(formula(after))
    d.subtract(_counts(formula(before)))
    parts = [f"{e}{v:+d}" for e, v in sorted(d.items()) if v]
    return "".join(parts) or "0"


def step_kind(by_rxn, r, src, tgt) -> str:
    """What one step does: the ring change, and the species in and out."""
    d = by_rxn[r]
    ins = sorted(formula(s) for s, c in d.items() if c < 0 and s != src for _ in range(-c))
    outs = sorted(formula(s) for s, c in d.items() if c > 0 and s != tgt for _ in range(c))
    core = delta(src, tgt)
    return f"{core}|in:{','.join(ins) or '-'}|out:{','.join(outs) or '-'}"


def motif(by_rxn, species, rxns) -> tuple[str, ...]:
    n = len(species)
    return tuple(sorted(
        step_kind(by_rxn, rxns[i], species[i], species[(i + 1) % n]) for i in range(n)
    ))


def carbons(smiles: str) -> int:
    m = Chem.MolFromSmiles(smiles)
    return sum(1 for a in m.GetAtoms() if a.GetSymbol() == "C") if m else 0


def coarse_step(by_rxn, r, src, tgt) -> str:
    """The step as carbon bookkeeping: what the ring gains or loses, and to whom.

    This is the mechanism level. An aldol addition of formaldehyde onto a triose
    and onto a hexose read the same, because they are the same transformation.
    """
    d = by_rxn[r]
    ins = sorted(carbons(s) for s, c in d.items() if c < 0 and s != src for _ in range(-c))
    outs = sorted(carbons(s) for s, c in d.items() if c > 0 and s != tgt for _ in range(c))
    n = carbons(tgt) - carbons(src)
    return (f"C{n:+d}|in:{','.join(f'C{i}' for i in ins) or '-'}"
            f"|out:{','.join(f'C{i}' for i in outs) or '-'}")


def coarse_motif(by_rxn, species, rxns) -> tuple[str, ...]:
    n = len(species)
    return tuple(sorted(
        coarse_step(by_rxn, rxns[i], species[i], species[(i + 1) % n]) for i in range(n)
    ))

"""Stoichiometric network analysis of one cycle, after Clarke.

Clarke (Cell Biophys. 12, 237-253, 1988) writes the steady state of a mechanism as
nu.v = 0 with v >= 0. The solutions form a convex polyhedral cone whose edges are the
extreme currents, and every steady state is a non-negative sum of them. Read on a single
drawn cycle this asks two questions `verify` does not: whether the intermediates on the
ring actually balance, and whether the drawing is one extreme current or two cycles
sharing a picture.

The ring and its shunt are one current. A fused sub-cycle is a separate current sharing a
molecule, so it is not folded in here; `autocycle verify` says so when one is present.
"""

from __future__ import annotations

import itertools
from dataclasses import dataclass

import numpy as np

from autocycle.spec import Cycle, SpecError


@dataclass(frozen=True)
class Current:
    """The cycle as a reaction current: one turn fires every step once."""

    species: list[str]
    matrix: np.ndarray     # nu, species x steps
    flux: np.ndarray       # v, one entry per step
    internal: list[str]    # held at steady state: ring and shunt, minus the seed
    seed: str | None

    @property
    def net(self) -> dict[str, float]:
        return dict(zip(self.species, self.matrix @ self.flux, strict=True))

    @property
    def imbalanced(self) -> list[str]:
        """Intermediates the drawing does not return, so nu.v != 0."""
        net = self.net
        return [s for s in self.internal if abs(net[s]) > 1e-9]

    @property
    def extreme(self) -> bool:
        """An edge of the cone: no sub-current of it is also a steady state.

        For a current with full support the rank criterion is exact, so this needs no
        cone enumeration.
        """
        support = [j for j, f in enumerate(self.flux) if f > 1e-9]
        if not self.internal:
            return len(support) == 1
        rows = [self.species.index(s) for s in self.internal]
        m = self.matrix[np.ix_(rows, support)]
        return int(np.linalg.matrix_rank(m)) == len(support) - 1

    @property
    def cone_dim(self) -> int:
        """Independent currents through these steps. One means a single extreme current."""
        rows = [self.species.index(s) for s in self.internal]
        if not rows:
            return self.matrix.shape[1]
        return self.matrix.shape[1] - int(np.linalg.matrix_rank(self.matrix[rows, :]))

    def decompose(self) -> list[np.ndarray]:
        """The extreme currents this one is made of, as flux vectors.

        Schmitz's method: hold `steps - rank - 1` components at zero, set one to 1 to
        avoid the trivial solution, solve the square system, and keep the solutions with
        no negative component.

        The number held is `cone_dim - 1`, so the enumeration is exhaustive over a
        handful of choices for any cycle that can be drawn as one.
        """
        rows = [self.species.index(s) for s in self.internal]
        nu = self.matrix[rows, :] if rows else self.matrix[:0, :]
        n = nu.shape[1]
        held = n - int(np.linalg.matrix_rank(nu)) - 1
        if held < 0:
            return []
        out: list[np.ndarray] = []
        for zeros in itertools.combinations(range(n), held):
            for one in range(n):
                if one in zeros:
                    continue
                rows_ = [nu[i] for i in range(nu.shape[0])]
                for z in zeros:
                    e = np.zeros(n); e[z] = 1.0; rows_.append(e)
                e = np.zeros(n); e[one] = 1.0; rows_.append(e)
                a = np.array(rows_)
                if a.shape[0] != n:
                    continue
                b = np.zeros(n); b[-1] = 1.0
                try:
                    v = np.linalg.solve(a, b)
                except np.linalg.LinAlgError:
                    continue
                if (v < -1e-9).any() or v.sum() <= 1e-9:
                    continue
                v = v / v[v > 1e-9].min()
                if not any(np.allclose(v, o, atol=1e-6) for o in out):
                    out.append(v)
        return out

    @property
    def overall(self) -> dict[str, int]:
        """The overall reaction left once the intermediates are eliminated."""
        return {
            s: c for s, c in self.net.items() if s not in self.internal and abs(c) > 1e-9
        }

    @property
    def atom_residual(self) -> dict[str, int] | None:
        """Atoms of the overall reaction, products minus reactants.

        All zero means the drawing states a closed balance. A drawing that suppresses
        water, protons or redox partners shows a residual in H and O; a residual in C
        usually means a carboxylation was left out. `*` counts pseudo-atom stubs.
        None where a species is a name rather than a structure, so there is nothing
        to count.
        """
        from rdkit import Chem

        out: dict[str, float] = {}
        for smi, n in self.overall.items():
            mol = Chem.MolFromSmiles(smi)
            if mol is None:
                return None
            for atom in Chem.AddHs(mol).GetAtoms():
                out[atom.GetSymbol()] = out.get(atom.GetSymbol(), 0.0) + n
        return {k: int(round(v)) for k, v in out.items() if abs(v) > 1e-9}

    @property
    def seed_yield(self) -> float | None:
        return None if self.seed is None else 1.0 + self.net[self.seed]


def current(cycle: Cycle) -> Current:
    steps = list(cycle.steps)
    n = len(cycle.nodes)
    edges = [(cycle.nodes[i].smiles, cycle.nodes[(i + 1) % n].smiles) for i in range(n)]

    if cycle.shunt is not None:
        if cycle.seed is None:
            raise SpecError("a shunt needs a seed to return to")
        chain = (
            [cycle.nodes[cycle.shunt.from_node].smiles]
            + [m.smiles for m in cycle.shunt.nodes]
            + [cycle.nodes[cycle.seed].smiles]
        )
        if len(cycle.shunt.steps) != len(chain) - 1:
            raise SpecError(
                f"shunt has {len(cycle.shunt.steps)} steps for {len(chain) - 1} edges"
            )
        steps += list(cycle.shunt.steps)
        edges += list(zip(chain, chain[1:], strict=False))

    species: list[str] = []
    for a, b in edges:
        for s in (a, b):
            if s not in species:
                species.append(s)
    for st in steps:
        for sp in list(st.consumes) + list(st.produces):
            if sp.smiles not in species:
                species.append(sp.smiles)

    nu = np.zeros((len(species), len(steps)), dtype=float)
    at = {s: i for i, s in enumerate(species)}
    for j, (st, (a, b)) in enumerate(zip(steps, edges, strict=True)):
        nu[at[a], j] -= 1
        nu[at[b], j] += 1
        for sp in st.consumes:
            nu[at[sp.smiles], j] -= sp.count
        for sp in st.produces:
            nu[at[sp.smiles], j] += sp.count

    seed = cycle.nodes[cycle.seed].smiles if cycle.seed is not None else None
    ring = [m.smiles for m in cycle.nodes]
    ring += [m.smiles for m in (cycle.shunt.nodes if cycle.shunt else [])]
    internal = [s for s in dict.fromkeys(ring) if s != seed]
    flux = np.array([st.flux for st in steps], dtype=float)
    return Current(
        species=species, matrix=nu, flux=flux, internal=internal, seed=seed
    )

"""Autocatalytic core test, after Blokhuis, Lacoste and Nghe (PNAS 2020, 117, 25230).

A core is n species and n reactions whose stoichiometric submatrix is square and
invertible and which admits a strictly positive flux making every core species
accumulate: v > 0 with nu @ v > 0. The test is a linear program, since the system
is homogeneous: maximise t subject to v >= t, nu v >= t, sum v = 1. A core exists
iff the optimum is positive.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import linprog


def positive_flux(nu: np.ndarray, tol: float = 1e-7) -> float:
    """The largest t with v >= t, nu v >= t, sum v = 1. Positive means autocatalytic."""
    n, m = nu.shape
    # variables [v (m), t (1)]; maximise t
    c = np.zeros(m + 1)
    c[-1] = -1.0
    A = np.zeros((m + n, m + 1))
    A[:m, :m] = -np.eye(m)          # -v + t <= 0
    A[:m, -1] = 1.0
    A[m:, :m] = -nu                 # -nu v + t <= 0
    A[m:, -1] = 1.0
    res = linprog(c, A_ub=A, b_ub=np.zeros(m + n),
                  A_eq=np.concatenate([np.ones(m), [0.0]])[None, :], b_eq=[1.0],
                  bounds=[(0, None)] * m + [(None, None)], method="highs")
    return float(res.x[-1]) if res.success else -1.0


def is_core(nu: np.ndarray, tol: float = 1e-7) -> bool:
    n, m = nu.shape
    if n != m or abs(np.linalg.det(nu)) < 1e-9:
        return False
    return positive_flux(nu) > tol


def core_type(nu: np.ndarray) -> str:
    """Blokhuis's types, by the number of graph cycles the core contains.

    PNAS 2020, 117, 25230: "The five types differ in their number of graph cycles
    and the way these cycles overlap. Type I consists of a single graph cycle that
    is weight asymmetric... Types II and III comprise two distinct but overlapping
    graph cycles, Type IV comprises three, and Type V more than three."

    The graph has the core species as nodes and an edge from each reaction's core
    reactant to each of its core products. II and III both have two cycles and are
    reported together, since they are told apart by how the cycles overlap rather
    than by how many there are.
    """
    import networkx as nx

    g = nx.DiGraph()
    n = nu.shape[0]
    g.add_nodes_from(range(n))
    for j in range(nu.shape[1]):
        src = [i for i in range(n) if nu[i, j] < 0]
        for s in src:
            for t in [i for i in range(n) if nu[i, j] > 0]:
                g.add_edge(s, t)
    k = sum(1 for _ in nx.simple_cycles(g))
    return {1: "I", 2: "II/III", 3: "IV"}.get(k, "V" if k > 3 else "none")

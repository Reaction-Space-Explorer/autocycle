"""Blokhuis's second criterion for the type, and whether it agrees with the first.

The paper gives two independent routes to a core's type. One counts graph cycles:
"Type I consists of a single graph cycle that is weight asymmetric... Types II and
III comprise two distinct but overlapping graph cycles, Type IV comprises three,
and Type V more than three." The other counts forks, a fork being "a reaction with
a single reactant and two products": "In type I cores, the fork ends with two
copies of the same compound, whereas, in types II to V, forks end with different
compounds. It is obligatory for types II and III to contain one fork with two
distinct products, for type IV to contain two such forks, and for type V to
contain three."

cores.core_type implements the first. This implements the second. They are stated
as equivalent, so disagreement is a finding either way.
"""
from __future__ import annotations

import numpy as np


def forks(nu: np.ndarray) -> tuple[int, int]:
    """(forks onto two copies of one species, forks onto two distinct species)."""
    same = distinct = 0
    for j in range(nu.shape[1]):
        col = nu[:, j]
        neg = np.flatnonzero(col < 0)
        pos = np.flatnonzero(col > 0)
        if len(neg) != 1 or col[neg[0]] != -1:
            continue
        if len(pos) == 1 and col[pos[0]] == 2:
            same += 1
        elif len(pos) == 2 and col[pos].sum() == 2:
            distinct += 1
    return same, distinct


def fork_type(nu: np.ndarray) -> str:
    same, dist = forks(nu)
    if dist == 0:
        return "I" if same else "none"
    return {1: "II/III", 2: "IV"}.get(dist, "V" if dist >= 3 else "none")

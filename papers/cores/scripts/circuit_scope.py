"""What the circuit search covers, and what it does not.

The search builds one closed ring of n species and n reactions and tests its
submatrix. A core need not have that shape. Golnik and co-workers show that cores
are MR-chordless circuits and their superpositions, so a core built from two or
more circuits is outside this search by construction, and a circuit in which one
reaction consumes two core species fails their unique-reactant condition and is
not a core at all.

This compares the circuit search against an exhaustive one that assumes neither:
every connected species set of size n, every assignment of a distinct consuming
reaction to each species, subject to that reaction consuming exactly one species
of the set. Strong connectivity is a property of cores, so the connectivity prune
keeps the comparison exhaustive.
"""
import collections
import itertools

import numpy as np
from _common import FOOD

from autocycle.cores.anchored import enumerate_cores
from autocycle.cores.cores import is_core
from autocycle.cores.enumerate_cores import load
from autocycle.cores.paths import RELS

NETS = [("glucose G2", "Glucose/GlucoseRels_2.tsv"),
        ("glucose+ammonia G2", "GlucoseAmm/GlucoseAmmRels_2.tsv"),
        ("glucose G3", "Glucose/GlucoseRels_3.tsv"),
        ("formose+ammonia G3", "FormoseAmm/FormoseAmmRels_3.tsv")]


def every_core(by, n):
    sp = sorted({s for d in by.values() for s in d if s not in FOOD})
    cons = {s: [r for r, d in by.items() if d.get(s, 0) < 0] for s in sp}
    adj = collections.defaultdict(set)
    for d in by.values():
        a = [s for s in d if d[s] < 0 and s not in FOOD]
        b = [s for s in d if d[s] > 0 and s not in FOOD]
        for x in a:
            for y in b:
                if x != y:
                    adj[x].add(y)
                    adj[y].add(x)
    subs, seen = set(), set()

    def grow(cur, frontier):
        if len(cur) == n:
            subs.add(frozenset(cur))
            return
        for y in sorted(frontier):
            key = (frozenset(cur), y)
            if key in seen:
                continue
            seen.add(key)
            grow(cur | {y}, (frontier | adj[y]) - cur - {y})

    for s in sp:
        grow({s}, adj[s])
    out = set()
    for X in subs:
        X = sorted(X)
        Xs = set(X)
        opts = []
        for x in X:
            o = [r for r in cons[x]
                 if sum(1 for y in Xs if by[r].get(y, 0) < 0) == 1
                 and any(by[r].get(y, 0) > 0 for y in Xs)]
            if not o:
                break
            opts.append(o)
        if len(opts) != n:
            continue
        for R in itertools.product(*opts):
            if len(set(R)) != n:
                continue
            nu = np.array([[by[r].get(s, 0) for r in R] for s in X], float)
            if is_core(nu):
                out.add((frozenset(X), frozenset(R)))
    return out


def chorded(by, X, R):
    X, R = sorted(X), sorted(R)
    nu = np.array([[by[r].get(s, 0) for r in R] for s in X], float)
    return any((nu[:, j] < 0).sum() > 1 for j in range(nu.shape[1]))


print(f"  {'network':<20} {'n':>2} {'all cores':>10} {'circuits':>9} "
      f"{'missed':>7} {'extra':>6} {'of which chorded':>17}")
for name, rel in NETS:
    by = load(RELS / rel)
    for n in (3, 4):
        ex = every_core(by, n)
        ours = {(frozenset(s), frozenset(r)) for s, r, _ in enumerate_cores(by, n, food=FOOD)[0]}
        extra = ours - ex
        ch = sum(chorded(by, X, R) for X, R in extra)
        print(f"  {name:<20} {n:>2} {len(ex):>10} {len(ours):>9} "
              f"{len(ex - ours):>7} {len(extra):>6} {ch:>17}", flush=True)

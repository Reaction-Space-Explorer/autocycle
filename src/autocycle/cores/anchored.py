"""Anchored core enumeration, streamed.

A cycle whose every reaction consumes one core species and produces one has
column sums of zero over the core, so the all-ones row annihilates it and the
matrix is singular. Every autocatalytic core therefore contains at least one
*amplifying* reaction, one with a net gain in core species, and those are a fifth
of a formose network. Anchoring on them and meeting in the middle replaces
walking every cycle.

Each cycle is emitted only from the smallest amplifying reaction it contains, so
it is generated once and no candidate list is held: cycles are verified in
chunks, a batched determinant rejects the singular ones, and the linear program
runs only on what survives.
"""

from __future__ import annotations

import collections
import sys
import time

import numpy as np

from autocycle.cores.cores import core_type, is_core
from autocycle.cores.enumerate_cores import FOOD, load

CHUNK = 200_000


def graph(by_rxn, food=FOOD):
    succ, pred = collections.defaultdict(list), collections.defaultdict(list)
    for r, d in by_rxn.items():
        cons = [s for s, c in d.items() if c < 0 and s not in food]
        prod = [s for s, c in d.items() if c > 0 and s not in food]
        for s in cons:
            for p in prod:
                if s != p:
                    succ[s].append((p, r))
                    pred[p].append((s, r))
    return succ, pred


def anchors(by_rxn, food=FOOD):
    """Reactions that could amplify inside some core.

    Lemma 2's positive column sum is over the core species, not over all non-food
    species, and a reaction can amplify inside a core while losing that sum to
    non-food species outside it, so the wider test drops cores. A reaction in a
    core consumes a core species at -1 or less, so its core sum is at most P - 1
    where P is what it produces, and amplifying implies P >= 2.
    """
    return [r for r, d in by_rxn.items()
            if sum(c for s, c in d.items() if s not in food and c > 0) >= 2]


def from_anchor(r0, by_rxn, n, food, succ, pred, amp):
    """Yield the n-cycles anchored at r0, once each.

    Split out of candidates() so the sharded enumerator runs the same unrolled
    walk rather than the general one, which is several times slower at this size.
    """
    def first(rxns):
        """Keep the cycle only at its smallest amplifying reaction."""
        return min(r for r in rxns if r in amp) == r0

    d = by_rxn[r0]
    starts = [s for s, c in d.items() if c > 0 and s not in food]
    ends = [s for s, c in d.items() if c < 0 and s not in food]
    for v1 in starts:
        for vn in ends:
            if n == 2:
                if v1 == vn:
                    continue
                for p, r in succ.get(v1, ()):
                    if p == vn and r != r0 and first((r, r0)):
                        yield (v1, vn), (r, r0)
            elif n == 3:
                tail = collections.defaultdict(list)
                for q, r in pred.get(vn, ()):
                    tail[q].append(r)
                for v2, r1 in succ.get(v1, ()):
                    if v2 in (v1, vn) or r1 == r0:
                        continue
                    for r2 in tail.get(v2, ()):
                        if r2 not in (r0, r1) and first((r1, r2, r0)):
                            yield (v1, v2, vn), (r1, r2, r0)
            elif n == 4:
                mid = collections.defaultdict(list)
                for q, r in pred.get(vn, ()):
                    mid[q].append(r)
                for v2, r1 in succ.get(v1, ()):
                    if v2 in (v1, vn) or r1 == r0:
                        continue
                    for v3, r2 in succ.get(v2, ()):
                        if v3 in (v1, v2, vn) or r2 in (r0, r1):
                            continue
                        for r3 in mid.get(v3, ()):
                            if r3 in (r0, r1, r2):
                                continue
                            if first((r1, r2, r3, r0)):
                                yield (v1, v2, v3, vn), (r1, r2, r3, r0)
            else:
                raise ValueError(f"cycle length {n} not implemented")


def candidates(by_rxn, n, food=FOOD):
    """Yield each n-species cycle containing an amplifying reaction, once."""
    succ, pred = graph(by_rxn, food)
    amp = set(anchors(by_rxn, food))
    for r0 in sorted(amp):
        yield from from_anchor(r0, by_rxn, n, food, succ, pred, amp)

def _check(by_rxn, batch):
    """Batched determinant first, then the program only for what survives."""
    if not batch:
        return []
    n = len(batch[0][0])
    nus = np.empty((len(batch), n, n))
    for i, (species, rxns) in enumerate(batch):
        for a, s in enumerate(species):
            for b, r in enumerate(rxns):
                nus[i, a, b] = by_rxn[r].get(s, 0)
    out = []
    for i in np.flatnonzero(np.abs(np.linalg.det(nus)) > 1e-9):
        if is_core(nus[i]):
            out.append((*batch[i], core_type(nus[i])))
    return out


def enumerate_cores(by_rxn, n, *, food=FOOD, chunk=CHUNK):
    """Stream the candidates and keep only the cores. Memory is O(chunk)."""
    found, batch, seen = [], [], 0
    for cand in candidates(by_rxn, n, food):
        batch.append(cand)
        if len(batch) == chunk:
            seen += len(batch)
            found += _check(by_rxn, batch)
            batch = []
    seen += len(batch)
    found += _check(by_rxn, batch)
    return found, seen


def distinct(by_rxn, found):
    sig = {r: tuple(sorted(d.items())) for r, d in by_rxn.items()}
    return {(frozenset(s), frozenset(sig[r] for r in rx)) for s, rx, _ in found}


if __name__ == "__main__":
    path, n = sys.argv[1], int(sys.argv[2])
    t0 = time.time()
    by_rxn = load(path)
    found, seen = enumerate_cores(by_rxn, n)
    print(f"  food {sorted(FOOD)}")      # a count without its food set is arbitrary
    print(f"  {len(by_rxn)} reactions | {len(anchors(by_rxn))} amplifying")
    print(f"  candidates streamed: {seen}")
    print(f"  cores {len(found)} | distinct {len(distinct(by_rxn, found))}  "
          f"total {time.time() - t0:.1f}s")
    print(f"  types {dict(sorted(collections.Counter(t for _, _, t in found).items()))}")

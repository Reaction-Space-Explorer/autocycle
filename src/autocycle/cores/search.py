"""General-n core search: meet in the middle, prune on the canonical anchor.

anchored.py hand-unrolls n = 2, 3 and 4 and refuses anything longer, which is why
the RAF check had to carry its own walker. The same idea works for any n. A cycle
is an anchor reaction r0 plus a simple path of n-1 edges from one of its products
to one of its reactants, so the path is built from both ends and joined in the
middle: ceil((n-1)/2) edges forward, the rest backward.

Two prunes. A cycle is emitted only at its smallest amplifying reaction, so a
partial path that already carries an amplifying reaction below the anchor is
abandoned where it is found rather than at the end. And the forward half is
indexed by its endpoint, so the backward half only visits midpoints that are
reachable at all.
"""

from __future__ import annotations

import collections

from autocycle.cores.anchored import _check, anchors, graph
from autocycle.cores.enumerate_cores import FOOD

CHUNK = 200_000


def _half(adj, start, k, amp, r0, block):
    """Simple paths of k edges from start. Yields (endpoint, species, rxns)."""
    out = collections.defaultdict(list)

    def walk(node, path, rxns):
        if len(rxns) == k:
            out[node].append((tuple(path), tuple(rxns)))
            return
        for nxt, r in adj.get(node, ()):
            if r == r0 or r in rxns or nxt in path or nxt in block:
                continue
            if r in amp and r < r0:            # not the canonical anchor
                continue
            walk(nxt, path + [nxt], rxns + [r])

    walk(start, [start], [])
    return out


def candidates(by_rxn, n, food=FOOD):
    if n < 2:
        raise ValueError("a core needs at least two species")
    succ, pred = graph(by_rxn, food)
    amp = set(anchors(by_rxn, food))
    fwd_k = n // 2
    bwd_k = (n - 1) - fwd_k

    for r0 in sorted(amp):
        d = by_rxn[r0]
        starts = [s for s, c in d.items() if c > 0 and s not in food]
        ends = [s for s, c in d.items() if c < 0 and s not in food]
        back = {vn: _half(pred, vn, bwd_k, amp, r0, frozenset()) for vn in ends}
        for v1 in starts:
            fwd = _half(succ, v1, fwd_k, amp, r0, frozenset())
            if not fwd:
                continue
            for vn in ends:
                if vn == v1:
                    continue
                for mid, tails in back[vn].items():
                    for fs, fr in fwd.get(mid, ()):
                        for bs, br in tails:
                            tail = tuple(reversed(bs))[1:]
                            if set(br) & set(fr) or set(tail) & set(fs):
                                continue
                            spec = fs + tail
                            if len(set(spec)) != n:
                                continue
                            yield spec, fr + tuple(reversed(br)) + (r0,)


def enumerate_cores(by_rxn, n, chunk=CHUNK, food=FOOD):
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

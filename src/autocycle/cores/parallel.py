"""Shard the enumeration over anchor reactions.

Each cycle is emitted only at its smallest amplifying reaction, so the anchors
partition the search exactly: shard by anchor and no cycle is found twice and
none is missed. The shards share nothing, so this is a process pool with no
communication beyond the cores each one returns.
"""

from __future__ import annotations

import collections
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

from autocycle.cores.anchored import _check, anchors, distinct, from_anchor, graph
from autocycle.cores.anchored import enumerate_cores as single
from autocycle.cores.enumerate_cores import FOOD, load
from autocycle.cores.search import _half

_STATE = {}


def _init(by_rxn, n, food):
    succ, pred = graph(by_rxn, food)
    _STATE.update(by=by_rxn, n=n, food=food, succ=succ, pred=pred,
                  amp=set(anchors(by_rxn, food)))


def _one(r0):
    """One shard: the cycles anchored at r0, verified.

    Uses the unrolled walk up to four species and the general meet-in-the-middle
    above it. The two agree core for core; the unrolled one is several times
    faster at the sizes the papers report, and sharding used to give that up.
    """
    by, n, food = _STATE["by"], _STATE["n"], _STATE["food"]
    succ, pred, amp = _STATE["succ"], _STATE["pred"], _STATE["amp"]
    if n <= 4:
        batch = list(from_anchor(r0, by, n, food, succ, pred, amp))
        return _check(by, batch), len(batch)
    fwd_k = n // 2
    bwd_k = (n - 1) - fwd_k
    d = by[r0]
    starts = [s for s, c in d.items() if c > 0 and s not in food]
    ends = [s for s, c in d.items() if c < 0 and s not in food]
    back = {vn: _half(pred, vn, bwd_k, amp, r0, frozenset()) for vn in ends}
    batch = []
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
                        batch.append((spec, fr + tuple(reversed(br)) + (r0,)))
    return _check(by, batch), len(batch)

def enumerate_cores(by_rxn, n, food=FOOD, workers=None):
    workers = workers or os.cpu_count()
    amp = sorted(anchors(by_rxn, food))
    found, seen = [], 0
    with Pool(workers, initializer=_init, initargs=(by_rxn, n, food)) as pool:
        for hits, k in pool.imap_unordered(_one, amp, chunksize=64):
            found += hits
            seen += k
    return found, seen


if __name__ == "__main__":
    path, n = sys.argv[1], int(sys.argv[2])
    w = int(sys.argv[3]) if len(sys.argv) > 3 else None
    food = {"O", "C=O", "C(=O)=O", "N"}
    by = load(path)
    t0 = time.time()
    found, seen = enumerate_cores(by, n, food, w)
    t = collections.Counter(x[2] for x in found)
    print(f"  {Path(path).stem}  n={n}  workers={w or os.cpu_count()}")
    print(f"    reactions {len(by)}  anchors {len(anchors(by, food))}  "
          f"candidates {seen}")
    print(f"    cores {len(found)}  distinct {len(distinct(by, found))}  "
          f"types {dict(sorted(t.items()))}  {time.time()-t0:.1f}s")


SHARD_ABOVE = 20_000     # anchors; see the note in auto()


def auto(by_rxn, n, food=FOOD, workers=None):
    """Enumerate, sharding only when the network is big enough to pay for it.

    Starting a pool copies the reaction table to every worker, which costs more
    than the search itself on a small network: glucose G5 at three species takes
    2.3 s in one process and 4.2 s across ten. It pays on the networks that take
    minutes. The threshold is a rule of thumb from those two points, not a
    measured optimum, and a caller who knows better should call the enumerator it
    wants directly.
    """
    if workers == 1 or len(anchors(by_rxn, food)) < SHARD_ABOVE:
        return single(by_rxn, n, food=food)   # third positional is chunk, not food
    return enumerate_cores(by_rxn, n, food, workers)

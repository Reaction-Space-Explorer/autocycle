"""The ladder at five species: counts only, never the cores themselves.

Three and a half million cores will not fit in a parent process, and do not need
to. Each shard returns only what the ladder needs -- its count and the sets of
distinct-stoichiometry and mechanism keys it saw -- and the parent takes unions.
"""

from __future__ import annotations

import collections
import os
import sys
import time
from multiprocessing import Pool
from pathlib import Path

from _common import FOOD

from autocycle.cores.anchored import anchors
from autocycle.cores.enumerate_cores import load
from autocycle.cores.motifs import coarse_motif, motif
from autocycle.cores.parallel import _STATE, _init, _one


def _shard(r0):
    hits, _ = _one(r0)
    by = _STATE["by"]
    sig = _STATE["sig"]
    dis, mot, mech = set(), set(), set()
    ty = collections.Counter()
    for sp, rx, t in hits:
        dis.add((frozenset(sp), frozenset(sig[r] for r in rx)))
        mot.add(motif(by, sp, rx))
        mech.add(coarse_motif(by, sp, rx))
        ty[t] += 1
    return len(hits), dis, mot, mech, ty


def _init2(by_rxn, n, food):
    _init(by_rxn, n, food)
    _STATE["sig"] = {r: tuple(sorted(d.items())) for r, d in by_rxn.items()}


if __name__ == "__main__":
    path, n = sys.argv[1], int(sys.argv[2])
    w = int(sys.argv[3]) if len(sys.argv) > 3 else os.cpu_count()
    by = load(path)
    amp = sorted(anchors(by, FOOD))
    t0 = time.time()
    cores, dis, mot, mech = 0, set(), set(), set()
    ty = collections.Counter()
    with Pool(w, initializer=_init2, initargs=(by, n, FOOD)) as pool:
        for c, d, m, k, t in pool.imap_unordered(_shard, amp, chunksize=16):
            cores += c
            dis |= d
            mot |= m
            mech |= k
            ty += t
    print(f"  {Path(path).stem}  n={n}")
    print(f"    cores {cores}  distinct {len(dis)}  motifs {len(mot)}  "
          f"mechanisms {len(mech)}")
    print(f"    cores per mechanism {cores/max(len(mech),1):.0f}:1   "
          f"types {dict(sorted(ty.items()))}   {time.time()-t0:.0f}s")

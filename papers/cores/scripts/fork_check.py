"""Blokhuis's two criteria for the type, run against each other.

The paper gives both a graph-cycle count and a fork count and states them as
equivalent, so disagreement is a finding either way.
"""
import collections

import numpy as np
from _common import FOOD

from autocycle.cores.anchored import enumerate_cores
from autocycle.cores.cores import core_type
from autocycle.cores.enumerate_cores import load
from autocycle.cores.forks import fork_type
from autocycle.cores.paths import RELS

NETS = [("formose G3", "Formose/FormoseRels_3.tsv"), ("glucose G5", "Glucose/GlucoseRels_5.tsv")]

total = agree = 0
for name, rel in NETS:
    by = load(RELS / rel)
    found, _ = enumerate_cores(by, 3, food=FOOD)
    c = collections.Counter()
    for sp, rx, _ in found:
        nu = np.array([[by[r].get(s, 0) for r in rx] for s in sp], dtype=float)
        c[(core_type(nu), fork_type(nu))] += 1
    hit = sum(v for k, v in c.items() if k[0] == k[1])
    total += len(found)
    agree += hit
    print(f"  {name:11s} {len(found):5d} cores, criteria agree on {hit:5d}")
print(f"  total {total} cores, agree {agree}, disagree {total-agree}")

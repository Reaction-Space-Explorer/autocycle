"""Are the type V cores genuinely type V, or a catch-all?

core_type counts simple cycles in the graph on the core species and maps one to I,
two to II/III, three to IV, and anything above three to V. Golnik and co-workers
report that irreversible networks hold cores outside Blokhuis's five types, so a
bin labelled "more than three" is the place such a structure would land silently.

Two checks. Blokhuis give a second, independent criterion: a type I core's fork
ends in two copies of one compound, types II to V in different ones, with one such
fork for II and III, two for IV and three for V. If the fork count agrees with the
graph-cycle count then V is V. And a core of n species cannot hold more than a
bounded number of simple cycles, so a bin far above that bound would be suspect.
"""
import collections

import numpy as np
from _common import FOOD, NETS

from autocycle.cores.cores import core_type
from autocycle.cores.enumerate_cores import load
from autocycle.cores.forks import fork_type
from autocycle.cores.parallel import auto as enumerate_cores
from autocycle.cores.paths import RELS

print(f"  {'network':<20} {'n':>2} {'cores':>7}  types (graph cycles)   forks agree")
for name, (rel, _) in NETS.items():
    by = load(RELS / rel)
    for n in (3, 4):
        found, _ = enumerate_cores(by, n, food=FOOD)
        if not found:
            continue
        kinds = collections.Counter()
        cycles = collections.Counter()
        agree = disagree = 0
        for sp, rx, _t in found:
            sp, rx = sorted(sp), sorted(rx)
            nu = np.array([[by[r].get(s, 0) for r in rx] for s in sp], float)
            t = core_type(nu)
            kinds[t] += 1
            try:
                f = fork_type(nu)
            except Exception:
                f = None
            if f is not None:
                agree += (f == t)
                disagree += (f != t)
        shown = "  ".join(f"{k}:{v}" for k, v in sorted(kinds.items()))
        print(f"  {name:<20} {n:>2} {len(found):>7,}  {shown:<22} "
              f"{agree} agree / {disagree} disagree", flush=True)

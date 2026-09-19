"""How long integer programming takes to *enumerate*, not just to find one.

Finding a single autocatalytic flow is not the expensive part. The claim worth
testing is about enumeration, so this solves repeatedly, each time forbidding the
supports already found with a no-good cut, and reports the time per solution.

    sum_{r in S} v_r <= |S| - 1     for each support S already returned

The cut forbids a superset of the previous support, which is what the flow
formulation returns, so distinct solutions here are distinct supports rather than
distinct cores; that is the comparison the literature's own timings describe.
"""
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import LinearConstraint, milp
from scipy.sparse import csc_matrix, vstack

from autocycle.cores.enumerate_cores import load
from autocycle.cores.flow_ilp import build

if __name__ == "__main__":
    path, target = sys.argv[1], sys.argv[2]
    k = int(sys.argv[3]) if len(sys.argv) > 3 else 20
    limit = float(sys.argv[4]) if len(sys.argv) > 4 else 300.0
    by = load(path)
    nu, species, rxns = build(by)
    n_s, n_r = nu.shape
    lo = np.zeros(n_s)
    lo[species.index(target)] = 1.0
    cuts, cut_lo, cut_hi = [], [], []
    print(f"  {Path(path).stem}: {n_r} reactions, {n_s} non-food species, "
          f"target {target}")
    t_all = time.time()
    for i in range(k):
        A = vstack([nu] + cuts).tocsc() if cuts else nu
        L = np.concatenate([lo] + [np.array(cut_lo)]) if cuts else lo
        H = np.concatenate([np.full(n_s, np.inf), np.array(cut_hi)]) if cuts \
            else np.full(n_s, np.inf)
        t0 = time.time()
        res = milp(c=np.ones(n_r),
                   constraints=LinearConstraint(A, L, H),
                   integrality=np.ones(n_r), bounds=(0, 10),
                   options={"time_limit": limit, "presolve": True})
        dt = time.time() - t0
        if res.x is None:
            print(f"  solution {i+1}: none ({res.message}) after {dt:.1f}s")
            break
        sup = np.flatnonzero(res.x > 0.5)
        row = np.zeros(n_r); row[sup] = 1.0
        cuts.append(csc_matrix(row))
        cut_lo.append(-np.inf); cut_hi.append(len(sup) - 1.0)
        print(f"  solution {i+1}: {len(sup)} reactions, cost {res.fun:.0f}, "
              f"{dt:6.1f}s   (cumulative {time.time()-t_all:6.1f}s)", flush=True)

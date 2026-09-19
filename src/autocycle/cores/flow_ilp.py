"""The integer-programming alternative, timed on the same network.

The comparison the paper needs. Flow search asks an integer program for a
hyperflow that produces a target species while every other non-food species is
at least at steady state; the smallest such flow is an autocatalytic subnetwork
containing that species. This is the formulation the literature reports as
becoming unusable at about a thousand molecules, and it has been quoted here from
someone else's timings on someone else's network. Here it runs on ours.

    minimise  sum v        subject to  (nu v)_x >= 1
                                       (nu v)_s >= 0   for every other non-food s
                                       v integer, 0 <= v <= cap

The two methods are not asked the same question and the difference is stated
rather than hidden: the program searches cores of any size and returns one, the
enumeration is bounded in size and returns all of them.
"""
import sys
import time
from pathlib import Path

import numpy as np
from scipy.optimize import LinearConstraint, milp
from scipy.sparse import csc_matrix

from autocycle.cores.enumerate_cores import load

FOOD = {"O", "C=O", "C(=O)=O", "N"}


def build(by_rxn, food=FOOD):
    rxns = sorted(by_rxn)
    species = sorted({s for d in by_rxn.values() for s in d} - food)
    ix = {s: i for i, s in enumerate(species)}
    rows, cols, vals = [], [], []
    for j, r in enumerate(rxns):
        for s, c in by_rxn[r].items():
            if s in ix:
                rows.append(ix[s]); cols.append(j); vals.append(float(c))
    nu = csc_matrix((vals, (rows, cols)), shape=(len(species), len(rxns)))
    return nu, species, rxns


def solve(nu, species, target, cap=10, limit=600.0):
    n_s, n_r = nu.shape
    lo = np.zeros(n_s)
    lo[species.index(target)] = 1.0
    t0 = time.time()
    res = milp(c=np.ones(n_r),
               constraints=LinearConstraint(nu, lo, np.full(n_s, np.inf)),
               integrality=np.ones(n_r),
               bounds=(0, cap),
               options={"time_limit": limit, "presolve": True})
    return res, time.time() - t0


if __name__ == "__main__":
    path, target = sys.argv[1], sys.argv[2]
    limit = float(sys.argv[3]) if len(sys.argv) > 3 else 600.0
    by = load(path)
    nu, species, rxns = build(by)
    print(f"  {Path(path).stem}: {len(rxns)} reactions, {len(species)} non-food species")
    print(f"  target {target}   time limit {limit:.0f}s", flush=True)
    if target not in species:
        print("  target is not a non-food species of this network")
        raise SystemExit(1)
    res, dt = solve(nu, species, target, limit=limit)
    print(f"  status {res.status}: {res.message}")
    print(f"  elapsed {dt:.1f}s")
    if res.x is not None:
        used = [(rxns[j], int(round(v))) for j, v in enumerate(res.x) if v > 0.5]
        print(f"  flow over {len(used)} reactions, total {sum(v for _, v in used)}")

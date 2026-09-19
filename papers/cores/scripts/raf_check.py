"""Theorem 11 of Golnik et al. (J. Theor. Biol. 2026, 635, 112580), constructively.

The theorem says the stoichiometric matrix of any RAF that is not itself a CAF
contains a square, well-formed, semipositive submatrix over non-food species
whose species and reactions alternate around an elementary circuit. That is the
object this project enumerates, so the theorem is checkable by exhibition.

One thing has to change first. The circuit in T2 runs on x being a *reactant* of
one reaction and a *net-product* of the next, and a catalyst is a reactant. A
catalyst contributes zero to the stoichiometric matrix, so a search that reads
its circuit off the sign pattern of that matrix -- which is what the networks in
this paper allow, having no explicit catalysts -- cannot see these cores at all.
Golnik and co-workers say as much: with catalysts present, "an inspection of the
network stoichiometry fully encoded by R, with its catalysts, and not only of its
net-stoichiometry, encoded by the stoichiometric matrix S alone, is needed."

So the circuit is built here from reactants together with catalysts, while the
matrix tested stays net-stoichiometric. The test itself is unchanged: cores.is_core.
"""
import time
from pathlib import Path

import numpy as np

from autocycle.cores.cores import is_core
from autocycle.io_crs import read_crs
from autocycle.raf import max_caf, max_raf, stoichiometry
from autocycle.spec import SpecError

CRS = Path(__file__).resolve().parents[1] / "data" / "crs"
MAX_N = 6


def circuits(adj, n):
    """Elementary circuits of n species and n distinct reactions."""
    order = {s: i for i, s in enumerate(sorted(adj))}
    out = []
    def walk(start, node, path, rxns):
        for nxt, r in adj.get(node, ()):
            if r in rxns:
                continue
            if nxt == start and len(path) == n:
                out.append((tuple(path), tuple(rxns + [r])))
            elif len(path) < n and nxt in order and order[nxt] > order[start] \
                    and nxt not in path:
                walk(start, nxt, path + [nxt], rxns + [r])
    for s in adj:
        walk(s, s, [s], [])
    return out


def search(rxns, raf, food, by):
    """Smallest core in the RAF, circuit taken over reactants and catalysts."""
    adj = {}
    for r in raf:
        react, prod, cats = rxns[r]
        sources = {s for s in react if s not in food}
        sources |= {c for g in cats for c in g if c not in food}
        targets = [s for s, v in by.get(r, {}).items() if v > 0]
        for s in sources:
            adj.setdefault(s, []).extend((p, r) for p in targets if p != s)
    for n in range(2, MAX_N + 1):
        for spec, rx in circuits(adj, n):
            nu = np.array([[by.get(r, {}).get(s, 0) for r in rx] for s in spec],
                          dtype=float)
            if is_core(nu):
                return n, spec, rx
    return 0, None, None


def stated(path):
    with open(path, encoding="utf-8", errors="replace") as _fh:
        _lines = _fh.readlines()
    for line in _lines:
        if line.startswith("#") and "maxRAF" in line:
            return line.lstrip("# ").strip()
    return ""


print(f"  {'system':32s} {'rxns':>5} {'RAF':>4} {'CAF':>4} {'RAF=CAF':>8} "
      f"{'core n':>7} {'time':>7}")
rows = []
for path in sorted(CRS.glob("*.crs")):
    try:
        system = read_crs(path)
    except SpecError as e:          # inhibition changes what a RAF is; not analysed here
        print(f"  {path.stem:32s} {str(e).split(':')[-1].strip()}")
        continue
    rxns = {r.name: (r.reactants, r.products, r.catalyst_groups)
            for r in system.reactions}
    food = system.food
    raf, caf = max_raf(system), max_caf(system)
    if not raf:
        print(f"  {path.stem:32s} {len(rxns):5d} {0:4d} {len(caf):4d} "
              f"{'-':>8} {'no RAF':>7}")
        continue
    is_a_caf = max_caf(system, raf) == raf
    by = stoichiometry(system, raf)
    t0 = time.time()
    n, spec, rx = search(rxns, raf, food, by)
    print(f"  {path.stem:32s} {len(rxns):5d} {len(raf):4d} {len(caf):4d} "
          f"{str(is_a_caf):>8} {n if n else '-':>7} {time.time()-t0:6.1f}s",
          flush=True)
    rows.append((path.stem, is_a_caf, n, spec, rx))

tested = [r for r in rows if not r[1]]
hit = [r for r in tested if r[2]]
print(f"\n  RAFs that are not CAFs: {len(tested)}")
print(f"    core exhibited:       {len(hit)}")
miss = [r[0] for r in tested if not r[2]]
print(f"    none found to n={MAX_N}:  {len(miss)}  {miss if miss else ''}")
if hit:
    name, _, n, spec, rx = hit[0]
    print(f"\n  smallest example, {name}: {n} species {list(spec)}")
    print(f"    reactions {list(rx)}")

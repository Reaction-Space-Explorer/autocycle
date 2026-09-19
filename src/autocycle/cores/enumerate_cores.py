"""Spike: enumerate autocatalytic cores in a Neo4j-import reaction network."""
import collections
import sys
import time

import numpy as np

from autocycle.cores.cores import core_type, is_core

FOOD = {"O", "C=O"}          # water and formaldehyde: the food, not core members

def load(path):
    nu = collections.defaultdict(int)
    with open(path, encoding="utf-8", errors="replace") as fh:
        rows = [line.rstrip("\n").split("\t") for line in fh]
    for p in rows:
        if len(p) < 4: continue
        try:                          # these tables may carry a header row
            coef = int(p[2])
        except ValueError:
            continue
        nu[(p[0], p[1])] += coef
    by_rxn = collections.defaultdict(dict)
    for (r, s), c in nu.items():
        if c: by_rxn[r][s] = c
    return by_rxn

def cycles(by_rxn, depth, food=FOOD):
    adj = collections.defaultdict(list)
    for r, d in by_rxn.items():
        cons = [s for s, c in d.items() if c < 0 and s not in food]
        prod = [s for s, c in d.items() if c > 0 and s not in food]
        for s in cons:
            for p in prod:
                if s != p: adj[s].append((p, r))
    order = {s: i for i, s in enumerate(sorted(adj))}
    out = []
    def walk(start, node, path, rxns):
        for nxt, r in adj.get(node, ()):
            if r in rxns: continue
            if nxt == start and len(path) >= 2:
                out.append((tuple(path), tuple(rxns + [r])))
            elif len(path) < depth and nxt in order and order[nxt] > order[start] and nxt not in path:
                walk(start, nxt, path + [nxt], rxns + [r])
    for s in adj:
        walk(s, s, [s], [])
    return out, len(adj)

def keep_cores(by_rxn, cand):
    """Square submatrix over the cycle's species and reactions; keep the cores."""
    found = []
    for species, rxns in cand:
        if len(set(rxns)) != len(species): continue
        nu = np.array([[by_rxn[r].get(s, 0) for r in rxns] for s in species], float)
        if (nu > 0).sum(axis=0).max() < 2 and abs(np.linalg.det(nu)) < 1e-9:
            continue
        if is_core(nu):
            found.append((species, rxns, core_type(nu)))
    return found

if __name__ == "__main__":
    path, depth = sys.argv[1], int(sys.argv[2])
    t0 = time.time(); by_rxn = load(path); t1 = time.time()
    cand, n = cycles(by_rxn, depth); t2 = time.time()
    found = keep_cores(by_rxn, cand); t3 = time.time()
    print(f"  network {len(by_rxn)} reactions | cycle graph {n} species")
    print(f"  cycles up to {depth} species: {len(cand)}   ({t2-t1:.1f}s)")
    print(f"  autocatalytic cores: {len(found)}   ({t3-t2:.1f}s LP)   total {t3-t0:.1f}s")
    types = collections.Counter(t for _, _, t in found)
    print(f"  by Blokhuis type: {dict(sorted(types.items()))}")
    for species, rxns, t in found[:4]:
        print(f"    type {t}: {' -> '.join(species)} via {', '.join(rxns)}")

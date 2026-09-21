"""Cores against network size, generation by generation.

Random-network theory says the chance of an autocatalytic set appearing rises
sharply once edge density passes a threshold. Generated chemistry is not random;
this measures what it does instead.
"""
import time

from _common import FOOD

from autocycle.cores.anchored import distinct, enumerate_cores
from autocycle.cores.enumerate_cores import load
from autocycle.cores.motifs import coarse_motif
from autocycle.cores.paths import RELS

R = RELS
NETS = [("Formose", "Formose/FormoseRels_%d.tsv", 6),
        ("Glucose", "Glucose/GlucoseRels_%d.tsv", 5),
        ("PyruvicAcid", "PyruvicAcid/PyruvicAcidRels_%d.tsv", 6)]

print(f"  {'network':12s} {'gen':>3} {'species':>8} {'rxns':>8} {'cores':>7} "
      f"{'distinct':>9} {'mech':>5} {'per 1k rxn':>11} {'time':>6}")
for name, pat, top in NETS:
    for g in range(1, top + 1):
        path = R / (pat % g)
        if not path.exists(): continue
        t0 = time.time()
        by = load(path)
        spec = {s for d in by.values() for s in d} - FOOD
        found, _ = enumerate_cores(by, 3, food=FOOD)
        mech = {coarse_motif(by, sp, rx) for sp, rx, _ in found} if found else set()
        nd = len(distinct(by, found)) if found else 0
        print(f"  {name:12s} {g:3d} {len(spec):8d} {len(by):8d} {len(found):7d} "
              f"{nd:9d} {len(mech):5d} {1000*len(found)/max(len(by),1):11.1f} "
              f"{time.time()-t0:5.0f}s", flush=True)
    print(flush=True)

"""The five study networks at their deepest generation."""
import collections
import time

from _common import FOOD

from autocycle.cores.anchored import distinct, enumerate_cores
from autocycle.cores.enumerate_cores import load
from autocycle.cores.motifs import coarse_motif
from autocycle.cores.paths import RELS
from autocycle.cores.triage import verdict

R = RELS
NETS = [("PyruvicAcid", "PyruvicAcid/PyruvicAcidRels_6.tsv"),
        ("Formose", "Formose/FormoseRels_6.tsv"),
        ("FormoseAmm", "FormoseAmm/FormoseAmmRels_4.tsv"),
        ("Glucose", "Glucose/GlucoseRels_5.tsv"),
        ("GlucoseAmm", "GlucoseAmm/GlucoseAmmRels_4.tsv")]

print(f"  {'network':13s} {'rxns':>8} {'cores':>7} {'distinct':>9} {'mech':>5} "
      f"{'serious':>8} {'cond':>6} {'artefact':>9} {'time':>6}")
for name, rel in NETS:
    t0 = time.time()
    path = R / rel
    by = load(path)
    rules = {}
    with open(path, encoding="utf-8", errors="replace") as _fh:
        _lines = _fh.readlines()
    for line in _lines:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4 and p[0] != "Index":
            rules.setdefault(p[0], p[3])
    found, _ = enumerate_cores(by, 3, food=FOOD)
    if not found:
        print(f"  {name:13s} {len(by):8d} {0:7d}", flush=True); continue
    mech = {coarse_motif(by, sp, rx) for sp, rx, _ in found}
    v = collections.Counter(verdict(tuple(sorted(rules.get(r, "") for r in rx)))[0]
                            for sp, rx, _ in found)
    print(f"  {name:13s} {len(by):8d} {len(found):7d} {len(distinct(by, found)):9d} "
          f"{len(mech):5d} {v['serious']:8d} {v['conditional']:6d} {v['artefact']:9d} "
          f"{time.time()-t0:5.0f}s", flush=True)

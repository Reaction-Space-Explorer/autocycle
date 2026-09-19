"""Which cores survive when a transformation rule is taken out of the network.

The mechanism classes are defined over the rule names the generator recorded. If
those classes are an artefact of how finely the rule library was split, removing
one rule should scatter the counts. If they track chemistry, removing a rule
should kill the mechanisms that need it and leave the rest alone.
"""
import collections
import time

from autocycle.cores.anchored import enumerate_cores
from autocycle.cores.enumerate_cores import load
from autocycle.cores.motifs import coarse_motif
from autocycle.cores.paths import RELS
from autocycle.cores.triage import verdict

FOOD = {"O", "C=O", "C(=O)=O", "N"}
NET = RELS / "Formose/FormoseRels_6.tsv"
DROP = ["Keto-enol migration twice", "Cannizarro 2, HCHO (reduction)",
        "Michael Addition 0,2, ", "Retro Aldol", "Aldol Condensation"]


def run(by, rules, label):
    found, _ = enumerate_cores(by, 3, food=FOOD)
    mech = {coarse_motif(by, sp, rx) for sp, rx, _ in found}
    v = collections.Counter(verdict(tuple(sorted(rules[r] for r in rx)))[0]
                            for sp, rx, _ in found)
    print(f"  {label:34s} {len(by):7d} {len(found):6d} {len(mech):5d} "
          f"{v['serious']:8d} {v['conditional']:6d} {v['artefact']:9d}", flush=True)
    return {coarse_motif(by, sp, rx) for sp, rx, _ in found}


rules = {}
with open(NET, encoding="utf-8", errors="replace") as _fh:
    _lines = _fh.readlines()
for line in _lines:
    p = line.rstrip("\n").split("\t")
    if len(p) >= 4:
        rules.setdefault(p[0], p[3])

full = load(NET)
print(f"  {'removed':34s} {'rxns':>7} {'cores':>6} {'mech':>5} "
      f"{'serious':>8} {'cond':>6} {'artefact':>9}")
base = run(full, rules, "nothing (reference)")

for name in DROP:
    t0 = time.time()
    kept = {r: d for r, d in full.items() if rules.get(r) != name}
    if len(kept) == len(full):
        print(f"  {name:34s} rule not present")
        continue
    mech = run(kept, rules, name)
    died, new = base - mech, mech - base
    print(f"  {'':34s} mechanisms lost {len(died):3d}, unchanged {len(base & mech):3d}, "
          f"new {len(new):3d}   ({time.time()-t0:.0f}s)", flush=True)

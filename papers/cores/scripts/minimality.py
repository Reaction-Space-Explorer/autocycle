"""Minimal four-species cores: the count after the passengers are removed.

A core is minimal by definition, so a four-species core that contains a
three-species core is not a core but a core with a passenger. The test is
containment of both the species and the reactions of an enumerated three-species
core -- enumerated, not any square submatrix: almost any 2x2 slice of a core is
invertible and admits a positive flux while corresponding to no cycle at all,
which would report every core as non-minimal.
"""
import collections

from _common import FOOD

from autocycle.cores.enumerate_cores import load
from autocycle.cores.motifs import coarse_motif
from autocycle.cores.parallel import auto as enumerate_cores
from autocycle.cores.paths import RELS
from autocycle.cores.triage import verdict

NETS = [("glucose G5", "Glucose/GlucoseRels_5.tsv"),
        ("formose G6", "Formose/FormoseRels_6.tsv"),
        ("pyruvic acid G6", "PyruvicAcid/PyruvicAcidRels_6.tsv")]

print(f"  {'network':18s} {'raw':>7} {'minimal':>8} {'dropped':>8} "
      f"{'mech':>5} {'serious':>8} {'cond':>6} {'artefact':>9} "
      f"{'I':>6} {'II/III':>7} {'IV':>5} {'V':>4}")
for name, rel in NETS:
    path = RELS / rel
    by = load(path)
    rules = {}
    with open(path, encoding="utf-8", errors="replace") as fh:
        for line in fh:
            cols = line.rstrip("\n").split("\t")
            if len(cols) >= 4:
                rules.setdefault(cols[0], cols[3])
    small, _ = enumerate_cores(by, 3, food=FOOD)
    big, _ = enumerate_cores(by, 4, food=FOOD)
    index = collections.defaultdict(list)
    for sp, rx, _ in small:
        index[frozenset(rx)].append(frozenset(sp))
    keys = list(index)
    minimal = []
    for sp, rx, t in big:
        R, S = frozenset(rx), frozenset(sp)
        if any(k <= R and any(v <= S for v in index[k]) for k in keys if k <= R):
            continue
        minimal.append((sp, rx, t))
    mech = {coarse_motif(by, sp, rx) for sp, rx, _ in minimal}
    v = collections.Counter(verdict(tuple(sorted(rules[r] for r in rx)))[0]
                            for sp, rx, _ in minimal)
    ty = collections.Counter(t for _, _, t in minimal)
    print(f"  {name:18s} {len(big):7d} {len(minimal):8d} "
          f"{len(big)-len(minimal):8d} {len(mech):5d} {v['serious']:8d} "
          f"{v['conditional']:6d} {v['artefact']:9d} {ty['I']:6d} "
          f"{ty['II/III']:7d} {ty['IV']:5d} {ty['V']:4d}", flush=True)

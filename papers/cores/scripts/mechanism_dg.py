"""Does the mechanism level predict free energy, which it was never built from?

A mechanism is a multiset of step descriptions, each the formula change a step
makes to the ring plus what it takes in and lets out. No energy enters that
definition and no structure beyond the formula, so if cores sharing a mechanism
also share a cycle free energy, the coarsening has caught chemistry it never saw.

It is not circular: a free energy depends on which isomer a formula denotes, and
the networks are full of isomers, so the same formula change carries different
energies on different substrates. It is not a strong test either, since steps that
make the same formula change tend to break similar bonds. What it rules out is the
coarsening being an arbitrary relabelling.

The cycle free energy is the sum over the ring, one turn, and only cores whose
every step carries a usable estimate are counted.
"""
import collections
import csv
import random
import statistics

from _common import FOOD, NETS

from autocycle.cores.enumerate_cores import load
from autocycle.cores.motifs import coarse_motif
from autocycle.cores.parallel import auto as enumerate_cores
from autocycle.cores.paths import ENERGIES, RELS

SENTINEL = 1e4
TRIALS = 2000


def pooled_within(groups):
    """Standard deviation within mechanism, pooled over mechanisms of size >= 2."""
    num = den = 0.0
    for vals in groups:
        if len(vals) < 2:
            continue
        m = statistics.fmean(vals)
        num += sum((v - m) ** 2 for v in vals)
        den += len(vals) - 1
    return (num / den) ** 0.5 if den else float("nan")


print(f"  {'network':<20} {'cores':>7} {'scored':>7} {'mech':>5} "
      f"{'within':>8} {'overall':>8} {'ratio':>6} {'null ratio':>11} {'p':>7}")
for name, (rel, stem) in NETS.items():
    by = load(RELS / rel)
    with open(ENERGIES / f"{stem}_energies_pH7.4.csv") as fh:
        en = {r["Index"]: r for r in csv.DictReader(fh)}

    def dg(rx, en=en):
        total = 0.0
        for r in rx:
            e = en.get(r)
            if e is None or not e["dG_prime_kJ_mol"].strip():
                return None
            if float(e["sigma_kJ_mol"]) >= SENTINEL:
                return None
            total += float(e["dG_prime_kJ_mol"])
        return total

    found, _ = enumerate_cores(by, 3, food=FOOD)
    pairs = []
    for sp, rx, _ in found:
        g = dg(rx)
        if g is not None:
            pairs.append((coarse_motif(by, sp, rx), g))
    if len(pairs) < 20:
        print(f"  {name:<20} {len(found):>7,} {len(pairs):>7,}  too few scored cores")
        continue

    byname = collections.defaultdict(list)
    for m, g in pairs:
        byname[m].append(g)
    within = pooled_within(byname.values())
    overall = statistics.stdev([g for _, g in pairs])

    labels = [m for m, _ in pairs]
    vals = [g for _, g in pairs]
    worse = 0
    for _ in range(TRIALS):
        random.shuffle(labels)
        sh = collections.defaultdict(list)
        for m, g in zip(labels, vals, strict=True):
            sh[m].append(g)
        if pooled_within(sh.values()) <= within:
            worse += 1
    print(f"  {name:<20} {len(found):>7,} {len(pairs):>7,} {len(byname):>5} "
          f"{within:>8.1f} {overall:>8.1f} {within / overall:>6.2f} "
          f"{'shuffled':>11} {(worse + 1) / (TRIALS + 1):>7.4f}", flush=True)

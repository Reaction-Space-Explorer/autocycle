"""How much of the reaction table is the same transformation listed twice.

The generator writes one row per rule that derives a hyperedge, so a reaction two
rules can produce appears twice. Counting distinct stoichiometries, over reactions
whose net stoichiometry is not zero, gives the generator's own edge count, and the
gap to the row count is aliasing the generator introduced rather than a judgement
made here. Also reports how many reactions the anchor test admits under each food
set, since that count is quoted alongside.
"""
import collections

from _common import BASE, FOOD

from autocycle.cores.anchored import anchors
from autocycle.cores.enumerate_cores import load
from autocycle.cores.paths import RELS

NETS = [("formose G6", "Formose/FormoseRels_6.tsv"),
        ("glucose G5", "Glucose/GlucoseRels_5.tsv")]

print(f"  {'network':<12} {'rows':>9} {'distinct edges':>15} {'inflation':>10} "
      f"{'amplifying, water+CH2O':>23} {'amplifying, +CO2,NH3':>21}")
for name, rel in NETS:
    by = load(RELS / rel)
    sigs = collections.Counter()
    for d in by.values():
        net = tuple(sorted((s, c) for s, c in d.items() if c))
        if net:                       # a reaction with zero net stoichiometry is not an edge
            sigs[net] += 1
    rows, edges = len(by), len(sigs)
    print(f"  {name:<12} {rows:>9,} {edges:>15,} {100 * rows / edges - 100:>9.0f}% "
          f"{len(anchors(by, BASE)):>23,} {len(anchors(by, FOOD)):>21,}", flush=True)

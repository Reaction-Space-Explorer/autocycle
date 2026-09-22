"""What survives if the generator had been guided by spontaneity.

The deposition carries every generation separately and the generations are
cumulative, so the expansion can be replayed under a filter without regenerating
anything: at each generation keep the reactions whose reactants are all food or
species that survived the generation before, and which the free energy calls
downhill; their products are the feedstock for the next.

Three classes, not two. A reaction whose free energy carries the sentinel variance
cannot be called either way, and those are 18% of formose and 46% of formose with
ammonia, so the replay is run twice: once admitting the undecidable reactions and
once refusing them. What survives both is what the filter actually supports.

The question this is for: methanol is made by a Cannizzaro, and the artefact class
needs methanol. If the pruned expansion never reaches it, that class depends on a
reaction a thermodynamically guided search would not have taken.
"""
import csv

from _common import BASE, FOOD

from autocycle.cores.enumerate_cores import load
from autocycle.cores.paths import ENERGIES, RELS

SENTINEL = 1e4
METHANOL = "CO"
NETS = [("formose", "Formose/FormoseRels_%d.tsv", "Formose_G6", 6),
        ("glucose", "Glucose/GlucoseRels_%d.tsv", "Glucose_G5", 5)]


def energies(stem):
    with open(ENERGIES / f"{stem}_energies_pH7.4.csv") as fh:
        return {r["Index"]: r for r in csv.DictReader(fh)}


def verdict(e):
    """downhill, uphill, or undecidable."""
    if e is None or not e["dG_prime_kJ_mol"].strip():
        return "undecidable"
    if float(e["sigma_kJ_mol"]) >= SENTINEL:
        return "undecidable"
    return "downhill" if float(e["dG_prime_kJ_mol"]) < 0 else "uphill"


def replay(pattern, stem, top, admit_undecidable):
    en = energies(stem)
    have = set(FOOD) | set(BASE)
    kept = 0
    for g in range(1, top + 1):
        by = load(RELS / (pattern % g))
        grew = True
        while grew:                       # a generation can unlock itself
            grew = False
            for r, d in by.items():
                need = [s for s, c in d.items() if c < 0]
                if not all(s in have for s in need):
                    continue
                v = verdict(en.get(r))
                if v == "uphill" or (v == "undecidable" and not admit_undecidable):
                    continue
                made = [s for s, c in d.items() if c > 0 and s not in have]
                if made:
                    have.update(made)
                    grew = True
        kept = len(have)
    return have, kept


print(f"  {'network':<10} {'undecidable':<12} {'species reached':>15} {'of full':>9} "
      f"{'methanol':>9}")
for name, pattern, stem, top in NETS:
    full = load(RELS / (pattern % top))
    allsp = {s for d in full.values() for s in d}
    for admit, label in ((True, "admitted"), (False, "refused")):
        have, _ = replay(pattern, stem, top, admit)
        reached = have & allsp
        print(f"  {name:<10} {label:<12} {len(reached):>15,} "
              f"{100 * len(reached) / len(allsp):>8.1f}% "
              f"{'yes' if METHANOL in have else 'NO':>9}", flush=True)

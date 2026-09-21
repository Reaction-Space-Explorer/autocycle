"""Anchored enumeration against a brute-force walk over every cycle.

Proposition 3 says the anchors partition the search; it does not say the walker
finds every cycle through an anchor, and that half is checked rather than proved.
The check only means something on a network where the anchor test and the
core-column condition come apart. They agree on every formose reaction, so
formose alone cannot catch a wrong anchor test. Glucose and the ammonia-fed
networks can, and are included for that reason.

Brute force walks every cycle, so it is run at the generations where that is
affordable rather than at the depths the paper reports.
"""
from autocycle.cores import anchored
from autocycle.cores.enumerate_cores import cycles, keep_cores, load
from autocycle.cores.paths import RELS

FOOD = {"O", "C=O", "C(=O)=O", "N"}
NETS = [("formose G3", "Formose/FormoseRels_3.tsv"),
        ("glucose G2", "Glucose/GlucoseRels_2.tsv"),
        ("glucose G3", "Glucose/GlucoseRels_3.tsv"),
        ("glucose+ammonia G2", "GlucoseAmm/GlucoseAmmRels_2.tsv"),
        ("glucose+ammonia G3", "GlucoseAmm/GlucoseAmmRels_3.tsv"),
        ("formose+ammonia G3", "FormoseAmm/FormoseAmmRels_3.tsv"),
        ("pyruvic acid G3", "PyruvicAcid/PyruvicAcidRels_3.tsv")]


def keys(found):
    return {(frozenset(sp), frozenset(rx)) for sp, rx, _ in found}


def apart(by_rxn):
    """Reactions where the two conditions can differ: P >= 2 but P - N <= 0."""
    k = 0
    for d in by_rxn.values():
        p = sum(c for s, c in d.items() if s not in FOOD and c > 0)
        n = -sum(c for s, c in d.items() if s not in FOOD and c < 0)
        k += p >= 2 and n >= p
    return k


print(f"  {'network':<20} {'reactions':>9} {'apart':>6}  {'n':>2} "
      f"{'brute':>6} {'anchored':>8}  agree")
for name, rel in NETS:
    path = RELS / rel
    if not path.exists():
        print(f"  {name:<20} absent")
        continue
    by = load(path)
    for n in (2, 3, 4):
        cand, _ = cycles(by, n, food=FOOD)
        brute = keys(keep_cores(by, [c for c in cand if len(c[0]) == n]))
        anch = keys(anchored.enumerate_cores(by, n, food=FOOD)[0])
        print(f"  {name:<20} {len(by):>9,} {apart(by):>6}  {n:>2} "
              f"{len(brute):>6} {len(anch):>8}  {brute == anch}")

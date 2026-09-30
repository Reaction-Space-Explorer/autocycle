"""Which enumerated core pairs have the structure a bistable pair needs?

Ivan and co-workers build bistability from two autocatalytic cycles drawing on one
resource, an annihilation that consumes a species of each and returns neither, and
a weak interconversion between them. Sakref and Rivoire, and Gagrani and co-workers,
likewise begin from a pair of cycles taken as given. The census here is of cores
rather than of pairs, so the question none of them was placed to ask is how many
such pairs a generated network actually offers.

Their R + X -> 2X, R + Y -> 2Y, X + Y -> P, X <-> Y is four conditions, and all
four are stoichiometric. The cores must be disjoint or they are not alternative
states. They must draw on a common input or they do not compete for it. Some
reaction must consume a species of each and return neither. And some reaction must
carry a species of one into the other, which is the mutation term, counted here in
its strict form of a reaction with one non-food reactant, since a reaction needing
a second substrate is a different process however it moves the species.

Competition is not free. Most cores here consume no food at all: they run on
intermediates that food produces further upstream, so the resource two of them
share is a species of the network rather than a member of the food set, and it has
to be looked for.

What this cannot say is whether any pair found is bistable. That needs rate
constants. A core is necessary and not sufficient, and a pair of them is no more.
"""
import collections
import itertools

from _common import FOOD, NETS

from autocycle.cores.enumerate_cores import load
from autocycle.cores.parallel import auto as enumerate_cores
from autocycle.cores.paths import RELS

N = 3


def sides(d):
    cons = [s for s, c in d.items() if c < 0 and s not in FOOD]
    prod = [s for s, c in d.items() if c > 0 and s not in FOOD]
    return cons, prod


def compositions(by, found):
    """Each distinct core species set, with what it draws on.

    Two cores over the same species by different reactions are one composition, and
    it is the composition a compartment is dominated by. What it draws on is
    everything its reactions consume and it is not made of, food included, which is
    the R the two cycles compete for.
    """
    draws = collections.defaultdict(set)
    for sp, rx, _ in found:
        k = frozenset(sp)
        for r in rx:
            draws[k] |= {s for s, c in by[r].items() if c < 0 and s not in k}
    return sorted(draws, key=sorted), draws


def index(by):
    """Reactions by the species pair they consume, and strict one-step conversion."""
    cross = collections.defaultdict(set)
    step = collections.defaultdict(set)
    for d in by.values():
        cons, prod = sides(d)
        prod = frozenset(prod)
        for a, b in itertools.combinations(sorted(set(cons)), 2):
            cross[(a, b)].add(prod)
        if len(set(cons)) == 1:
            step[cons[0]] |= prod
    return cross, step


def screen(by, sets, draws):
    """The funnel, each stage a condition of their model added to the one before.

    Crossing the cores that hold each reactant of each reaction is the obvious way
    and the wrong one: a species sits in many cores, so that product runs to
    millions where the answer is tens of thousands. Indexing once by species pair
    and asking each pair of compositions its own question turns it around.
    """
    cross, step = index(by)
    n = disjoint = competing = exclusive = motif = 0
    for i, j in itertools.combinations(range(len(sets)), 2):
        si, sj = sets[i], sets[j]
        n += 1
        if si & sj:
            continue
        disjoint += 1
        if not draws[si] & draws[sj]:
            continue
        competing += 1
        union = si | sj
        if not any(p.isdisjoint(union)
                   for a in si for b in sj
                   for p in cross.get((a, b) if a < b else (b, a), ())):
            continue
        exclusive += 1
        if any(step[a] & sj for a in si) and any(step[b] & si for b in sj):
            motif += 1
    return n, disjoint, competing, exclusive, motif


def main():
    print(f"  {'network':<20} {'cores':>7} {'states':>7} {'pairs':>11} "
          f"{'disjoint':>11} {'competing':>11} {'exclusive':>11} {'motif':>9}")
    for name, (rel, _) in NETS.items():
        by = load(RELS / rel)
        found, _ = enumerate_cores(by, N, food=FOOD)
        sets, draws = compositions(by, found)
        n, dis, comp, ex, mot = screen(by, sets, draws)
        print(f"  {name:<20} {len(found):>7,} {len(sets):>7,} {n:>11,} "
              f"{dis:>11,} {comp:>11,} {ex:>11,} {mot:>9,}", flush=True)


if __name__ == "__main__":
    main()

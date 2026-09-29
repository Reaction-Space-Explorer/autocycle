"""Check the enumeration against autogatito, the authors' own implementation.

Golnik and co-workers publish autogatito, which enumerates cores as MR-chordless
circuits and their superpositions. It is not a dependency of this package: it
wants SBML rather than the relation tables, and it is run separately. This script
takes what it wrote and answers the two questions that matter here. How many of
its cores at each size are single circuits, which is what the search in this paper
covers, and do those counts match.

Run autogatito first, with the dump its analysis does not write by default:

    python partitionNetwork.py -i net.xml -t 4 -o name
    AUTOGATITO_DUMP=cores.pkl python partitionAnalysis.py -x net.xml \
        -i PickleFiles/name/partitionTree0.pkl -b 1000 -n True -s name \
        -t 4 -e 1000 -c -p -o cycleData

then point this at the partition tree and the dump.
"""
import collections
import pickle
import sys


def classify(net, core):
    """True when the core is one closed ring rather than several joined.

    autogatito stores a core as the edges pairing each species with the reaction
    consuming it. The core is a single circuit when those species can be ordered
    so that each reaction produces the next one around the ring.
    """
    cons = {}
    for u, v in core:
        s, r = (u, v) if net.nodes[u]["Type"] == "Species" else (v, u)
        cons[s] = r
    species = set(cons)
    succ = {s: {t for t in net.successors(cons[s]) if t in species and t != s}
            for s in species}
    start = next(iter(species))

    def walk(cur, seen):
        if len(seen) == len(species):
            return start in succ[cur]
        return any(walk(t, seen | {t}) for t in succ[cur] - seen)

    return walk(start, {start})


def main(tree, dump):
    with open(tree, "rb") as fh:
        net = pickle.load(fh)[5]
    with open(dump, "rb") as fh:
        cores = pickle.load(fh)
    tally = collections.defaultdict(lambda: [0, 0])
    for c in cores:
        tally[len(c)][0 if classify(net, c) else 1] += 1
    print(f"  {'n':>3} {'cores':>8} {'circuits':>9} {'superpositions':>15}")
    for n in sorted(tally):
        circ, sup = tally[n]
        print(f"  {n:>3} {circ + sup:>8,} {circ:>9,} {sup:>15,}")


if __name__ == "__main__":
    main(*sys.argv[1:3])

"""Do the five chemistries draw on one vocabulary of mechanisms, or five?

A mechanism is a multiset of step descriptions, each the formula change a step
makes to the ring plus what it takes in and lets out. That description carries no
rule name and no substrate, so the same label in two networks means the same
transformation reached by whatever chemistry each network had available. The
labels are therefore comparable across networks, which the counts in Table 1 are
not.

Reports how many mechanisms each network has, how many it shares with each other
network, and how many are common to all five.
"""
import collections

from _common import FOOD, NETS

from autocycle.cores.anchored import enumerate_cores
from autocycle.cores.enumerate_cores import load
from autocycle.cores.motifs import coarse_motif
from autocycle.cores.paths import RELS

# one network takes an hour, and the machine has gone down mid-run before, so each
# is cached as it finishes and a rerun picks up where the last one stopped
CACHE = Path(__file__).resolve().parents[1] / "results" / "shared_mechanisms_cache.json"
cache = json.loads(CACHE.read_text()) if CACHE.exists() else {}

mech, cores = {}, {}
for name, (rel, _) in NETS.items():
    if name in cache:
        counts = collections.Counter({tuple(k.split("\t")): v for k, v in cache[name]["m"].items()})
        mech[name], cores[name] = counts, cache[name]["n"]
        print(f"  {name:<20} {cores[name]:>6,} cores  {len(counts):>4} mechanisms  (cached)", flush=True)
        continue
    by = load(RELS / rel)
    found, _ = enumerate_cores(by, 3, food=FOOD)
    counts = collections.Counter(coarse_motif(by, sp, rx) for sp, rx, _ in found)
    mech[name], cores[name] = counts, len(found)
    cache[name] = {"n": len(found), "m": {"\t".join(k): v for k, v in counts.items()}}
    CACHE.write_text(json.dumps(cache))
    print(f"  {name:<20} {len(found):>6,} cores  {len(counts):>4} mechanisms", flush=True)

names = list(mech)
print(f"\n  shared mechanisms, pairwise ({'|'.join(n.split()[0] for n in names)}):")
print(f"  {'':<20} " + " ".join(f"{n.split()[0][:8]:>8}" for n in names))
for a in names:
    row = "".join(f"{len(set(mech[a]) & set(mech[b])):>9,}" for b in names)
    print(f"  {a:<20}{row}")

common = set.intersection(*(set(m) for m in mech.values()))
union = set().union(*(set(m) for m in mech.values()))
print(f"\n  {len(union)} distinct mechanisms across all five networks")
print(f"  {len(common)} common to all five "
      f"({100 * len(common) / len(union):.1f}% of the union)")
for k in range(5, 0, -1):
    n = sum(1 for m in union if sum(m in mech[x] for x in names) == k)
    share = sum(c for x in names for m, c in mech[x].items()
                if sum(m in mech[y] for y in names) == k)
    print(f"    in {k} network{'s' if k > 1 else ' '}: {n:>5} mechanisms, "
          f"{share:>7,} cores ({100 * share / sum(cores.values()):.1f}%)")

print("\n  the mechanisms common to all five, by total cores carried:")
tot = {m: sum(mech[x].get(m, 0) for x in names) for m in common}
for m, c in sorted(tot.items(), key=lambda kv: -kv[1])[:10]:
    print(f"    {c:>6,}  {'  '.join(m)}")

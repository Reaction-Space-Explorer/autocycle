"""What survives if every unresolved core went the worst way, and the best.

A core whose cycle free energy has no usable estimate is not a core that failed
the test; it is one the test cannot reach, since an unbounded variance can never
satisfy dG + 1.96 sigma < 0. The true spontaneous count therefore lies between
the resolved count and the resolved count plus every unresolved core. Where that
interval is narrow the conclusion holds whatever the blind cores are; where it is
wide the honest report is the interval.

Arithmetic over the counts in results/paper_numbers.txt. Regenerate that first if
the enumeration or the food set changes.
"""
import re
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
TEXT = (ROOT / "results" / "paper_numbers.txt").read_text()

nets = {}
for block in TEXT.split("=== ")[1:]:
    name = block.split(" ===")[0]
    counts = {}
    for verdict, state, n in re.findall(
            r"^\s+(serious|conditional|artefact)\s+(spontaneous|not|no estimate)\s+(\d+)$",
            block, re.M):
        counts[(verdict, state)] = int(n)
    nets[name] = counts

print(f"  {'network':<20} {'class':<12} {'cores':>6} {'blind':>6} "
      f"{'low':>6} {'high':>6}   spontaneous fraction")
for name, c in nets.items():
    for v in ("serious", "conditional", "artefact"):
        spont, nots = c.get((v, "spontaneous"), 0), c.get((v, "not"), 0)
        blind = c.get((v, "no estimate"), 0)
        tot = spont + nots + blind
        if not tot:
            continue
        lo, hi = spont / tot, (spont + blind) / tot
        flag = "  vacuous" if hi - lo > 0.5 else ""
        print(f"  {name:<20} {v:<12} {tot:>6,} {blind:>6,} "
              f"{spont:>6,} {spont + blind:>6,}   {lo:5.1%} to {hi:5.1%}{flag}")
    print()

print("  blind rate by class, the artefact column is the one to look at:")
for name, c in nets.items():
    row = []
    for v in ("serious", "conditional", "artefact"):
        tot = sum(c.get((v, s), 0) for s in ("spontaneous", "not", "no estimate"))
        if tot:
            row.append(f"{v} {c.get((v, 'no estimate'), 0) / tot:.1%}")
    print(f"    {name:<20} {'  '.join(row)}")

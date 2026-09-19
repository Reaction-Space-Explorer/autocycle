"""The Calvin and reverse Krebs cycles, from KEGG, against Blokhuis's claim.

Blokhuis, Lacoste and Nghe state that "autocatalytic cores of types I and III can
be found in the Calvin cycle and reverse Krebs cycle". The modules are KEGG
M00165 and M00173; equations come from the reaction records, with coefficients,
so nothing here depends on reading a figure.

Food is the cofactor pool the cell buffers rather than the cycle making it:
water, carbon dioxide, protons, the adenine and nicotinamide nucleotides,
phosphate and pyrophosphate, coenzyme A, ferredoxin, quinones. That choice is a
result, not a parameter, so it is varied at the end.
"""
import collections
import re
from pathlib import Path

import numpy as np

from autocycle.cores.cores import core_type, is_core
from autocycle.cores.enumerate_cores import cycles
from autocycle.cores.forks import fork_type, forks

KEGG = Path(__file__).resolve().parents[1] / "data" / "kegg"

FOOD = {
    "C00001",  # H2O
    "C00011",  # CO2
    "C00080",  # H+
    "C00002", "C00008", "C00020",           # ATP, ADP, AMP
    "C00003", "C00004", "C00005", "C00006",  # NAD+, NADH, NADPH, NADP+
    "C00009", "C00013",                      # phosphate, pyrophosphate
    "C00010",                                # CoA
    "C00138", "C00139",                      # reduced / oxidised ferredoxin
    "C00399", "C00390",                      # ubiquinone / ubiquinol
    "C00007", "C00282",                      # O2, H2
}
TERM = re.compile(r"^(?:(\d+)\s+)?(C\d{5})$")


def equations(path=KEGG / "equations.txt"):
    out, entry = {}, None
    with open(path, encoding="utf-8") as _fh:
        _lines = _fh.readlines()
    for line in _lines:
        if line.startswith("ENTRY"):
            entry = line.split()[1]
        elif line.startswith("EQUATION") and entry:
            out[entry] = line.split(None, 1)[1].strip()
    return out


def parse_equation(eq):
    """Net stoichiometry, and whether KEGG writes the reaction as reversible."""
    left, right = re.split(r"<=>|=>|<=", eq, maxsplit=1)
    col = collections.Counter()
    for side, sign in ((left, -1), (right, +1)):
        for term in side.split(" + "):
            m = TERM.match(term.strip())
            if m:
                col[m.group(2)] += sign * int(m.group(1) or 1)
    return {s: c for s, c in col.items() if c}, "<=>" in eq


def network(module, eqs, food, reversible=True):
    """Both directions of a reversible reaction, since the direction a cycle needs
    is not the direction KEGG happens to write. A reaction and its reverse give
    columns that are negatives of each other, so no pair of them is invertible and
    the core test rejects that degenerate case on its own."""
    with open(KEGG / f"{module}.rxns") as fh:
        ids = {ln.strip() for ln in fh}
    by = {}
    for r in ids:
        if r not in eqs:
            continue
        net, rev = parse_equation(eqs[r])
        col = {s: c for s, c in net.items() if s not in food}
        if col:
            by[r] = col
            if rev and reversible:
                by[r + "_rev"] = {s: -c for s, c in col.items()}
    return by


def search(by, food, upto=12):
    """Cores up to `upto` species. These cycles are long: the Calvin core runs
    from ribulose bisphosphate through glycerate and back, which is more species
    than any core in the generated networks."""
    hits, seen = [], set()
    cand, _ = cycles(by, upto, food=food)
    for spec, rx in cand:
        if len(set(rx)) != len(spec):
            continue
        key = (frozenset(spec), frozenset(rx))
        if key in seen:
            continue
        seen.add(key)
        nu = np.array([[by[r].get(s, 0) for r in rx] for s in spec], dtype=float)
        if is_core(nu):
            hits.append((spec, rx, core_type(nu), fork_type(nu), forks(nu)))
    return hits


eqs = equations()
NAMES = {"M00165": "Calvin cycle", "M00173": "reverse Krebs cycle"}
for module, label in NAMES.items():
    by = network(module, eqs, FOOD)
    hits = search(by, FOOD)
    types = collections.Counter(h[2] for h in hits)
    print(f"  {label:22s} {len(by):3d} reactions, "
          f"{len({s for d in by.values() for s in d}):3d} non-food species")
    print(f"    cores {len(hits):3d}   types " +
          "  ".join(f"{k}:{v}" for k, v in sorted(types.items())) or "    none")
    bad = [h for h in hits if h[2] != h[3]]
    print(f"    cycle and fork criteria agree on {len(hits)-len(bad)} of {len(hits)}")
    for t in ("I", "II/III", "IV", "V"):
        ex = [h for h in hits if h[2] == t]
        if ex:
            spec, rx, _, ft, (same, dist) = min(ex, key=lambda h: len(h[0]))
            print(f"      type {t:6s} {len(ex):3d} cores, smallest {len(spec)} species; "
                  f"forks: {same} onto one species, {dist} onto two -> {ft}")
    print(flush=True)

print("  food-set sensitivity, Calvin cycle:")
for drop, note in [(set(), "as above"), ({"C00009"}, "phosphate not food"),
                   ({"C00002", "C00008"}, "ATP/ADP not food"),
                   ({"C00011"}, "CO2 not food")]:
    f = FOOD - drop
    by = network("M00165", eqs, f)
    hits = search(by, f)
    t = collections.Counter(x[2] for x in hits)
    print(f"    {note:22s} cores {len(hits):3d}   " +
          "  ".join(f"{k}:{v}" for k, v in sorted(t.items())))

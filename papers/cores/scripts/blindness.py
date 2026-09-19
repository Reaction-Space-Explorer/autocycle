"""Why the thermodynamic filter cannot see a core.

Component contribution reports three different failures as one "unestimable":
an unbounded variance (the sentinel, the value beside it arbitrary), a null
0 +/- 0 returned where reagents and products share a decomposition, and no
value at all. A core inherits the worst reason among its reactions.
"""
import collections
import csv

from autocycle.cores.anchored import enumerate_cores
from autocycle.cores.enumerate_cores import load
from autocycle.cores.paths import ENERGIES, RELS
from autocycle.cores.triage import verdict


def _rows(path):
    """csv.DictReader over a file that is closed when the read finishes."""
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


SENTINEL, NULL_DG = 1e4, 1e-3
FOOD = {"O", "C=O", "C(=O)=O", "N"}
NETS = [("glucose G5", "Glucose/GlucoseRels_5.tsv", "Glucose_G5"),
        ("glucose+ammonia G4", "GlucoseAmm/GlucoseAmmRels_4.tsv", "GlucoseAmm_G4"),
        ("formose G6", "Formose/FormoseRels_6.tsv", "Formose_G6"),
        ("formose+ammonia G4", "FormoseAmm/FormoseAmmRels_4.tsv", "FormoseAmm_G4"),
        ("pyruvic acid G6", "PyruvicAcid/PyruvicAcidRels_6.tsv", "PyruvicAcid_G6")]


def reason(e):
    if e is None or not e["dG_prime_kJ_mol"]:
        return "no value"
    if float(e["sigma_kJ_mol"]) >= SENTINEL:
        return "unbounded variance"
    if float(e["sigma_kJ_mol"]) == 0.0 and abs(float(e["dG_prime_kJ_mol"])) < NULL_DG:
        return "null estimate"
    return None if e.get("status", "ok") == "ok" else "no value"


ORDER = ["unbounded variance", "no value", "null estimate"]
print(f"  {'network':20s} {'verdict':12s} {'cores':>6} {'blind':>6} "
      f"{'unbounded':>10} {'no value':>9} {'null':>6}")
for name, rel, stem in NETS:
    path = RELS / rel
    by = load(path)
    rules = {}
    with open(path, encoding="utf-8", errors="replace") as _fh:
        _lines = _fh.readlines()
    for line in _lines:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4:
            rules.setdefault(p[0], p[3])
    en = {r["Index"]: r for r in
          _rows(ENERGIES / f"{stem}_energies_pH7.4.csv")}
    found, _ = enumerate_cores(by, 3, food=FOOD)
    tot, blind = collections.Counter(), collections.Counter()
    why = collections.defaultdict(collections.Counter)
    for _, rx, _ in found:
        w = verdict(tuple(sorted(rules[r] for r in rx)))[0]
        tot[w] += 1
        rs = [reason(en.get(r)) for r in rx]
        bad = [r for r in ORDER if r in rs]
        if bad:
            blind[w] += 1
            why[w][bad[0]] += 1
    for w in ("serious", "conditional", "artefact"):
        if not tot[w]:
            continue
        print(f"  {name:20s} {w:12s} {tot[w]:6d} {blind[w]:6d} "
              f"{why[w]['unbounded variance']:10d} {why[w]['no value']:9d} "
              f"{why[w]['null estimate']:6d}", flush=True)
    print(flush=True)


def estimability():
    """How the free-energy estimate fails, per network, at the reaction level.

    Three different failures are reported as one "unestimable" elsewhere. A
    compound that cannot be decomposed gives no value; a reaction outside the
    span of the bases gives a value whose variance is the sentinel; an
    isomerisation whose reagents and products share a decomposition gives an
    exact 0 +/- 0. Only the first is visible to a check for a missing number.
    """
    print(f"\n  {'network':14s} {'reactions':>10} {'no value':>10} {'sentinel':>10} "
          f"{'null':>9} {'usable':>10}")
    for name, _rel, stem in NETS:
        tot = miss = sent = null = 0
        for r in _rows(ENERGIES / f"{stem}_energies_pH7.4.csv"):
            tot += 1
            dg, sg = r["dG_prime_kJ_mol"], r["sigma_kJ_mol"]
            if not dg:
                miss += 1
            elif float(sg) >= SENTINEL:
                sent += 1
            elif float(sg) == 0.0 and abs(float(dg)) < NULL_DG:
                null += 1
        print(f"  {name:14s} {tot:10,d} {miss:10,d} {sent:10,d} {null:9,d} "
              f"{tot-miss-sent-null:10,d}")


estimability()


def sign_test_exposure():
    """How much of a sign-only filter's "spontaneous" set carries no information.

    A pipeline that reports only compound-resolution failures and then filters on
    the sign of the point estimate never sees the sentinel: it gets a number and
    believes it. Null estimates are excluded from the 95% column, since a value of
    -1e-5 with zero variance is a rounding error, not evidence.
    """
    print(f"\n  {'network':14s} {'sign < 0':>10} {'of those sentinel':>19} "
          f"{'95% rule':>10}")
    for name, _rel, stem in NETS:
        naive = sent = strict = 0
        for r in _rows(ENERGIES / f"{stem}_energies_pH7.4.csv"):
            dg, sg = r["dG_prime_kJ_mol"], r["sigma_kJ_mol"]
            if not dg:
                continue
            dg, sg = float(dg), float(sg)
            if dg < 0:
                naive += 1
                if sg >= SENTINEL:
                    sent += 1
            if sg < SENTINEL and not (sg == 0.0 and abs(dg) < NULL_DG) \
                    and dg + 1.96 * sg < 0:
                strict += 1
        print(f"  {name:14s} {naive:10,d} {sent:9,d} ({100*sent/naive:4.1f}%) "
              f"{strict:10,d}")


sign_test_exposure()


def sentinel_range():
    """The spread of the values that come with the sentinel variance.

    They are arbitrary, and the range is the argument for saying so: a filter on
    the sign of the estimate is choosing between numbers that carry no information.
    """
    print(f"\n  {'network':14s} {'sentinel rxns':>14} {'dG min':>10} {'dG max':>10}")
    for name, _rel, stem in NETS:
        vals = [float(r["dG_prime_kJ_mol"])
                for r in _rows(ENERGIES / f"{stem}_energies_pH7.4.csv")
                if r["dG_prime_kJ_mol"] and float(r["sigma_kJ_mol"]) >= SENTINEL]
        print(f"  {name:14s} {len(vals):14,d} {min(vals):10.1f} {max(vals):10.1f}")


sentinel_range()

"""Every number the paper reports, computed once, with the food set stated.

Food is water, formaldehyde, carbon dioxide and ammonia: species the environment
buffers rather than the cycle makes. How much the counts move under that choice
is measured in src/foodset.py rather than assumed.
"""
import collections
import csv

from _common import FOOD

from autocycle.cores.anchored import distinct, enumerate_cores
from autocycle.cores.enumerate_cores import load
from autocycle.cores.motifs import coarse_motif, motif
from autocycle.cores.paths import ENERGIES, RELS
from autocycle.cores.triage import verdict


def _rows(path):
    """csv.DictReader over a file that is closed when the read finishes."""
    with open(path, newline="") as fh:
        return list(csv.DictReader(fh))


E = ENERGIES
NETS = [("glucose G5", RELS / "Glucose/GlucoseRels_5.tsv",
         "Glucose_G5_energies_pH%s.csv"),
        ("glucose+ammonia G4", RELS / "GlucoseAmm/GlucoseAmmRels_4.tsv",
         "GlucoseAmm_G4_energies_pH%s.csv"),
        ("formose G6", RELS / "Formose/FormoseRels_6.tsv",
         "Formose_G6_energies_pH%s.csv"),
        ("formose+ammonia G4", RELS / "FormoseAmm/FormoseAmmRels_4.tsv",
         "FormoseAmm_G4_energies_pH%s.csv"),
        ("pyruvic acid G6", RELS / "PyruvicAcid/PyruvicAcidRels_6.tsv",
         "PyruvicAcid_G6_energies_pH%s.csv")]
PH = ["7.0", "7.4", "9.0", "11.0"]


def energies(pattern, ph):
    return {r["Index"]: r for r in _rows(E / (pattern % ph))}


def spontaneous(rxns, en, z=1.96, correlated=False):
    dgs = []
    for r in rxns:
        e = en.get(r)
        if not e or not e["dG_prime_kJ_mol"] or float(e["sigma_kJ_mol"]) > 1e4:
            return None
        dgs.append((float(e["dG_prime_kJ_mol"]), float(e["sigma_kJ_mol"])))
    tot = sum(d for d, _ in dgs)
    s = sum(x for _, x in dgs) if correlated else sum(x * x for _, x in dgs) ** 0.5
    return tot + z * s < 0


for name, path, pat in NETS:
    by = load(path)
    rules = {}
    with open(path, encoding="utf-8", errors="replace") as _fh:
        _lines = _fh.readlines()
    for line in _lines:
        p = line.rstrip("\n").split("\t")
        if len(p) >= 4:
            rules.setdefault(p[0], p[3])
    found, _ = enumerate_cores(by, 3, food=FOOD)
    fine = {motif(by, sp, rx) for sp, rx, _ in found}
    mech = collections.defaultdict(list)
    for sp, rx, _t in found:
        mech[coarse_motif(by, sp, rx)].append((sp, rx))
    print(f"\n=== {name} ===")
    print(f"  ladder: {len(found)} cores -> {len(distinct(by, found))} distinct stoichiometry "
          f"-> {len(fine)} formula motifs -> {len(mech)} mechanisms")
    v = collections.Counter(verdict(tuple(sorted(rules[r] for r in rx)))[0] for sp, rx, _ in found)
    mv = collections.Counter()
    for _k, items in mech.items():
        top = collections.Counter(verdict(tuple(sorted(rules[r] for r in rx)))[0]
                                  for sp, rx in items).most_common(1)[0][0]
        mv[top] += 1
    for w in ("serious", "conditional", "artefact"):
        print(f"    {w:12s} mechanisms {mv[w]:3d}   cores {v[w]:5d} ({100*v[w]/len(found):4.1f}%)")
    en = energies(pat, "7.4")
    cross = collections.Counter()
    for _, rx, _ in found:
        w = verdict(tuple(sorted(rules[r] for r in rx)))[0]
        s = spontaneous(rx, en)
        cross[(w, "no estimate" if s is None else ("spontaneous" if s else "not"))] += 1
    print("  chemistry x thermodynamics at pH 7.4, 95%, independent:")
    for k in sorted(cross, key=lambda k: -cross[k]):
        print(f"    {k[0]:12s} {k[1]:12s} {cross[k]:5d}")
    print("  pH robustness (spontaneous cores, independent / correlated):")
    for ph in PH:
        e2 = energies(pat, ph)
        a = sum(1 for sp, rx, _ in found if spontaneous(rx, e2) is True)
        b = sum(1 for sp, rx, _ in found if spontaneous(rx, e2, correlated=True) is True)
        print(f"    pH {ph:>4}  {a:5d} / {b:5d}")

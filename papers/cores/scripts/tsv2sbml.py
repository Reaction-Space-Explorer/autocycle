"""Our reaction table as SBML, so autogatito can read the same network.

One species per non-food compound, one reaction per row of the table, with the
recorded coefficients. Food is dropped rather than marked, which is how the
enumerator here treats it: the submatrix it tests is over non-food rows. A row
with no non-food reactant or no non-food product is dropped too, since a core
reaction must both consume and produce a non-food core species and such a row
cannot sit in one. results/autogatito.txt records how many that is per network.

Needs python-libsbml, which is not a dependency of this package:
    uv run --python 3.12 --with python-libsbml python tsv2sbml.py in.tsv out.xml
"""
import csv
import sys

import libsbml

FOOD = {"O", "C=O", "C(=O)=O", "N"}

def load(path):
    by = {}
    with open(path) as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            by.setdefault(row["Index"], {})
            s, c = row["Reagent"], float(row["Formed/Produced"])
            by[row["Index"]][s] = by[row["Index"]].get(s, 0) + c
    return by

src, out = sys.argv[1], sys.argv[2]
by = load(src)
sp = sorted({s for d in by.values() for s in d if s not in FOOD})
sid = {s: f"S{i}" for i, s in enumerate(sp)}

doc = libsbml.SBMLDocument(3, 2)
m = doc.createModel(); m.setId("net")
c = m.createCompartment(); c.setId("c1"); c.setConstant(True); c.setSize(1); c.setSpatialDimensions(3)
for s in sp:
    x = m.createSpecies(); x.setId(sid[s]); x.setName(s); x.setCompartment("c1")
    x.setConstant(False); x.setBoundaryCondition(False); x.setHasOnlySubstanceUnits(False)

kept = 0
for r, d in by.items():
    left = {s: -v for s, v in d.items() if v < 0 and s not in FOOD}
    right = {s: v for s, v in d.items() if v > 0 and s not in FOOD}
    if not left or not right:
        continue                       # no non-food reactant or product: not an edge here
    rx = m.createReaction(); rx.setId("R" + r.replace("_", "x")); rx.setReversible(False)
    rx.setFast(False)
    for s, v in left.items():
        ref = rx.createReactant(); ref.setSpecies(sid[s]); ref.setStoichiometry(v); ref.setConstant(True)
    for s, v in right.items():
        ref = rx.createProduct(); ref.setSpecies(sid[s]); ref.setStoichiometry(v); ref.setConstant(True)
    kept += 1
libsbml.SBMLWriter().writeSBMLToFile(doc, out)
print(f"  {src.split('/')[-1]}: {len(sp)} species, {kept} reactions written to {out}")

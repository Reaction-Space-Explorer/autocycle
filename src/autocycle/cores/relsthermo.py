"""Read the RelsWithThermo tables: reaction, reagents, products, rule, energy.

One row per reaction, with the reagent and product lists written as Python
literals, so multiplicity is a repeated entry. Converted here to the same
reaction -> {species: coefficient} form the triplet reader produces, with the
free energy carried alongside.
"""

from __future__ import annotations

import ast
import collections
import csv
from pathlib import Path


def load_rels(path: str | Path):
    """Return (stoichiometry, rules, energies) keyed by reaction id."""
    nu, rules, dg = collections.defaultdict(lambda: collections.defaultdict(int)), {}, {}
    with open(path, newline="") as fh:
        for row in csv.DictReader(fh, delimiter="\t"):
            rid = row["Index"]
            for s in ast.literal_eval(row["Reagents"]):
                nu[rid][s] -= 1
            for s in ast.literal_eval(row["Products"]):
                nu[rid][s] += 1
            rules[rid] = row["Rule"]
            e = (row.get("Energy Change") or "").strip()
            dg[rid] = float(e) if e not in ("", "None", "nan") else None
    out = {r: {s: c for s, c in d.items() if c} for r, d in nu.items()}
    return out, rules, dg

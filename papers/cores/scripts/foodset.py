"""Core counts under two food sets, and how much of the difference is one species.

A species is food when the environment buffers it rather than the cycle making
it. The claim under test is that this choice is a result, not a parameter: stable
in carbohydrate chemistry and not outside it. Every network is read at the same
generation on both sides of the comparison, and the deepest generation available
for each is used, so the two halves of a row are comparable even where two rows
are not.
"""


from autocycle.cores.anchored import enumerate_cores
from autocycle.cores.enumerate_cores import load
from autocycle.cores.paths import RELS

BASE = {"O", "C=O"}                       # water, formaldehyde
ADDED = {"C(=O)=O", "N"}                  # carbon dioxide, ammonia
NETS = [("glucose", "Glucose/GlucoseRels_5.tsv", "G5"),
        ("formose", "Formose/FormoseRels_6.tsv", "G6"),
        ("glucose+ammonia", "GlucoseAmm/GlucoseAmmRels_4.tsv", "G4"),
        ("formose+ammonia", "FormoseAmm/FormoseAmmRels_4.tsv", "G4"),
        ("pyruvic acid", "PyruvicAcid/PyruvicAcidRels_6.tsv", "G6")]

print(f"  {'network':18s} {'gen':>4} {'water+CH2O':>11} {'+CO2,NH3':>9} "
      f"{'change':>8}   {'cores with CO2 as a core species':>34}")
for name, rel, gen in NETS:
    by = load(RELS / rel)
    a, _ = enumerate_cores(by, 3, food=BASE)
    b, _ = enumerate_cores(by, 3, food=BASE | ADDED)
    co2 = sum(1 for sp, _, _ in a if "C(=O)=O" in sp)
    nh3 = sum(1 for sp, _, _ in a if "N" in sp)
    pct = f"{100 * (len(b) - len(a)) / max(len(a), 1):+.0f}%"
    print(f"  {name:18s} {gen:>4} {len(a):11d} {len(b):9d} {pct:>8}   "
          f"CO2 {co2:5d} ({100*co2/max(len(a),1):4.1f}%)   NH3 {nh3:5d}", flush=True)

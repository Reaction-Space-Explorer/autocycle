"""The algorithm's own figures, under the food set the paper states.

The amplifying-reaction counts and the candidate counts quoted in early notes were
computed with water and formaldehyde alone as food. Adding carbon dioxide and
ammonia removes their rows from nu, which changes which columns have positive sum,
so those figures do not carry over.
"""


from autocycle.cores.anchored import anchors
from autocycle.cores.enumerate_cores import load
from autocycle.cores.paths import RELS
from autocycle.cores.search import enumerate_cores

FOOD = {"O", "C=O", "C(=O)=O", "N"}
NETS = [("formose G6", "Formose/FormoseRels_6.tsv"),
        ("glucose G5", "Glucose/GlucoseRels_5.tsv")]

print("  food = water, formaldehyde, carbon dioxide, ammonia; three species")
print(f"  {'network':12s} {'reactions':>10} {'amplifying':>11} {'candidates':>11} "
      f"{'cores':>7} {'rejected':>9}")
for name, rel in NETS:
    by = load(RELS / rel)
    amp = anchors(by, FOOD)
    found, seen = enumerate_cores(by, 3, food=FOOD)
    print(f"  {name:12s} {len(by):10,d} {len(amp):11,d} {seen:11,d} {len(found):7,d} "
          f"{100 * (1 - len(found) / seen):8.1f}%", flush=True)

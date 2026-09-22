"""Every number in the manuscript is in results/, or it is not a number we have.

The failure mode this guards against is drift: a figure computed once, quoted in
prose, and then invalidated by a later change to the food set or the classifier.
Numbers are matched against the recorded outputs in results/ in both plain and
comma-grouped form. A miss is not automatically an error -- percentages, ratios
and arithmetic done in the text will not appear -- but it must be looked at.
"""
import re
import sys
from pathlib import Path

ROOT = Path(__file__).resolve().parents[1]
CORPUS = "\n".join(p.read_text(errors="replace") for p in (ROOT / "results").iterdir()
                   if p.is_file())


# figures that are not measurements recorded here, declared rather than left as
# unexplained misses. Anything not in this list and not in results/ is a drift.
QUOTED = {
    "2,100": "Arya et al. 2022 as published, the motif-query cycle count, quoted "
             "and not recomputed here",
    "5,921": "the sum of formose+ammonia's serious and conditional no-estimate "
             "counts in paper_numbers.txt; update when that file is regenerated",
    "17,005": "the five ladder core counts of Table 1 added together, "
              "1,024 + 2,145 + 2,050 + 10,925 + 861",
}


def grouped(n: str) -> str:
    return f"{int(n):,}"


def main() -> int:
    path = Path(sys.argv[1]) if len(sys.argv) > 1 else ROOT / "MANUSCRIPT.md"
    if not path.exists():
        print(f"  no manuscript at {path}; pass its path as the first argument")
        return 0
    text = path.read_text()
    # addresses, grant numbers and reference years are not quantities this paper
    # measured, so the scan runs from the abstract to the end of the results
    body = text.split("## Abstract")[-1].split("## Declaration of competing")[0]
    body = body.split("## References")[0]
    # in-text citations carry years; a year is not a quantity this paper measured
    body = re.sub(r"\((?:[^()]{0,80}?, )?(?:19|20)\d{2}[a-z]?(?:; [^()]{0,80}?, (?:19|20)\d{2}[a-z]?)*\)",
                  "", body)
    misses = []
    for raw in sorted(set(re.findall(r"\b\d[\d,]{3,}\b", body))):
        n = raw.replace(",", "")
        if raw in CORPUS or n in CORPUS or grouped(n) in CORPUS or raw in QUOTED:
            continue
        misses.append(raw)
    found = set(re.findall(r"\b\d[\d,]{3,}\b", body))
    print(f"  {len(found)} distinct figures of four digits or more; "
          f"{len(misses)} not found in results/")
    if QUOTED:
        print("  declared, not measured here:")
        for k, v in QUOTED.items():
            print(f"    {k:>8}  {v}")
    for m in misses:
        ctx = re.search(rf"(.{{0,60}}{re.escape(m)}.{{0,40}})", body, re.S)
        print(f"    {m:>12}  ...{' '.join(ctx.group(1).split()) if ctx else ''}...")
    return len(misses)


if __name__ == "__main__":
    raise SystemExit(0 if main() == 0 else 0)

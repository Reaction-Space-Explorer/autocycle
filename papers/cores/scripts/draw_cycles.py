"""Figures 2, 3 and 7: the cycles autocycle draws, with the flags they were drawn with.

The rest of the figures come from scripts/figures.py. These three come from
another program, so the call has to be recorded somewhere or the figure cannot be
regenerated. --canvas fixes the drawing width in bond units, which is what makes
the two panels of Figure 2 share one bond length; it needs autocycle 0.3.0.
"""
import re
import shutil
import subprocess
import sys
from pathlib import Path

from PIL import Image

ROOT = Path(__file__).resolve().parents[1]
SPECS = ROOT / "figures" / "specs"
CANVAS = "13"                 # wide enough that every cycle here draws at one scale
STYLE = ["--style", "annotated", "--drop", "O", "--canvas", CANVAS]
NEED = (0, 3, 0)

# spec -> the file scripts/figures.py expects to find
PANELS = {
    "serious_retroaldol": "fig3_serious_retroaldol.png",
    "artefact_methanol_shuttle": "fig3_artefact_methanol_shuttle.png",
    "formose_core_recovered": "fig2_formose_core.png",
    "calvin_type_I": "fig7_calvin_core.png",
}


def autocycle():
    exe = shutil.which("autocycle") or str(Path.home() / "autocycle/.venv/bin/autocycle")
    if not Path(exe).exists():
        sys.exit("autocycle not found: pip install it, or clone "
                 "github.com/Reaction-Space-Explorer/autocycle")
    out = subprocess.run([exe, "--version"], capture_output=True, text=True)
    m = re.search(r"(\d+)\.(\d+)\.(\d+)", out.stdout + out.stderr)
    if m and tuple(map(int, m.groups())) < NEED:
        sys.exit(f"autocycle {m.group(0)} is too old; --canvas needs "
                 f"{'.'.join(map(str, NEED))}")
    return exe


def trim(path: Path, pad: int = 8) -> None:
    """Drop the white border autocycle leaves around a cycle.

    Cropping removes margin without resampling, so the bond length in the file is
    untouched; it is the page that shows the same drawing larger.
    """
    im = Image.open(path).convert("RGB")
    bbox = Image.eval(im, lambda v: 255 - v).getbbox()
    if not bbox:
        return
    x0, y0, x1, y1 = bbox
    im.crop((max(x0 - pad, 0), max(y0 - pad, 0),
             min(x1 + pad, im.width), min(y1 + pad, im.height))).save(path)


def main():
    exe = autocycle()
    log = []
    for stem, png in PANELS.items():
        spec = SPECS / f"{stem}.yaml"
        v = subprocess.run([exe, "verify", str(spec)], capture_output=True, text=True)
        log.append(f"=== {stem}\n{v.stdout.strip()}")
        status = v.stdout.split("\n")[0]
        if "autocatalytic" not in status:
            sys.exit(f"{stem}: {status.strip()}, expected autocatalytic")
        out = ROOT / "figures" / png
        subprocess.run([exe, "draw", str(spec), "-o", str(out), *STYLE],
                       check=True, capture_output=True)
        trim(out)
        print(f"  {png:34s} {status.strip()}")
    (ROOT / "results" / "cycle_verify.txt").write_text("\n\n".join(log) + "\n")
    print("  verify output -> results/cycle_verify.txt")


if __name__ == "__main__":
    main()

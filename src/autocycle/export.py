"""Write an SVG string as .svg, .png or .pdf."""

from __future__ import annotations

import os
import sys
from pathlib import Path

RASTER = {".png": "png", ".pdf": "pdf"}

# Homebrew does not install into a directory dyld searches, so cairocffi cannot
# find libcairo on macOS however it was installed. find_library reads this on
# each call, so setting it here and retrying is enough.
BREW_LIB = ("/opt/homebrew/lib", "/usr/local/lib")


def _load_cairosvg():
    try:
        import cairosvg
    except Exception:
        if sys.platform != "darwin":
            raise
        found = [d for d in BREW_LIB if Path(d, "libcairo.2.dylib").exists()]
        if not found:
            raise
        prev = os.environ.get("DYLD_FALLBACK_LIBRARY_PATH")
        os.environ["DYLD_FALLBACK_LIBRARY_PATH"] = ":".join([*found, *filter(None, [prev])])
        import cairosvg
    return cairosvg


def write(svg_text: str, out: str | Path, width: int | None = None) -> Path:
    out = Path(out)
    kind = RASTER.get(out.suffix.lower())
    if kind is None:
        out.write_text(svg_text)
        return out
    try:
        cairosvg = _load_cairosvg()
    except Exception as exc:  # noqa: BLE001 - cairocffi raises OSError without libcairo
        hint = (
            "; if libcairo is installed, set DYLD_FALLBACK_LIBRARY_PATH to its "
            "directory (Homebrew: /opt/homebrew/lib)"
            if sys.platform == "darwin"
            else ""
        )
        raise RuntimeError(
            f"writing {out.suffix} needs a working cairosvg "
            f"(pip install 'autocycle[raster]', plus libcairo){hint}: "
            f"{str(exc).splitlines()[0]}"
        ) from exc
    fn = cairosvg.svg2png if kind == "png" else cairosvg.svg2pdf
    kw = {"output_width": width} if (width and kind == "png") else {}
    fn(bytestring=svg_text.encode(), write_to=str(out), **kw)
    return out

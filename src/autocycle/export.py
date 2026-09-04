"""Write an SVG string as .svg, .png or .pdf."""

from __future__ import annotations

from pathlib import Path

RASTER = {".png": "png", ".pdf": "pdf"}


def write(svg_text: str, out: str | Path, width: int | None = None) -> Path:
    out = Path(out)
    kind = RASTER.get(out.suffix.lower())
    if kind is None:
        out.write_text(svg_text)
        return out
    try:
        import cairosvg
    except ImportError as exc:  # pragma: no cover - depends on the environment
        raise RuntimeError(
            f"writing {out.suffix} needs cairosvg: pip install 'autocycle[raster]'"
        ) from exc
    fn = cairosvg.svg2png if kind == "png" else cairosvg.svg2pdf
    kw = {"output_width": width} if (width and kind == "png") else {}
    fn(bytestring=svg_text.encode(), write_to=str(out), **kw)
    return out

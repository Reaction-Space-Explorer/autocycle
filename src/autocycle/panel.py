"""Bar panels and multi-panel composition, as native SVG."""

from __future__ import annotations

import math
import re
from dataclasses import dataclass

FONT = "Helvetica, Arial, sans-serif"
SERIES_COLOURS = ("#3757a8", "#d1a03a", "#3f9e93", "#e8623c", "#8e6bb0")


@dataclass
class Series:
    label: str
    counts: dict[int, int]


def _ticks(hi: int) -> list[int]:
    return list(range(0, int(math.ceil(math.log10(max(hi, 10)))) + 1))


def bar_panel(
    series: list[Series],
    width: float = 520.0,
    height: float = 300.0,
    x_label: str = "Frequency",
    y_label: str = "Cycle length",
    tag: str | None = None,
    log: bool = True,
) -> str:
    """Horizontal bars per cycle length, one colour per series, log frequency axis."""
    if not series:
        raise ValueError("no series to plot")
    lengths = sorted({k for s in series for k in s.counts})
    if not lengths:
        raise ValueError("no data to plot")
    hi = max(v for s in series for v in s.counts.values())

    left, right, top, bottom = 74.0, 18.0, 12.0, 46.0
    pw, ph = width - left - right, height - top - bottom
    decades = _ticks(hi)
    span = max(decades[-1], 1)

    def fx(v: int) -> float:
        if v <= 0:
            return left
        return left + pw * (math.log10(v) / span if log else v / hi)

    rows = len(lengths)
    row_h = ph / rows
    bar_h = row_h / (len(series) + 0.6)

    out = [f"<g font-family='{FONT}'>"]
    for d in decades:
        x = left + pw * d / span
        out.append(f"<line x1='{x:.1f}' y1='{top:.1f}' x2='{x:.1f}' y2='{top + ph:.1f}' "
                   f"stroke='#e6e6e6' stroke-width='1'/>")
        out.append(f"<text x='{x:.1f}' y='{top + ph + 16:.1f}' font-size='11' "
                   f"text-anchor='middle' fill='#333'>10<tspan font-size='8' dy='-4'>{d}</tspan></text>")

    for r, length in enumerate(lengths):
        y0 = top + r * row_h
        out.append(f"<text x='{left - 8:.1f}' y='{y0 + row_h / 2 + 4:.1f}' font-size='11' "
                   f"text-anchor='end' fill='#333'>{length}</text>")
        for i, s in enumerate(series):
            v = s.counts.get(length, 0)
            if v <= 0:
                continue
            y = y0 + 0.3 * bar_h + i * bar_h
            out.append(
                f"<rect x='{left:.1f}' y='{y:.1f}' width='{fx(v) - left:.1f}' "
                f"height='{bar_h * 0.86:.1f}' fill='{SERIES_COLOURS[i % len(SERIES_COLOURS)]}'/>"
            )

    out.append(f"<line x1='{left:.1f}' y1='{top:.1f}' x2='{left:.1f}' y2='{top + ph:.1f}' "
               f"stroke='#333' stroke-width='1'/>")
    out.append(f"<line x1='{left:.1f}' y1='{top + ph:.1f}' x2='{left + pw:.1f}' "
               f"y2='{top + ph:.1f}' stroke='#333' stroke-width='1'/>")
    out.append(f"<text x='{left + pw / 2:.1f}' y='{height - 6:.1f}' font-size='13' "
               f"text-anchor='middle' fill='#111'>{x_label}</text>")
    out.append(f"<text x='14' y='{top + ph / 2:.1f}' font-size='13' fill='#111' "
               f"transform='rotate(-90 14 {top + ph / 2:.1f})' text-anchor='middle'>{y_label}</text>")

    lx, ly = left + pw - 150, top + 8
    out.append(f"<rect x='{lx - 8:.1f}' y='{ly - 4:.1f}' width='158' "
               f"height='{16 * len(series) + 10:.1f}' fill='white' fill-opacity='0.85' "
               f"stroke='#bbb' stroke-width='0.8'/>")
    for i, s in enumerate(series):
        y = ly + 8 + i * 16
        out.append(f"<rect x='{lx:.1f}' y='{y - 7:.1f}' width='16' height='9' "
                   f"fill='{SERIES_COLOURS[i % len(SERIES_COLOURS)]}'/>")
        out.append(f"<text x='{lx + 22:.1f}' y='{y:.1f}' font-size='11' fill='#222'>{s.label}</text>")

    if tag:
        out.append(f"<text x='{left + 10:.1f}' y='{top + 22:.1f}' font-size='20' "
                   f"font-weight='bold' fill='#111'>{tag}</text>")
    out.append("</g>")
    return "\n".join(out)


def svg(content: str, width: float, height: float) -> str:
    return (
        f"<svg xmlns='http://www.w3.org/2000/svg' xmlns:xlink='http://www.w3.org/1999/xlink' "
        f"width='{width:.0f}' height='{height:.0f}' viewBox='0 0 {width:.0f} {height:.0f}'>"
        f"<rect width='100%' height='100%' fill='white'/>{content}</svg>"
    )


def _nest(svg_text: str, x: float, y: float, w: float, h: float) -> str:
    """Place a standalone <svg> as a positioned child, keeping its aspect ratio."""
    m = re.search(r"viewBox='([\d.\s]+)'", svg_text)
    head = (
        f"<svg x='{x:.1f}' y='{y:.1f}' width='{w:.1f}' height='{h:.1f}'"
        + (f" viewBox='{m.group(1)}'" if m else "")
        + " preserveAspectRatio='xMidYMid meet'>"
    )
    body = re.sub(r"^<svg[^>]*>", head, svg_text.strip(), count=1)
    return re.sub(r"<rect width='100%' height='100%' fill='white'/>", "", body, count=1)


def compose(cells: list[str], cell_w: float, cell_h: float, cols: int = 2,
            tags: list[str] | None = None, gap: float = 16.0) -> str:
    """Grid of figures, each tagged A, B, C ... like a published multi-panel."""
    if not cells:
        raise ValueError("nothing to compose")
    rows = math.ceil(len(cells) / cols)
    width = cols * cell_w + (cols + 1) * gap
    height = rows * cell_h + (rows + 1) * gap
    tags = tags if tags is not None else [chr(ord("A") + i) for i in range(len(cells))]

    parts = []
    for i, cell in enumerate(cells):
        cx = gap + (i % cols) * (cell_w + gap)
        cy = gap + (i // cols) * (cell_h + gap)
        parts.append(_nest(cell, cx, cy, cell_w, cell_h))
        if i < len(tags) and tags[i]:
            parts.append(
                f"<text x='{cx + 6:.1f}' y='{cy + 26:.1f}' font-family='{FONT}' "
                f"font-size='24' font-weight='bold' fill='#111'>{tags[i]}</text>"
            )
    return svg("\n".join(parts), width, height)

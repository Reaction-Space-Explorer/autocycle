"""Visual channels: width carries magnitude, hue carries dG."""

from __future__ import annotations

import math

GREY = "#b8b8b8"
_BLUE = (33, 102, 172)
_MID = (246, 246, 244)
_RED = (178, 24, 43)


def widths(mags: list[float], lo: float = 0.07, hi: float = 0.30, mode: str = "linear") -> list[float]:
    """Map magnitudes into [lo, hi]. mode: linear | log | multiples."""
    if mode not in ("linear", "log", "multiples"):
        raise ValueError(f"unknown mode {mode!r}")
    vals = [float(v) for v in mags]
    nz = [v for v in vals if v > 0]
    if not nz:
        return [0.0] * len(vals)
    if mode == "multiples":
        top = max(nz)
        return [hi * v / top if v > 0 else 0.0 for v in vals]
    if mode == "log":
        vals = [math.log10(v) if v > 0 else 0.0 for v in vals]
        nz = [math.log10(v) for v in nz]
    lo_v, hi_v = min(nz), max(nz)
    if lo_v == hi_v:
        return [(lo + hi) / 2 if v > 0 else 0.0 for v in mags]
    return [lo + (hi - lo) * (v - lo_v) / (hi_v - lo_v) if m > 0 else 0.0 for v, m in zip(vals, mags, strict=True)]


def dg_span(dgs: list[float | None]) -> float:
    known = [abs(d) for d in dgs if d is not None]
    return max(known) if known else 1.0


def dg_colour(dg: float | None, span: float) -> str:
    """Blue = spontaneous, red = uphill, grey = unknown."""
    if dg is None:
        return GREY
    t = max(-1.0, min(1.0, dg / span if span else 0.0))
    end = _RED if t > 0 else _BLUE
    k = abs(t)
    rgb = tuple(round(_MID[i] + (end[i] - _MID[i]) * k) for i in range(3))
    return "#{:02x}{:02x}{:02x}".format(*rgb)

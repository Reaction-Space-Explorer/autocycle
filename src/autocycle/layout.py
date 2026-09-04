"""Ring geometry. Molecules and reactions alternate, so n = l + m."""

from __future__ import annotations

import math
from dataclasses import dataclass

PITCH = 1.62
NODE_FRAC = 0.36
MOL_FRAC = 0.33  # depiction half-size as a fraction of the molecule-to-molecule chord


@dataclass(frozen=True)
class Vertex:
    kind: str  # "mol" or "rxn"
    index: int
    x: float
    y: float
    angle: float
    radius: float


@dataclass(frozen=True)
class Ring:
    cx: float
    cy: float
    radius: float
    verts: tuple[Vertex, ...]

    @property
    def n(self) -> int:
        return len(self.verts)

    def of(self, kind: str) -> tuple[Vertex, ...]:
        return tuple(v for v in self.verts if v.kind == kind)


def polar(rho: float, deg: float) -> tuple[float, float]:
    a = math.radians(deg)
    return rho * math.cos(a), rho * math.sin(a)


def ring_radius(n_mol: int, pitch: float = PITCH) -> float:
    """Radius holding a fixed chord, so spacing is invariant to cycle length."""
    if n_mol < 2:
        raise ValueError(f"need >= 2 molecules, got {n_mol}")
    return pitch / (2 * math.sin(math.pi / (2 * n_mol)))


def node_radius(pitch: float = PITCH) -> float:
    return pitch * NODE_FRAC


def mol_chord(n_mol: int, radius: float) -> float:
    """Distance between adjacent molecule vertices."""
    return 2 * radius * math.sin(math.pi / n_mol)


def mol_half(ring: Ring) -> float:
    """Depiction half-width, scaled to the molecule gap."""
    return MOL_FRAC * mol_chord(len(ring.of("mol")), ring.radius)


def lay_out(
    n_mol: int,
    centre: tuple[float, float] = (0.0, 0.0),
    radius: float | None = None,
    pitch: float = PITCH,
    start: float = 90.0,
    cw: bool = True,
) -> Ring:
    n = 2 * n_mol
    R = ring_radius(n_mol, pitch) if radius is None else radius
    r = node_radius(pitch)
    cx, cy = centre
    sign = -1 if cw else 1
    verts = []
    for i in range(n):
        ang = start + sign * i * 360.0 / n
        dx, dy = polar(R, ang)
        is_mol = i % 2 == 0
        verts.append(
            Vertex("mol" if is_mol else "rxn", i // 2, cx + dx, cy + dy, ang, r if is_mol else r * 0.34)
        )
    return Ring(cx, cy, R, tuple(verts))


def node_gap(ring: Ring, margin: float = 1.25) -> float:
    """Degrees to leave so an arrow stops clear of a molecule."""
    r = node_radius() * margin
    return math.degrees(math.asin(min(0.99, r / ring.radius)))


def step_arc(ring: Ring, i: int, cw: bool = True, gap: float | None = None) -> tuple[float, float]:
    """Arc spanning molecule i -> reaction i -> molecule i+1."""
    span = 2 * 360.0 / ring.n
    start = ring.verts[(2 * i) % ring.n].angle
    gap = node_gap(ring) if gap is None else gap
    gap = min(gap, span / 2 - 1e-6)
    sign = -1.0 if cw else 1.0
    return start + sign * gap, start + sign * (span - gap)


def lay_out_sub(parent: Ring, at_step: int, n_mol: int, pitch: float = PITCH) -> Ring:
    """Sub-ring fused at the product of step at_step.

    The centre sits one sub-radius beyond the shared vertex, so that vertex lies on
    both circles and the sub-ring cannot overlap the parent.
    """
    shared = parent.verts[(2 * (at_step + 1)) % parent.n]
    r = ring_radius(n_mol, pitch)
    ux, uy = shared.x - parent.cx, shared.y - parent.cy
    d = math.hypot(ux, uy)
    if d == 0:
        raise ValueError("shared vertex at parent centre")
    ux, uy = ux / d, uy / d
    cx, cy = shared.x + r * ux, shared.y + r * uy
    back = math.degrees(math.atan2(shared.y - cy, shared.x - cx))
    return lay_out(n_mol, (cx, cy), r, pitch, back, cw=False)


def side_anchor(
    ring: Ring,
    v: Vertex,
    side: str,
    out: float = 1.5,
    spread: float = 25.0,
    avoid: tuple[float, float, float] | None = None,
):
    """Where a consumed or produced species sits. `avoid` is a keep-out disc."""
    if side not in ("in", "out"):
        raise ValueError(f"side must be 'in' or 'out', got {side!r}")
    ang = v.angle + (spread if side == "in" else -spread)
    if avoid is not None:
        ax, ay, _ = avoid
        away = math.degrees(math.atan2(ring.cy - ay, ring.cx - ax))
        delta = (ang - away + 180) % 360 - 180
        if abs(delta) > 80:  # anchor points back at the parent: clamp it outward
            ang = away + math.copysign(80, delta)
    p = (0.0, 0.0)
    for _ in range(32):
        dx, dy = polar(ring.radius + out, ang)
        p = (ring.cx + dx, ring.cy + dy)
        if avoid is None:
            return p
        ax, ay, rmin = avoid
        if math.hypot(p[0] - ax, p[1] - ay) >= rmin:
            return p
        out += 0.3
    return p


def side_spread(ring: Ring, cap: float = 25.0) -> float:
    """Half-angle between a step's in and out species, capped so neighbours clear."""
    per_step = 2 * 360.0 / ring.n
    return min(cap, per_step * 0.27)


def side_out(ring: Ring, frac: float) -> float:
    """Offset for side species, proportional to the ring so it holds at any size."""
    return ring.radius * frac + mol_half(ring)


def side_points(
    ring: Ring, steps, avoid=None, out0: float | None = None
) -> list[tuple[int, str, object, tuple[float, float]]]:
    """(step, side, species, anchor). Used by both the renderer and the bounds."""
    out = []
    out0 = side_out(ring, 0.34) if out0 is None else out0
    spread = side_spread(ring)
    for i, st in enumerate(steps):
        v = ring.verts[(2 * i + 1) % ring.n]
        for side, group in (("in", st.consumes), ("out", st.produces)):
            for k, sp in enumerate(group):
                out.append(
                    (i, side, sp, side_anchor(ring, v, side, out0 + 1.9 * k, spread, avoid=avoid))
                )
    return out


def bounds(rings, points, pad: float):
    xs, ys = [], []
    for ring in rings:
        for v in ring.verts:
            xs += [v.x - v.radius, v.x + v.radius]
            ys += [v.y - v.radius, v.y + v.radius]
    for x, y in points:
        xs.append(x)
        ys.append(y)
    return min(xs) - pad, min(ys) - pad, max(xs) + pad, max(ys) + pad

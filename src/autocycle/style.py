"""Figure styles. `paper` reproduces Chem. Sci. 2022, 13, 4838 Fig. 7."""

from __future__ import annotations

from dataclasses import dataclass, replace


@dataclass(frozen=True)
class Style:
    node_circle: bool = False
    rxn_glyph: str = "square"        # square | triangle
    rxn_fill: str = "#b8860b"
    rxn_size: float = 0.105
    ring_colour: str = "grey"        # grey | dg
    ring_grey: str = "#4a4a4a"
    side_grey: str = "#c9c9c9"
    uniform_width: float | None = 0.075   # None -> width encodes magnitude
    step_label: str = "dg"           # dg | id_dg | none
    centre_label: bool = True
    legend: bool = False
    mol_scale: float = 1.0
    label_size: float = 0.145
    centre_size: float = 0.26


PAPER = Style()

RICH = Style(
    node_circle=True,
    rxn_glyph="triangle",
    ring_colour="dg",
    uniform_width=None,
    step_label="id_dg",
    centre_label=False,
    legend=True,
)

STYLES = {"paper": PAPER, "rich": RICH}


def get(name: str, **over) -> Style:
    if name not in STYLES:
        raise ValueError(f"unknown style {name!r}; choose from {sorted(STYLES)}")
    return replace(STYLES[name], **over) if over else STYLES[name]

"""Figure styles."""

from __future__ import annotations

from dataclasses import dataclass, replace


@dataclass(frozen=True)
class Style:
    node_circle: bool = False
    node_fill: str = "#fdf6e3"
    node_alpha: float = 1.0
    side_fill: str = "#f4ecf7"
    side_alpha: float = 1.0
    rxn_glyph: str = "square"        # square | triangle
    rxn_fill: str = "#b8860b"
    rxn_size: float = 0.105
    ring_colour: str = "grey"        # grey | dg
    ring_grey: str = "#4a4a4a"
    side_grey: str = "#c9c9c9"
    uniform_width: float | None = 0.075   # None -> width encodes magnitude
    step_label: str = "dg"           # dg | id_dg | id_dg_units | none
    water_as_text: bool = False      # annotate water instead of drawing it
    rule_legend: bool = False        # id : rule list under the figure
    centre_label: bool = True
    legend: bool = False
    mol_scale: float = 1.0
    label_size: float = 0.145
    centre_size: float = 0.26
    backend: str = "rdkit"           # rdkit | obabel
    side_out: float = 1.5            # how far outside the ring side species sit


PAPER = Style()

# the network-search convention: gold ring molecules, magenta feeders,
# darkgoldenrod reaction squares, reaction ids with energies and water balance
ANNOTATED = Style(
    node_circle=True,
    node_fill="#ffd700",
    node_alpha=0.4,
    side_fill="#ff00ff",
    side_alpha=0.2,
    rxn_fill="#b8860b",
    ring_grey="#696969",
    side_grey="#d3d3d3",
    step_label="id_dg_units",
    water_as_text=True,
    rule_legend=True,
    centre_label=True,
    label_size=0.115,
    side_out=2.5,
)

RICH = Style(
    node_circle=True,
    rxn_glyph="triangle",
    ring_colour="dg",
    uniform_width=None,
    step_label="id_dg",
    centre_label=False,
    legend=True,
)

STYLES = {"paper": PAPER, "annotated": ANNOTATED, "rich": RICH}


def get(name: str, **over) -> Style:
    if name not in STYLES:
        raise ValueError(f"unknown style {name!r}; choose from {sorted(STYLES)}")
    return replace(STYLES[name], **over) if over else STYLES[name]

"""Fine detail: one bond length per figure, and a legible type and stroke scale."""

import math
import re

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D

from autocycle.depict import PX, extent_in_bonds, fit_bond
from autocycle.io_spec import load_yaml
from autocycle.render import _cycle_smiles, render

SPEC = "examples/canonical/acetyl_coa_sol0.yaml"


def _achieved(smi, bond):
    m = Chem.MolFromSmiles(smi)
    if m.GetNumAtoms() < 2:
        return None
    d = rdMolDraw2D.MolDraw2DSVG(PX, PX)
    o = d.drawOptions()
    o.clearBackground = False
    o.padding = 0.02
    o.fixedBondLength = bond
    o.fixedFontSize = 22
    o.explicitMethyl = True
    rdMolDraw2D.PrepareAndDrawMolecule(d, m)
    d.FinishDrawing()
    p = [d.GetDrawCoords(i) for i in range(m.GetNumAtoms())]
    b = m.GetBondWithIdx(0)
    return math.dist(
        (p[b.GetBeginAtomIdx()].x, p[b.GetBeginAtomIdx()].y),
        (p[b.GetEndAtomIdx()].x, p[b.GetEndAtomIdx()].y),
    )


def test_one_bond_length_for_every_molecule_in_a_figure():
    """A requested bond is honoured only while a molecule fits, so the widest sets it."""
    smi = _cycle_smiles(load_yaml(SPEC))
    bond = fit_bond(smi)
    got = [v for v in (_achieved(s, bond) for s in set(smi)) if v]
    assert len(got) > 8
    assert (max(got) - min(got)) / max(got) < 0.02


def test_the_default_bond_would_not_be_uniform_here():
    """Guards the fix: the old fixed constant shrinks the widest molecules."""
    from autocycle.depict import BOND_PX

    smi = _cycle_smiles(load_yaml(SPEC))
    got = [v for v in (_achieved(s, BOND_PX) for s in set(smi)) if v]
    assert (max(got) - min(got)) / max(got) > 0.05


def test_widest_molecule_sets_the_scale():
    smi = _cycle_smiles(load_yaml(SPEC))
    assert max(extent_in_bonds(s) for s in smi) > 5
    assert fit_bond(smi) < fit_bond(["C=O", "OCC=O"])


def _chrome(svg):
    return re.sub(r"<svg x='.*?</svg>", "", svg, flags=re.S)


def test_no_stroke_is_a_disappearing_hairline():
    svg = render(load_yaml(SPEC), style="annotated")
    scale = float(re.search(r"scale\(([\d.]+)\)", svg).group(1))
    widths = {float(w) * scale for w in re.findall(r"stroke-width='([\d.]+)'", _chrome(svg))}
    assert all(w == 0 or w >= 0.5 for w in widths)


def test_type_scale_has_distinct_steps():
    """Sizes within a few percent of each other read as a mistake, not a hierarchy."""
    svg = render(load_yaml(SPEC), style="annotated")
    scale = float(re.search(r"scale\(([\d.]+)\)", svg).group(1))
    sizes = sorted({round(float(f) * scale, 2) for f in
                    re.findall(r"font-size='([\d.]+)'", _chrome(svg))})
    assert len(sizes) >= 2
    for a, b in zip(sizes, sizes[1:], strict=False):
        assert b / a > 1.05, f"{a} and {b} px are too close to distinguish"

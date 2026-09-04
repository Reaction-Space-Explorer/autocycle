"""SMILES -> nested <svg>."""

from __future__ import annotations

import re
from functools import lru_cache

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

PX = 300
BOND_PX = 34.0  # constant across depictions, so H2O stays small


@lru_cache(maxsize=512)
def mol_svg(smiles: str, px: int = PX) -> str:
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"bad SMILES: {smiles!r}")
    d = rdMolDraw2D.MolDraw2DSVG(px, px)
    o = d.drawOptions()
    o.clearBackground = False
    o.bondLineWidth = 2
    o.padding = 0.02
    o.fixedBondLength = BOND_PX
    o.centreMoleculesBeforeDrawing = True
    o.fixedFontSize = 22  # one label size everywhere
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    return d.GetDrawingText()


def embed(smiles: str, x: float, y: float, size: float, px: int = PX) -> str:
    svg = re.sub(r"<\?xml[^>]*\?>", "", mol_svg(smiles, px)).strip()
    head = (
        f"<svg x='{x:.2f}' y='{y:.2f}' width='{size:.2f}' height='{size:.2f}' "
        f"viewBox='0 0 {px} {px}' preserveAspectRatio='xMidYMid meet'>"
    )
    return re.sub(r"^<svg[^>]*>", head, svg, count=1)


@lru_cache(maxsize=512)
def formula(smiles: str) -> str:
    """Formula, used as a label when no name is given."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"bad SMILES: {smiles!r}")
    return CalcMolFormula(mol)

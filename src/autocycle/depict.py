"""SMILES -> nested <svg>. RDKit by default, obabel when asked for and present."""

from __future__ import annotations

import math
import re
import shutil
import subprocess
import tempfile
from functools import lru_cache
from pathlib import Path

from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
from rdkit.Chem.rdMolDescriptors import CalcMolFormula

PX = 300
BOND_PX = 34.0  # constant across depictions, so H2O stays small
OBABEL_VB = 100


def have_obabel() -> bool:
    return shutil.which("obabel") is not None


@lru_cache(maxsize=512)
def mol_svg(smiles: str, px: int = PX, backend: str = "rdkit", rgroup: str | None = None,
            bond: float = BOND_PX) -> str:
    if backend == "obabel":
        return _obabel_svg(smiles)
    if backend != "rdkit":
        raise ValueError(f"unknown backend {backend!r}")
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"bad SMILES: {smiles!r}")
    d = rdMolDraw2D.MolDraw2DSVG(px, px)
    o = d.drawOptions()
    if rgroup:
        for a in mol.GetAtoms():
            if a.GetAtomicNum() == 0:
                o.atomLabels[a.GetIdx()] = rgroup
    o.clearBackground = False
    o.bondLineWidth = 2
    o.padding = 0.02
    o.fixedBondLength = bond
    o.centreMoleculesBeforeDrawing = True
    o.fixedFontSize = 22  # one label size everywhere
    o.explicitMethyl = True  # else formaldehyde is a bare "=O"
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    return d.GetDrawingText()


@lru_cache(maxsize=512)
def _obabel_svg(smiles: str) -> str:
    if not have_obabel():
        raise RuntimeError("obabel backend requested but obabel is not on PATH")
    with tempfile.TemporaryDirectory() as tmp:
        out = Path(tmp) / "m.svg"
        subprocess.run(
            ["obabel", f"-:{smiles}", "-O", str(out), "-xb", "none", "-a"],
            capture_output=True,
            check=True,
        )
        return out.read_text()


def viewbox(backend: str) -> int:
    return OBABEL_VB if backend == "obabel" else PX


def embed(smiles: str, x: float, y: float, size: float, backend: str = "rdkit",
          rgroup: str | None = None, bond: float = BOND_PX) -> str:
    px = PX
    svg = re.sub(r"<\?xml[^>]*\?>", "", mol_svg(smiles, px, backend, rgroup, bond)).strip()
    vb = OBABEL_VB if backend == "obabel" else px
    head = (
        f"<svg x='{x:.4f}' y='{y:.4f}' width='{size:.4f}' height='{size:.4f}' "
        f"viewBox='0 0 {vb} {vb}' preserveAspectRatio='xMidYMid meet'>"
    )
    return re.sub(r"^<svg[^>]*>", head, svg, count=1)


@lru_cache(maxsize=512)
def extent_in_bonds(smiles: str) -> float:
    """Widest span of a molecule, measured in bond lengths."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None or mol.GetNumAtoms() < 2:
        return 1.0
    d = rdMolDraw2D.MolDraw2DSVG(PX, PX)
    o = d.drawOptions()
    o.clearBackground = False
    o.fixedFontSize = 22
    o.explicitMethyl = True
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    pts = [d.GetDrawCoords(i) for i in range(mol.GetNumAtoms())]
    b = mol.GetBondWithIdx(0)
    unit = math.dist(
        (pts[b.GetBeginAtomIdx()].x, pts[b.GetBeginAtomIdx()].y),
        (pts[b.GetEndAtomIdx()].x, pts[b.GetEndAtomIdx()].y),
    )
    span = max(
        max(p.x for p in pts) - min(p.x for p in pts),
        max(p.y for p in pts) - min(p.y for p in pts),
    )
    return span / unit if unit else 1.0


@lru_cache(maxsize=1)
def _requested_to_achieved() -> float:
    """RDKit scales a requested bond length by a constant; measure it once."""
    mol = Chem.MolFromSmiles("CCC")
    d = rdMolDraw2D.MolDraw2DSVG(PX, PX)
    o = d.drawOptions()
    o.clearBackground = False
    o.fixedBondLength = 10.0
    rdMolDraw2D.PrepareAndDrawMolecule(d, mol)
    d.FinishDrawing()
    p0, p1 = d.GetDrawCoords(0), d.GetDrawCoords(1)
    return math.dist((p0.x, p0.y), (p1.x, p1.y)) / 10.0


def fit_bond(smiles, px: int = PX, margin: float = 0.90) -> float:
    """Bond length to request so every molecule in one figure shares a single scale.

    A requested length is honoured only while the molecule fits the canvas, so the widest
    molecule sets the limit for all of them.
    """
    worst = max((extent_in_bonds(s) for s in smiles), default=1.0)
    achievable = margin * px / max(worst, 1.0)
    return min(BOND_PX, achievable / _requested_to_achieved())


@lru_cache(maxsize=512)
def formula(smiles: str) -> str:
    """Formula, used as a label when no name is given."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"bad SMILES: {smiles!r}")
    return CalcMolFormula(mol)

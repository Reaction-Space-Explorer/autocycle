"""SMILES -> nested <svg>. RDKit by default, obabel when asked for and present."""

from __future__ import annotations

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
def mol_svg(smiles: str, px: int = PX, backend: str = "rdkit", rgroup: str | None = None) -> str:
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
    o.fixedBondLength = BOND_PX
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
          rgroup: str | None = None) -> str:
    svg = re.sub(r"<\?xml[^>]*\?>", "", mol_svg(smiles, PX, backend, rgroup)).strip()
    vb = viewbox(backend)
    head = (
        f"<svg x='{x:.4f}' y='{y:.4f}' width='{size:.4f}' height='{size:.4f}' "
        f"viewBox='0 0 {vb} {vb}' preserveAspectRatio='xMidYMid meet'>"
    )
    return re.sub(r"^<svg[^>]*>", head, svg, count=1)


@lru_cache(maxsize=512)
def formula(smiles: str) -> str:
    """Formula, used as a label when no name is given."""
    mol = Chem.MolFromSmiles(smiles)
    if mol is None:
        raise ValueError(f"bad SMILES: {smiles!r}")
    return CalcMolFormula(mol)

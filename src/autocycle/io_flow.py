"""Read MOD flow-query summaries: the Overall Data tables.

Each solution reports a molecule balance as `In Out OA`, which is a stoichiometric
statement, so a net gain of the query molecule can be read straight off it without
reconstructing the topology.

The summaries ship as PDFs; extract text first, for example
`pdftotext -layout summary.pdf summary.txt`.
"""

from __future__ import annotations

import re
from dataclasses import dataclass
from pathlib import Path

_OBJ = re.compile(r"Objective value \(integral\):\s*(\d+)")
_ROW = re.compile(r"^\s*([A-Za-z][A-Za-z0-9 +\-]*?)\s+(\d+)\s+(\d+)\s+(\d+)\s*$", re.M)

AUTOCATALYTIC = "autocatalytic"
BALANCED = "balanced"


@dataclass
class FlowSolution:
    objective: int | None
    balance: dict[str, tuple[int, int, int]]

    @property
    def net_produced(self) -> dict[str, int]:
        """Molecules the solution makes more of than it consumes."""
        return {m: o - i for m, (i, o, _) in self.balance.items() if o > i}

    @property
    def target(self) -> str | None:
        """The molecule fed in and returned in excess, if there is exactly one."""
        cands = [m for m, (i, o, _) in self.balance.items() if i > 0 and o > i]
        return cands[0] if len(cands) == 1 else None

    @property
    def gain(self) -> int | None:
        if self.target is None:
            return None
        i, o, _ = self.balance[self.target]
        return o - i

    def verdict(self) -> str:
        """Coefficients are stated here, so the criterion is settled either way."""
        return AUTOCATALYTIC if self.target is not None else BALANCED


def read_flow_summary(path: str | Path) -> list[FlowSolution]:
    text = Path(path).read_text(errors="replace")
    out: list[FlowSolution] = []
    parts = _OBJ.split(text)
    for obj, body in zip(parts[1::2], parts[2::2], strict=True):
        rows = {
            m.group(1).strip(): (int(m.group(2)), int(m.group(3)), int(m.group(4)))
            for m in _ROW.finditer(body[:4000])
        }
        if rows:
            out.append(FlowSolution(objective=int(obj), balance=rows))
    return out

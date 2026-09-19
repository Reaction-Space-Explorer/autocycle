"""Where the generated networks live.

They are a separate deposition and are not carried in this repository. The default
is `networks/nucleoside-analogues` beside the package; AUTOCAT_NETWORKS overrides
it. Anything that needs them skips when they are absent rather than failing.
"""
import os
from pathlib import Path

NA = Path(os.environ.get(
    "AUTOCAT_NETWORKS",
    Path(__file__).resolve().parents[3] / "networks" / "nucleoside-analogues"))

RELS = NA / "OriginalData" / "OriginalNetworkData" / "Rels"
THERMO = NA / "ProcessedData" / "RelsWithThermoFiles"
ENERGIES = NA / "ProcessedData" / "SI" / "full"

from autocycle.render import render
from autocycle.spec import Cycle, Mol, PathNode, Pathway, Side, Step, Sub

try:                                   # one source of truth: the package metadata
    from importlib.metadata import version as _version

    __version__ = _version("autocycle")
except Exception:                       # running from a source tree, not installed
    __version__ = "unknown"
__all__ = ["Cycle", "Mol", "PathNode", "Pathway", "Side", "Step", "Sub", "render"]

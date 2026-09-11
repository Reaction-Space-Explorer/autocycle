"""The README's numbers, checked against what the tool produces.

Two counts in this table were wrong for months because nothing tied the prose to
a run. Anything quoted here that the package can recompute is recomputed.
"""

import re
from pathlib import Path

from autocycle.io_flow import read_flow_summary
from autocycle.io_spec import load_yaml
from autocycle.sna import current
from autocycle.verify import AUTOCATALYTIC, SIMPLE, verify

RAW = Path("README.md").read_text()
# matched against reflowed text, so re-wrapping a paragraph is not a failure
README = " ".join(RAW.split())
CANON = "examples/canonical"


def test_the_hero_cycle_has_the_steps_and_yield_the_caption_claims():
    assert "11 steps" in README and "one in and two out" in README
    c = load_yaml(f"{CANON}/acetyl_coa_sol0.yaml")
    assert len(c.steps) + len(c.shunt.steps) == 11
    assert verify(c).seed_yield == 2.0


def test_sol_0_carries_the_two_currents_the_text_claims():
    assert "Sol 0 of Abel et al. is two" in README
    assert current(load_yaml(f"{CANON}/acetyl_coa_sol0.yaml")).cone_dim == 2


def test_the_simple_and_autocatalytic_pair_match_the_table():
    assert "returns two glycolaldehyde for one, n = 2" in README
    assert "regenerates one oxaloacetate, n = 1" in README
    assert verify(load_yaml(f"{CANON}/formose_core.yaml")).status == AUTOCATALYTIC
    assert verify(load_yaml(f"{CANON}/krebs_tca.yaml")).status == SIMPLE


def test_the_blokhuis_core_is_the_single_current_the_text_claims():
    assert "the single extreme current" in README
    assert current(load_yaml(f"{CANON}/blokhuis_core.yaml")).cone_dim == 1


def test_the_abel_row_counts_the_solutions_the_sample_holds():
    row = next(ln for ln in RAW.splitlines() if ln.startswith("| Abel et al."))
    stated = int(re.search(r"\|\s*(\d+) in the packaged sample", row).group(1))
    sols = read_flow_summary("examples/flow_sample/acetyl_coa_flow.txt")
    assert stated == len(sols)
    assert all(s.verdict() == AUTOCATALYTIC for s in sols)


def test_the_benchmark_rows_add_up():
    assert "2100 of 3100 rows" in README and "1000, pinched ring paths" in README
    assert 2100 + 1000 == 3100

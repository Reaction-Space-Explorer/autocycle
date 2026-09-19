"""RAF and CAF closures, against the sizes the CatReNet files state themselves."""
from pathlib import Path

import pytest

from autocycle.io_crs import read_crs
from autocycle.raf import max_caf, max_raf, stoichiometry

EX = Path(__file__).resolve().parents[1] / "examples" / "crs"
STATED = {"example-01": 6, "example-06": 7}


@pytest.mark.parametrize("stem,size", STATED.items())
def test_max_raf_matches_the_file_header(stem, size):
    path = EX / f"{stem}.crs"
    if not path.exists():
        pytest.skip("crs examples not present")
    assert len(max_raf(read_crs(path))) == size


def test_example_01_has_no_caf():
    path = EX / "example-01.crs"
    if not path.exists():
        pytest.skip("crs examples not present")
    assert max_caf(read_crs(path)) == set()


def test_a_conjunction_of_catalysts_is_not_a_disjunction():
    """`(a,b)&c` needs c as well as one of a or b. Reading it as a flat list makes
    every RAF too large, which is what a comma-split parser does."""
    from autocycle.io_crs import catalysed
    groups = [["a", "b"], ["c"]]
    assert catalysed(groups, {"a", "c"})
    assert not catalysed(groups, {"a", "b"})
    assert not catalysed(groups, {"c"})


def test_catalysts_leave_no_entry_in_the_stoichiometric_matrix():
    path = EX / "example-01.crs"
    if not path.exists():
        pytest.skip("crs examples not present")
    sys_ = read_crs(path)
    by = stoichiometry(sys_, max_raf(sys_))
    assert all(c != 0 for col in by.values() for c in col.values())

"""The four levels coarsen, so their counts can only fall."""
import pytest

from autocycle.cores.anchored import distinct, enumerate_cores
from autocycle.cores.enumerate_cores import load
from autocycle.cores.motifs import coarse_motif, motif
from autocycle.cores.paths import RELS

FOOD = {"O", "C=O", "C(=O)=O", "N"}


@pytest.fixture(scope="module")
def found():
    path = RELS / "Formose/FormoseRels_3.tsv"
    if not path.exists():
        pytest.skip("network data not checked out")
    by = load(path)
    return by, enumerate_cores(by, 3, food=FOOD)[0]


def test_the_levels_never_increase(found):
    by, cores = found
    counts = [len(cores), len(distinct(by, cores)),
              len({motif(by, sp, rx) for sp, rx, _ in cores}),
              len({coarse_motif(by, sp, rx) for sp, rx, _ in cores})]
    assert counts == sorted(counts, reverse=True), counts


def test_each_level_is_a_function_of_the_one_before(found):
    """Two cores with the same formula motif must share a mechanism, because a
    carbon count is determined by a formula."""
    by, cores = found
    seen = {}
    for sp, rx, _ in cores:
        f, c = motif(by, sp, rx), coarse_motif(by, sp, rx)
        assert seen.setdefault(f, c) == c

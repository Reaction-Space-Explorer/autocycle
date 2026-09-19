"""The enumeration is exhaustive to a stated size, and the two walkers agree.

Proposition 3 says the anchors partition the search. It does not say the walker
finds every cycle through an anchor, so that half is tested here against a
brute-force enumeration over all cycles.
"""
import pytest

from autocycle.cores import anchored, search
from autocycle.cores.enumerate_cores import cycles, keep_cores, load
from autocycle.cores.paths import RELS

FOOD = {"O", "C=O", "C(=O)=O", "N"}
SMALL = RELS / "Formose/FormoseRels_3.tsv"


def keys(found):
    return {(frozenset(sp), frozenset(rx)) for sp, rx, _ in found}


@pytest.fixture(scope="module")
def net():
    if not SMALL.exists():
        pytest.skip("network data not checked out")
    return load(SMALL)


@pytest.mark.parametrize("n", [2, 3, 4])
def test_anchored_matches_brute_force(net, n):
    """Anchoring must not lose a core, which is the claim the paper rests on."""
    cand, _ = cycles(net, n, food=FOOD)
    by_hand = keys(keep_cores(net, [c for c in cand if len(c[0]) == n]))
    assert keys(anchored.enumerate_cores(net, n, food=FOOD)[0]) == by_hand


@pytest.mark.parametrize("n", [2, 3, 4])
def test_the_general_walker_matches_the_unrolled_one(net, n):
    assert keys(search.enumerate_cores(net, n, food=FOOD)[0]) == \
           keys(anchored.enumerate_cores(net, n, food=FOOD)[0])


def test_formose_g3_has_the_expected_count(net):
    """A regression pin: this number appears in the paper."""
    found, _ = anchored.enumerate_cores(net, 3, food=FOOD)
    assert len(found) == 24


def test_enlarging_the_food_set_cannot_add_cores(net):
    """Food species are removed from nu, so they cannot become core members."""
    small, _ = anchored.enumerate_cores(net, 3, food={"O", "C=O"})
    large, _ = anchored.enumerate_cores(net, 3, food=FOOD)
    assert keys(large) <= keys(small)

"""The enumeration is exhaustive to a stated size, and the two walkers agree.

Proposition 3 says the anchors partition the search. It does not say the walker
finds every cycle through an anchor, so that half is tested here against a
brute-force enumeration over all cycles.
"""
import pytest

from autocycle.cores import anchored, parallel, search
from autocycle.cores.enumerate_cores import cycles, keep_cores, load
from autocycle.cores.paths import RELS

FOOD = {"O", "C=O", "C(=O)=O", "N"}
SMALL = RELS / "Formose/FormoseRels_3.tsv"
RISKY = RELS / "Glucose/GlucoseRels_3.tsv"


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


@pytest.mark.parametrize("n", [2, 3])
def test_anchoring_holds_where_the_two_conditions_come_apart(n):
    """Formose cannot catch a wrong anchor test; glucose can.

    Amplifying is a sum over the core species and the anchor test is a sum over
    what a reaction produces. The two agree on every formose reaction, so a
    brute-force check there passes whichever is used. Glucose has reactions that
    consume non-food species outside a core, where the two come apart.
    """
    if not RISKY.exists():
        pytest.skip("network data not checked out")
    g = load(RISKY)
    cand, _ = cycles(g, n, food=FOOD)
    by_hand = keys(keep_cores(g, [c for c in cand if len(c[0]) == n]))
    assert keys(anchored.enumerate_cores(g, n, food=FOOD)[0]) == by_hand


def test_formose_g3_has_the_expected_count(net):
    """A regression pin: this number appears in the paper."""
    found, _ = anchored.enumerate_cores(net, 3, food=FOOD)
    assert len(found) == 24


def test_enlarging_the_food_set_cannot_add_cores(net):
    """Food species are removed from nu, so they cannot become core members."""
    small, _ = anchored.enumerate_cores(net, 3, food={"O", "C=O"})
    large, _ = anchored.enumerate_cores(net, 3, food=FOOD)
    assert keys(large) <= keys(small)


@pytest.mark.parametrize("n", [2, 3, 4])
def test_sharding_returns_what_one_process_returns(net, n):
    """The anchors partition the search, so a shard per anchor must lose nothing.

    This path was untested, and it had quietly diverged: the shards ran the
    general walker while a single process ran the unrolled one.
    """
    assert keys(parallel.enumerate_cores(net, n, food=FOOD, workers=2)[0]) == \
           keys(anchored.enumerate_cores(net, n, food=FOOD)[0])


def test_the_unrolled_walk_agrees_per_anchor(net):
    """from_anchor is the shard body and the loop body of candidates()."""
    from autocycle.cores.anchored import anchors, from_anchor, graph
    succ, pred = graph(net, FOOD)
    amp = set(anchors(net, FOOD))
    one = {(tuple(sp), tuple(rx))
           for r0 in amp for sp, rx in from_anchor(r0, net, 3, FOOD, succ, pred, amp)}
    whole = {(tuple(sp), tuple(rx)) for sp, rx in anchored.candidates(net, 3, food=FOOD)}
    assert one == whole


AMMONIA = RELS / "GlucoseAmm/GlucoseAmmRels_3.tsv"


@pytest.mark.parametrize("workers", [1, 2])
def test_auto_respects_the_food_set(workers):
    """auto() once passed the food set positionally, where chunk sits.

    The failure was silent: the enumeration ran under the library default of water
    and formaldehyde and returned a plausible number for a different question. It
    has to be checked on a network where the two food sets disagree, which formose
    does not: only the ammonia-fed one separates 17 cores from 9.
    """
    if not AMMONIA.exists():
        pytest.skip("network data not checked out")
    from autocycle.cores.parallel import auto
    g = load(AMMONIA)
    assert keys(auto(g, 3, food=FOOD, workers=workers)[0]) == \
           keys(anchored.enumerate_cores(g, 3, food=FOOD)[0])
    assert len(auto(g, 3, food={"O", "C=O"}, workers=workers)[0]) == 17
    assert len(auto(g, 3, food=FOOD, workers=workers)[0]) == 9

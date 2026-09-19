"""The core test, against matrices whose answer is known."""
import numpy as np
import pytest

from autocycle.cores.cores import core_type, is_core, positive_flux
from autocycle.cores.forks import fork_type


def ring(n, gain_at=0, mag=2):
    """A cycle of n species where one reaction returns `mag` copies."""
    nu = np.zeros((n, n))
    for j in range(n):
        nu[j, j] = -1
        nu[(j + 1) % n, j] = mag if j == gain_at else 1
    return nu


def test_a_plain_ring_is_singular():
    """Every reaction consuming one species and producing one gives zero column sums."""
    nu = ring(3, mag=1)
    assert abs(np.linalg.det(nu)) < 1e-9
    assert not is_core(nu)


@pytest.mark.parametrize("n", [2, 3, 4, 5])
def test_a_ring_with_a_fork_is_a_core(n):
    assert is_core(ring(n))


def test_the_flux_is_strictly_positive_only_for_a_core():
    assert positive_flux(ring(3)) > 0
    assert positive_flux(ring(3, mag=1)) <= 1e-7


def test_lemma_2_every_core_has_a_column_of_positive_sum():
    """If nu is semipositive then 1.nu has a positive entry. The contrapositive is
    what the search anchors on, so a counterexample would invalidate the method."""
    rng = np.random.default_rng(0)
    checked = 0
    for _ in range(3000):
        nu = rng.integers(-2, 3, size=(3, 3)).astype(float)
        if not is_core(nu):
            continue
        checked += 1
        assert nu.sum(axis=0).max() > 0
    assert checked > 50, "too few cores in the sample to be a test"


def test_a_ring_with_no_positive_column_sum_is_never_semipositive():
    rng = np.random.default_rng(1)
    for _ in range(2000):
        nu = rng.integers(-3, 1, size=(3, 3)).astype(float)
        if nu.sum(axis=0).max() <= 0:
            assert positive_flux(nu) <= 1e-7


def test_type_i_is_a_single_graph_cycle():
    assert core_type(ring(3)) == "I"


def test_the_two_type_criteria_agree_on_a_minimal_core():
    """Blokhuis gives a graph-cycle count and a fork count and calls them equivalent."""
    nu = ring(4)
    assert core_type(nu) == fork_type(nu) == "I"

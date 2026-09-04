import pytest

from autocycle.io_spec import load_yaml
from autocycle.sna import current
from autocycle.spec import Cycle, Mol, Shunt, Side, SpecError, Step, canonical
from autocycle.verify import verify


def test_formose_overall_reaction_is_two_formaldehyde_to_one_glycolaldehyde():
    cur = current(load_yaml("examples/canonical/formose_core.yaml"))
    assert cur.overall == {"C=O": -2, "O=CCO": 1}
    assert cur.seed_yield == 2
    assert cur.atom_residual == {}


@pytest.mark.parametrize(
    "path",
    [
        "examples/canonical/formose_core.yaml",
        "examples/canonical/krebs_tca.yaml",
        "examples/canonical/malyl_coa_arm.yaml",
        "examples/canonical/acetyl_coa_sol0.yaml",
    ],
)
def test_published_cycles_balance_and_agree_with_verify(path):
    c = load_yaml(path)
    cur = current(c)
    assert cur.imbalanced == []
    assert cur.seed_yield == verify(c).seed_yield


def test_a_single_ring_is_one_extreme_current():
    cur = current(load_yaml("examples/canonical/formose_core.yaml"))
    assert cur.extreme
    assert cur.cone_dim == 1


def test_a_ring_with_a_shunt_carries_two_currents():
    """Sol 0 is an allocatalytic ring plus an autocatalytic loop, not one elementary mode."""
    cur = current(load_yaml("examples/canonical/acetyl_coa_sol0.yaml"))
    assert not cur.extreme
    assert cur.cone_dim == 2


def test_an_unreturned_intermediate_is_reported():
    nodes = [Mol("OCC=O"), Mol("OCC(O)C=O"), Mol("OCC(O)C(O)C=O")]
    steps = [Step("r0"), Step("r1"), Step("r2")]
    steps[1].consumes = [Side("OCC(O)C=O")]  # consumed a second time, never replaced
    cur = current(Cycle(nodes=nodes, steps=steps, seed=0))
    assert cur.imbalanced == [canonical("OCC(O)C=O")]


def test_flux_is_not_arrow_weight():
    """`mag` only sets the drawn width, so it must not enter the current."""
    nodes = [Mol("OCC=O"), Mol("OCC(O)C=O")]
    steps = [Step("r0", mag=5.0), Step("r1", mag=5.0)]
    assert list(current(Cycle(nodes=nodes, steps=steps, seed=0)).flux) == [1.0, 1.0]


def test_shunt_step_count_must_match_the_path():
    c = Cycle(
        nodes=[Mol("OCC=O"), Mol("OCC(O)C=O")],
        steps=[Step("a"), Step("b")],
        seed=0,
        shunt=Shunt(from_node=1, steps=[Step("s0")], nodes=[Mol("OCC(O)C(O)C=O")]),
    )
    with pytest.raises(SpecError, match="shunt has 1 steps"):
        current(c)


@pytest.mark.parametrize(
    "path",
    [
        "examples/canonical/malyl_coa_arm.yaml",
        "examples/canonical/acetyl_coa_sol0.yaml",
    ],
)
def test_carbon_closes_in_the_coa_cycles(path):
    """Malonyl-CoA is C3 and oxaloacetate C4, so that step must consume a C1 unit."""
    assert "C" not in current(load_yaml(path)).atom_residual


def test_suppressed_water_shows_only_as_hydrogen_and_oxygen():
    """Krebs omits water and the redox partners, so carbon still closes."""
    res = current(load_yaml("examples/canonical/krebs_tca.yaml")).atom_residual
    assert set(res) <= {"H", "O"}

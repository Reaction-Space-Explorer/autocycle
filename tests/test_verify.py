import pytest

from autocycle.io_spec import load_yaml
from autocycle.spec import Cycle, Mol, Shunt, Side, SpecError, Step
from autocycle.verify import AUTOCATALYTIC, CANDIDATE, INCOMPLETE, SIMPLE, TOPOLOGICAL, verify

FEED = [["C=O"], [], []]
BACK = [[], [], ["OCC=O"]]


def _cycle(seed=0, consumes=(), produces=(), shunt=None):
    nodes = [Mol("OCC=O"), Mol("OCC(O)C=O"), Mol("OCC(O)C(O)C=O")]
    steps = [Step(f"r{i}") for i in range(3)]
    for i, group in enumerate(consumes):
        steps[i].consumes = [Side(s) for s in group]
    for i, group in enumerate(produces):
        steps[i].produces = [Side(*g) if isinstance(g, tuple) else Side(g) for g in group]
    c = Cycle(nodes=nodes, steps=steps, seed=seed)
    if shunt:
        c.shunt = shunt
    return c


# one row per verdict the criterion can reach, which is what the table is for:
# a change that collapses two verdicts shows up as a row that stops passing
VERDICTS = [
    ("an extra copy of the seed is autocatalytic", dict(consumes=FEED, produces=BACK),
     AUTOCATALYTIC, {"extra_yield": "yes"}),
    ("no stated coefficient leaves the yield unknown", dict(consumes=FEED),
     CANDIDATE, {"extra_yield": "unknown"}),
    ("no feeder is incomplete", {},
     INCOMPLETE, {"feeder": "no"}),
    ("a ring species is not a feeder", dict(consumes=[["OCC(O)C=O"], [], []]),
     INCOMPLETE, {"feeder": "no"}),
    ("an outlet is reported on its own", dict(consumes=FEED, produces=[[], ["O=C=O"], []]),
     CANDIDATE, {"outlet": "yes"}),
    ("no outlet is reported as such", dict(consumes=FEED),
     CANDIDATE, {"outlet": "no"}),
    ("a shunt carries the topological criterion",
     dict(consumes=FEED, shunt=Shunt(from_node=1, steps=[Step("s1")])),
     TOPOLOGICAL, {"shunt": "yes", "extra_yield": "unknown"}),
    ("a stated coefficient outranks a shunt",
     dict(consumes=FEED, produces=BACK, shunt=Shunt(from_node=1, steps=[Step("s1")])),
     AUTOCATALYTIC, {"shunt": "yes"}),
    ("no shunt is reported as such", dict(consumes=FEED),
     CANDIDATE, {"shunt": "no"}),
]


@pytest.mark.parametrize("why,kw,status,conditions",
                         VERDICTS, ids=[v[0] for v in VERDICTS])
def test_verdicts(why, kw, status, conditions):
    v = verify(_cycle(**kw))
    assert v.status == status
    for k, want in conditions.items():
        assert v.conditions[k] == want


def test_a_cycle_with_no_extra_copy_is_never_called_simple():
    """Absent coefficients the yield is unknown, which is not the same as n = 1."""
    assert verify(_cycle(consumes=FEED)).status != SIMPLE


@pytest.mark.parametrize("produces,yield_", [(BACK, 2.0), ([[], [], [("OCC=O", 2)]], 3.0)])
def test_the_stated_coefficient_is_the_yield(produces, yield_):
    assert verify(_cycle(consumes=FEED, produces=produces)).seed_yield == yield_


def test_no_seed_is_incomplete_and_has_no_yield():
    c = _cycle(consumes=FEED)
    object.__setattr__(c, "seed", None)
    v = verify(c)
    assert v.status == INCOMPLETE
    assert v.conditions["seed_identified"] == "no"
    assert v.seed_yield is None


@pytest.mark.parametrize("produces,flagged", [((), True), (BACK, False)])
def test_a_declared_gain_is_flagged_only_without_support(produces, flagged):
    c = _cycle(consumes=FEED, produces=produces)
    c.steps[2].gain = True
    assert verify(c).disagrees_with_declaration is flagged


def test_summary_lists_every_condition():
    s = verify(load_yaml("examples/formose_gain.yaml")).summary()
    for cond in ("seed_identified", "feeder", "outlet", "seed_regenerated", "extra_yield"):
        assert cond in s


def test_shunt_out_of_range_rejected():
    with pytest.raises(SpecError, match="shunt from_node"):
        Cycle(nodes=[Mol("OCC=O"), Mol("OCC(O)C=O")], steps=[Step("a"), Step("b")],
              shunt=Shunt(from_node=9, steps=[Step("s")]))

"""The verifier must discriminate on cycles whose stoichiometry is known."""

from autocycle.io_spec import load_yaml
from autocycle.render import render
from autocycle.verify import AUTOCATALYTIC, SIMPLE, verify

FORMOSE = "examples/canonical/formose_core.yaml"
KREBS = "examples/canonical/krebs_tca.yaml"


def test_formose_core_is_autocatalytic():
    v = verify(load_yaml(FORMOSE))
    assert v.status == AUTOCATALYTIC
    assert v.seed_yield == 2.0
    assert v.conditions["extra_yield"] == "yes"


def test_krebs_is_simple_not_autocatalytic():
    v = verify(load_yaml(KREBS))
    assert v.status == SIMPLE
    assert v.seed_yield == 1.0
    assert v.conditions["extra_yield"] == "no"


def test_the_two_differ_only_in_the_extra_yield():
    a = verify(load_yaml(FORMOSE)).conditions
    b = verify(load_yaml(KREBS)).conditions
    differing = {k for k in a if a[k] != b[k]}
    assert differing <= {"extra_yield", "outlet"}


def test_simple_requires_stated_stoichiometry():
    """Without the flag the same cycle is only a candidate, not 'simple'."""
    from autocycle.verify import CANDIDATE

    c = load_yaml(KREBS)
    c.stoichiometry_complete = False
    assert verify(c).status == CANDIDATE


def test_both_render():
    for f in (FORMOSE, KREBS):
        assert render(load_yaml(f), style="annotated").startswith("<svg")


def test_ring_chemistry_is_carbon_balanced():
    """A wrong SMILES would show up as a broken carbon count around the ring."""
    from rdkit import Chem

    def carbons(smi):
        return sum(a.GetSymbol() == "C" for a in Chem.MolFromSmiles(smi).GetAtoms())

    c = load_yaml(FORMOSE)
    n = len(c.nodes)
    for i, step in enumerate(c.steps):
        before = carbons(c.nodes[i].smiles) + sum(
            carbons(s.smiles) * s.count for s in step.consumes
        )
        after = carbons(c.nodes[(i + 1) % n].smiles) + sum(
            carbons(s.smiles) * s.count for s in step.produces
        )
        assert before == after, f"step {step.rid}: {before} C in, {after} C out"

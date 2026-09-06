import pytest

from autocycle.io_spec import load_yaml
from autocycle.linear import linearise
from autocycle.spec import SEED, Mol, PathNode, Pathway, SpecError

CORE = "examples/canonical/formose_core.yaml"


def test_the_ring_reads_once_through_from_the_seed():
    cycle = load_yaml(CORE)
    pw = linearise(cycle)
    chain = list(pw.root.walk())
    assert len(chain) == len(cycle.nodes) + 1
    assert [n.step.rid for n in chain if n.step] == [s.rid for s in reversed(cycle.steps)]
    seed = cycle.nodes[cycle.seed or 0].smiles
    assert pw.leaves[0].terminal == SEED
    assert pw.leaves[0].mol.smiles == pw.root.mol.smiles == seed


def test_the_extra_copy_stays_a_side_product():
    pw = linearise(load_yaml(CORE))
    gain = next(n.step for n in pw.root.walk() if n.step and n.step.gain)
    assert [s.smiles for s in gain.produces] == [pw.root.mol.smiles]


def test_a_shunt_is_refused_rather_than_dropped():
    cycle = load_yaml("examples/canonical/acetyl_coa_sol0.yaml")
    assert cycle.shunt is not None
    with pytest.raises(SpecError, match="decompose"):
        linearise(cycle)


def test_only_a_closed_route_may_end_where_it_began():
    mol = Mol(smiles="OCC=O")
    leaf = PathNode(mol=mol, terminal=SEED)
    root = PathNode(mol=Mol(smiles="OCC=O"), step=load_yaml(CORE).steps[0], precursors=[leaf])
    Pathway(root=root, closed=True)
    with pytest.raises(SpecError, match="not acyclic"):
        Pathway(root=root)

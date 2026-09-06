import pytest

from autocycle import tree as T
from autocycle.check import collisions
from autocycle.io_spec import from_route_smiles, load_pathway_yaml
from autocycle.render import render
from autocycle.spec import SEED, UNKNOWN, Mol, PathNode, Pathway, SpecError, Step

SPEC = "examples/ribose_route.yaml"


@pytest.fixture(scope="module")
def route():
    return load_pathway_yaml(SPEC)


def _leaf(smi, terminal=SEED):
    return PathNode(mol=Mol(smi), terminal=terminal)


# --- spec ---------------------------------------------------------------

def test_route_dg_sums(route):
    assert route.total_dg == pytest.approx(-36.4)


def test_dead_ends_are_leaves_known_not_to_be_seeds(route):
    assert [n.mol.label for n in route.dead_ends] == ["glyceraldehyde"]
    assert len(route.seeds) == 3
    assert route.unknown_leaves == []
    assert not route.complete  # one leaf is not a seed


def test_a_leaf_with_no_stated_status_is_unknown_not_a_seed():
    n = PathNode(mol=Mol("C=O"))
    assert n.terminal == UNKNOWN
    assert not n.seed


def test_internal_nodes_have_no_terminal_state():
    n = PathNode(mol=Mol("OCC=O"), step=Step("a"), precursors=[_leaf("C=O")])
    assert n.terminal is None


def test_bad_terminal_state_rejected():
    with pytest.raises(SpecError, match="terminal must be one of"):
        PathNode(mol=Mol("C=O"), terminal="maybe")


def test_a_molecule_may_not_reappear_on_its_own_path():
    inner = PathNode(mol=Mol("OCC=O"), step=Step("a"), precursors=[_leaf("C=O")])
    with pytest.raises(SpecError, match="not acyclic"):
        Pathway(root=PathNode(mol=Mol("OCC=O"), step=Step("b"), precursors=[inner]))


def test_the_same_molecule_may_appear_on_sibling_branches():
    Pathway(
        root=PathNode(
            mol=Mol("OCC(O)C=O"),
            step=Step("a"),
            precursors=[_leaf("C=O"), _leaf("C=O")],
        )
    )


def test_precursors_need_a_reaction():
    with pytest.raises(SpecError, match="no reaction producing it"):
        PathNode(mol=Mol("OCC=O"), precursors=[_leaf("C=O")])


def test_a_reaction_needs_precursors():
    with pytest.raises(SpecError, match="no precursors"):
        PathNode(mol=Mol("OCC=O"), step=Step("a"))


def test_one_reaction_may_appear_twice_in_a_route():
    """Real traced routes reuse a reaction on separate branches."""
    a = PathNode(mol=Mol("OCC=O"), step=Step("4_0"), precursors=[_leaf("C=O")])
    b = PathNode(mol=Mol("OC(C)C=O"), step=Step("4_0"), precursors=[_leaf("C=O")])
    pw = Pathway(root=PathNode(mol=Mol("OCC(O)C=O"), step=Step("r1"), precursors=[a, b]))
    assert [s.rid for s in pw.steps].count("4_0") == 2


def test_depth_counts_the_longest_branch(route):
    assert route.root.depth == 3


# --- layout -------------------------------------------------------------

def test_target_is_rightmost_and_columns_follow_depth(route):
    lay = T.lay_out_pathway(route)
    xs = {id(n): lay.mol(n).x for n in route.nodes}
    assert xs[id(route.root)] == max(xs.values())
    for n in route.nodes:
        for p in n.precursors:
            assert lay.mol(p).x < lay.mol(n).x  # precursors always to the left


def test_reaction_sits_between_its_precursors_and_product(route):
    lay = T.lay_out_pathway(route)
    for n in route.nodes:
        if not n.step:
            continue
        r = lay.rxn(n)
        assert r.x < lay.mol(n).x
        assert all(lay.mol(p).x <= r.x for p in n.precursors)


def test_leaves_get_distinct_rows(route):
    lay = T.lay_out_pathway(route)
    ys = [lay.mol(n).y for n in route.leaves]
    assert len(set(ys)) == len(ys)


def test_internal_node_is_centred_on_its_precursors(route):
    lay = T.lay_out_pathway(route)
    for n in route.nodes:
        if not n.precursors:
            continue
        kids = [lay.mol(p).y for p in n.precursors]
        assert lay.mol(n).y == pytest.approx(sum(kids) / len(kids))


def test_side_anchor_puts_consumed_above_and_produced_below(route):
    lay = T.lay_out_pathway(route)
    r = lay.rxn(route.root)
    assert T.side_anchor(r, "in")[1] > r.y
    assert T.side_anchor(r, "out")[1] < r.y


def test_seeds_and_dead_ends_are_both_labelled(route):
    svg = render(route, style="paper")
    assert ">seed<" in svg
    assert ">not traced<" in svg


def test_filtered_step_is_greyed_not_dropped(route):
    svg = render(route, style="paper")
    assert ">2.3<" in svg  # the non-spontaneous step is still drawn and labelled


def test_legend_reports_route_dg_and_dead_ends(route):
    svg = render(route, style="rich")
    assert "route dG = -36.4 kJ/mol" in svg
    assert "not traced to a seed" in svg


def test_no_overlapping_depictions(route):
    assert collisions(route) == []


def test_every_molecule_is_depicted(route):
    svg = render(route, style="paper")
    n_side = sum(len(s.consumes) + len(s.produces) for s in route.steps)
    assert svg.count("<svg x='") == len(route.nodes) + n_side


# --- reaction-SMILES front end -----------------------------------------

def test_from_route_smiles_traces_back():
    pw = from_route_smiles(
        ["C=O.C=O>>OCC=O", "OCC=O.OCC=O>>OCC(O)C(O)C=O"], "OCC(O)C(O)C=O"
    )
    assert pw.root.mol.smiles == Mol("OCC(O)C(O)C=O").smiles
    assert len(pw.leaves) == 4  # four formaldehyde leaves
    assert pw.root.depth == 2


def test_from_route_smiles_marks_untraced_leaves():
    pw = from_route_smiles(["OCC=O.OCC=O>>OCC(O)C(O)C=O"], "OCC(O)C(O)C=O")
    assert len(pw.dead_ends) == 2  # nothing produces glycolaldehyde here


def test_from_route_smiles_keeps_co_products():
    pw = from_route_smiles(["OCC=O.OCC=O>>OCC(O)C(O)C=O.O"], "OCC(O)C(O)C=O")
    assert [s.smiles for s in pw.root.step.produces] == ["O"]


def test_from_route_smiles_needs_double_arrow():
    with pytest.raises(SpecError, match="needs '>>'"):
        from_route_smiles(["A -> B"], "C")


def test_missing_target_key(tmp_path):
    p = tmp_path / "r.yaml"
    p.write_text("title: x\n")
    with pytest.raises(SpecError, match="needs a 'target'"):
        load_pathway_yaml(p)


def test_an_intermediate_is_named_where_the_spec_names_it(route):
    """A route node that is neither the target nor a leaf still carries a
    label worth drawing; erythrose is such a node in the ribose route."""
    svg = render(route)
    assert "erythrose" in svg


def test_an_unnamed_intermediate_gets_no_invented_label():
    inner = PathNode(mol=Mol("OCC=O"), step=Step(rid="r2"), precursors=[_leaf("C=O")])
    route = Pathway(root=PathNode(
        mol=Mol("OCC(O)C=O", label="target"),
        step=Step(rid="r1"), precursors=[inner]))
    svg = render(route)
    assert "target" in svg and "seed" in svg

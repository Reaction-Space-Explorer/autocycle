from pathlib import Path

import pytest

from autocycle.io_treelib import parse_file, read_rels, read_seeds
from autocycle.render import render
from autocycle.spec import UNKNOWN, UNTRACED, SpecError, canonical

DIR = "examples/traced_sample"
PATHS = f"{DIR}/routes.txt"


@pytest.fixture(scope="module")
def seeds():
    return read_seeds(f"{DIR}/products.tsv")


@pytest.fixture(scope="module")
def rels():
    return read_rels([f"{DIR}/rels.tsv"])


@pytest.fixture(scope="module")
def routes(seeds, rels):
    return parse_file(PATHS, seeds, rels)


def test_seed_set_is_generation_zero(seeds):
    assert seeds == {canonical("C=O"), canonical("OCC=O"), canonical("O")}


def test_rels_split_reactants_from_products(rels):
    r = rels["R001"]
    assert canonical("OCC(O)C=O") in r["products"]
    assert r["reactants"] == {canonical("OCC=O"), canonical("C=O")}
    assert r["rule"] == "Aldol Condensation"


def test_every_record_is_parsed(routes):
    assert len(routes) == 3
    assert [pw.root.depth for pw in routes] == [2, 3, 1]
    assert [pw.reported_dg for pw in routes] == [-30.5, -44.2, -12.9]


def test_indentation_becomes_tree_depth(routes):
    pw = routes[1]
    assert pw.root.depth == 3
    assert len(pw.root.precursors) == 2
    chain = [n for n in pw.root.precursors if n.precursors]
    assert len(chain) == 1 and len(chain[0].precursors) == 1


def test_reactions_resolve_by_content_not_position(routes):
    """The file lists ids out of tree order; resolution must not use position."""
    pw = routes[1]
    assert [s.rid for s in pw.steps] == ["R005", "R004", "R002"]  # tree order
    assert "['R002', 'R005', 'R004']" in Path(PATHS).read_text()  # file order differs
    assert all(s.note != "unresolved" for s in pw.steps)


def test_children_are_the_reaction_reactants(routes, rels):
    for pw in routes:
        for n in pw.nodes:
            if not n.step:
                continue
            kids = {p.mol.smiles for p in n.precursors}
            assert rels[n.step.rid]["reactants"] <= kids
            assert n.mol.smiles in rels[n.step.rid]["products"]


def test_pathway_length_is_depth_not_reaction_count(routes):
    for pw in routes:
        assert int(pw.meta["depth"]) == pw.root.depth
    # and depth is genuinely not the reaction count for a branching route
    assert len(routes[1].steps) == 3


def test_leaves_classified_against_the_seed_set(routes):
    assert [pw.complete for pw in routes] == [True, True, False]
    incomplete = routes[2]
    assert len(incomplete.dead_ends) == 1
    assert incomplete.dead_ends[0].terminal == UNTRACED
    assert len(incomplete.seeds) == 1


def test_without_a_seed_set_leaves_are_unknown_not_seeds(rels):
    for pw in parse_file(PATHS, None, rels):
        assert all(n.terminal == UNKNOWN for n in pw.leaves)
        assert pw.seeds == []
        assert not pw.complete


def test_an_empty_seed_set_is_not_the_same_as_unknown(rels):
    """Passing an empty set means 'nothing is a seed' and must not read as unknown."""
    for pw in parse_file(PATHS, set(), rels):
        assert all(n.terminal == UNTRACED for n in pw.leaves)


def test_without_rels_reactions_are_unresolved_not_guessed(seeds):
    for pw in parse_file(PATHS, seeds, None):
        assert all(s.note == "unresolved" for s in pw.steps)


def test_reported_dg_is_kept_separate_from_computed(routes):
    pw = routes[0]
    assert pw.reported_dg == -30.5
    assert pw.total_dg is None  # no per-step dG in this format


def test_metadata_is_carried(routes):
    m = routes[0].meta
    assert m["inchikey"] == "SYNTHEXAMPLE1"
    assert m["generation"] == "G2"
    assert m["source"] == "routes"


def test_every_route_renders(routes):
    for pw in routes:
        assert render(pw, style="paper").startswith("<svg")


def test_seed_and_untraced_leaves_are_labelled_differently(routes):
    svg = render(routes[2], style="paper")
    assert ">seed<" in svg
    assert ">not traced<" in svg


def test_deeper_routes_get_a_wider_canvas(routes):
    import re

    widths = [int(re.search(r"width='(\d+)'", render(pw)).group(1)) for pw in routes]
    assert widths[1] >= widths[2]


def test_g3_style_block_headers_are_accepted(tmp_path):
    """The other published shape: one header, several 'Pathway N' blocks."""
    p = tmp_path / "g3.txt"
    p.write_text(
        "Network Smiles: OCC(O)C=O\nGeneration: G1\n\n Pathway 0 \n"
        "OCC(O)C=O\n├── C=O\n└── OCC=O\n"
        "Energy Change: -1.0\nReaction IDs: ['R001']\n"
        " Pathway 1 \nOCC(O)C=O\n├── C=O\n└── OCC=O\n"
        "Energy Change: -2.0\nReaction IDs: ['R001']\n"
    )
    pws = parse_file(p)
    assert len(pws) == 2
    assert pws[0].meta["depth"] is None  # this shape states no length
    assert [pw.reported_dg for pw in pws] == [-1.0, -2.0]


def test_malformed_tree_line_is_named(tmp_path):
    p = tmp_path / "bad.txt"
    p.write_text("Network Smiles: C=O\n\n Pathway 0 \nC=O\n  ~~~ nonsense\n")
    with pytest.raises(SpecError, match="cannot parse tree line"):
        parse_file(p)


def test_reaction_id_count_mismatch_is_named(tmp_path):
    p = tmp_path / "bad.txt"
    p.write_text(
        "Network Smiles: OCC=O\n\n Pathway 0 \nOCC=O\n"
        "├── C=O\n└── C=O\nReaction IDs: ['a', 'b']\n"
    )
    with pytest.raises(SpecError, match="2 reaction ids but 1 internal nodes"):
        parse_file(p)

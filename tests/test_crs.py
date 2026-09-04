import pytest

from autocycle.io_crs import find_cycles, read_crs, to_cycle, to_graph
from autocycle.render import render
from autocycle.spec import SpecError

EX1 = "examples/crs/example-01.crs"
EX6 = "examples/crs/example-06.crs"


def test_parses_reactions_food_and_catalysts():
    s = read_crs(EX1)
    assert len(s.reactions) == 6
    r = s.reactions[0]
    assert r.name == "r1"
    assert r.reactants == ["a1", "b1"]
    assert r.products == ["c1"]
    assert r.catalysts == ["c3"]
    assert "a1" in s.food and "b1" in s.food
    assert "c1" not in s.food


def test_species_covers_every_name():
    s = read_crs(EX6)
    assert s.food == {"f"}
    assert {"x1", "x7"} <= s.species


def test_comments_and_blank_lines_ignored():
    s = read_crs(EX6)
    assert all(not r.name.startswith("#") for r in s.reactions)


def test_quoted_reaction_names_survive():
    assert "r'1" in {r.name for r in read_crs(EX1).reactions}


def test_the_loop_is_catalytic_not_mass_flow():
    """In `r: a+b [c] -> d` the loop runs through catalysis; flow finds nothing here."""
    s = read_crs(EX1)
    found = find_cycles(s, "catalysis")
    assert len(found) == 1
    assert set(found[0]) == {"c1", "c2", "c3"}   # rotation is not guaranteed
    assert find_cycles(s, "flow") == []


def test_example_06_has_several_catalytic_cycles():
    assert len(find_cycles(read_crs(EX6), "catalysis")) == 4


def test_unknown_graph_mode_rejected():
    with pytest.raises(ValueError, match="unknown mode"):
        to_graph(read_crs(EX1), "sideways")


def test_cycle_species_are_abstract_and_render_as_names():
    s = read_crs(EX1)
    c = to_cycle(s, find_cycles(s)[0])
    assert [m.structure for m in c.nodes] == [False] * 3
    svg = render(c, style="paper")
    assert ">c1<" in svg
    assert "<svg x='" not in svg  # no structure depictions were embedded


def test_food_only_appears_as_a_feeder():
    s = read_crs(EX1)
    c = to_cycle(s, find_cycles(s)[0])
    fed = {sp.smiles for st in c.steps for sp in st.consumes}
    assert fed <= s.food
    assert fed


def test_a_ring_that_is_not_a_cycle_is_named():
    s = read_crs(EX1)
    with pytest.raises(SpecError, match="not a cycle"):
        to_cycle(s, ["c1", "a1"])


def test_reverse_arrow_is_normalised(tmp_path):
    p = tmp_path / "r.crs"
    p.write_text("r1: y [c] <- x\nFood: x\n")
    r = read_crs(p).reactions[0]
    assert r.reactants == ["x"]
    assert r.products == ["y"]


def test_reversible_arrow_is_flagged(tmp_path):
    p = tmp_path / "r.crs"
    p.write_text("r1: x [c] <-> y\nFood: x\n")
    assert read_crs(p).reactions[0].reversible


def test_unparseable_line_is_named(tmp_path):
    p = tmp_path / "bad.crs"
    p.write_text("this is not a reaction\n")
    with pytest.raises(SpecError, match="cannot parse line"):
        read_crs(p)


def test_no_reactions_is_named(tmp_path):
    p = tmp_path / "empty.crs"
    p.write_text("# only a comment\nFood: a\n")
    with pytest.raises(SpecError, match="no reactions"):
        read_crs(p)

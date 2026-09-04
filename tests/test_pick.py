from autocycle.io_cypher import from_cypher_row, read_rows
from autocycle.io_spec import load_yaml
from autocycle.pick import distance, farthest_first, tokens


def _cycles():
    import contextlib

    from autocycle.spec import SpecError

    out = []
    for row in read_rows("examples/cypher_sample.csv"):
        with contextlib.suppress(SpecError):
            out.append(from_cypher_row(row))
    return out


def test_identical_cycles_are_distance_zero():
    c = load_yaml("examples/formose_gain.yaml")
    assert distance(tokens(c), tokens(c)) == 0.0


def test_disjoint_cycles_are_distance_one():
    assert distance(frozenset({"a"}), frozenset({"b"})) == 1.0
    assert distance(frozenset(), frozenset()) == 0.0


def test_tokens_capture_length_rules_and_sides():
    t = tokens(load_yaml("examples/formose_gain.yaml"))
    assert "len:4" in t
    assert any(x.startswith("rule:") for x in t)
    assert any(x.startswith("side:") for x in t)
    assert any(x.startswith("seed:") for x in t)
    assert "subs:1" in t


def test_picks_are_unique_and_respect_k():
    cycles = _cycles()
    items = [(tokens(c), c) for c in cycles]
    picked = farthest_first(items, 3)
    assert len(picked) == 3
    assert len({id(p) for p in picked}) == 3


def test_asking_for_more_than_available_returns_everything():
    items = [(tokens(c), c) for c in _cycles()]
    assert len(farthest_first(items, 99)) == len(items)


def test_second_pick_is_the_least_similar_candidate():
    a = (frozenset({"x", "y", "z"}), "rich")
    near = (frozenset({"x", "y"}), "near")
    far = (frozenset({"q"}), "far")
    assert farthest_first([a, near, far], 2) == ["rich", "far"]

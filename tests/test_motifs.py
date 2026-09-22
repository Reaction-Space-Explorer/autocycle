"""The mechanism level, which is the claim the paper's headline rests on.

A mechanism is meant to be substrate-agnostic and alias-proof: the same
transformation on a triose and on a hexose must read the same, and a rule renamed
must not create a second mechanism. Those two properties are what makes the
collapse from cores to mechanisms a measurement rather than a relabelling, so they
are tested directly rather than through an enumeration.
"""
from autocycle.cores.motifs import carbons, coarse_motif, coarse_step, delta, formula, motif

GLYC, TRIOSE, TETROSE = "OCC=O", "OCC(O)C=O", "OCC(O)C(O)C=O"
HCHO = "C=O"


def _rxn(reactants, products):
    d = {}
    for s in reactants:
        d[s] = d.get(s, 0) - 1
    for s in products:
        d[s] = d.get(s, 0) + 1
    return d


def test_formula_and_carbons_read_the_molecule():
    assert formula(GLYC) == "C2H4O2"
    assert carbons(GLYC) == 2 and carbons(TETROSE) == 4


def test_delta_is_the_signed_change_along_the_ring():
    assert delta(GLYC, TRIOSE) == "C+1H+2O+1"
    assert delta(GLYC, GLYC) == "0"


def test_an_aldol_reads_the_same_on_a_triose_and_on_a_tetrose():
    """The substrate-agnostic claim: same transformation, same mechanism step."""
    small = {"r": _rxn([GLYC, HCHO], [TRIOSE])}
    large = {"r": _rxn([TRIOSE, HCHO], [TETROSE])}
    assert coarse_step(small, "r", GLYC, TRIOSE) == coarse_step(large, "r", TRIOSE, TETROSE)


def test_the_formula_level_separates_what_the_mechanism_level_joins():
    """The two rungs are different, or the ladder would have nothing to collapse."""
    small = {"r": _rxn([GLYC, HCHO], [TRIOSE])}
    large = {"r": _rxn([TRIOSE, HCHO], [TETROSE])}
    assert motif(small, [GLYC, TRIOSE], ["r", "r"]) != motif(large, [TRIOSE, TETROSE], ["r", "r"])


def test_a_renamed_rule_does_not_make_a_second_mechanism():
    """The alias-proof claim: the label is built from the chemistry, not the name."""
    by_a = {"Aldol Condensation": _rxn([GLYC, HCHO], [TRIOSE])}
    by_b = {"Knoevenagel C": _rxn([GLYC, HCHO], [TRIOSE])}
    assert (coarse_motif(by_a, [GLYC, TRIOSE], ["Aldol Condensation"] * 2)
            == coarse_motif(by_b, [GLYC, TRIOSE], ["Knoevenagel C"] * 2))


def test_the_label_does_not_depend_on_where_the_ring_is_entered():
    """A cycle has no first step, so the mechanism must be order-free."""
    by = {"a": _rxn([GLYC, HCHO], [TRIOSE]), "b": _rxn([TRIOSE], [GLYC, HCHO])}
    assert (coarse_motif(by, [GLYC, TRIOSE], ["a", "b"])
            == coarse_motif(by, [TRIOSE, GLYC], ["b", "a"]))

"""The triage is a stated rule, so it can be checked rather than argued with."""
import pytest

from autocycle.cores.triage import verdict

CASES = [
    ("ordinary sugar chemistry is serious",
     ("Aldol Condensation", "Retro Aldol", "Keto-enol migration twice"), "serious"),
    ("a Cannizzaro alone is conditional",
     ("Aldol Condensation", "Cannizarro 2, HCHO (oxidation)"), "conditional"),
    ("methanol as donor with its own supply is the artefact",
     ("Michael Addition 0,2, ", "Cannizarro 2, HCHO (reduction)"), "artefact"),
    ("a reverse benzilic rearrangement is an artefact",
     ("Benzilic Acid Rearrangement (inverse)",), "artefact"),
]


@pytest.mark.parametrize("why,rules,want", CASES, ids=[c[0] for c in CASES])
def test_the_rule_gives_the_stated_verdict(why, rules, want):
    assert verdict(rules)[0] == want


def test_the_artefact_is_named_by_its_chemistry():
    _, why = verdict(("Michael Addition 0,2, ", "Cannizarro 2, HCHO (reduction)"))
    assert "methanol" in why


def test_the_rule_does_not_depend_on_the_order_of_the_steps():
    a = ("Michael Addition 0,2, ", "Cannizarro 2, HCHO (reduction)")
    assert verdict(a) == verdict(tuple(reversed(a)))

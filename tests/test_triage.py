"""The triage is a stated rule, so it can be checked rather than argued with."""
from autocycle.cores.triage import verdict


def test_ordinary_sugar_chemistry_is_serious():
    assert verdict(("Aldol Condensation", "Retro Aldol",
                    "Keto-enol migration twice"))[0] == "serious"


def test_a_cannizzaro_alone_is_conditional():
    assert verdict(("Aldol Condensation", "Cannizarro 2, HCHO (oxidation)"))[0] == "conditional"


def test_methanol_as_donor_with_its_own_supply_is_the_artefact():
    v, why = verdict(("Michael Addition 0,2, ", "Cannizarro 2, HCHO (reduction)"))
    assert v == "artefact" and "methanol" in why


def test_a_reverse_benzilic_rearrangement_is_an_artefact():
    assert verdict(("Benzilic Acid Rearrangement (inverse)",))[0] == "artefact"


def test_the_rule_does_not_depend_on_the_order_of_the_steps():
    a = ("Michael Addition 0,2, ", "Cannizarro 2, HCHO (reduction)")
    assert verdict(a) == verdict(tuple(reversed(a)))

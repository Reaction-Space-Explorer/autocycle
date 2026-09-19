"""Chemical triage of the mechanisms, by what each step needs to run.

The verdicts are a stated rule, not a hand list, so they can be re-derived and
argued with. A mechanism is judged by the transformations it contains:

serious      every step is ordinary aqueous sugar chemistry: aldol addition or
             retro-aldol between sugars or with formaldehyde, keto-enol
             isomerisation, hydration or dehydration, decarboxylation of an
             alpha-keto acid. These run at neutral pH without help.

conditional  a step needs a condition the network does not model but that can be
             argued prebiotically: a Cannizzaro, which needs strong base, or a
             conjugate (Michael) addition or its reverse, which competes with
             water as nucleophile, or fixation of CO2 onto an enol, which needs
             activation.

artefact     the amplification depends on chemistry that does not run in the
             direction required: methanol acting as the Michael donor while
             Cannizzaro supplies it, in an aqueous formaldehyde soup where water
             outcompetes methanol by orders of magnitude; or a benzilic acid
             rearrangement run backwards, which is effectively irreversible.
"""

from __future__ import annotations

CANNIZZARO = "Cannizarro"
MICHAEL = "Michael Addition"
BENZILIC = "Benzilic Acid Rearrangement (inverse)"
METHANOL_MAKER = "Cannizarro 2, HCHO (reduction)"


def verdict(rules: tuple[str, ...]) -> tuple[str, str]:
    joined = " ".join(rules)
    if BENZILIC in joined:
        return "artefact", "benzilic acid rearrangement run in reverse"
    if MICHAEL in joined and METHANOL_MAKER in joined:
        return "artefact", "methanol shuttle: methanol as Michael donor, Cannizzaro supplying it"
    if CANNIZZARO in joined and MICHAEL in joined:
        return "conditional", "Cannizzaro needs strong base; conjugate addition competes with water"
    if CANNIZZARO in joined:
        return "conditional", "Cannizzaro needs strong base"
    if MICHAEL in joined:
        return "conditional", "conjugate addition or its reverse competes with water"
    if "Decarboxylation" in joined and "C(=O)=O" in joined:
        return "conditional", "CO2 fixation onto an enol needs activation"
    return "serious", "ordinary aqueous sugar chemistry"

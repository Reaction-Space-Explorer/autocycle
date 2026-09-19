# Validating the core test against Blokhuis, Lacoste and Nghe

Source read directly: the paper (PNAS 2020, 117, 25230) and its SI, both in
~/Downloads. Their definition, verbatim:

> "The five types differ in their number of graph cycles and the way these cycles
> overlap. Type I consists of a single graph cycle that is weight asymmetric,
> defined as the product of the stoichiometric coefficients of its reaction
> products being different than that of its reactants... Types II and III comprise
> two distinct but overlapping graph cycles, Type IV comprises three, and Type V
> more than three."

> "every forward reaction of a core involves only one core species as a reactant"

## What was wrong

Our type classifier counted *forks*, reactions producing more than one core
species. That is not their definition, which counts graph cycles. It is now the
cycle count: nodes are core species, edges run from each reaction's core reactant
to each of its core products.

For three-species cores the two happen to agree, because one fork makes exactly
one extra cycle, so the reported numbers do not move: glucose 695 Type I and 347
Type II/III either way. At four species and beyond they would have diverged.

II and III are reported together. They have the same number of cycles and are
told apart by how those cycles overlap, which the SI defines in its ear
decomposition and which we do not yet implement.

## What passes

| check | result |
|---|---|
| their toy formose core | Type I, matching the paper |
| a single cycle that is weight symmetric | singular, rejected, as their definition requires |
| a constructed two-cycle core | Type II/III |
| every forward reaction consumes one core species | 0 violations in 1,042 glucose cores |
| every Type I core is weight asymmetric | 0 violations in 695 |
| their cavitand core, SI eq. S28 | core, two graph cycles, consistent with their Type III |

The one-reactant constraint holding with zero violations was not guaranteed by
our construction; it was measured.

## A second published type, from the text

SI eq. S28 gives the cavitand-amplification matrix in full, and the text names
the core: "we obtain an autocatalytic core of Type III consisting of the species
DCC, I1, DCU and Z". Encoded from that equation, our test returns a core with
det -1, a positive flux, and two graph cycles, which is the II/III bucket and so
consistent with their Type III.

The criterion therefore touches the literature at two published examples of two
types, both taken from equations rather than figures. That is how the autocycle
note did it too: its Abel spec is transcribed from the authors' own MOD flow
output rather than the printed figure, and its Blokhuis spec cites SI eq. S13.

## What is still not validated, and why

The SI's other worked examples are figures: Fig. S4a gives a Type II core for the
reverse Krebs cycle, Fig. S4b three Type I and four Type II cores for the Calvin
cycle, Fig. S5 a Type III for the cavitand amplification network. Their species
and reaction sets are not in the extractable text, so reproducing them needs a
person to transcribe the figures. Until that is done, the criterion touches the
literature at one published example, of one type, plus the constructed cases
above.

Rather than transcribe or OCR them, take the chemistry from a database. The
reverse Krebs and Calvin cycles are in KEGG and MetaCyc with exact
stoichiometry, and the SI states the answer to check against: Type II for rTCA,
three Type I and four Type II for Calvin. That validates against real chemistry
from a citable source, and it is the same move that produced the two validations
above. Optical extraction of a reaction scheme fails silently, which is the one
failure mode the elemental residual exists to catch.

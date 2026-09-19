# Inputs

Networks are Neo4j import tables from the CRNR generator: `rels_N.txt` is a
stoichiometric matrix in triplet form (reaction, species SMILES, signed
coefficient, rule), with multiplicity written as repeated rows. No MØD needed.

| network | archive | reactions | species | on this machine |
|---|---|---|---|---|
| Formose G6 | `Formose (6 rounds)-20260902T150858Z-1-001.zip` | 306,244 | 117,874 | yes, ~/Downloads/popvax |
| Glucose G5 | `Glucose (5 Rounds)-20260911T165210Z-1-001.zip` | 120,952 | 48,403 | yes, ~/Downloads |
| FormoseAmm G4 | — | 145,820 | 35,318 | energies only |
| GlucoseAmm G4 | — | 108,352 | — | energies only |
| Maillard G3 | — | 8,856 | 2,923 | no |
| PyruvicAcid G3 | — | — | — | energies only |

Energies: `Jim_NA/nucleoside-analogues/ProcessedData/SI/full/*_energies_pH*.csv`,
at pH 7.0, 7.4, 9.0 and 11.0. Columns are dG_prime_kJ_mol, sigma_kJ_mol,
estimable, status.

**Trap.** 55,721 of 306,244 Formose G6 rows (18.2%) carry a dG value but are
flagged `estimable=False` with sigma = 100000, the unbounded-variance sentinel,
and values such as +373 kJ/mol. Filtering on dG without checking sigma admits
them silently. Only 553 rows have no value at all. Report three columns:
survives the filter, fails it, has no usable estimate.

The 95% filter is `dG + 1.96*sigma < 0`; that reproduces the stored
`spontaneous_95` flag on 527/527 estimable rows of Formose G3, where both are
present. 1.645 reproduces 526/527, so the convention is inferred, not
documented. Confirm which was intended before it carries a result.

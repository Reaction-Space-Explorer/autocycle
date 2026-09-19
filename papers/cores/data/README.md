# Inputs

Networks are Neo4j import tables from the CRNR generator: `rels_N.txt` is a
stoichiometric matrix in triplet form (reaction, species SMILES, signed
coefficient, rule), with multiplicity written as repeated rows. Reading one does
not need MØD.

They are a separate deposition and are not carried here. Set `AUTOCAT_NETWORKS` to
a checkout of `Reaction-Space-Explorer/nucleoside-analogues`, where the reaction
tables are `OriginalData/OriginalNetworkData/Rels` and the free energies
`ProcessedData/SI/full`, at pH 7.0, 7.4, 9.0 and 11.0, with columns
dG_prime_kJ_mol, sigma_kJ_mol, estimable and status.

The five networks analysed here, counted as the loader parses them:

| network | reactions | species |
|---|---|---|
| Glucose G5 | 120,952 | 48,403 |
| Glucose + ammonia G4 | 108,334 | 26,126 |
| Formose G6 | 306,145 | 117,874 |
| Formose + ammonia G4 | 145,811 | 35,318 |
| Pyruvic acid G6 | 131,007 | 63,280 |

The same deposition carries HCN, Maillard and Urey-Miller tables, which are not
used here.

One trap in the energy tables. Of the 306,244 rows in the Formose G6 table at pH
7.4, 55,721 (18.2%) carry a dG value while flagged `estimable=False`, with sigma =
100000, the unbounded-variance sentinel, and values such as +373 kJ/mol. Filtering
on dG without also checking sigma admits them silently. Only 553 rows have no value
at all. Report three columns: survives the filter, fails it, has no usable estimate.

The 95% filter is `dG + 1.96*sigma < 0`, which reproduces the stored
`spontaneous_95` flag on 527 of 527 estimable rows of Formose G3, where both are
present. Using 1.645 reproduces 526, so the convention is inferred rather than
documented. Confirm which was intended before it carries a result.

# Synthetic traced-pathway sample

Hand-written example data in the format `autocycle.io_treelib` reads: a treelib-rendered
route tree per record, with a seed table (`products.tsv`, generation-0 molecules) and a
reaction table (`rels.tsv`, reactants at `-1` and products at `1`).

The chemistry is generic formose sugar chemistry — aldol condensation, keto-enol migration —
which has been textbook since Butlerow (1861). Reaction ids, energies and InChIKeys are
invented. Nothing here is derived from any particular network or study.

It is built to exercise the awkward parts of the format:

* **Reaction ids are deliberately out of tree order** (`['R002', 'R005', 'R004']`), because
  treelib sorts children alphabetically when printing while an id list is in insertion
  order. Reactions must be resolved by content, not position.
* **`Pathway Length` is the tree depth**, not the reaction count.
* **Record 3 does not bottom out**: one leaf is not a generation-0 seed, so it is `untraced`
  and the route is not `complete`.

# autocycle

Publication figures for reaction **cycles** and synthetic **routes**, with 2D structures
drawn inline at one constant bond length.

Two layouts behind one spec:

| | Layout | What the figure shows |
|---|---|---|
| `Cycle` | ring | the step where the autocatalyst returns, and the thermodynamic bottleneck |
| `Pathway` | layered, target on the right | where a route bottoms out in seeds, and where it merely stopped |

Cycle *finders* are well covered ([CatReNet](https://github.com/husonlab/catrenet), `autogato`,
MØD). `autocycle` is only the renderer, and is meant to sit behind any of them.

![cycle](examples/gallery/cycle_03.png)

## Install

```bash
pip install -e ".[dev]"
```

## Examples

A cycle from a YAML spec:

```bash
autocycle draw examples/formose_gain.yaml -o cycle.svg --drop O
```

A cycle found in a `source,target` edge list:

```bash
autocycle list examples/edges.csv                     # what cycles are in there
autocycle from-edges examples/edges.csv --cycle 0 --gain-at 3 -o cycle.svg
```

A synthetic route:

```bash
autocycle route examples/ribose_route.yaml -o route.svg
autocycle route-smiles \
  --reaction "C=O.C=O>>OCC=O" --reaction "OCC=O.OCC=O>>OCC(O)C(O)C=O" \
  --target "OCC(O)C(O)C=O" -o route.svg
```

![route](examples/routes/route_01.png)

From Python:

```python
from autocycle import Cycle, Mol, Side, Step, render

c = Cycle(
    nodes=[Mol("OCC=O"), Mol("OCC(O)C=O"), Mol("OCC(O)C(O)C=O")],
    steps=[
        Step("r1", "Aldol", dg=-16.4, consumes=[Side("C=O")]),
        Step("r2", "Aldol", dg=-13.4, consumes=[Side("C=O")]),
        Step("r3", "Retro-aldol", dg=-9.4, mag=2.0, gain=True, produces=[Side("OCC=O")]),
    ],
    seed=0,
)
open("cycle.svg", "w").write(render(c))
```

Whole corpora, with a report:

```bash
autocycle bench cycles_dir/ --sample 5                  # Cypher ring-query CSVs
autocycle bench-routes routes_dir/ --seeds products.tsv --rels rels.tsv --sample 5
```

## Styles

`paper` (default) reproduces *Chem. Sci.* 2022, **13**, 4838 Fig. 7: bare structures, ochre
reaction squares, dark-grey flow arrows, pale-grey side arrows, bare ΔG numbers.

`rich` turns the spare ink into data: arrow **width** = magnitude (`--mode linear|log|multiples`),
**hue** = ΔG, **concentric pair** = reversible, **green band** = the gain step.

## What it will not do

- A step failing a ΔG filter is **greyed, not deleted**.
- `total_dg` is `None`, never a partial sum, if any step lacks a ΔG.
- A route leaf is `seed`, `untraced`, or `unknown` — a leaf whose status the source never
  stated is never assumed to be a seed.
- Reactions in a treelib route file are matched by **content**, not by position in the
  `Reaction IDs` list: treelib sorts children when printing, so positions disagree.
  Unmatched reactions stay `unresolved` and are counted.
- `bench-routes` counts targets whose file has a header but no tree. No route found is a
  result, not something to skip.

## Credits

Arrow width encoding magnitude, with an explicit linear/log choice, and reversible steps as
concentric pairs, are from Catacycle — McFarlane, Henderson, Donnecke & McIndoe,
*Organometallics* 2019, **38**, 4051 ([code](https://github.com/brettrhenderson/Catacycle_Web), MIT).
The geometry here is an independent SVG implementation.

The cycle layout and thermodynamic annotation follow Arya, Ray, Sharma, Cruz Simbron,
Lozano, Smith, Andersen, Chen, Meringer & Cleaves, *An open source computational workflow
for the discovery of autocatalytic networks in abiotic reactions*, *Chem. Sci.* 2022, **13**, 4838.

MIT.

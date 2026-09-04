# autocycle

Publication figures for reaction **cycles** and synthetic **routes**, with 2D structures
drawn inline at one constant bond length.

Two layouts behind one spec:

| | Layout | What the figure shows |
|---|---|---|
| `Cycle` | ring | the step where the autocatalyst returns, and the thermodynamic bottleneck |
| `Pathway` | layered, target on the right | where a route bottoms out in seeds, and where it merely stopped |

![cycle](examples/gallery/cycle_03.png)

*A cycle found in a reaction network: gold ring molecules, magenta feeders and products,
reaction ids with ΔG, the autocatalyst ringed in green, rules listed underneath.*

## Install

Needs Python >= 3.10.

```bash
uv pip install git+https://github.com/Reaction-Space-Explorer/autocycle
```

From a checkout, or to run the tests:

```bash
uv venv --python 3.12 && uv pip install -e ".[dev]" && pytest -q
```

Or without installing anything:

```bash
uvx --python 3.12 --from git+https://github.com/Reaction-Space-Explorer/autocycle autocycle --help
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

## Choosing which cycle to show

A network search returns thousands of cycles, most of them not worth a figure. `list` and
`bench` report two things that decide it:

- **distinct feeder count** — a cycle fed by one molecule is a stronger result than one
  needing three. `from-edges --rank` orders by fewest feeders, then lightest feeder.
- **restricted molecules on the ring** — a cycle running *through* formaldehyde or methanol
  is usually not what you want to present: methanol in particular is hard to oxidise or
  reduce without a catalyst, so such a cycle is hard to interpret. Flagged, never dropped,
  and the set is configurable.

## Styles

`--style paper` (default) follows *Chem. Sci.* 2022, **13**, 4838 Fig. 7: bare structures,
ochre reaction squares, dark-grey flow arrows, pale-grey side arrows, bare ΔG numbers.

`--style annotated` is the network-search convention: gold ring molecules, magenta feeders,
reaction ids with energies, water balance annotated rather than drawn, and a rule list under
the figure. This is what the gallery above uses.

`--style rich` turns the spare ink into data: arrow **width** = magnitude
(`--mode linear|log|multiples`), **hue** = ΔG, **concentric pair** = reversible,
**green band** = the gain step.

Neither style is a pixel reproduction of a published figure: the layout and encodings match,
the structures are drawn by RDKit or obabel rather than by the original code.

### Structure depictions

RDKit by default, which needs nothing extra. `--backend obabel` uses Open Babel if it is on
your `PATH` (`brew install open-babel`, `apt install openbabel`), which draws small molecules
more legibly — formaldehyde as `O = CH₂` rather than a bare `=O` — and matches the published
figures more closely. Every depiction is drawn at one constant bond length either way, so a
hexose and a water molecule keep their relative size.

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

## Citing

The cycle layout, the feeder/consumer convention and the thermodynamic annotation follow:

> Arya, A.; Ray, J.; Sharma, S.; Cruz Simbron, R.; Lozano, A.; Smith, H. B.; Andersen, J. L.;
> Chen, H.; Meringer, M.; Cleaves, H. J. **An open source computational workflow for the
> discovery of autocatalytic networks in abiotic reactions.** *Chem. Sci.* **2022**, *13*,
> 4838–4853. [doi:10.1039/D2SC00256F](https://doi.org/10.1039/D2SC00256F)

> Cruz-Simbron, R.; Sharma, S.; Arya, A.; Ray, J.; Lozano, A.; Andersen, J. L.; Chen, H.;
> Cleaves, H. J. **Combined Network and High Resolution Mass Spectrometry Analysis of the
> Formose Reaction Reveals Mechanisms for Emergent Behaviors.** *ChemRxiv* **2024**.
> [doi:10.26434/chemrxiv-2024-nj0p6](https://doi.org/10.26434/chemrxiv-2024-nj0p6)

Arrow width encoding magnitude, with an explicit linear/log choice, and reversible steps as
concentric arrow pairs, are from Catacycle:

> McFarlane, J.; Henderson, B.; Donnecke, S.; McIndoe, J. S. **An Information-Rich Graphical
> Representation of Catalytic Cycles.** *Organometallics* **2019**, *38*, 4051–4053.
> [doi:10.1021/acs.organomet.9b00563](https://doi.org/10.1021/acs.organomet.9b00563)
> · [code](https://github.com/brettrhenderson/Catacycle_Web) (MIT)

The geometry here is an independent SVG implementation of that design.

MIT.

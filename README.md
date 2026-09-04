# autocycle

[![CI](https://github.com/Reaction-Space-Explorer/autocycle/actions/workflows/ci.yml/badge.svg)](https://github.com/Reaction-Space-Explorer/autocycle/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)

Publication figures for reaction **cycles** and synthetic **routes**, from a YAML spec or
straight out of a network search. Structures are drawn inline at one constant bond length.

<p align="center"><img src="examples/canonical/formose_core.png" width="560"></p>

The formose core cycle, from [examples/canonical/formose_core.yaml](examples/canonical/formose_core.yaml).
A disc marks a role: **gold** the autocatalyst, **magenta** the off-cycle feeders and
products, nothing the intermediates. One turn consumes one glycolaldehyde and two
formaldehyde and returns two glycolaldehyde, so the gain step is labelled and
`verify` calls it autocatalytic.

A whole search becomes one figure: frequency by cycle length stratified by feeder count,
composed with example cycles.

<p align="center"><img src="examples/figure_panel.png" width="880"></p>

## Install

Needs Python >= 3.10.

```bash
uv pip install git+https://github.com/Reaction-Space-Explorer/autocycle
```

From a checkout, or to run the tests:

```bash
uv venv --python 3.12 && uv pip install -e ".[dev]" && pytest -q
```

## Use

```bash
# a cycle or a route from YAML
autocycle draw examples/canonical/formose_core.yaml --style annotated -o cycle.svg
autocycle route examples/ribose_route.yaml -o route.pdf

# a cycle found in a source,target edge list
autocycle list examples/edges.csv
autocycle from-edges examples/edges.csv --cycle 0 --rank -o cycle.svg

# a whole corpus: report, then the multi-panel figure
autocycle bench cycles_dir/ --sample 5
autocycle panel cycles_dir/ --sample 3 -o figure.png
```

Output is `.svg`, or `.png` / `.pdf` with `pip install "autocycle[raster]"`.

```python
from autocycle import Cycle, Mol, Side, Step, render
from autocycle.verify import verify

c = Cycle(
    nodes=[Mol("OCC=O"), Mol("OCC(O)C=O"), Mol("OCC(O)C(O)C=O")],
    steps=[
        Step("r1", "Aldol", dg=-16.4, consumes=[Side("C=O")]),
        Step("r2", "Aldol", dg=-13.4, consumes=[Side("C=O")]),
        Step("r3", "Retro-aldol", dg=-9.4, gain=True, produces=[Side("OCC=O")]),
    ],
    seed=0,
    stoichiometry_complete=True,
)
verify(c).status          # 'autocatalytic'
open("cycle.svg", "w").write(render(c))
```

## Inputs and exports

| Read | |
|---|---|
| YAML spec | cycles and routes, the canonical form; a `shunt` block is supported |
| reaction SMILES | `from-smiles`, `route-smiles` |
| `source,target` CSV | `list`, `from-edges`, via networkx |
| Cypher ring-query CSV | `bench`, `panel` — ring, shunt, feeders, generations |
| treelib route files | `bench-routes`, with a seed table and reaction table |
| CatReNet `.crs` | `from-crs`; species are abstract names, so figures carry no chemistry |

| Write | |
|---|---|
| SVG, PNG, PDF | any command, via `-o` |
| SBGN-PD (SBGN-ML) | `autocycle sbgn spec.yaml` — lossy, see below |

## Styles

| Style | Cycle | Route |
|---|---|---|
| **`paper`** (default)<br>Bare structures, ochre reaction squares, bare ΔG. Follows *Chem. Sci.* 2022, **13**, 4838 Fig. 7. | <img src="examples/styles/cycle_paper.png" width="290"> | <img src="examples/styles/route_paper.png" width="290"> |
| **`annotated`**<br>Discs mark roles, reaction ids with energies inside the ring, water balance as text, rules listed underneath. | <img src="examples/styles/cycle_annotated.png" width="290"> | <img src="examples/styles/route_annotated.png" width="290"> |
| **`rich`**<br>Arrow **width** = magnitude (`--mode linear\|log\|multiples`), **hue** = ΔG, **concentric** = reversible, **green band** = gain. | <img src="examples/styles/cycle_rich.png" width="290"> | <img src="examples/styles/route_rich.png" width="290"> |

CoA thioesters and other R-groups work: a pseudo-atom such as the `[CoA]` stub used by
MØD parses directly, so SMILES can be pasted from a paper's supplementary data and the
stub is drawn as `CoA—S—`. See
[examples/canonical/malyl_coa_arm.yaml](examples/canonical/malyl_coa_arm.yaml).

No style is a pixel reproduction of a published figure: the layout and encodings match, the
structures are drawn by RDKit, or by Open Babel with `--backend obabel`, which renders small
molecules more legibly (`O = CH₂` rather than a bare `=O`).

## Does the cycle qualify?

Orgel's distinction, as stated by Peretó: a **simple** cycle regenerates one reactant
stoichiometrically, so the molecule produced replaces the one consumed; an
**autocatalytic** cycle "exhibits an additional yield of the feeder", n > 1. The
distinction is a stoichiometric coefficient, so `verify` reports conditions, not a boolean.

```python
verify(cycle).status    # autocatalytic | topological | simple | candidate | incomplete
verify(cycle).summary() # 'candidate (seed_identified=yes, ... extra_yield=unknown)'
```

- **`autocatalytic`** — an extra copy of the seed is stated. `examples/canonical/formose_core.yaml`.
- **`topological`** — a *shunt* bridges a ring molecule back to the seed. Side species on
  the shunt count toward the yield, which is usually where the extra copy is made:
  `examples/canonical/acetyl_coa_sol0.yaml` is a 5-molecule ring plus a 6-step shunt,
  11 steps in total, and comes out `autocatalytic` with a seed yield of 2. The published
  search treats that as the minimal criterion for n > 1, while noting it carries no flow
  constraint for mass balance. All 2100 cycles in the glucose corpus are this.
- **`simple`** — stoichiometry is stated and there is no extra copy, n = 1.
  `examples/canonical/krebs_tca.yaml`, which Peretó gives as a simple cycle.
- **`candidate`** — the structure holds but nothing settles the yield. Never asserted, never
  dismissed.

`verify` also flags a spec that declares a gain step its conditions cannot support.

## Choosing which cycle to show

A search returns thousands. `list` and `bench` report what decides it:

- **distinct feeder count** — one feeder is a stronger result than three. `--rank` orders by
  fewest feeders, then lightest, then earliest generation.
- **cycle centrality** — `select.cycle_centrality`, the lowest centrality on the ring, for
  when the feeder count cannot separate candidates.
- **restricted molecules on the ring** — methanol is hard to oxidise or reduce without a
  catalyst, so a cycle running *through* it is hard to interpret. Flagged, never dropped.

## What it will not do

- A step failing a ΔG filter is greyed, not deleted; `total_dg` is `None`, never a partial sum.
- A route leaf is `seed`, `untraced` or `unknown` — never assumed to be a seed.
- Reactions in a treelib file are matched by content, not by position in `Reaction IDs`:
  treelib sorts children when printing, so positions disagree. Unmatched stay `unresolved`.
- `bench-routes` counts targets whose file has a header but no tree. No route found is a result.
- SBGN-PD has no glyph for a shunt or for stoichiometric gain, so that export records them
  only as notes.

## Benchmark

| | cycles | routes |
|---|---|---|
| source | 31 Cypher result files | three treelib corpora |
| parsed | 2100 of 3100 rows | 8396 routes |
| rejected | 1000, pinched ring paths | 585 targets with no tree |
| reactions resolved by content | — | 54254, none unresolved |
| overlapping depictions | 0 | 0 |

The 1000 rejections visit a molecule twice, so `ringMols` has repeats and does not match
`ringRels`. They are named, not repaired. Re-run with `bench` and `bench-routes`;
`check.collisions` is the overlap gate the tests use, so layout regressions fail the suite.

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

The simple-versus-autocatalytic criterion is Orgel's, as stated in:

> Peretó, J. **Out of fuzzy chemistry: from prebiotic chemistry to metabolic networks.**
> *Chem. Soc. Rev.* **2012**, *41*, 5394. [doi:10.1039/C2CS35054H](https://doi.org/10.1039/C2CS35054H)

Arrow width encoding magnitude, with an explicit linear/log choice, and reversible steps as
concentric arrow pairs, are from Catacycle:

> McFarlane, J.; Henderson, B.; Donnecke, S.; McIndoe, J. S. **An Information-Rich Graphical
> Representation of Catalytic Cycles.** *Organometallics* **2019**, *38*, 4051–4053.
> [doi:10.1021/acs.organomet.9b00563](https://doi.org/10.1021/acs.organomet.9b00563)
> · [code](https://github.com/brettrhenderson/Catacycle_Web) (MIT)

MIT.

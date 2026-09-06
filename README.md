# autocycle

[![CI](https://github.com/Reaction-Space-Explorer/autocycle/actions/workflows/ci.yml/badge.svg)](https://github.com/Reaction-Space-Explorer/autocycle/actions/workflows/ci.yml)
[![License: MIT](https://img.shields.io/badge/License-MIT-yellow.svg)](https://opensource.org/licenses/MIT)
[![Python 3.10+](https://img.shields.io/badge/python-3.10%2B-blue.svg)](https://www.python.org/downloads/)

Figures for reaction **cycles** and synthetic **routes**, from a YAML spec or from a network
search. Structures are drawn inline, all at one bond length.

<p align="center"><img src="examples/canonical/acetyl_coa_sol0.png" width="620"></p>

The shortest autocatalytic acetyl-CoA cycle of
[Abel et al. 2026](https://doi.org/10.1038/s41540-025-00641-8), Fig. 3 Solution 0, from
[one spec](examples/canonical/acetyl_coa_sol0.yaml). Discs mark roles: gold the
autocatalyst, magenta off-cycle feeders and products, none for intermediates. Five
molecules on the inner ring, a six-step shunt on the outer arc with its own intermediates,
11 steps, reproducing that paper's objective value and its acetyl-CoA balance of one in and
two out.

## Simple or autocatalytic

| | |
|---|---|
| <img src="examples/canonical/formose_core.png" width="330"> | <img src="examples/canonical/krebs_tca.png" width="330"> |
| `autocatalytic`: the formose core returns two glycolaldehyde for one, n = 2, the overall reaction `HOCH2CHO + 2 H2CO -> 2 HOCH2CHO` as stated by Andersen et al. | `simple`: each turn of the TCA cycle regenerates one oxaloacetate, n = 1. |

Orgel's distinction is a stoichiometric coefficient. Andersen, Flamm, Merkle and Stadler
give the network form: a composite reaction is formally autocatalytic for x when it reads
`(A) + m x -> n x + (W)` with n > m. `verify` reports the conditions behind that rather
than a boolean:

```python
verify(cycle).status    # autocatalytic | topological | simple | candidate | incomplete
verify(cycle).summary() # 'candidate (seed_identified=yes, ... extra_yield=unknown)'
```

- `autocatalytic`: an extra copy of the seed is stated. Side species on a shunt count, which
  is usually where it is made.
- `topological`: a shunt bridges a ring molecule back to the seed, but no coefficient is
  recorded. The published search treats this as the minimal criterion for n > 1, noting it
  carries no flow constraint for mass balance. All 2100 cycles in the glucose corpus.
- `simple`: stoichiometry is stated and there is no extra copy, n = 1.
- `candidate`: the structure holds but nothing settles the yield. Never asserted, never
  dismissed.

`verify` also flags a spec that declares a gain step its conditions cannot support.

### Steady state

`verify` reads the coefficients a source states. `sna` reads the cycle as Clarke does, as a
stoichiometric matrix with a current `v`:

```
$ autocycle verify examples/canonical/acetyl_coa_sol0.yaml
status     autocatalytic
seed yield n = 2
overall    2 O=C(O)O  ->  1 *SC(C)=O
atoms      * +1, S +1, O -5, H -1
currents   2 (not a single extreme current)
```

- **intermediates balance**: a ring molecule the drawing never returns is named, exit 1.
- **atoms balance**: products minus reactants. Carbon is the load-bearing one; suppressing
  water and protons leaves a residual in H and O. This found a missing carboxylation in two
  of the specs here.
- **one extreme current, or a sum**: Sol 0 of Abel et al. is two, an allocatalytic ring plus
  a loop through the shunt that carries all the autocatalysis.

`flux` is that current; `mag` only sets arrow width.

## Whole searches

Frequency by cycle length, stratified by distinct feeder count, with example cycles.

<p align="center"><img src="examples/figure_panel.png" width="860"></p>

## Install

Python >= 3.10.

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

# check the stoichiometry without drawing anything
autocycle verify examples/canonical/formose_core.yaml

# a cycle found in a source,target edge list
autocycle list examples/edges.csv
autocycle from-edges examples/edges.csv --cycle 0 --rank -o cycle.svg

# a whole corpus: report, then the multi-panel figure
autocycle bench cycles_dir/ --sample 5
autocycle panel cycles_dir/ --sample 3 -o figure.png
```

Output is `.svg`, or `.png` / `.pdf` once `cairosvg` and libcairo are present.

The specs travel with an install, so the commands above run anywhere:

```bash
autocycle draw $(autocycle examples)/canonical/formose_core.yaml -o cycle.svg
```

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
| YAML spec | cycles and routes, the canonical form, including a `shunt` block |
| reaction SMILES | `from-smiles`, `route-smiles` |
| `source,target` CSV | `list`, `from-edges`, via networkx |
| Cypher ring-query CSV | `bench`, `panel`: ring, shunt, feeders, generations |
| treelib route files | `bench-routes`, with a seed table and reaction table |
| CatReNet `.crs` | `from-crs`. Species are abstract names, so figures carry no chemistry |

Unknown YAML keys are an error, so a typo cannot quietly change a verdict:

| | Keys |
|---|---|
| top level | `title`, `seed`, `stoichiometry_complete`, `nodes`, `steps`, `subcycles`, `shunt` |
| node | `smiles`, `label`, `generation`, `structure` |
| step | `id`, `rule`, `dg`, `mag`, `flux`, `rev_mag`, `consumes`, `produces`, `gain`, `filtered`, `note` |
| side species | `smiles`, `count`, `generation`, `structure` |
| `shunt` | `from_node`, `nodes`, `steps` |
| `subcycles` entry | `at_step`, `nodes`, `steps`, `label` |
| route node | `mol`, `from`, `reaction`, `terminal`, `generation` |

| Write | |
|---|---|
| SVG, PNG, PDF | any command, via `-o` |
| SBGN-PD (SBGN-ML) | `autocycle sbgn spec.yaml`. Lossy, see below |

## Styles

| Style | Cycle | Route |
|---|---|---|
| **`paper`** (default)<br>Bare structures, ochre reaction squares, bare ΔG. Follows *Chem. Sci.* 2022, **13**, 4838 Fig. 7. | <img src="examples/styles/cycle_paper.png" width="290"> | <img src="examples/styles/route_paper.png" width="290"> |
| **`annotated`**<br>Discs mark roles, reaction ids with energies inside the ring, water balance as text, rules listed underneath. | <img src="examples/styles/cycle_annotated.png" width="290"> | <img src="examples/styles/route_annotated.png" width="290"> |
| **`rich`**<br>Arrow width = magnitude (`--mode linear\|log\|multiples`), hue = ΔG, concentric = reversible, green band = gain. | <img src="examples/styles/cycle_rich.png" width="290"> | <img src="examples/styles/route_rich.png" width="290"> |

Structures come from RDKit, or from Open Babel with `--backend obabel`, which renders small
molecules more legibly (`O = CH₂` rather than a bare `=O`). Within a figure every depiction
uses one bond length: a requested length is honoured only while a molecule fits its canvas,
so the widest molecule sets the scale for all of them.

R-groups work. A pseudo-atom such as the `[CoA]` stub used by MØD parses directly, so SMILES
can be pasted from supplementary data, and the stub draws as `CoA-S-`.

No style is a pixel reproduction of a published figure. The layout and encodings match; the
structures are drawn by RDKit or Open Babel rather than by the original code.

## Choosing which cycle to show

A search returns thousands. `list` and `bench` report what decides it:

- **distinct feeder count**: one feeder is a stronger result than three. `--rank` orders by
  fewest feeders, then lightest, then earliest generation.
- **cycle centrality**: `select.cycle_centrality`, the lowest closeness centrality on the
  ring, which is Zubarev's definition, for when feeder count cannot separate candidates.
- **restricted molecules on the ring**: methanol is hard to oxidise or reduce without a
  catalyst, so a cycle running through it is hard to interpret. Flagged, never dropped.

Unedited search output, fused seven-membered rings included:

<p align="center"><img src="examples/gallery/cycle_03.png" width="520"></p>

## What it will not do

- A step failing a ΔG filter is greyed, not deleted. `total_dg` is `None`, never a partial sum.
- A route leaf is `seed`, `untraced` or `unknown`, never assumed to be a seed.
- Reactions in a treelib file are matched by content, not by position in `Reaction IDs`:
  treelib sorts children when printing, so positions disagree. Unmatched stay `unresolved`.
- `bench-routes` counts targets whose file has a header but no tree. No route found is a result.
- SBGN-PD has no glyph for a shunt or for stoichiometric gain, so that export records them
  as notes.
- `verify` reads a ring and its shunt as one current. A fused sub-cycle is a separate
  current, reported as a note rather than folded into the balance.

## Benchmark

| | cycles | routes |
|---|---|---|
| source | 31 Cypher result files | three treelib corpora |
| parsed | 2100 of 3100 rows | 8396 routes |
| rejected | 1000, pinched ring paths | 585 targets with no tree |
| reactions resolved by content | | 54254, none unresolved |
| overlapping depictions | 0 | 0 |

The 1000 rejections visit a molecule twice, so `ringMols` has repeats and does not match
`ringRels`. They are named, not repaired. `check.collisions` is the overlap gate the tests
use, so layout regressions fail the suite.

## Published cycles surveyed

`verify` applied to four families of published autocatalytic cycles, each found by a
different method, and a fifth quoted from its own paper:

| Source | Method | Cycles | Verdict |
|---|---|---|---|
| Arya et al. 2022, glucose corpus | Cypher motif query | 2100 | 2100 `topological`, 0 coefficient-confirmed |
| Abel et al. 2026, flow solutions | MØD with ILP | 15 | 15 `autocatalytic`, net gain 1 each |
| CatReNet examples | RAF sets | 498 | 0 coefficient-confirmed; the `.crs` format records no multiplicities |
| Blokhuis et al. 2020, toy formose | stoichiometric core enumeration | 1 transcribed | `autocatalytic`, n = 2, one extreme current, matching their Type I |
| Zubarev et al. 2015, rTCA supernetwork | combinatorial expansion | 1881 | 758 carry one branching point forming an autocatalytic loop. Their counts, not re-run here |

The difference is what each method records, not how good the cycles are. An ILP flow query
constrains the target's outflow to exceed its inflow, so its solutions state a coefficient:
`io_flow` reads the balance off the Overall Data tables, 11 steps for acetyl-CoA and 12 for
malate, matching that paper's Table 1. A motif query matches a shunt instead, which is
topological evidence without mass balance. A `.crs` file has nowhere to write a multiplicity
at all, which says nothing against RAF sets: Golnik et al. prove that under mild conditions
any RAF is stoichiometrically autocatalytic.

Cores separate two things that are easy to conflate. A core is defined by an invertible
stoichiometric matrix whose inverse is an elementary mode of autocatalysis, so that method
settles the criterion by construction, and it does so over carbon counts rather than
structures: what a RAF file lacks is multiplicities, not chemistry.
`examples/canonical/blokhuis_core.yaml` is their toy formose, and the single extreme current
`sna` finds there is their Type I reached by another route.

The last row is quoted rather than measured here. Zubarev, Rappoport and Aspuru-Guzik expand
a 175-molecule, 444-reaction supernetwork around the reverse TCA cycle and count branching
points forming an autocatalytic loop, the shunt under another name: 758 of 1881 carry one,
174 two, 20 three. Another chemistry and another search, the same topological signature at a
similar rate, equally without a stated coefficient.

Reproduce with `bench`, `bench-routes`, `io_flow.read_flow_summary` on a summary exported
by `pdftotext -layout`, and `from-crs` on a `.crs` file. Every figure here is built by
`make figures`; the two search figures need the full glucose corpus, which is not
redistributed, so they are `make gallery panel CORPUS=/path/to/csvs`.

## Citing

What each source is used for:

- **Cycle layout, feeder and consumer convention, ΔG annotation** — Arya et al., *An open
  source computational workflow for the discovery of autocatalytic networks in abiotic
  reactions*, Chem. Sci. **2022**, 13, 4838–4853.
  [doi](https://doi.org/10.1039/D2SC00256F) · Cruz-Simbron et al., ChemRxiv **2024**.
  [doi](https://doi.org/10.26434/chemrxiv-2024-nj0p6)
- **The criterion `(A) + m x -> n x + (W)` with n > m, their Definition 1** — Andersen,
  Flamm, Merkle, Stadler, J. Syst. Chem. **2020**, 8, 121–133.
  [arXiv](https://arxiv.org/abs/2107.03086)
- **Simple against autocatalytic, Orgel's distinction** — Peretó, Chem. Soc. Rev. **2012**,
  41, 5394. [doi](https://doi.org/10.1039/C2CS35054H)
- **`nu.v = 0` with `v >= 0`, and its cone of extreme currents** — Clarke, *Stoichiometric
  Network Analysis*, Cell Biophys. **1988**, 12, 237–253.
- **Listing those currents (`decompose`)** — Schmitz, Kolar-Anić, Anić, Čupić, *Stoichiometric
  Network Analysis and Associated Dimensionless Kinetic Equations*, J. Phys. Chem. A **2008**,
  112, 13452–13457. [doi](https://doi.org/10.1021/jp8056674)
- **Cycle centrality, and the rTCA counts in the survey** — Zubarev, Rappoport,
  Aspuru-Guzik, Sci. Rep. **2015**, 5, 8009. [doi](https://doi.org/10.1038/srep08009)
- **Autocatalytic cores, and `examples/canonical/blokhuis_core.yaml`** — Blokhuis, Lacoste,
  Nghe, PNAS **2020**, 117, 25230–25236. [doi](https://doi.org/10.1073/pnas.2013527117)
- **RAF sets against stoichiometric autocatalysis** — Golnik, Gatter, Hordijk, Stadler,
  Vassena, **2026**. [arXiv](https://arxiv.org/abs/2605.25523)
- **Arrow width as magnitude, reversible steps as concentric pairs** — Catacycle:
  McFarlane, Henderson, Donnecke, McIndoe, Organometallics **2019**, 38, 4051–4053.
  [doi](https://doi.org/10.1021/acs.organomet.9b00563) ·
  [code](https://github.com/brettrhenderson/Catacycle_Web)

Molecules are drawn by [RDKit](https://www.rdkit.org), or [Open Babel](https://openbabel.org)
with `--backend obabel`; cycles in an edge list are found with
[NetworkX](https://networkx.org).

Adjacent, and not dependencies: nesting one SVG inside another is
[skunk](https://github.com/whitead/skunk)'s trick, done there through matplotlib where this
emits SVG directly, which is what holds one bond length across a figure. Depiction itself is
solved by [SmilesDrawer](https://doi.org/10.1021/acs.jcim.7b00425) in the browser and by
[mols2grid](https://doi.org/10.5281/zenodo.6591473) over a grid of structures. Neither draws
a cycle.

MIT.

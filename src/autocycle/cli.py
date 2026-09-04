from __future__ import annotations

import argparse
import sys
from pathlib import Path

from autocycle import ingest
from autocycle.io_spec import (
    from_reaction_smiles,
    from_route_smiles,
    load_pathway_yaml,
    load_yaml,
)
from autocycle.render import render
from autocycle.spec import SpecError, drop_side
from autocycle.style import get

MODES = ("linear", "log", "multiples")


def _style(a):
    """Always a Style, so the collision check tests what is actually drawn."""
    return get(a.style, backend=a.backend) if getattr(a, "backend", None) else get(a.style)


def _write(cycle, a) -> None:
    if getattr(a, "drop", None):
        drop_side(cycle, a.drop)
    Path(a.out).write_text(render(cycle, mode=a.mode, style=_style(a), legend=a.legend))
    print(f"wrote {a.out}")


def main(argv: list[str] | None = None) -> int:
    ap = argparse.ArgumentParser(prog="autocycle", description="autocatalytic cycle figures")
    sub = ap.add_subparsers(dest="cmd", required=True)

    common = argparse.ArgumentParser(add_help=False)
    common.add_argument("-o", "--out", default="cycle.svg")
    common.add_argument("--mode", choices=MODES, default="linear", help="arrow width scaling")
    common.add_argument("--drop", action="append", default=[],
                        help="side species SMILES to omit, repeatable (e.g. --drop O)")
    common.add_argument("--style", choices=("paper", "annotated", "rich"), default="paper")
    common.add_argument("--backend", choices=("rdkit", "obabel"), default=None,
                        help="structure depiction engine (default rdkit)")
    common.add_argument("--legend", dest="legend", action="store_true", default=None)
    common.add_argument("--no-legend", dest="legend", action="store_false")

    d = sub.add_parser("draw", parents=[common], help="draw a YAML spec")
    d.add_argument("spec")

    s = sub.add_parser("from-smiles", parents=[common], help="draw from reaction SMILES")
    s.add_argument("--reaction", action="append", required=True)
    s.add_argument("--order", required=True, help="comma-separated ring molecule SMILES")
    s.add_argument("--title")

    rt = sub.add_parser("route", parents=[common], help="draw a synthetic route from a YAML spec")
    rt.add_argument("spec")

    rs = sub.add_parser(
        "route-smiles", parents=[common], help="trace a route back from a target"
    )
    rs.add_argument("--reaction", action="append", required=True)
    rs.add_argument("--target", required=True)
    rs.add_argument("--title")

    e = sub.add_parser("from-edges", parents=[common], help="draw a cycle found in an edge list")
    e.add_argument("csv")
    e.add_argument("--cycle", type=int, default=0, help="which cycle to draw (see 'list')")
    e.add_argument("--all", action="store_true", help="write every cycle, suffixed by index")
    e.add_argument("--min-len", type=int, default=3)
    e.add_argument("--max-len", type=int, default=12)
    e.add_argument("--cut", type=float, default=0.0, help="dg below this counts as spontaneous")
    e.add_argument("--gain-at", type=int, default=None, help="step index producing the gain")
    e.add_argument("--rank", action="store_true",
                   help="order cycles for presentation: fewest feeders, then lightest feeder")

    b = sub.add_parser("bench", help="run the renderer over a directory of Cypher result CSVs")
    b.add_argument("dir")
    b.add_argument("--limit", type=int, default=0, help="stop after N rows (0 = all)")
    b.add_argument("--sample", type=int, default=0, help="write this many example SVGs")
    b.add_argument("--sample-dir", default="bench_out")
    b.add_argument("--style", choices=("paper", "annotated", "rich"), default="paper")
    b.add_argument("--backend", choices=("rdkit", "obabel"), default=None)

    br = sub.add_parser(
        "bench-routes", help="run the renderer over a directory of treelib pathway files"
    )
    br.add_argument("dir")
    br.add_argument("--seeds", help="seed table (TSV: Generation, Smiles)")
    br.add_argument("--rels", action="append", default=[], help="reaction table TSV, repeatable")
    br.add_argument("--limit", type=int, default=0)
    br.add_argument("--sample", type=int, default=0)
    br.add_argument("--sample-dir", default="bench_routes_out")
    br.add_argument("--style", choices=("paper", "annotated", "rich"), default="paper")
    br.add_argument("--backend", choices=("rdkit", "obabel"), default=None)

    pn = sub.add_parser("panel", parents=[common],
                        help="one multi-panel figure: stratified histogram plus example cycles")
    pn.add_argument("dir")
    pn.add_argument("--sample", type=int, default=2, help="example cycles to include")
    pn.add_argument("--cols", type=int, default=2)
    pn.add_argument("--limit", type=int, default=0)

    ls = sub.add_parser("list", help="list cycles in an edge list, optionally to CSV")
    ls.add_argument("csv")
    ls.add_argument("--min-len", type=int, default=3)
    ls.add_argument("--max-len", type=int, default=12)
    ls.add_argument("--cut", type=float, default=0.0)
    ls.add_argument("--csv-out")

    a = ap.parse_args(argv)
    try:
        if a.cmd == "draw":
            _write(load_yaml(a.spec), a)
            return 0

        if a.cmd == "route":
            _write(load_pathway_yaml(a.spec), a)
            return 0

        if a.cmd == "route-smiles":
            _write(from_route_smiles(a.reaction, a.target, a.title), a)
            return 0

        if a.cmd == "from-smiles":
            order = [x.strip() for x in a.order.split(",") if x.strip()]
            _write(from_reaction_smiles(a.reaction, order, a.title), a)
            return 0

        if a.cmd == "bench":
            return _bench(a)

        if a.cmd == "bench-routes":
            return _bench_routes(a)

        if a.cmd == "panel":
            return _panel(a)

        g = ingest.spontaneous(ingest.read_edges(a.csv), a.cut)
        cycles = ingest.find_cycles(g, a.min_len, a.max_len)
        if not cycles:
            print(f"no cycles of length {a.min_len}-{a.max_len}", file=sys.stderr)
            return 1

        if a.cmd == "list":
            from autocycle.select import feeders, role_violations
            for i, c in enumerate(cycles):
                dgs = [
                    min(e.values(), key=lambda x: (x["dg"] is None, x["dg"] or 0.0))["dg"]
                    for e in (g.get_edge_data(c[j], c[(j + 1) % len(c)]) for j in range(len(c)))
                ]
                tot = "n/a" if any(d is None for d in dgs) else f"{sum(dgs):+.1f}"
                cyc = ingest.to_cycle(g, c)
                bad = role_violations(cyc)
                flag = f"  ring-role: {','.join(bad)}" if bad else ""
                print(
                    f"{i:4d}  len={len(c):2d}  feeders={len(feeders(cyc))}  "
                    f"total_dg={tot:>8}{flag}"
                )
            if a.csv_out:
                ingest.summary_csv(g, cycles, a.csv_out)
                print(f"wrote {a.csv_out}")
            return 0

        if getattr(a, "rank", False):
            from autocycle.select import rank_key
            cycles.sort(key=lambda ring: rank_key(ingest.to_cycle(g, ring)))
        targets = range(len(cycles)) if a.all else [a.cycle]
        for i in targets:
            if not 0 <= i < len(cycles):
                print(f"cycle {i} out of range (0-{len(cycles) - 1})", file=sys.stderr)
                return 1
            c = ingest.to_cycle(g, cycles[i], title=f"cycle {i}", gain_at=a.gain_at)
            if a.all:
                a.out = f"{Path(a.out).stem}_{i}{Path(a.out).suffix}"
            _write(c, a)
        return 0

    except SpecError as exc:
        print(f"error: {exc}", file=sys.stderr)
        return 2


def _bench_routes(a) -> int:
    from collections import Counter

    from autocycle.check import collisions
    from autocycle.io_treelib import parse_file, read_rels, read_seeds
    from autocycle.pick import farthest_first, tokens

    files = sorted(Path(a.dir).rglob("*.txt"))
    if not files:
        print(f"no .txt pathway files under {a.dir}", file=sys.stderr)
        return 1

    seeds = read_seeds(a.seeds) if a.seeds else None
    if seeds is not None and not seeds:
        # no G0 rows means unknown, not "nothing is a seed"
        print("  seed table has no generation-0 rows; treating seeds as unknown")
        seeds = None
    rels = read_rels(a.rels) if a.rels else None
    print(
        f"  {len(seeds) if seeds is not None else 'unknown'} seeds, "
        f"{len(rels) if rels else 0} reactions"
    )

    style = _style(a)
    routes = 0
    ok = 0
    failed: Counter[str] = Counter()
    depths: Counter[int] = Counter()
    leafstate: Counter[str] = Counter()
    resolved = unresolved = 0
    complete = 0
    collided = 0
    no_route = 0
    nets: Counter[str] = Counter()
    candidates: list = []

    print("loading tables")
    for path in files:
        try:
            pws = parse_file(path, seeds, rels)
        except SpecError as exc:
            failed[str(exc).split(":")[-1].strip()[:48]] += 1
            continue
        if not pws:
            # header but no tree: no route found. Counted, not skipped.
            no_route += 1
            continue
        for pw in pws:
            routes += 1
            if a.limit and routes > a.limit:
                break
            nets[path.parent.name] += 1
            depths[pw.root.depth] += 1
            for n in pw.leaves:
                leafstate[n.terminal] += 1
            for st in pw.steps:
                if st.note == "unresolved":
                    unresolved += 1
                else:
                    resolved += 1
            complete += pw.complete
            if collisions(pw, style):
                collided += 1
            try:
                svg = render(pw, style=style)
            except Exception as exc:  # noqa: BLE001 - a bench reports, never crashes
                failed[type(exc).__name__] += 1
                continue
            ok += 1
            if a.sample:
                candidates.append((tokens(pw), (pw, svg)))
        if a.limit and routes > a.limit:
            break

    print(f"files            {len(files)}")
    print(f"targets w/o route {no_route}  (header present, no pathway found)")
    print(f"routes           {routes}")
    print(f"rendered         {ok}")
    print(f"groups           {dict(nets.most_common(8))}")
    print(f"depths           {dict(sorted(depths.items()))}")
    print(f"reactions        {resolved} resolved by content, {unresolved} unresolved")
    print(f"leaf status      {dict(leafstate.most_common())}")
    print(f"complete routes  {complete} / {routes} (every leaf a known seed)")
    print(f"overlapping depictions  {collided} / {ok} figures")
    if failed:
        print("failures")
        for reason, n in failed.most_common(8):
            print(f"  {n:6d}  {reason}")

    if candidates:
        d = Path(a.sample_dir)
        d.mkdir(parents=True, exist_ok=True)
        picked = farthest_first(candidates, a.sample)
        print(f"most distinct {len(picked)} of {len(candidates)} routes:")
        for i, (pw, svg) in enumerate(picked):
            (d / f"route_{i:02d}.svg").write_text(svg)
            print(
                f"  route_{i:02d}  depth={pw.root.depth}  steps={len(pw.steps)}  "
                f"leaves={len(pw.leaves)}  complete={pw.complete}  {pw.meta.get('source', '')[:34]}"
            )
        print(f"wrote {len(picked)} SVGs to {d}/")
    return 0


def _panel(a) -> int:
    from collections import Counter

    from autocycle.io_cypher import from_cypher_row, read_rows
    from autocycle.panel import Series, bar_panel, compose, svg
    from autocycle.pick import farthest_first, tokens
    from autocycle.select import feeders

    files = sorted(Path(a.dir).rglob("*.csv"))
    total: Counter[int] = Counter()
    by_feeders: dict[int, Counter[int]] = {}
    candidates: list = []
    n = 0
    for path in files:
        batch = read_rows(path)
        if not batch or "ringMols" not in batch[0]:
            continue
        for row in batch:
            try:
                cycle = from_cypher_row(row)
            except SpecError:
                continue
            n += 1
            if a.limit and n > a.limit:
                break
            length = len(cycle.nodes)
            total[length] += 1
            by_feeders.setdefault(len(feeders(cycle)), Counter())[length] += 1
            candidates.append((tokens(cycle), cycle))
        if a.limit and n > a.limit:
            break

    if not total:
        print(f"no cycles found under {a.dir}", file=sys.stderr)
        return 1

    series = [Series("Total", dict(total))]
    for k in sorted(by_feeders, reverse=True):
        word = {1: "One", 2: "Two", 3: "Three"}.get(k, str(k))
        series.append(Series(f"{word} feeder" + ("" if k == 1 else "s"), dict(by_feeders[k])))

    cell_w, cell_h = 700.0, 620.0
    cells = [svg(bar_panel(series, cell_w, cell_h), cell_w, cell_h)]
    for cycle in farthest_first(candidates, a.sample):
        cells.append(render(cycle, mode=a.mode, style=_style(a)))

    Path(a.out).write_text(compose(cells, cell_w, cell_h, cols=a.cols))
    strata = {k: sum(v.values()) for k, v in sorted(by_feeders.items())}
    print(f"wrote {a.out}: {n} cycles, lengths {dict(sorted(total.items()))}, "
          f"{len(cells)} panels")
    print(f"  cycles by distinct feeder count: {strata}")
    if len(strata) == 1:
        print("  only one stratum: this source attaches a single feeder per cycle")
    return 0


def _bench(a) -> int:
    from collections import Counter

    from autocycle.check import collisions
    from autocycle.io_cypher import from_cypher_row, read_rows
    from autocycle.pick import farthest_first, tokens
    from autocycle.select import feeders, role_violations
    from autocycle.verify import verify

    files = sorted(Path(a.dir).rglob("*.csv"))
    if not files:
        print(f"no CSVs under {a.dir}", file=sys.stderr)
        return 1

    style = _style(a)
    rows = ok = skipped = 0
    rejected: Counter[str] = Counter()
    lengths: Counter[int] = Counter()
    render_fail: Counter[str] = Counter()
    collided = 0
    worst = []
    role_bad = 0
    feeder_hist: Counter[int] = Counter()
    verdicts: Counter[str] = Counter()
    unsupported = 0
    candidates: list = []

    for path in files:
        batch = read_rows(path)
        if not batch or "ringMols" not in batch[0]:
            skipped += 1
            continue
        for row in batch:
            rows += 1
            if a.limit and rows > a.limit:
                break
            try:
                cycle = from_cypher_row(row, title=path.stem)
            except SpecError as exc:
                rejected[str(exc).split(":")[0]] += 1
                continue
            lengths[len(cycle.nodes)] += 1
            feeder_hist[len(feeders(cycle))] += 1
            role_bad += bool(role_violations(cycle))
            v = verify(cycle)
            verdicts[v.status] += 1
            unsupported += v.disagrees_with_declaration
            hits = collisions(cycle, style)
            if hits:
                collided += 1
                worst.append((hits[0][2], path.name, hits[0][:2]))
            try:
                svg = render(cycle, style=style)
            except Exception as exc:  # noqa: BLE001 - a bench must report, not crash
                render_fail[type(exc).__name__] += 1
                continue
            ok += 1
            if a.sample:
                candidates.append((tokens(cycle), (cycle, svg)))
        if a.limit and rows > a.limit:
            break

    print(f"files            {len(files)} ({skipped} skipped: not Cypher results)")
    print(f"rows             {rows}")
    print(f"rendered         {ok}")
    print(f"rejected         {sum(rejected.values())}")
    for reason, n in rejected.most_common():
        print(f"  {n:6d}  {reason}")
    if render_fail:
        print("render failures")
        for name, n in render_fail.most_common():
            print(f"  {n:6d}  {name}")
    print(f"ring lengths     {dict(sorted(lengths.items()))}")
    print(f"distinct feeders {dict(sorted(feeder_hist.items()))}")
    print(f"restricted molecule on the ring  {role_bad} / {ok} cycles")
    print(f"autocatalysis     {dict(verdicts.most_common())}")
    if unsupported:
        print(f"  {unsupported} cycles declare a gain step the conditions do not support")
    print(f"overlapping depictions  {collided} / {ok} figures")
    for frac, name, pair in sorted(worst, reverse=True)[:5]:
        print(f"  {frac:.2f}  {name}  {pair[0][:28]} vs {pair[1][:28]}")

    if candidates:
        d = Path(a.sample_dir)
        d.mkdir(parents=True, exist_ok=True)
        picked = farthest_first(candidates, a.sample)
        print(f"most distinct {len(picked)} of {len(candidates)} cycles:")
        for i, (cycle, svg) in enumerate(picked):
            n_side = sum(len(s.consumes) + len(s.produces) for s in cycle.steps)
            name = f"cycle_{i:02d}"
            (d / f"{name}.svg").write_text(svg)
            rules = sorted({s.rule for s in cycle.steps if s.rule})
            print(
                f"  {name}  len={len(cycle.nodes)}  side={n_side}  "
                f"rules={len(rules)}  {'; '.join(rules)[:64]}"
            )
        print(f"wrote {len(picked)} SVGs to {d}/")
    return 0


if __name__ == "__main__":
    raise SystemExit(main())

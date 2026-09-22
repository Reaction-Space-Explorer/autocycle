"""The paper's figures, in the group's figures4papers house style."""
import re
from pathlib import Path

import matplotlib

matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from autocycle.cores import plot_style as PS

PS.apply()
plt.rcParams.update({"font.size": 9, "axes.titlesize": 10, "axes.labelsize": 9,
                     "xtick.labelsize": 8, "ytick.labelsize": 8,
                     "legend.fontsize": 8, "axes.titleweight": "bold",
                     "axes.linewidth": 1.1})
INK, MUTED, GRID = PS.INK, PS.MUTED, PS.GRID
BLUE, BLUE2 = PS.BLUE_MAIN, PS.BLUE_SECOND
GREEN, RED, NEUTRAL = PS.GREEN_3, PS.RED_STRONG, PS.NEUTRAL
TEAL, VIOLET = PS.TEAL, PS.VIOLET
# the four rungs, dark to light: each level holds fewer objects than the last
RUNG = [BLUE, BLUE2, TEAL, PS.GREEN_2]
# the three verdicts, fixed everywhere they appear
VERDICT = {"serious": PS.GREEN_3, "conditional": PS.BLUE_SECOND, "artefact": PS.RED_STRONG}
OUT = Path(__file__).resolve().parents[1] / "figures"

NETS = ["glucose\nG5", "glucose+NH$_3$\nG4", "formose\nG6", "formose+NH$_3$\nG4",
        "pyruvate\nG6"]


def _recorded():
    """Read the ladder and the triage out of results/paper_numbers.txt.

    They were held here as literals and went stale the first time the enumeration
    changed, which is the drift the manuscript audit exists to catch and cannot
    see from inside a figure. Blocks are in the order of NETS.
    """
    text = (OUT.parent / "results" / "paper_numbers.txt").read_text()
    ladder = {k: [] for k in ("cores", "distinct", "motifs", "mechanisms")}
    triage = {k: [] for k in ("serious", "conditional", "artefact")}
    for block in text.split("=== ")[1:]:
        c, d, m, x = (int(v) for v in re.search(
            r"ladder: (\d+) cores -> (\d+) distinct stoichiometry -> "
            r"(\d+) formula motifs -> (\d+) mechanisms", block).groups())
        for k, v in zip(ladder, (c, d, m, x), strict=True):
            ladder[k].append(v)
        got = dict(re.findall(r"^\s+(serious|conditional|artefact)\s+mechanisms"
                              r"\s+\d+\s+cores\s+(\d+)", block, re.M))
        for k in triage:
            triage[k].append(int(got.get(k, 0)))
    return ladder, triage


LADDER, TRIAGE = _recorded()


def fig1_ladder():
    fig, ax = plt.subplots(figsize=(6.6, 3.6))
    steps = list(LADDER)
    x = np.arange(len(NETS))
    w = 0.2
    for i, k in enumerate(steps):
        ax.bar(x + (i - 1.5) * w, LADDER[k], w, label=k, color=RUNG[i % len(RUNG)],
               edgecolor="white", linewidth=0.5)
    ax.set_yscale("log")
    ax.set_xticks(x); ax.set_xticklabels(NETS)
    ax.set_ylabel("count (log scale)")
    ax.set_title("Cores collapse to mechanisms by an order of magnitude or more")
    ax.legend(frameon=False, ncol=4, loc="upper center", bbox_to_anchor=(0.5, -0.13))
    ax.set_ylim(10, 4e4)
    for i, (c, m) in enumerate(zip(LADDER["cores"], LADDER["mechanisms"], strict=False)):
        ax.text(i, 1.9e4, f"{c/m:.0f}:1", ha="center", fontsize=8, color=MUTED)
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(OUT / "fig1_ladder.png", dpi=300)
    print("  fig1_ladder.png")


def fig4_triage_thermo():
    """Chemistry against thermodynamics, with the blind column resolved."""
    cross = {   # (verdict): [downhill, uphill, blind]  at pH 7.4
        "glucose G5":        {"serious": (165, 706, 83), "conditional": (2, 12, 54),
                              "artefact": (0, 2, 0)},
        "formose G6":        {"serious": (182, 1086, 149), "conditional": (43, 155, 25),
                              "artefact": (2, 2, 406)},
        "pyruvate G6":       {"serious": (53, 488, 116), "conditional": (33, 79, 17),
                              "artefact": (1, 1, 73)},
    }
    fig, axes = plt.subplots(1, 3, figsize=(8.4, 3.2), sharey=False)
    cols = [VERDICT["serious"], VERDICT["artefact"], NEUTRAL]
    labels = ["downhill", "uphill", "no usable estimate"]
    for ax, (net, d) in zip(axes, cross.items(), strict=True):
        ks = list(d)
        bottom = np.zeros(len(ks))
        for j, (c, lab) in enumerate(zip(cols, labels, strict=False)):
            vals = np.array([d[k][j] for k in ks], float)
            ax.bar(ks, vals, 0.6, bottom=bottom, color=c, label=lab,
                   edgecolor="white", linewidth=0.5)
            bottom += vals
        ax.set_title(net)
        ax.tick_params(axis="x", rotation=20)
        ax.spines[["top", "right"]].set_visible(False)
    axes[0].set_ylabel("cores")
    axes[1].legend(frameon=False, ncol=3, loc="upper center",
                   bbox_to_anchor=(0.5, -0.30))
    fig.suptitle("The thermodynamic filter cannot see the artefacts", y=1.0)
    fig.tight_layout(); fig.savefig(OUT / "fig4_triage_thermo.png", dpi=300,
                                    bbox_inches="tight")
    print("  fig4_triage_thermo.png")


def fig5_depth():
    """Density and collapse against depth. Points on fewer than 20 cores are
    hollow: a rate over 1 or 9 cores says nothing, and the early generations
    are where the counts are that small."""
    gens = {  # generation, cores/1000 rxn, cores/mechanism, core count
        "formose":  ([3, 4, 5, 6], [41.4, 23.3, 14.2, 6.7],
                     [12.0, 17.7, 40.7, 51.3], [24, 159, 774, 2050]),
        "glucose":  ([2, 3, 4, 5], [9.4, 4.0, 12.3, 8.5],
                     [2.0, 2.3, 13.1, 42.7], [2, 9, 236, 1024]),
        "pyruvate": ([3, 4, 5, 6], [4.2, 20.6, 9.5, 6.6],
                     [1.0, 6.8, 9.2, 21.0], [1, 41, 166, 861]),
    }
    fig, (a, b) = plt.subplots(1, 2, figsize=(7.6, 3.3))
    for i, (name, (g, dens, ratio, n)) in enumerate(gens.items()):
        solid = [k for k, c in enumerate(n) if c >= 20]
        for ax, y in ((a, dens), (b, ratio)):
            ax.plot(g, y, "-", color=RUNG[i % len(RUNG)], lw=1.4, alpha=0.5)
            ax.plot([g[k] for k in solid], [y[k] for k in solid], "o",
                    color=RUNG[i % len(RUNG)], ms=5, label=name if ax is a else None)
            hollow = [k for k in range(len(n)) if k not in solid]
            ax.plot([g[k] for k in hollow], [y[k] for k in hollow], "o",
                    mfc="white", mec=RUNG[i % len(RUNG)], ms=5, mew=1.2)
    a.set_xlabel("generation"); a.set_ylabel("cores per 1000 reactions")
    a.set_title("Past the peak, cores get rarer per reaction")
    b.set_xlabel("generation"); b.set_ylabel("cores per mechanism")
    b.set_title("and more of them are the same chemistry")
    a.plot([], [], "o", mfc="white", mec=MUTED, ms=5, label="fewer than 20 cores")
    for ax in (a, b):
        ax.spines[["top", "right"]].set_visible(False)
        ax.set_xticks([2, 3, 4, 5, 6])
    a.legend(frameon=False, fontsize=8)
    fig.tight_layout(); fig.savefig(OUT / "fig5_depth.png", dpi=300)
    print("  fig5_depth.png")




def fig6_rule_removal():
    """What each verdict class loses when one transformation is deleted."""
    rules = ["keto-enol\nmigration", "Cannizzaro\n(HCHO red.)", "Michael\naddition",
             "retro-aldol", "aldol\ncondensation"]
    base = {"serious": 1417, "conditional": 223, "artefact": 410}
    after = {"serious": [1407, 1417, 1417, 216, 813],
             "conditional": [223, 213, 166, 56, 108],
             "artefact": [410, 7, 7, 403, 410]}
    fig, ax = plt.subplots(figsize=(7.0, 3.5))
    x = np.arange(len(rules)); w = 0.26
    for i, (k, c) in enumerate(zip(["serious", "conditional", "artefact"],
                                   [VERDICT["serious"], VERDICT["conditional"], VERDICT["artefact"]], strict=False)):
        frac = [100 * v / base[k] for v in after[k]]
        ax.bar(x + (i - 1) * w, frac, w, color=c, label=k,
               edgecolor="white", linewidth=0.5)
    ax.axhline(100, color=MUTED, lw=0.8, ls=":")
    ax.set_xticks(x); ax.set_xticklabels(rules)
    ax.set_ylabel("percent of class surviving")
    ax.set_ylim(0, 118)
    ax.set_title("Deleting the rule that makes methanol removes the artefacts\n"
                 "and leaves the serious cores untouched")
    ax.legend(frameon=False, ncol=3, loc="upper center", bbox_to_anchor=(0.5, -0.17))
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(OUT / "fig6_rule_removal.png", dpi=300, bbox_inches="tight")
    print("  fig6_rule_removal.png")




def fig3_paired():
    """The artefact beside a serious core, drawn identically.

    Both panels come from autocycle with the same style flag, the same layout and
    the same role colours. autocycle sizes its canvas to the cycle, so the two
    come back at different heights; they are padded to a common canvas here
    rather than rescaled, which would change the bond lengths and break the one
    thing the pairing is for. Nothing in the rendering differs but the chemistry.
    """
    import matplotlib.image as mpimg
    panels = [(OUT / "fig3_serious_retroaldol.png", "a",
               "serious: retro-aldol cleavage"),
              (OUT / "fig3_artefact_methanol_shuttle.png", "b",
               "artefact: methanol shuttle")]
    def crop(i, pad=8):
        """Drop the white border. Cropping does not rescale, so bond lengths hold."""
        ink = (i[:, :, :3] < 0.99).any(axis=2)
        rows, cols = np.where(ink.any(axis=1))[0], np.where(ink.any(axis=0))[0]
        r0, r1 = max(rows[0] - pad, 0), min(rows[-1] + pad + 1, i.shape[0])
        c0, c1 = max(cols[0] - pad, 0), min(cols[-1] + pad + 1, i.shape[1])
        return i[r0:r1, c0:c1]

    imgs = [crop(mpimg.imread(p)) for p, _, _ in panels]
    h = max(i.shape[0] for i in imgs)
    w = max(i.shape[1] for i in imgs)
    padded = []
    for i in imgs:
        canvas = np.ones((h, w, i.shape[2]), dtype=i.dtype)
        top = (h - i.shape[0]) // 2
        left = (w - i.shape[1]) // 2
        canvas[top:top + i.shape[0], left:left + i.shape[1]] = i
        padded.append(canvas)
    fig, axes = plt.subplots(1, 2, figsize=(7.6, 7.6 * h / (2 * w) + 0.35))
    for ax, img, (_, tag, title) in zip(axes, padded, panels, strict=False):
        ax.imshow(img)
        ax.axis("off")
        # leading space keeps the title clear of the panel letter
        ax.set_title("    " + title, fontsize=9, loc="left")
        ax.text(0.0, 1.0, tag, transform=ax.transAxes, fontsize=11,
                fontweight="bold", va="bottom", ha="left")
    fig.subplots_adjust(wspace=0.06)
    fig.savefig(OUT / "fig3_paired.png", dpi=300, bbox_inches="tight")
    print("  fig3_paired.png")


def fig_bound():
    """What the free-energy filter can and cannot decide, per network.

    Table 5 as a picture: the resolved spontaneous fraction, the interval the
    unresolved cores leave it in, and the half-way mark the conclusion turns on.
    """
    import re
    rows = []
    for block in (OUT.parent / "results" / "paper_numbers.txt").read_text().split("=== ")[1:]:
        name = block.split(" ===")[0]
        th = {(k, st): int(v) for k, st, v in re.findall(
            r"^\s+(serious|conditional|artefact)\s+(spontaneous|not|no estimate)\s+(\d+)",
            block, re.M)}
        sp, no = th.get(("serious", "spontaneous"), 0), th.get(("serious", "not"), 0)
        bl = th.get(("serious", "no estimate"), 0)
        rows.append((name.replace("+", " + "), sp, bl, sp + no + bl))

    fig, ax = plt.subplots(figsize=(6.6, 2.9))
    y = np.arange(len(rows))[::-1]
    for i, (_, sp, bl, tot) in zip(y, rows, strict=True):
        lo, hi = 100 * sp / tot, 100 * (sp + bl) / tot
        decided = hi < 50
        ax.plot([lo, hi], [i, i], lw=7, solid_capstyle="butt",
                color=PS.GREEN_2 if decided else PS.NEUTRAL, zorder=2)
        ax.plot([lo], [i], "o", ms=7, color=BLUE, zorder=3)
        ax.text(hi + 1.5, i, f"{lo:.0f} to {hi:.0f}%", va="center", fontsize=8,
                color=INK if decided else MUTED)
    ax.axvline(50, color=PS.RED_STRONG, lw=1.2, ls="--", zorder=1)
    ax.text(51.5, len(rows) - 0.42, "half", ha="left", fontsize=8, color=PS.RED_STRONG)
    ax.set_yticks(y); ax.set_yticklabels([r[0] for r in rows], fontsize=9)
    ax.set_xlabel("serious cores that are spontaneous (%)")
    ax.set_xlim(0, 104); ax.set_ylim(-0.6, len(rows) - 0.05)
    ax.set_title("Two networks the filter cannot decide either way")
    ax.spines[["top", "right", "left"]].set_visible(False)
    ax.tick_params(axis="y", length=0)
    fig.tight_layout(); fig.savefig(OUT / "fig_bound.png", dpi=300)
    print("  fig_bound.png")


def fig_coresize():
    """The ladder against the number of species admitted to a core."""
    import re
    txt = (OUT.parent / "results" / "ladder.txt").read_text()
    got = {}
    for m in re.finditer(r"(\w+)Rels_\d+\s+n=(\d+)\s+cores (\d+)\s+distinct (\d+)"
                         r"\s+motifs (\d+)\s+mechanisms (\d+)", txt.replace("\n", " ")):
        got[(m.group(1), int(m.group(2)))] = tuple(int(g) for g in m.groups()[2:])
    nets = ["Glucose", "Formose", "PyruvicAcid"]
    label = {"Glucose": "glucose G5", "Formose": "formose G6", "PyruvicAcid": "pyruvate G6"}
    fig, ax = plt.subplots(figsize=(6.0, 3.3))
    for net, colour, mark in zip(nets, (BLUE, TEAL, PS.RED_STRONG), ("o", "s", "^"), strict=True):
        ns = [n for n in (3, 4, 5) if (net, n) in got]
        ratio = [got[(net, n)][0] / got[(net, n)][3] for n in ns]
        ax.plot(ns, ratio, marker=mark, color=colour, lw=2, ms=6, label=label[net])
        ax.annotate(f"{ratio[-1]:.0f}:1", (ns[-1], ratio[-1]), textcoords="offset points",
                    xytext=(6, -2), fontsize=8, color=colour)
    ax.set_yscale("log")
    ax.set_xticks([3, 4, 5]); ax.set_xlim(2.85, 5.4)
    ax.set_xlabel("species admitted to a core")
    ax.set_ylabel("cores per mechanism (log)")
    ax.set_title("The collapse widens with core size")
    ax.legend(frameon=False, fontsize=8, loc="upper left")
    ax.spines[["top", "right"]].set_visible(False)
    fig.tight_layout(); fig.savefig(OUT / "fig_coresize.png", dpi=300)
    print("  fig_coresize.png")


if __name__ == "__main__":
    fig1_ladder(); fig3_paired(); fig4_triage_thermo(); fig5_depth(); fig6_rule_removal()
    fig_bound(); fig_coresize()

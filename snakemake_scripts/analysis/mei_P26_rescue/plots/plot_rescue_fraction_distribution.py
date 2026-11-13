#!/usr/bin/env python
"""Figure 5: how much of the defect wMel reverses, as a distribution.

Plots results/rescue/summary/rescue_genes.tsv against the permutation null in
results/rescue/rescue_null/rescue_null_permutations.tsv.

An alternative to fig5_rescue_categories, which shows the same quantity as
stacked category counts and a bar of the median. Three things are wrong with
that form, and this one fixes all three.

  The categories are thresholds, not measurements. "full" is rescue_fraction
  >= 0.75 and no longer significant, "partial" >= 0.25, "exacerbated" <=
  -0.25. Counting genes into those bins throws away the shape of the
  distribution and makes the result look like a property of the data when it
  is partly a property of where the lines were drawn. Here the thresholds
  become annotations on a continuous axis, and the reader can see how much
  would move if they were set elsewhere.

  A bar of the median has no uncertainty on it. rescue_fraction is a ratio of
  two fitted effects, taken over genes selected for significance in its own
  denominator, so regression to the mean pushes it positive before any
  biology. A median plotted alone invites exactly the over-reading this needs
  to avoid.

  The null is the whole point and the old panel does not show it. Shuffling
  the labels produces median rescue fractions of 0.10 to 0.51 in this dataset,
  and four of nine cell types sit inside that range. Those 16 balanced
  relabellings are drawn here as points on each row, so whether the observed
  median clears them is visible rather than asserted in the caption.

  A  per-gene rescue fraction per cell type, as a violin, with the observed
     median and all 16 null medians on the same row. Filled when the observed
     median exceeds every null draw.
  B  the denominator: how many genes were DE in the mutant at all. A high
     rescue fraction over 30 genes is not the same result as one over 3,000.

Usage:
  python plot_rescue_fraction_distribution.py \
      --genes results/rescue/summary/rescue_genes.tsv \
      --null results/rescue/rescue_null/rescue_null_permutations.tsv \
      --outdir results/rescue/figures
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (empty_figure, COL, LOW_CONFIDENCE, WIDE, apply_style,
                      ct_short, note, panel_letter, save)
import matplotlib.pyplot as plt

NON_OVARIAN = ("adult tracheal cell", "oviduct")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--genes", required=True, help="summary/rescue_genes.tsv")
    p.add_argument("--null", default=None,
                   help="rescue_null/rescue_null_permutations.tsv; without it "
                        "the panel loses the only thing that makes the "
                        "rescue fraction interpretable, so pass it")
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig5_rescue_fraction")
    p.add_argument("--clip", type=float, default=2.0,
                   help="rescue fraction axis limit. The tail runs to about "
                        "+-4.6 because the ratio explodes on small "
                        "denominators, and those few genes would set the axis "
                        "for everything else")
    p.add_argument("--min_de", type=int, default=50,
                   help="cell types below this many mutant-DE genes are drawn "
                        "but marked, since a fraction over a handful of genes "
                        "is not an effect size")
    p.add_argument("--drop_labels",
                   default=",".join(tuple(LOW_CONFIDENCE) + NON_OVARIAN))
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def build(g, nl, args, suffix, note_txt):
    n_de = g.groupby("tissue").size()
    med = g.groupby("tissue")["rescue_fraction"].median()
    order = med.sort_values().index.tolist()
    y = np.arange(len(order))
    print(f"\n[{suffix}] {len(order)} cell types; {note_txt}")

    fig = plt.figure(figsize=(WIDE, max(2.6, 0.42 * len(order) + 1.3)))
    gs = fig.add_gridspec(1, 2, wspace=0.14, width_ratios=[2.9, 1.0])

    # ---- A: the distribution, the observed median, and the null ----
    ax = fig.add_subplot(gs[0, 0])
    ax.axvline(0, color=COL["black"], lw=0.8, zorder=2)
    ax.axvline(1, color=COL["grey"], lw=0.8, ls=":", zorder=2)
    ax.text(1, len(order) - 0.32, " full rescue", fontsize=5.5,
            color=COL["grey"], va="top")
    ax.text(0, len(order) - 0.32, " no change", fontsize=5.5,
            color=COL["grey"], va="top")

    data = [np.clip(g.loc[g["tissue"] == t, "rescue_fraction"],
                    -args.clip, args.clip) for t in order]
    # vert= was deprecated for orientation= in matplotlib 3.11. Support both,
    # because the cluster env is not necessarily on the same version as here.
    import matplotlib
    horiz = ({"orientation": "horizontal"}
             if tuple(int(x) for x in matplotlib.__version__.split(".")[:2]) >= (3, 11)
             else {"vert": False})
    vp = ax.violinplot(data, positions=y, widths=0.82, showextrema=False,
                       showmedians=False, **horiz)
    for b in vp["bodies"]:
        b.set_facecolor(COL["lightgrey"])
        b.set_edgecolor(COL["grey"])
        b.set_linewidth(0.5)
        b.set_alpha(0.85)

    beats = {}
    for i, t in enumerate(order):
        obs = med[t]
        if nl is not None and t in nl.index:
            draws = np.asarray(nl.loc[[t], "median_rf"], float)
            ax.scatter(draws, np.full(len(draws), i), s=7, marker="|",
                       color=COL["vermillion"], zorder=4, linewidths=0.9)
            beats[t] = bool(obs > np.nanmax(draws))
        else:
            beats[t] = False
        ax.scatter([obs], [i], s=40, marker="D", zorder=5,
                   color=COL["green"] if beats[t] else "white",
                   edgecolors=COL["black"], linewidths=0.6)

    ax.set_yticks(y, [ct_short([t])[0] + ("" if n_de[t] >= args.min_de else " *")
                      for t in order], fontsize=6)
    ax.set_xlim(-args.clip - 0.05, args.clip + 0.05)
    ax.set_ylim(-0.7, len(order) - 0.25)
    ax.set_xlabel("per-gene rescue fraction")
    ax.set_title("Effect size against the shuffled null", pad=4,
                 fontsize=7.5)
    cap = ("violin: all mutant-DE genes, clipped to "
           f"±{args.clip:g}; diamond: observed median, filled when it "
           "exceeds every null draw")
    if nl is not None:
        cap += "\nred ticks: the 16 balanced relabellings"
    if (n_de < args.min_de).any():
        cap += f"\n* fewer than {args.min_de} mutant-DE genes"
    note(ax, cap, y=-0.13)
    panel_letter(ax, "A", dx=-0.26)

    # ---- B: the denominator ----
    ax = fig.add_subplot(gs[0, 1])
    ax.barh(y, [n_de[t] for t in order], height=0.68, zorder=3,
            color=[COL["blue"] if n_de[t] >= args.min_de else COL["lightgrey"]
                   for t in order],
            edgecolor=COL["black"], linewidth=0.4)
    ax.set_yticks(y, [])
    ax.set_ylim(-0.7, len(order) - 0.25)
    ax.set_xlabel("mutant-DE genes")
    ax.set_title("Genes to rescue", pad=4, fontsize=7.5)
    panel_letter(ax, "B", dx=-0.12)

    save(fig, args.outdir, f"{args.stem}_{suffix}")

    out = pd.DataFrame({"n_mutant_DE": n_de.reindex(order),
                        "median_rf": med.reindex(order).round(3),
                        "frac_below_0": [float((g.loc[g["tissue"] == t,
                                                      "rescue_fraction"] < 0).mean())
                                         for t in order],
                        "beats_null": [beats[t] for t in order]})
    out["frac_below_0"] = out["frac_below_0"].round(3)
    print(out.to_string())
    won = [t for t in order if beats[t]]
    if nl is not None:
        print(f"\n{len(won)} of {len(order)} cell types have an observed median "
              f"rescue fraction above every one of the 16 relabellings: "
              f"{', '.join(ct_short(won)) or 'none'}. The rest are inside the "
              "range shuffling produces and must not be reported as rescue.")


def main():
    args = parse_args()
    apply_style()
    g = pd.read_csv(args.genes, sep="\t")
    g["rescue_fraction"] = pd.to_numeric(g["rescue_fraction"], errors="coerce")
    g = g.dropna(subset=["rescue_fraction"])
    drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}
    g = g[~g["tissue"].isin(drop)]
    print(f"{g['tissue'].nunique()} cell types, {len(g):,} mutant-DE gene by "
          f"cell-type pairs; dropped {sorted(drop) or 'none'}")

    nl = None
    if args.null and os.path.exists(args.null):
        nl = pd.read_csv(args.null, sep="\t")
        nl["median_rf"] = pd.to_numeric(nl["median_rf"], errors="coerce")
        nl = nl[~nl["tissue"].isin(drop)].set_index("tissue")
        print(f"null: {len(nl)} permutation rows over "
              f"{nl.index.nunique()} cell types")
    else:
        print("No --null given. The panel will show the distribution without "
              "the comparison that makes it interpretable; pass the "
              "permutations file.")

    n_de = g.groupby("tissue").size()
    for v in (["robust", "all"] if args.version == "both" else [args.version]):
        if v == "robust":
            keep = set(n_de[n_de >= args.min_de].index)
            sub = g[g["tissue"].isin(keep)]
            txt = f"cell types with at least {args.min_de} mutant-DE genes"
        else:
            sub, txt = g, "every cell type"
        if not sub["tissue"].nunique():
            empty_figure(args.outdir, f"{args.stem}_{v}",
                         f"No cell type survives the {v} filter.\n"
                         "An empty set here is a result, not a failure.")
            continue
        build(sub, nl, args, v, txt)


if __name__ == "__main__":
    main()

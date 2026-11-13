#!/usr/bin/env python
"""Figure 5: how much of the mei-P26 defect wMel reverses, as diverging bars.

Plots results/rescue/summary/rescue_category_counts.tsv beside the
permutation null in results/rescue/rescue_null/rescue_null_permutations.tsv.

An alternative to fig5_rescue_categories. The categories are the same; what
changes is that the bars diverge from zero and the null sits next to them.

Why diverging. In a left-aligned stacked bar the exacerbated class starts at a
different x in every row, so the one comparison a reader most wants, how much
goes the wrong way against how much goes the right way, cannot be made by eye.
Anchoring the scale at zero puts every exacerbated segment on the same
baseline and makes the balance readable down the column. The neutral class
straddles zero, half to each side, which is the usual convention for an
ordered scale with a middle category: exacerbated < none < partial < full.

Why the null panel. rescue_fraction is a ratio of two fitted effects taken
over genes selected for significance in its own denominator, so regression to
the mean pushes it positive before any biology. Shuffled labels give median
rescue fractions of 0.10 to 0.51 in this dataset. Panel B draws all 16
balanced relabellings against the observed median on an axis scaled to where
those numbers actually are, so a cell type inside the null range is visible
rather than asserted in a caption.

The denominator is folded into the y labels rather than given a panel. A
rescue fraction over 287 genes is a different object from one over 3,526, and
the number has to be in front of the reader, but it does not need a third of
the figure's width.

  A  category composition per cell type, diverging from zero
  B  observed median rescue fraction against the 16 relabellings

Usage:
  python plot_rescue_diverging.py \
      --input results/rescue/summary/rescue_category_counts.tsv \
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

CATS = {"full": COL["green"], "partial": COL["skyblue"],
        "none": COL["lightgrey"], "exacerbated": COL["vermillion"]}
NON_OVARIAN = ("adult tracheal cell", "oviduct")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True,
                   help="summary/rescue_category_counts.tsv")
    p.add_argument("--null", default=None,
                   help="rescue_null/rescue_null_permutations.tsv; without it "
                        "panel B is dropped and the figure loses the only "
                        "thing that makes a rescue fraction interpretable")
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig5_rescue_diverging")
    p.add_argument("--min_de", type=int, default=50,
                   help="cell types below this many mutant-DE genes are drawn "
                        "but marked")
    p.add_argument("--drop_labels",
                   default=",".join(tuple(LOW_CONFIDENCE) + NON_OVARIAN))
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def build(d, nl, args, suffix, note_txt):
    cats = [c for c in ("exacerbated", "none", "partial", "full")
            if c in d.columns]
    n = d[cats].sum(axis=1)
    frac = d[cats].div(n, axis=0)
    # Order by what the figure is about: how much of the defect is reversed.
    rank = (frac.get("full", 0) + frac.get("partial", 0)).sort_values()
    order = rank.index.tolist()
    y = np.arange(len(order))
    print(f"\n[{suffix}] {len(order)} cell types; {note_txt}")

    ncol = 2 if nl is not None else 1
    fig = plt.figure(figsize=(WIDE, max(2.5, 0.40 * len(order) + 1.3)))
    gs = fig.add_gridspec(1, ncol, wspace=0.10,
                          width_ratios=[2.2, 1.0][:ncol])

    # ---- A: diverging stacked bars ----
    ax = fig.add_subplot(gs[0, 0])
    # Everything to the left of zero is the defect made worse; the neutral
    # class straddles zero so the two sides stay comparable.
    left0 = -(frac.get("exacerbated", 0) + frac.get("none", 0) / 2)
    for t_i, t in enumerate(order):
        x = left0[t]
        for c in cats:
            w = frac.loc[t, c]
            ax.barh(t_i, w, left=x, height=0.70, color=CATS[c], zorder=3,
                    edgecolor="white", linewidth=0.4)
            x += w
    ax.axvline(0, color=COL["black"], lw=0.9, zorder=4)
    ax.set_yticks(y, [f"{ct_short([t])[0]} ({int(n[t]):,})"
                      + ("" if n[t] >= args.min_de else " *") for t in order],
                  fontsize=6)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.set_xlabel("fraction of mutant-DE genes")
    ax.xaxis.set_major_formatter(
        plt.FuncFormatter(lambda v, _: f"{abs(v):.0%}"))
    ax.set_title("What wMel does to the defect", pad=4, fontsize=7.5)
    h = [plt.Rectangle((0, 0), 1, 1, color=CATS[c]) for c in cats]
    ax.legend(h, cats, loc="lower left", bbox_to_anchor=(0, 1.06), ncol=4,
              fontsize=5.5, frameon=False, handlelength=1.1,
              columnspacing=1.0)
    cap = ("left of zero the defect is worse, right of zero reversed;\n"
           "the neutral class straddles zero. (n) = mutant-DE genes")
    if (n < args.min_de).any():
        cap += f"\n* fewer than {args.min_de} genes, so the fractions are coarse"
    note(ax, cap, y=-0.13)
    panel_letter(ax, "A", dx=-0.30)

    beats = {}
    if nl is not None:
        # ---- B: observed median against the relabellings ----
        ax = fig.add_subplot(gs[0, 1])
        obs = d["median_rescue_fraction"]
        for i, t in enumerate(order):
            if t in nl.index:
                draws = np.asarray(nl.loc[[t], "median_rf"], float)
                ax.scatter(draws, np.full(len(draws), i), s=9, marker="|",
                           color=COL["grey"], zorder=3, linewidths=0.9)
                beats[t] = bool(obs[t] > np.nanmax(draws))
            else:
                beats[t] = False
            ax.scatter([obs[t]], [i], s=34, marker="D", zorder=5,
                       color=COL["green"] if beats[t] else "white",
                       edgecolors=COL["black"], linewidths=0.6)
        ax.axvline(0, color=COL["black"], lw=0.7)
        ax.set_yticks(y, [])
        ax.set_ylim(-0.6, len(order) - 0.4)
        ax.set_xlabel("median rescue fraction")
        ax.set_title("Against shuffled labels", pad=4, fontsize=7.5)
        note(ax, "ticks: the 16 balanced relabellings\ndiamond filled when the "
                 "observed median clears all of them", y=-0.13)
        panel_letter(ax, "B", dx=-0.08)

    save(fig, args.outdir, f"{args.stem}_{suffix}")

    out = frac.reindex(order).round(3)
    out["n_mutant_DE"] = n.reindex(order).astype(int)
    if nl is not None:
        out["median_rf"] = d["median_rescue_fraction"].reindex(order).round(3)
        out["beats_null"] = [beats[t] for t in order]
    print(out.to_string())
    if nl is not None:
        won = [t for t in order if beats[t]]
        print(f"\n{len(won)} of {len(order)} cell types clear every "
              f"relabelling: {', '.join(ct_short(won)) or 'none'}. The rest sit "
              "inside the range shuffling produces, whatever their bars look "
              "like, and must not be reported as rescue.")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t", index_col=0)
    drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}
    d = d[~d.index.isin(drop)]
    for c in d.columns:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    print(f"{len(d)} cell types; dropped {sorted(drop) or 'none'}")

    nl = None
    if args.null and os.path.exists(args.null):
        nl = pd.read_csv(args.null, sep="\t")
        nl["median_rf"] = pd.to_numeric(nl["median_rf"], errors="coerce")
        nl = nl[~nl["tissue"].isin(drop)].set_index("tissue")
        print(f"null: {len(nl)} permutation rows")
    else:
        print("No --null given, so panel B is dropped. Pass the permutations "
              "file; the categories alone cannot say whether shuffling "
              "produces the same thing.")

    cats = [c for c in ("exacerbated", "none", "partial", "full")
            if c in d.columns]
    n = d[cats].sum(axis=1)
    for v in (["robust", "all"] if args.version == "both" else [args.version]):
        sub = d[n >= args.min_de] if v == "robust" else d
        txt = (f"cell types with at least {args.min_de} mutant-DE genes"
               if v == "robust" else "every cell type")
        if not len(sub):
            empty_figure(args.outdir, f"{args.stem}_{v}",
                         f"No cell type survives the {v} filter.\n"
                         "An empty set here is a result, not a failure.")
            continue
        build(sub, nl, args, v, txt)


if __name__ == "__main__":
    main()

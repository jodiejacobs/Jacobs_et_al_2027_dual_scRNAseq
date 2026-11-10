#!/usr/bin/env python
"""Figure 5: how much of the mei-P26 transcriptional defect wMel reverses.

Plots results/rescue/summary/rescue_category_counts.tsv.

This is the paper's central quantitative claim, and the categories behind it
are definitions rather than measurements, so the panel shows the composition
rather than a single headline percentage. A gene counts as rescued when its
mutant effect shrinks by at least 75% and is no longer significant with wMel;
partial at 25%; exacerbated when wMel makes the defect worse by more than 25%.
Those thresholds are choices, which is why the exacerbated class is plotted
alongside rather than dropped: a rescue claim that ignores genes moving the
wrong way is not a rescue claim.

The denominator matters as much as the fractions. A cell type with 30 mutant
DE genes can show 80% rescue and mean almost nothing, so the counts are drawn
beside the proportions and cell types below --min_de are greyed rather than
silently included.

  A  category composition per cell type, as fractions
  B  the denominator: how many genes are DE in the mutant at all
  C  median rescue fraction, the continuous version of the same quantity

Two versions are written. "robust" keeps cell types with at least --min_de
mutant DE genes; "all" keeps every cell type.

Usage:
  python plot_rescue_category_counts.py \
      --input results/rescue/summary/rescue_category_counts.tsv \
      --outdir figures/fig5
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (empty_figure, COL, LOW_CONFIDENCE, WIDE, apply_style, ct_short, note,
                      panel_letter, save)
import matplotlib.pyplot as plt

# Matches CAT_COLORS in 05_rescue_summary.py, recoloured to the Okabe-Ito set.
CATS = {"full": COL["green"], "partial": COL["skyblue"],
        "none": COL["lightgrey"], "exacerbated": COL["vermillion"]}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig5_rescue_categories")
    p.add_argument("--min_de", type=int, default=50,
                   help="mutant DE genes a cell type needs before its fractions "
                        "are worth reading")
    p.add_argument("--drop_labels", default=",".join(LOW_CONFIDENCE))
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def build(d, args, suffix, note_txt):
    cats = [c for c in CATS if c in d.columns]
    d = d.copy()
    tot = d[cats].sum(axis=1).replace(0, np.nan)
    frac = d[cats].div(tot, axis=0)
    rescued = frac[[c for c in ("full", "partial") if c in frac.columns]].sum(axis=1)
    d = d.loc[rescued.sort_values().index]
    frac, rescued = frac.loc[d.index], rescued.loc[d.index]
    print(f"\n[{suffix}] {len(d)} cell types; {note_txt}")

    fig = plt.figure(figsize=(WIDE, 2.9))
    gs = fig.add_gridspec(1, 3, wspace=0.16, width_ratios=[1.45, 1.0, 1.0])
    y = np.arange(len(d))
    thin = d["n_mutant_DE"] < args.min_de if "n_mutant_DE" in d.columns else \
        pd.Series(False, index=d.index)

    # ---- A: category composition ----
    ax = fig.add_subplot(gs[0, 0])
    left = np.zeros(len(d))
    for c in cats:
        v = frac[c].fillna(0).to_numpy()
        ax.barh(y, v, left=left, height=0.70, color=CATS[c], label=c, zorder=3)
        left += v
    for yi, t in zip(y, thin):
        if t:
            ax.text(1.01, yi, "✻", fontsize=7, va="center", color=COL["grey"])
    ax.set_yticks(y, ct_short(d.index), fontsize=6)
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of mutant DE genes")
    ax.set_title("What wMel does to the defect", pad=4)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.10), ncol=4, fontsize=6)
    note(ax, f"✻ fewer than {args.min_de} DE genes, fractions unstable", y=-0.17)
    panel_letter(ax, "A", dx=-0.52, dy=1.24)

    # ---- B: the denominator ----
    ax = fig.add_subplot(gs[0, 1])
    if "n_mutant_DE" in d.columns:
        ax.barh(y, d["n_mutant_DE"], height=0.70, zorder=3,
                color=[COL["lightgrey"] if t else COL["blue"] for t in thin])
        ax.axvline(args.min_de, color=COL["vermillion"], lw=0.9, ls="--", zorder=4)
        ax.set_xscale("log")
    ax.set_yticks(y, [])
    ax.set_xlabel("genes DE in the mutant")
    ax.set_title("How much there was to rescue", pad=4)
    panel_letter(ax, "B", dx=-0.10)

    # ---- C: the continuous version ----
    ax = fig.add_subplot(gs[0, 2])
    if "median_rescue_fraction" in d.columns:
        v = pd.to_numeric(d["median_rescue_fraction"], errors="coerce")
        ax.barh(y, v, height=0.70, zorder=3,
                color=[COL["lightgrey"] if t else
                       (COL["green"] if x > 0 else COL["vermillion"])
                       for t, x in zip(thin, v)])
        ax.axvline(0, color=COL["black"], lw=0.7)
        ax.axvline(1, color=COL["grey"], lw=0.6, ls="--")
        ax.text(1, len(d) - 0.3, " full", fontsize=6, color=COL["grey"], va="top")
    ax.set_yticks(y, [])
    ax.set_xlabel("median rescue fraction")
    ax.set_title("Effect size, not categories", pad=4)
    panel_letter(ax, "C", dx=-0.10)

    fig.suptitle(note_txt, y=1.14, fontsize=6.5, color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")

    out = d.copy()
    out["frac_rescued"] = rescued
    cols = [c for c in (cats + ["n_mutant_DE", "median_rescue_fraction",
                                "frac_rescued"]) if c in out.columns]
    print(out[cols].round(3).to_string())
    if thin.any():
        print(f"{int(thin.sum())} cell types have fewer than {args.min_de} mutant "
              "DE genes. Their rescue fractions rest on too few genes to quote.")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t", index_col=0)
    for c in d.columns:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}
    d = d[~d.index.isin(drop)]
    print(f"Labels excluded as low-confidence: {sorted(drop) or 'none'}")
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        if v == "robust" and "n_mutant_DE" in d.columns:
            sub = d[d["n_mutant_DE"] >= args.min_de]
            txt = f"cell types with at least {args.min_de} mutant DE genes"
        else:
            sub = d
            txt = "every cell type, however few genes it has"
        if not len(sub):
            empty_figure(args.outdir, f"{args.stem}_{v}",
                         f"No rows survive the {v} filter.\n"
                         f"An empty set here is a result, not a failure:\n"
                         f"say so in the text rather than citing the\n"
                         f"unfiltered figure as though it passed.")
            continue
        build(sub, args, v, txt)


if __name__ == "__main__":
    main()

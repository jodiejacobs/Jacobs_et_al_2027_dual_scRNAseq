#!/usr/bin/env python
"""Figure 5: does wMel rescue genes that go up or genes that go down in mei-P26?

Plots results/rescue/summary/rescue_genes.tsv.

Each gene that is DE in the mutant is split by the direction of its mutant
effect, and each half is classified the same way as fig5_rescue_categories
(full, partial, none, exacerbated). Panel A is the genes up in mei-P26,
panel B the genes down. If rescue were direction-blind the two panels would
match.

Read the comparison with two caveats in mind. wMel itself pushes more host
genes down than up (about 57% of its significant genes), so genes that are up
in the mutant have more room to be pulled back by infection alone. And
rescue_fraction is a ratio taken over genes selected for significance in its
own denominator, so regression to the mean inflates rescue before any biology
(see fig5_rescue_null). The permutation null was not run per direction, so a
difference between the panels is descriptive, not a test.

Two versions are written. "robust" keeps cell types with at least --min_de
mutant DE genes in total; "all" keeps every cell type. Both drop the labels in
--drop_labels. The bottom row pools the cell types shown.

Usage:
  python plot_rescue_direction.py \
      --input results/rescue/summary/rescue_genes.tsv \
      --outdir results/rescue/figures
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (empty_figure, COL, LOW_CONFIDENCE, WIDE, apply_style, ct_short,
                      panel_letter, save)
import matplotlib.pyplot as plt

# Same colours and order as fig5_rescue_categories.
CATS = {"full": COL["green"], "partial": COL["skyblue"],
        "none": COL["lightgrey"], "exacerbated": COL["vermillion"]}
PANELS = [("up", "Up in mei-P26"), ("down", "Down in mei-P26")]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig5_rescue_direction")
    p.add_argument("--min_de", type=int, default=50,
                   help="mutant DE genes a cell type needs before it is shown "
                        "in the robust version")
    p.add_argument("--drop_labels", default=",".join(LOW_CONFIDENCE))
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def tally(d):
    """Gene counts per category for one direction, indexed by cell type."""
    return (d.groupby(["tissue", "category"]).size().unstack(fill_value=0)
            .reindex(columns=list(CATS), fill_value=0))


def build(d, args, suffix, note_txt):
    tissues = d["tissue"].unique()
    counts = {k: tally(d[d["direction_in_mutant"] == k]).reindex(tissues, fill_value=0)
              for k, _ in PANELS}
    pooled = {k: v.sum().rename("All shown") for k, v in counts.items()}
    rescued = lambda c: (c["full"] + c["partial"]) / c.sum(axis=1).replace(0, np.nan)
    both = counts["up"] + counts["down"]
    order = list(rescued(both).sort_values().index)          # lowest rescue at the bottom
    rows = order + ["All shown"]
    print(f"\n[{suffix}] {len(order)} cell types; {note_txt}")

    fig, axes = plt.subplots(1, 2, figsize=(WIDE, 0.22 * len(rows) + 1.0), sharey=True,
                             gridspec_kw={"wspace": 0.55})
    y = np.arange(len(rows))
    for ax, (k, title), letter in zip(axes, PANELS, "AB"):
        c = pd.concat([counts[k], pooled[k].to_frame().T]).loc[rows]
        n = c.sum(axis=1)
        frac = c.div(n.replace(0, np.nan), axis=0).fillna(0)
        left = np.zeros(len(rows))
        for cat, col in CATS.items():
            ax.barh(y, frac[cat], left=left, height=0.72, color=col, label=cat, zorder=3)
            left += frac[cat].to_numpy()
        for yi, r in zip(y, rows):
            pct = 100 * (c.loc[r, "full"] + c.loc[r, "partial"]) / max(n[r], 1)
            ax.text(1.02, yi, f"{pct:.0f}%  n={int(n[r]):,}", fontsize=6,
                    va="center", color=COL["black"] if r == "All shown" else COL["grey"])
        ax.axhline(len(rows) - 1.5, color=COL["black"], lw=0.5)
        ax.set_xlim(0, 1)
        ax.set_ylim(-0.6, len(rows) - 0.4)
        ax.set_xlabel("fraction of genes")
        ax.set_title(title, pad=4)
        panel_letter(ax, letter, dx=-0.52 if letter == "A" else -0.08, dy=1.10)
        print(f"  {k}: " + ", ".join(f"{r} {100*(c.loc[r,'full']+c.loc[r,'partial'])/max(n[r],1):.0f}%"
                                    for r in rows[-1:]))
    axes[0].set_yticks(y, ct_short(rows), fontsize=6)
    axes[0].get_yticklabels()[-1].set_fontweight("bold")
    axes[0].legend(loc="lower left", bbox_to_anchor=(0, 1.12), ncol=4, fontsize=6)
    axes[1].text(1.0, 1.12, "label: % rescued (full + partial)", transform=axes[1].transAxes,
                 fontsize=6, color=COL["grey"], ha="right", va="bottom")
    fig.suptitle(note_txt, y=1.06, fontsize=6.5, color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")

    for k, _ in PANELS:
        print(f"\n{k} in mutant, genes per category")
        print(pd.concat([counts[k], pooled[k].to_frame().T]).loc[rows].to_string())


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}
    d = d[~d["tissue"].isin(drop)]
    print(f"Labels excluded as low-confidence: {sorted(drop) or 'none'}")
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        if v == "robust":
            n = d.groupby("tissue").size()
            sub = d[d["tissue"].isin(n[n >= args.min_de].index)]
            txt = f"cell types with at least {args.min_de} mutant DE genes"
        else:
            sub = d
            txt = "every cell type, however few genes it has"
        if not len(sub):
            empty_figure(args.outdir, f"{args.stem}_{v}",
                         f"No rows survive the {v} filter.")
            continue
        build(sub, args, v, txt)


if __name__ == "__main__":
    main()

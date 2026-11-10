#!/usr/bin/env python
"""Figure 2: whether per-cell titer separates infected from uninfected cells.

Plots results/rescue/uninfected_floor/titer_separation.tsv.

Every per-cell titer claim in the paper assumes titer measures infection in
that cell. This panel tests the assumption directly: within each cell type, how
well does titer alone tell an infected cell from an uninfected one? An AUROC
near 1 means the measurement works there. An AUROC near 0.5 means it does not,
and no titer claim in that cell type is supported at the per-cell level
whatever the group means look like.

The companion axis is the fraction of infected cells clearing the floor, which
answers a different question: even where separation is good, a low fraction
means most infected cells carry no usable signal and per-cell analyses are
running on a minority.

Two versions are written. "robust" keeps cell types that clear --min_auroc,
which are the ones a per-cell claim can rest on; "all" keeps every cell type
with enough cells, so the weak ones stay visible.

Usage:
  python plot_titer_separation.py \
      --input results/rescue/uninfected_floor/titer_separation.tsv \
      --outdir figures/fig2
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


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig2_titer_separation")
    p.add_argument("--min_auroc", type=float, default=0.7,
                   help="separation below this does not support a per-cell claim")
    p.add_argument("--drop_labels", default=",".join(LOW_CONFIDENCE))
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def build(d, args, suffix, note_txt):
    d = d.sort_values("auroc")
    y = np.arange(len(d))
    weak = d["auroc"] < args.min_auroc
    print(f"\n[{suffix}] {len(d)} cell types, {int(weak.sum())} below AUROC "
          f"{args.min_auroc}; {note_txt}")

    fig = plt.figure(figsize=(WIDE, 2.8))
    gs = fig.add_gridspec(1, 3, wspace=0.16, width_ratios=[1.25, 1.0, 1.0])

    # ---- A: separation ----
    ax = fig.add_subplot(gs[0, 0])
    ax.barh(y, d["auroc"], height=0.68, zorder=3,
            color=[COL["lightgrey"] if w else COL["blue"] for w in weak])
    ax.axvline(0.5, color=COL["black"], lw=0.7, zorder=4)
    ax.axvline(args.min_auroc, color=COL["vermillion"], lw=0.9, ls="--", zorder=4)
    ax.set_yticks(y, ct_short(d["cell_type"]), fontsize=6)
    ax.set_xlim(0.4, 1.02)
    ax.set_xlabel("AUROC, titer vs infection")
    ax.set_title("Does titer identify infection?", pad=4)
    note(ax, f"red line: {args.min_auroc}; 0.5 is chance", y=-0.20)
    panel_letter(ax, "A", dx=-0.52)

    # ---- B: fraction of cells clearing the floor ----
    ax = fig.add_subplot(gs[0, 1])
    ax.barh(y - 0.18, d["frac_infected_above_floor"], height=0.34,
            color=COL["green"], label="infected", zorder=3)
    if "frac_uninf_above_floor" in d.columns:
        ax.barh(y + 0.18, d["frac_uninf_above_floor"], height=0.34,
                color=COL["grey"], label="uninfected", zorder=3)
    ax.set_yticks(y, [])
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction above the floor")
    ax.set_title("How many cells carry signal", pad=4)
    ax.legend(loc="lower right", fontsize=6)
    panel_letter(ax, "B", dx=-0.10)

    # ---- C: median titer, infected against uninfected ----
    ax = fig.add_subplot(gs[0, 2])
    ax.barh(y - 0.18, d["median_titer_infected"], height=0.34,
            color=COL["green"], zorder=3)
    if "median_titer_uninf" in d.columns:
        ax.barh(y + 0.18, d["median_titer_uninf"], height=0.34,
                color=COL["grey"], zorder=3)
    ax.set_yticks(y, [])
    ax.set_xscale("symlog", linthresh=1e-6)
    ax.set_xlabel("median titer")
    ax.set_title("Median signal", pad=4)
    panel_letter(ax, "C", dx=-0.10)

    fig.suptitle(note_txt, y=1.04, fontsize=6.5, color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")

    cols = [c for c in ("cell_type", "n_infected", "n_uninf", "auroc", "pvalue",
                        "frac_infected_above_floor", "frac_uninf_above_floor")
            if c in d.columns]
    print(d[cols].round(4).to_string(index=False))
    if weak.any():
        print(f"Titer does not separate infection status in "
              f"{', '.join(ct_short(d.loc[weak, 'cell_type']))}. Do not make "
              "per-cell titer claims in these cell types.")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    for c in ("auroc", "pvalue", "frac_infected_above_floor",
              "frac_uninf_above_floor", "median_titer_infected",
              "median_titer_uninf", "n_infected", "n_uninf"):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}
    d = d[~d["cell_type"].isin(drop)]
    print(f"Labels excluded as low-confidence: {sorted(drop) or 'none'}")
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        sub = d[d["auroc"] >= args.min_auroc] if v == "robust" else d
        txt = (f"cell types where titer separates infection at AUROC ≥ {args.min_auroc}"
               if v == "robust" else "every cell type with enough cells of both kinds")
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

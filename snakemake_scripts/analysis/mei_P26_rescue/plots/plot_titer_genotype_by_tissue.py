#!/usr/bin/env python
"""Figure 3: does wMel titer differ between mei-P26 and OreR, by cell type?

Plots results/rescue/titer/titer_genotype_by_tissue.tsv.

If wMel rescues the mutant, an obvious first question is whether there is
simply more of it there. The model is quasi-binomial on wMel fraction with a
tissue-by-genotype term, so each cell type gets a titer ratio with a confidence
interval rather than a difference of means over four samples.

The uninfected pseudobulk titer is plotted alongside, because it is the floor
for the same cell type measured in animals with no Wolbachia. A ratio in a cell
type whose infected titer is close to its own noise floor is not interpretable
however tight its interval, and the panel marks those rather than leaving the
reader to cross-reference another figure.

Two versions are written. "robust" keeps cell types whose infected titer clears
the uninfected floor by --min_fold; "all" keeps every modelled cell type.

Usage:
  python plot_titer_genotype_by_tissue.py \
      --input results/rescue/titer/titer_genotype_by_tissue.tsv \
      --outdir figures/fig3
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (empty_figure, COL, LOW_CONFIDENCE, WIDE, apply_style, ct_short,
                      fdr_stars, note, panel_letter, save)
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig3_titer_genotype")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--min_fold", type=float, default=5.0,
                   help="infected titer must exceed the uninfected floor by this "
                        "factor for the ratio to be interpretable")
    p.add_argument("--drop_labels", default=",".join(LOW_CONFIDENCE))
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def mean_cols(d):
    return [c for c in d.columns if c.startswith("mean_pb_titer_")]


def build(d, args, suffix, note_txt):
    d = d.sort_values("titer_ratio")
    y = np.arange(len(d))
    sig = d["fdr"] < args.fdr
    print(f"\n[{suffix}] {len(d)} cell types, {int(sig.sum())} at FDR < {args.fdr}; "
          f"{note_txt}")

    fig = plt.figure(figsize=(WIDE, 2.9))
    gs = fig.add_gridspec(1, 2, wspace=0.12, width_ratios=[1.25, 1.0])

    # ---- A: the titer ratio with its interval ----
    ax = fig.add_subplot(gs[0, 0])
    for yi, (_, r) in zip(y, d.iterrows()):
        ax.plot([r["ci_low"], r["ci_high"]], [yi, yi], color=COL["grey"],
                lw=1.0, zorder=2)
    ax.scatter(d["titer_ratio"], y, s=26, zorder=3,
               color=[COL["vermillion"] if s else COL["lightgrey"] for s in sig],
               edgecolors=COL["black"], linewidths=0.4)
    for yi, (_, r) in zip(y, d.iterrows()):
        s = fdr_stars(r["fdr"])
        if s:
            ax.text(r["ci_high"] * 1.04, yi, s, fontsize=6, va="center")
    ax.axvline(1, color=COL["black"], lw=0.7)
    ax.set_xscale("log")
    ax.set_yticks(y, ct_short(d["tissue"]), fontsize=6)
    ax.set_xlabel("titer ratio, mei-P26 / OreR")
    ax.set_title("Is there more wMel in the mutant?", pad=4)
    note(ax, "bars are 95% CI; 1 is no difference", y=-0.16)
    panel_letter(ax, "A", dx=-0.52)

    # ---- B: the levels, against the uninfected floor ----
    ax = fig.add_subplot(gs[0, 1])
    mc = mean_cols(d)
    for i, c in enumerate(mc):
        ax.scatter(d[c], y + (i - (len(mc) - 1) / 2) * 0.22, s=22, zorder=3,
                   color=COL["vermillion"] if "Mei" in c else COL["blue"],
                   edgecolors=COL["black"], linewidths=0.3,
                   label=c.replace("mean_pb_titer_", ""))
    if "uninf_pb_titer" in d.columns:
        ax.scatter(d["uninf_pb_titer"], y, s=30, marker="|", color=COL["black"],
                   zorder=4, label="uninfected floor")
    ax.set_xscale("log")
    ax.set_yticks(y, [])
    ax.set_xlabel("pseudobulk wMel fraction")
    ax.set_title("Signal against its own floor", pad=4)
    ax.legend(loc="lower right", fontsize=5.5)
    panel_letter(ax, "B", dx=-0.06)

    fig.suptitle(note_txt, y=1.04, fontsize=6.5, color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")

    cols = ["tissue", "titer_ratio", "ci_low", "ci_high", "pvalue", "fdr",
            "uninf_pb_titer"] + mc
    print(d[[c for c in cols if c in d.columns]].round(5).to_string(index=False))
    if sig.any():
        print(f"\nTiter differs by genotype in "
              f"{', '.join(ct_short(d.loc[sig, 'tissue']))}. A titer difference is "
              "not itself rescue: it says how much wMel is present, not what it "
              "does.")
    else:
        print("\nNo cell type shows a genotype difference in titer. The rescue is "
              "therefore not explained by the mutant simply carrying more wMel, "
              "which is worth stating explicitly.")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    for c in d.columns:
        if c not in ("tissue", "contrast"):
            d[c] = pd.to_numeric(d[c], errors="coerce")
    drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}
    d = d[~d["tissue"].isin(drop)]
    print(f"Labels excluded as low-confidence: {sorted(drop) or 'none'}")
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        if v == "robust" and "uninf_pb_titer" in d.columns and mean_cols(d):
            lowest = d[mean_cols(d)].min(axis=1)
            sub = d[lowest > args.min_fold * d["uninf_pb_titer"].fillna(0)]
            txt = (f"cell types whose infected titer exceeds their uninfected floor "
                   f"by more than {args.min_fold:g}x")
        else:
            sub = d
            txt = "every modelled cell type, floor cleared or not"
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

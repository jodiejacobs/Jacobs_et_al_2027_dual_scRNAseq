#!/usr/bin/env python
"""Figure 5: the rescue signal against a permutation null.

Plots results/rescue/rescue_null/rescue_null_summary.tsv.

The rescue fraction is a ratio of two fitted effects, and ratios of noisy
quantities look impressive by construction: shuffling the labels and
recomputing will produce a respectable median rescue fraction in a cell type
where nothing is happening. This panel is the answer to that objection, and it
is the single most important control for the paper's central claim. For each
cell type it plots the observed statistic against the distribution obtained by
permuting the labels, so the reader can see how much of the rescue is more
than the null produces.

Two statistics are shown because they fail differently. The median rescue
fraction is the effect size and can be inflated by a few genes with small
denominators. The count of interaction-significant genes is discrete and
cannot be, but it is less sensitive. A cell type that beats the null on both
is the strongest case the design supports.

Two versions are written. "robust" keeps cell types beating the null on both
statistics; "all" keeps every cell type, so the ones that do not are visible.

Usage:
  python plot_rescue_null_summary.py \
      --input results/rescue/rescue_null/rescue_null_summary.tsv \
      --outdir figures/fig5
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
    p.add_argument("--stem", default="fig5_rescue_null")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--drop_labels", default=",".join(LOW_CONFIDENCE))
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def panel(ax, d, obs, nmean, nmax, pcol, xlabel, title, args):
    d = d.sort_values(obs)
    y = np.arange(len(d))
    sig = d[pcol] < args.alpha
    for yi, (_, r) in zip(y, d.iterrows()):
        ax.plot([r[nmean], r[nmax]], [yi, yi], color=COL["lightgrey"], lw=3.2,
                solid_capstyle="butt", zorder=2)
    ax.scatter(d[nmean], y, s=10, marker="|", color=COL["grey"], zorder=3)
    ax.scatter(d[obs], y, s=28, zorder=4,
               color=[COL["green"] if s else COL["lightgrey"] for s in sig],
               edgecolors=COL["black"], linewidths=0.4)
    for yi, (_, r) in zip(y, d.iterrows()):
        s = fdr_stars(r[pcol])
        if s:
            ax.text(r[obs], yi + 0.26, s, ha="center", fontsize=6)
    ax.set_yticks(y, ct_short(d["tissue"]), fontsize=6)
    ax.set_xlabel(xlabel)
    ax.set_title(title, pad=4)
    ax.margins(x=0.14)
    return d


def build(d, args, suffix, note_txt):
    print(f"\n[{suffix}] {len(d)} cell types; {note_txt}")
    fig = plt.figure(figsize=(WIDE, 3.0))
    gs = fig.add_gridspec(1, 2, wspace=0.10)

    ax = fig.add_subplot(gs[0, 0])
    panel(ax, d, "obs_median_rf", "null_median_rf_mean", "null_median_rf_max",
          "p_median_rf", "median rescue fraction",
          "Effect size vs permuted labels", args)
    note(ax, "grey bar: null mean to null max; point: observed", y=-0.15)
    panel_letter(ax, "A", dx=-0.56)

    ax = fig.add_subplot(gs[0, 1])
    d2 = panel(ax, d, "obs_n_interaction", "null_n_interaction_mean",
               "null_n_interaction_max", "p_n_interaction",
               "interaction-significant genes",
               "Gene count vs permuted labels", args)
    ax.set_yticks(np.arange(len(d2)), [])
    panel_letter(ax, "B", dx=-0.06)

    fig.suptitle(note_txt, y=1.03, fontsize=6.5, color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")

    cols = ["tissue", "obs_median_rf", "null_median_rf_mean", "null_median_rf_max",
            "p_median_rf", "obs_n_interaction", "null_n_interaction_mean",
            "null_n_interaction_max", "p_n_interaction"]
    print(d[[c for c in cols if c in d.columns]].round(4).to_string(index=False))
    # A permutation p cannot go below 1/(n_perm + 1). When every p sits on that
    # floor, "nothing is significant" says only that too few permutations were
    # run, and reading it as evidence against rescue inverts the result.
    pmin = float(np.nanmin(np.r_[d["p_median_rf"], d["p_n_interaction"]]))
    implied = int(round(1 / pmin - 1)) if pmin > 0 else None
    if pmin > args.alpha:
        n_at = int(np.isclose(np.r_[d["p_median_rf"], d["p_n_interaction"]],
                              pmin).sum())
        print(f"\nThe smallest p-value in this table is {pmin:.4f} = "
              f"1/({implied}+1), the floor for {implied} permutations, and "
              f"{n_at} of {2 * len(d)} tests sit exactly on it. "
              f"No cell type can reach p < {args.alpha} "
              "however large its effect, so the absence of significance here is a "
              "property of the permutation set, not of the data.\n"
              "This is not a count you can simply raise: with 2 replicates per "
              "group the relabellings are exhaustible, so the floor is set by how "
              "many distinct ones the scheme admits. Using every within-genotype "
              "relabelling, choose(4,2) = 6 per genotype, gives 35 null draws and "
              f"a floor of {1/36:.4f}, which clears {args.alpha}. A floor of "
              "0.0588 means the run used only the 16 mixed assignments; "
              "regenerate it with the full set before quoting any p-value here.")
    else:
        both = d[(d["p_median_rf"] < args.alpha) & (d["p_n_interaction"] < args.alpha)]
        print(f"\n{len(both)} cell types beat the null on both statistics: "
              f"{', '.join(ct_short(both['tissue'])) or 'none'}. These are the only "
              "ones where the rescue claim is supported against permuted labels.")

    # Separation from the null is readable whatever the permutation count.
    clear = d[d["obs_n_interaction"] > d["null_n_interaction_max"]]
    if len(clear):
        print(f"\nSeparation from the null, which the p-value floor does not "
              f"affect: {len(clear)} of {len(d)} cell types have more "
              "interaction-significant genes than any permutation produced "
              f"(observed {int(clear['obs_n_interaction'].min())} to "
              f"{int(clear['obs_n_interaction'].max())} against a null maximum of "
              f"{int(d['null_n_interaction_max'].max())}). That gap, not the "
              "p-value, is the evidence to report until the null is rerun.")
    inflated = d[(d["obs_median_rf"] <= d["null_median_rf_max"])
                 & (d["obs_median_rf"] > d["null_median_rf_mean"])]
    if len(inflated):
        print(f"{', '.join(ct_short(inflated['tissue']))} have an observed rescue "
              "fraction inside the null range. A respectable-looking number there "
              "is what shuffling the labels also produces.")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    for c in d.columns:
        if c != "tissue":
            d[c] = pd.to_numeric(d[c], errors="coerce")
    drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}
    d = d[~d["tissue"].isin(drop)]
    print(f"Labels excluded as low-confidence: {sorted(drop) or 'none'}")
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        if v == "robust":
            sub = d[(d["p_median_rf"] < args.alpha)
                    & (d["p_n_interaction"] < args.alpha)]
            txt = "cell types beating the permutation null on both statistics"
        else:
            sub = d
            txt = "every cell type, including those the null explains"
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

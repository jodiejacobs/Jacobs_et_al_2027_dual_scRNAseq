#!/usr/bin/env python
"""Eggshell figure: more cells secreting, or each cell secreting more?

Plots results/rescue/eggshell/eggshell_by_sample.tsv with the contrasts in
eggshell_tests.tsv, both from 18_eggshell.py.

A pseudobulk increase in eggshell transcripts has two possible sources, and
they are different claims. If more follicle cells reach the secreting stages,
the fraction of cells with the program on rises. If each secreting cell makes
more, the score among on cells rises. This figure puts the two side by side,
one column per stage class, one point per library.

  top row     fraction of follicle cells above their library's ambient floor
  bottom row  mean score among the cells that are on

"On" means above the 99th percentile of the same score in that library's
marker-defined germline cells, which cannot make eggshell, so the floor is the
library's own ambient level. Each panel carries the interaction estimate and
its BH FDR from the 4-df OLS over the eight libraries.

Usage:
  python plot_eggshell_score.py \
      --input results/rescue/eggshell/eggshell_by_sample.tsv \
      --tests results/rescue/eggshell/eggshell_tests.tsv \
      --outdir results/rescue/figures
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (COL, GROUP_BASE, GROUP_ORDER, WIDE, apply_style,
                      ct_short, note, panel_letter, save)
from gene_sets import EGGSHELL_ORDER
import matplotlib.pyplot as plt

NICE = {"OreR_uninf": "OreR", "OreR_wMel": "OreR\n+ wMel",
        "MeiP26_uninf": "mei-P26", "MeiP26_wMel": "mei-P26\n+ wMel"}
METRICS = [("frac_on", "fraction of cells on"),
           ("score_on", "score among on cells")]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="eggshell_by_sample.tsv")
    p.add_argument("--tests", required=True, help="eggshell_tests.tsv")
    p.add_argument("--tissue", default="all follicle cells",
                   help="a follicle-cell type, or the pooled population")
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig_eggshell_score")
    return p.parse_args()


def main():
    args = parse_args()
    apply_style()
    b = pd.read_csv(args.input, sep="\t")
    t = pd.read_csv(args.tests, sep="\t")
    b = b[b["tissue"] == args.tissue]
    t = t[(t["tissue"] == args.tissue) & (t["contrast"] == "interaction")]
    if not len(b):
        raise SystemExit(f"no rows for tissue '{args.tissue}'")
    classes = [c for c in EGGSHELL_ORDER if c in set(b["class"])]

    fig, axes = plt.subplots(2, len(classes), figsize=(WIDE, 4.2),
                             squeeze=False, sharex=True)
    rng = np.random.default_rng(0)
    for j, c in enumerate(classes):
        for i, (m, lab) in enumerate(METRICS):
            ax = axes[i][j]
            s = b[b["class"] == c]
            for k, g in enumerate(GROUP_ORDER):
                v = s.loc[s["group"] == g, m].to_numpy(float)
                if not len(v):
                    continue
                x = k + rng.uniform(-0.12, 0.12, len(v))
                ax.scatter(x, v, s=26, color=GROUP_BASE[g], zorder=3,
                           edgecolors=COL["black"], linewidths=0.4)
                ax.plot([k - 0.25, k + 0.25], [np.nanmean(v)] * 2,
                        color=COL["black"], lw=1.0, zorder=4)
            r = t[(t["class"] == c) & (t["metric"] == m)]
            if len(r):
                est, q = float(r["estimate"].iloc[0]), r["fdr"].iloc[0]
                qs = "not tested" if pd.isna(q) else f"FDR {q:.2g}"
                unit = " (logit)" if m == "frac_on" else ""
                ax.text(0.02, 0.97, f"interaction {est:+.2f}{unit}\n{qs}",
                        transform=ax.transAxes, va="top", fontsize=5.5,
                        color=COL["black"] if (pd.notna(q) and q < 0.05)
                        else COL["grey"])
            if m == "frac_on":
                ax.set_ylim(0, 1.0)
            else:
                # Headroom above the data, so the statistics never sit on a
                # point.
                v = s[m].to_numpy(float)
                lo, hi = np.nanmin(v), np.nanmax(v)
                r_ = (hi - lo) or 1.0
                ax.set_ylim(lo - 0.10 * r_, hi + 0.45 * r_)
            ax.set_xticks(range(len(GROUP_ORDER)),
                          [NICE[g] for g in GROUP_ORDER], fontsize=5.5)
            if j == 0:
                ax.set_ylabel(lab, fontsize=6.5)
            if i == 0:
                ax.set_title(c, pad=4, fontsize=7.5)
            panel_letter(ax, "ABCDEF"[i * len(classes) + j], dx=-0.22, dy=1.05)
    note(axes[1][0], f"{ct_short([args.tissue])[0]}; one point per library, "
                     "bar = group mean; on = above the library's germline "
                     "ambient floor (q0.99)", y=-0.30)
    fig.tight_layout(h_pad=1.4, w_pad=1.2)
    save(fig, args.outdir, args.stem)

    piv = b.groupby(["class", "group"])[["frac_on", "score_on"]].mean() \
        .unstack("group").round(3)
    print(piv.to_string())
    print()
    print(t[["class", "metric", "estimate", "pvalue", "fdr"]].round(4)
          .to_string(index=False))


if __name__ == "__main__":
    main()

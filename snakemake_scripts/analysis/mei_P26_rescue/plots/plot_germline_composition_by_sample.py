#!/usr/bin/env python
"""Figure 4: germline content of each sample, under three denominators.

Plots results/rescue/germline/germline_composition_by_sample.tsv.

mei-P26 ovaries lose germline cells, and whether wMel restores them is the
most direct cellular form of the rescue claim. The difficulty is the
denominator. Germline fraction of all barcodes is sensitive to how much debris
each library carries; fraction of non-artefact cells removes that but still
moves with somatic yield; fraction of germline plus Traffic-jam soma is the
most dissection-independent of the three. They can disagree, and a rescue that
appears under one denominator and not the others is a dissection artifact
rather than a result, so all three are plotted side by side rather than one
being chosen.

Each point is a sample, with four samples per genotype-infection group, so the
panel also shows how much of the apparent effect is replicate spread.

Two versions are written. "robust" keeps only the dissection-independent
denominator, which is the one to quote; "all" shows all three so the reader can
see they agree.

Usage:
  python plot_germline_composition_by_sample.py \
      --input results/rescue/germline/germline_composition_by_sample.tsv \
      --outdir figures/fig4
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import empty_figure, COL, WIDE, apply_style, note, panel_letter, save
import matplotlib.pyplot as plt

GROUPS = ["OreR_uninf", "OreR_wMel", "MeiP26_uninf", "MeiP26_wMel"]
GCOL = {"OreR_uninf": COL["skyblue"], "OreR_wMel": COL["blue"],
        "MeiP26_uninf": COL["orange"], "MeiP26_wMel": COL["vermillion"]}
METRICS = {"frac_of_all": "of all barcodes",
           "frac_of_not_artefact": "of non-artefact cells",
           "frac_of_germ_plus_tj_soma": "of germline + Tj soma"}
PREFERRED = "frac_of_germ_plus_tj_soma"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig4_germline_composition")
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def build(d, args, metrics, suffix, note_txt):
    groups = [g for g in GROUPS if g in set(d["group"])]
    print(f"\n[{suffix}] {len(d)} samples, {len(groups)} groups, "
          f"{len(metrics)} denominator(s); {note_txt}")

    fig = plt.figure(figsize=(WIDE if len(metrics) > 1 else WIDE * 0.55, 2.6))
    gs = fig.add_gridspec(1, len(metrics) + 1, wspace=0.38,
                          width_ratios=[1.0] * len(metrics) + [0.85])

    for k, m in enumerate(metrics):
        ax = fig.add_subplot(gs[0, k])
        for i, g in enumerate(groups):
            v = d.loc[d["group"] == g, m].to_numpy(float)
            j = (np.random.default_rng(i).random(len(v)) - 0.5) * 0.22
            ax.scatter(np.full(len(v), i) + j, v, s=26, color=GCOL.get(g, COL["grey"]),
                       edgecolors=COL["black"], linewidths=0.3, zorder=3)
            ax.plot([i - 0.26, i + 0.26], [v.mean()] * 2, color=COL["black"],
                    lw=1.2, zorder=4)
        ax.set_xticks(range(len(groups)),
                      [g.replace("_", "\n") for g in groups], fontsize=6)
        ax.set_ylabel("germline fraction" if k == 0 else "")
        ax.set_title(METRICS.get(m, m), pad=4, fontsize=7.5)
        ax.margins(y=0.18)
        panel_letter(ax, "ABC"[k], dx=-0.30)

    # ---- the counts behind the fractions ----
    ax = fig.add_subplot(gs[0, len(metrics)])
    for i, g in enumerate(groups):
        v = d.loc[d["group"] == g, "n_germline"].to_numpy(float)
        j = (np.random.default_rng(i).random(len(v)) - 0.5) * 0.22
        ax.scatter(np.full(len(v), i) + j, v, s=26, color=GCOL.get(g, COL["grey"]),
                   edgecolors=COL["black"], linewidths=0.3, zorder=3)
        ax.plot([i - 0.26, i + 0.26], [v.mean()] * 2, color=COL["black"],
                lw=1.2, zorder=4)
    ax.set_xticks(range(len(groups)), [g.replace("_", "\n") for g in groups],
                  fontsize=6)
    ax.set_ylabel("germline cells")
    ax.set_title("Raw counts", pad=4, fontsize=7.5)
    ax.margins(y=0.18)
    note(ax, "a fraction can rise while the count falls", y=-0.30)
    panel_letter(ax, "ABCD"[len(metrics)], dx=-0.34)

    fig.suptitle(note_txt, y=1.05, fontsize=6.5, color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")

    summ = d.groupby("group")[list(metrics) + ["n_germline", "n_cells"]].mean()
    print(summ.reindex(groups).round(4).to_string())
    if len(metrics) > 1:
        rank = {m: summ[m].reindex(groups).rank().tolist() for m in metrics}
        if len({tuple(v) for v in rank.values()}) > 1:
            print("The denominators disagree on the group ordering. Quote the "
                  "dissection-independent one and say the others differ, rather "
                  "than picking whichever supports the claim.")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    for c in list(METRICS) + ["n_cells", "n_germline", "n_not_artefact", "n_tj_soma"]:
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    have = [m for m in METRICS if m in d.columns]
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        if v == "robust":
            ms = [PREFERRED] if PREFERRED in have else have[:1]
            txt = ("germline fraction of germline plus Traffic-jam soma, the "
                   "denominator least sensitive to dissection and debris")
        else:
            ms = have
            txt = "all three denominators, so disagreement between them is visible"
        if not ms:
            empty_figure(args.outdir, f"{args.stem}_{v}",
                         f"No rows survive the {v} filter.\n"
                         f"An empty set here is a result, not a failure:\n"
                         f"say so in the text rather than citing the\n"
                         f"unfiltered figure as though it passed.")
            continue
        build(d, args, ms, v, txt)


if __name__ == "__main__":
    main()

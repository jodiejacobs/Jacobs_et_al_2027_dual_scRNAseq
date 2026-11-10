#!/usr/bin/env python
"""Figure 4: germline subclusters along the differentiation trajectory.

Plots results/rescue/germline/germline_subclusters.tsv.

If wMel rescues mei-P26 by letting cysts past a differentiation block, the
signature is a group shift in where germ cells sit along the trajectory, not a
change in how many clusters exist. The panel therefore orders subclusters by
pseudotime and plots group composition within each, so the reader can see
whether mutant cells pile up early and whether wMel moves them.

The early and late marker scores are plotted beside the pseudotime because
pseudotime has no units and its direction is set by a chosen root cell. The
marker scores are independent of that choice, so when the two agree the
ordering is trustworthy and when they do not the pseudotime should not be
read as differentiation.

Two versions are written. "robust" keeps subclusters with at least --min_cells
cells and contributions from more than one sample group; "all" keeps every
subcluster.

Usage:
  python plot_germline_subclusters.py \
      --input results/rescue/germline/germline_subclusters.tsv \
      --outdir figures/fig4
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import empty_figure, COL, WIDE, apply_style, note, panel_letter, save, shorten
import matplotlib.pyplot as plt

GROUPS = ["OreR_uninf", "OreR_wMel", "MeiP26_uninf", "MeiP26_wMel"]
GCOL = {"OreR_uninf": COL["skyblue"], "OreR_wMel": COL["blue"],
        "MeiP26_uninf": COL["orange"], "MeiP26_wMel": COL["vermillion"]}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig4_germline_subclusters")
    p.add_argument("--min_cells", type=int, default=50)
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def build(d, args, suffix, note_txt):
    fr = [f"frac_{g}" for g in GROUPS if f"frac_{g}" in d.columns]
    d = d.sort_values("median_pseudotime")
    y = np.arange(len(d))
    print(f"\n[{suffix}] {len(d)} subclusters, "
          f"{int(d['n_cells'].sum()):,} germ cells; {note_txt}")

    fig = plt.figure(figsize=(WIDE, 2.8))
    gs = fig.add_gridspec(1, 3, wspace=0.20, width_ratios=[1.0, 1.0, 1.25])

    # ---- A: group composition per subcluster ----
    ax = fig.add_subplot(gs[0, 0])
    left = np.zeros(len(d))
    for c in fr:
        v = pd.to_numeric(d[c], errors="coerce").fillna(0).to_numpy()
        g = c.replace("frac_", "")
        ax.barh(y, v, left=left, height=0.70, color=GCOL.get(g, COL["grey"]),
                label=g.replace("_", " "), zorder=3)
        left += v
    ax.set_yticks(y, [f"c{c}" for c in d["germ_cluster"]], fontsize=6)
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of the subcluster")
    ax.set_title("Who is in each subcluster", pad=4)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.08), ncol=2, fontsize=5.5)
    note(ax, "ordered by pseudotime, earliest at the bottom", y=-0.18)
    panel_letter(ax, "A", dx=-0.26, dy=1.22)

    # ---- B: pseudotime and size ----
    ax = fig.add_subplot(gs[0, 1])
    ax.barh(y, d["median_pseudotime"], height=0.70, color=COL["grey"], zorder=3)
    for yi, n in zip(y, d["n_cells"]):
        ax.text(0.01, yi, f" n={int(n)}", fontsize=5.5, va="center",
                transform=ax.get_yaxis_transform())
    ax.set_yticks(y, [])
    ax.set_xlabel("median pseudotime")
    ax.set_title("Position on the trajectory", pad=4)
    panel_letter(ax, "B", dx=-0.08)

    # ---- C: marker scores, independent of the root choice ----
    ax = fig.add_subplot(gs[0, 2])
    for c, col, lab in (("mean_early_score", COL["green"], "early (bam, bgcn, Dad)"),
                        ("mean_late_score", COL["purple"], "late (orb, osk, grk…)")):
        if c in d.columns:
            ax.plot(pd.to_numeric(d[c], errors="coerce"), y, "o-", color=col,
                    markersize=3.5, lw=1.0, label=lab, zorder=3)
    ax.axvline(0, color=COL["grey"], lw=0.6, ls="--")
    ax.set_yticks(y, [])
    ax.set_xlabel("mean marker score")
    ax.set_title("Does pseudotime track differentiation?", pad=4, fontsize=7.5)
    ax.legend(loc="lower right", fontsize=5.5)
    panel_letter(ax, "C", dx=-0.08)

    fig.suptitle(note_txt, y=1.14, fontsize=6.5, color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")

    cols = ["germ_cluster", "n_cells", "median_pseudotime", "mean_early_score",
            "mean_late_score"] + fr
    print(d[[c for c in cols if c in d.columns]].round(3).to_string(index=False))
    if "top_markers" in d.columns:
        for _, r in d.iterrows():
            print(f"  c{r['germ_cluster']}: {shorten([r['top_markers']], 90)[0]}")
    if {"mean_early_score", "mean_late_score"} <= set(d.columns):
        e = pd.to_numeric(d["mean_early_score"], errors="coerce")
        l = pd.to_numeric(d["mean_late_score"], errors="coerce")
        if len(d) > 2 and not (e.is_monotonic_decreasing or l.is_monotonic_increasing):
            print("Neither marker score moves monotonically with pseudotime, so "
                  "the ordering is not clearly differentiation. Do not describe "
                  "it as a differentiation axis without saying so.")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    for c in d.columns:
        if c not in ("germ_cluster", "top_markers"):
            d[c] = pd.to_numeric(d[c], errors="coerce")
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        if v == "robust":
            fr = [f"frac_{g}" for g in GROUPS if f"frac_{g}" in d.columns]
            mixed = (d[fr] > 0.02).sum(axis=1) > 1 if fr else True
            sub = d[(d["n_cells"] >= args.min_cells) & mixed]
            txt = (f"subclusters with ≥{args.min_cells} cells drawn from more than "
                   "one group, where a composition shift can be read")
        else:
            sub = d
            txt = "every subcluster, however small or group-specific"
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

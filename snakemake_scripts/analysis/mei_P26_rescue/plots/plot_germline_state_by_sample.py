#!/usr/bin/env python
"""Figure 4: germline differentiation state in each sample.

Plots results/rescue/germline/germline_state_by_sample.tsv.

The rescue hypothesis at the cellular level is that mei-P26 germ cells stall
early and wMel lets them through. This table is the per-sample test of that:
median pseudotime, the fraction of cells before the early cutoff, the early and
late marker scores, and the proportion of each germline subcluster.

Each metric is plotted with every sample shown, because the test has four
samples per group and a shift visible only in group means is not a shift. The
marker scores are plotted beside the pseudotime for the reason given in the
subcluster panel: pseudotime direction depends on a chosen root, and the marker
scores do not, so agreement between them is what makes the axis trustworthy.

Two versions are written. "robust" shows the metrics that do not depend on the
pseudotime root (marker scores and subcluster proportions); "all" adds the
pseudotime-derived ones.

Usage:
  python plot_germline_state_by_sample.py \
      --input results/rescue/germline/germline_state_by_sample.tsv \
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
NICE = {"median_pt": "median pseudotime", "early_by_pt": "fraction early",
        "early_score": "early marker score", "late_score": "late marker score",
        "titer": "wMel titer"}
# Metrics whose value does not depend on which cell was chosen as the root.
ROOT_FREE = ("early_score", "late_score", "titer")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig4_germline_state")
    p.add_argument("--max_panels", type=int, default=6)
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def build(d, args, metrics, suffix, note_txt):
    groups = [g for g in GROUPS if g in set(d["group"])] or sorted(set(d["group"]))
    metrics = metrics[:args.max_panels]
    n = len(metrics)
    print(f"\n[{suffix}] {len(d)} samples, {n} metric(s): "
          f"{', '.join(metrics)}; {note_txt}")

    ncol = min(n, 3)
    nrow = int(np.ceil(n / ncol))
    fig = plt.figure(figsize=(WIDE * min(1.0, ncol / 3 + 0.12), 2.3 * nrow))
    gs = fig.add_gridspec(nrow, ncol, hspace=0.70, wspace=0.46)

    for k, m in enumerate(metrics):
        ax = fig.add_subplot(gs[k // ncol, k % ncol])
        for i, g in enumerate(groups):
            v = pd.to_numeric(d.loc[d["group"] == g, m], errors="coerce").dropna()
            j = (np.random.default_rng(i).random(len(v)) - 0.5) * 0.22
            ax.scatter(np.full(len(v), i) + j, v, s=24, zorder=3,
                       color=GCOL.get(g, COL["grey"]),
                       edgecolors=COL["black"], linewidths=0.3)
            if len(v):
                ax.plot([i - 0.26, i + 0.26], [v.mean()] * 2, color=COL["black"],
                        lw=1.2, zorder=4)
        ax.set_xticks(range(len(groups)),
                      [g.replace("_", "\n") for g in groups], fontsize=5.5)
        ax.set_title(NICE.get(m, m.replace("_", " ")), pad=4, fontsize=7.5)
        ax.margins(y=0.20)
        if m not in ROOT_FREE and not m.startswith("prop_c"):
            ax.text(0.98, 0.02, "root-dependent", transform=ax.transAxes,
                    fontsize=5, color=COL["grey"], ha="right")
        panel_letter(ax, "ABCDEF"[k], dx=-0.32)

    fig.suptitle(note_txt, y=1.0 + 0.03 / nrow, fontsize=6.5, color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")

    summ = d.groupby("group")[metrics].mean().reindex(groups)
    print(summ.round(4).to_string())
    if "MeiP26_uninf" in summ.index and "MeiP26_wMel" in summ.index:
        for m in metrics:
            a, b = summ.at["MeiP26_uninf", m], summ.at["MeiP26_wMel", m]
            if np.isfinite(a) and np.isfinite(b):
                print(f"  {NICE.get(m, m)}: mei-P26 {a:.4f} -> with wMel {b:.4f}")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    props = [c for c in d.columns if c.startswith("prop_c")]
    base = [c for c in ("median_pt", "early_by_pt", "early_score", "late_score",
                        "titer") if c in d.columns]
    for c in base + props:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    if "group" not in d.columns:
        raise SystemExit(f"expected a group column, got {list(d.columns)}")
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        if v == "robust":
            ms = [c for c in base if c in ROOT_FREE] + props
            txt = ("metrics independent of the pseudotime root: marker scores and "
                   "subcluster proportions")
        else:
            ms = base + props
            txt = "every metric, including those that depend on the chosen root"
        ms = [m for m in ms if m in d.columns]
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

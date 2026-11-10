#!/usr/bin/env python
"""Figure 3: cell-type proportions in every sample.

Plots results/rescue/composition/composition_proportions.tsv.

The composition tests report contrasts. This shows the data under them, which
is the only way to see whether a significant contrast rests on a consistent
shift across four samples or on one outlying library. At this sample size that
distinction decides whether a composition result is worth reporting, so the
panel plots every sample rather than group means with error bars.

Proportions are compositional: they sum to one, so one cell type rising forces
the others down. A coordinated shift across many cell types is usually one real
change plus its arithmetic shadow, and the panel orders cell types by abundance
so the dominant populations driving that are obvious.

Two versions are written. "robust" keeps cell types above --min_prop, where a
proportion rests on enough cells to be stable; "all" keeps every cell type.

Usage:
  python plot_composition_proportions.py \
      --input results/rescue/composition/composition_proportions.tsv \
      --outdir figures/fig3
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

META = ("sample", "genotype", "infection_status", "group", "replicate", "condition")
GCOL = {"OreR_uninf": COL["skyblue"], "OreR_wMel": COL["blue"],
        "MeiP26_uninf": COL["orange"], "MeiP26_wMel": COL["vermillion"]}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig3_composition_proportions")
    p.add_argument("--min_prop", type=float, default=0.01,
                   help="mean proportion a cell type needs for the robust version")
    p.add_argument("--drop_labels", default=",".join(LOW_CONFIDENCE))
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def group_of(row):
    g = str(row.get("genotype", "")).replace("-", "")
    i = str(row.get("infection_status", ""))
    return f"{g}_{i}" if g and i else str(row.get("group", "all"))


def build(P, meta, args, suffix, note_txt):
    order = P.mean(axis=0).sort_values().index
    P = P[order]
    groups = meta["group"].tolist()
    uniq = [g for g in ("OreR_uninf", "OreR_wMel", "MeiP26_uninf", "MeiP26_wMel")
            if g in set(groups)] or sorted(set(groups))
    y = np.arange(len(order))
    print(f"\n[{suffix}] {P.shape[1]} cell types x {P.shape[0]} samples; {note_txt}")

    fig = plt.figure(figsize=(WIDE, max(2.6, 0.24 * len(order) + 1.4)))
    gs = fig.add_gridspec(1, 2, wspace=0.10, width_ratios=[1.4, 1.0])

    # ---- A: every sample as a point ----
    ax = fig.add_subplot(gs[0, 0])
    for gi, g in enumerate(uniq):
        rows = [i for i, gg in enumerate(groups) if gg == g]
        off = (gi - (len(uniq) - 1) / 2) * 0.17
        for r in rows:
            ax.scatter(P.iloc[r].to_numpy(), y + off, s=14,
                       color=GCOL.get(g, COL["grey"]), edgecolors="none",
                       zorder=3, label=g.replace("_", " ") if r == rows[0] else None)
    ax.set_xscale("symlog", linthresh=1e-3)
    ax.set_yticks(y, ct_short(order), fontsize=6)
    ax.set_xlabel("proportion of cells")
    ax.set_title("Every sample, every cell type", pad=4)
    ax.legend(loc="lower right", fontsize=5.5)
    note(ax, "four samples per group; spread between them is the real error bar",
         y=-0.09)
    panel_letter(ax, "A", dx=-0.50)

    # ---- B: group means, as a difference from OreR uninfected ----
    ax = fig.add_subplot(gs[0, 1])
    gm = P.groupby(pd.Index(groups, name="group")).mean().reindex(uniq)
    ref = "OreR_uninf" if "OreR_uninf" in gm.index else gm.index[0]
    lfc = np.log2((gm + 1e-4).div(gm.loc[ref] + 1e-4, axis=1))
    lim = float(np.nanmax(np.abs(lfc.to_numpy()))) or 1.0
    im = ax.imshow(lfc.T.to_numpy(), cmap="RdBu_r", vmin=-lim, vmax=lim,
                   aspect="auto")
    ax.set_xticks(range(len(gm.index)),
                  [g.replace("_", "\n") for g in gm.index], fontsize=5.5)
    ax.set_yticks(y, [])
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.015)
    cb.set_label(f"log$_2$ vs {ref.replace('_', ' ')}", fontsize=6)
    cb.ax.tick_params(labelsize=5.5)
    ax.set_title("Group means", pad=4)
    panel_letter(ax, "B", dx=-0.06)

    fig.suptitle(note_txt, y=1.02, fontsize=6.5, color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")
    print(gm.T.round(4).to_string())


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t", index_col=0)
    mcols = [c for c in d.columns if c in META]
    meta = d[mcols].copy()
    P = d.drop(columns=mcols).apply(pd.to_numeric, errors="coerce")
    drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}
    P = P[[c for c in P.columns if c not in drop]]
    meta["group"] = meta.apply(group_of, axis=1)
    print(f"Labels excluded as low-confidence: {sorted(drop) or 'none'}")
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        if v == "robust":
            keep = P.columns[P.mean(axis=0) >= args.min_prop]
            sub, txt = P[keep], (f"cell types averaging at least "
                                 f"{args.min_prop:.1%} of cells")
        else:
            sub, txt = P, "every cell type, however rare"
        if not sub.shape[1]:
            empty_figure(args.outdir, f"{args.stem}_{v}",
                         f"No rows survive the {v} filter.\n"
                         f"An empty set here is a result, not a failure:\n"
                         f"say so in the text rather than citing the\n"
                         f"unfiltered figure as though it passed.")
            continue
        build(sub, meta, args, v, txt)


if __name__ == "__main__":
    main()

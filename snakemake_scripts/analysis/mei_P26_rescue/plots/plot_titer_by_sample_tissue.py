#!/usr/bin/env python
"""Figure 3: wMel titer in every sample and cell type.

Plots results/rescue/titer/titer_by_sample_tissue.tsv.

The genotype model summarises this table into one ratio per cell type. Here is
what it summarises, which matters because the model has four infected
pseudobulks per cell type and a ratio from four points can be set by one of
them. Plotting every pseudobulk lets a reader see whether a cell type's ratio
is a consistent offset or a single outlier.

Uninfected samples are drawn in the same panel rather than dropped, because
they are the measured floor for each cell type and a comparison against them is
the whole reason to trust an infected value. The fraction of cells above the
floor is plotted beside the pseudobulk titer, since a high pseudobulk value can
come either from many lightly infected cells or from a few heavily infected
ones, and those mean different things.

Two versions are written. "robust" keeps the pseudobulks that passed the
minimum cell count, which are the ones the model used; "all" keeps every row.

Usage:
  python plot_titer_by_sample_tissue.py \
      --input results/rescue/titer/titer_by_sample_tissue.tsv \
      --outdir figures/fig3
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (empty_figure, COL, LOW_CONFIDENCE, WIDE, apply_style, ct_short, note,
                      panel_letter, save, truthy)
import matplotlib.pyplot as plt

GCOL = {"OreR_uninf": COL["skyblue"], "OreR_wMel": COL["blue"],
        "MeiP26_uninf": COL["orange"], "MeiP26_wMel": COL["vermillion"]}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig3_titer_by_sample_tissue")
    p.add_argument("--drop_labels", default=",".join(LOW_CONFIDENCE))
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def build(d, args, suffix, note_txt):
    order = (d[d["infection_status"] == "wMel"].groupby("tissue")["pb_titer"]
             .median().sort_values().index.tolist())
    order = order or sorted(d["tissue"].unique())
    y = np.arange(len(order))
    groups = [g for g in GCOL if g in set(d["group"])] or sorted(set(d["group"]))
    print(f"\n[{suffix}] {len(d)} pseudobulks, {len(order)} cell types; {note_txt}")

    fig = plt.figure(figsize=(WIDE, max(2.6, 0.24 * len(order) + 1.4)))
    gs = fig.add_gridspec(1, 2, wspace=0.10, width_ratios=[1.35, 1.0])

    # ---- A: pseudobulk titer ----
    ax = fig.add_subplot(gs[0, 0])
    for gi, g in enumerate(groups):
        s = d[d["group"] == g]
        off = (gi - (len(groups) - 1) / 2) * 0.17
        xs = [order.index(t) for t in s["tissue"] if t in order]
        vs = s.loc[s["tissue"].isin(order), "pb_titer"].to_numpy()
        ax.scatter(np.clip(vs, 1e-7, None), np.array(xs) + off, s=14, zorder=3,
                   color=GCOL.get(g, COL["grey"]), edgecolors="none",
                   label=g.replace("_", " "))
    ax.set_xscale("log")
    ax.set_yticks(y, ct_short(order), fontsize=6)
    ax.set_xlabel("pseudobulk wMel fraction")
    ax.set_title("Titer in every sample", pad=4)
    ax.legend(loc="lower right", fontsize=5.5)
    note(ax, "uninfected points are this cell type's own noise floor", y=-0.09)
    panel_letter(ax, "A", dx=-0.50)

    # ---- B: how many cells carry the signal ----
    ax = fig.add_subplot(gs[0, 1])
    if "frac_above_floor" in d.columns:
        for gi, g in enumerate(groups):
            s = d[d["group"] == g]
            off = (gi - (len(groups) - 1) / 2) * 0.17
            xs = [order.index(t) for t in s["tissue"] if t in order]
            vs = s.loc[s["tissue"].isin(order), "frac_above_floor"].to_numpy()
            ax.scatter(vs, np.array(xs) + off, s=14, zorder=3,
                       color=GCOL.get(g, COL["grey"]), edgecolors="none")
        ax.set_xlim(-0.02, 1.02)
    ax.set_yticks(y, [])
    ax.set_xlabel("fraction of cells above the floor")
    ax.set_title("Many light or few heavy?", pad=4)
    panel_letter(ax, "B", dx=-0.06)

    fig.suptitle(note_txt, y=1.02, fontsize=6.5, color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")

    summ = (d.groupby(["tissue", "group"])["pb_titer"].median().unstack()
            .reindex(order))
    print(summ.round(5).to_string())


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    for c in ("pb_titer", "frac_above_floor", "n_cells", "wmel_umis", "host_umis"):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    if "pb_titer" not in d.columns and {"wmel_umis", "host_umis"} <= set(d.columns):
        d["pb_titer"] = d["wmel_umis"] / (d["wmel_umis"] + d["host_umis"])
    for c in ("genotype", "infection_status", "tissue", "sample"):
        d[c] = d[c].astype(str)
    d["group"] = d["genotype"].str.replace("-", "", regex=False) + "_" + \
        d["infection_status"]
    drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}
    d = d[~d["tissue"].isin(drop) & d["pb_titer"].notna()]
    print(f"Labels excluded as low-confidence: {sorted(drop) or 'none'}")
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        if v == "robust" and "passes_min_cells" in d.columns:
            sub = d[truthy(d["passes_min_cells"]).to_numpy()]
            txt = "pseudobulks that passed the minimum cell count, as modelled"
        else:
            sub = d
            txt = "every pseudobulk, cell-count filter or not"
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

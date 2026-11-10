#!/usr/bin/env python
"""Figure 2: the wMel signal floor measured in uninfected cells.

Plots results/rescue/uninfected_floor/uninfected_floor.tsv.

Uninfected ovaries contain no Wolbachia, so any wMel reads assigned to their
cells are the assay's noise: index hopping, ambient RNA, misassignment. That
makes the uninfected samples a measured floor rather than an assumed one, and
every per-cell titer claim in the paper has to clear it. The panel plots the
whole quantile ladder rather than one chosen threshold, because which quantile
is used as the cutoff is a judgement and the reader should see what the choice
costs.

The second panel breaks the floor down by cell type, which matters because a
pooled floor hides cell types where the noise is much higher, and those are
exactly the cell types where a titer claim would be unsafe.

There is no robust/all split: this table is the noise measurement itself, and
filtering it would defeat its purpose.

Usage:
  python plot_uninfected_floor.py \
      --input results/rescue/uninfected_floor/uninfected_floor.tsv \
      --outdir figures/fig2
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import COL, WIDE, apply_style, ct_short, note, panel_letter, save
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig2_uninfected_floor")
    p.add_argument("--floor_quantile", type=float, default=0.99,
                   help="the quantile used as the floor elsewhere in the pipeline")
    p.add_argument("--pooled_label", default="pooled_uninfected")
    return p.parse_args()


def qcols(d):
    out = []
    for c in d.columns:
        if str(c).startswith("q"):
            try:
                out.append((float(str(c)[1:]), c))
            except ValueError:
                pass
    return [c for _, c in sorted(out)], [q for q, _ in sorted(out)]


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    cols, qs = qcols(d)
    if not cols:
        raise SystemExit("no quantile columns (q0.5, q0.9, ...) in this table")
    for c in cols + ["n_cells"]:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    pooled = d[d["group"] == args.pooled_label]
    per_ct = d[d["group"] != args.pooled_label].copy()
    chosen = f"q{args.floor_quantile:g}"
    print(f"{len(per_ct)} cell types plus the pooled row; quantiles: "
          f"{', '.join(cols)}")
    if chosen in d.columns and len(pooled):
        print(f"Pooled floor at {chosen} = {pooled[chosen].iloc[0]:.3e}")

    fig = plt.figure(figsize=(WIDE, 2.6))
    gs = fig.add_gridspec(1, 2, wspace=0.30, width_ratios=[1.0, 1.25])

    # ---- A: the quantile ladder ----
    ax = fig.add_subplot(gs[0, 0])
    for _, r in per_ct.iterrows():
        ax.plot(qs, r[cols].to_numpy(float), color=COL["lightgrey"], lw=0.7, zorder=2)
    if len(pooled):
        ax.plot(qs, pooled[cols].to_numpy(float).ravel(), color=COL["black"],
                lw=1.6, zorder=4, label="pooled uninfected")
    if chosen in d.columns:
        ax.axvline(args.floor_quantile, color=COL["vermillion"], lw=0.9, ls="--",
                   zorder=3, label=f"floor used ({chosen})")
    ax.set_yscale("symlog", linthresh=1e-6)
    ax.set_xlabel("quantile of uninfected cells")
    ax.set_ylabel("wMel titer")
    ax.set_title("Where the noise sits", pad=4)
    ax.legend(loc="upper left", fontsize=6)
    note(ax, "grey lines are individual cell types", y=-0.26)
    panel_letter(ax, "A", dx=-0.24)

    # ---- B: the chosen floor, per cell type ----
    ax = fig.add_subplot(gs[0, 1])
    if chosen in per_ct.columns and len(per_ct):
        s = per_ct.sort_values(chosen)
        y = np.arange(len(s))
        ax.barh(y, s[chosen], height=0.68, color=COL["blue"], zorder=3)
        if len(pooled):
            pv = float(pooled[chosen].iloc[0])
            ax.axvline(pv, color=COL["vermillion"], lw=1.0, ls="--", zorder=4)
            ax.text(pv, len(s) - 0.3, " pooled", color=COL["vermillion"],
                    fontsize=6, va="top")
            worse = s[s[chosen] > pv]
            if len(worse):
                print(f"{len(worse)} cell types have a higher floor than the pooled "
                      f"value: {', '.join(ct_short(worse['group']))}. A per-cell "
                      "titer claim in these is weaker than the pooled floor implies.")
        ax.set_yticks(y, ct_short(s["group"]), fontsize=6)
        ax.set_xlabel(f"wMel titer at {chosen} of uninfected cells")
        ax.set_title("Floor by cell type", pad=4)
        ax.set_xscale("symlog", linthresh=1e-6)
    panel_letter(ax, "B", dx=-0.46)

    save(fig, args.outdir, args.stem)
    print(d[["group", "n_cells"] + cols].to_string(index=False))


if __name__ == "__main__":
    main()

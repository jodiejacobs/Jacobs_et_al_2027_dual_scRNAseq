#!/usr/bin/env python
"""Supplementary: cells kept and removed as debris, by sample and cell type.

Plots results/rescue/label_qc/cell_counts_by_sample_tissue.tsv.

Debris removal is the single filtering decision with the largest downstream
effect in this analysis: it took ambient wMel from 1.21% to 0.39% of the
infected signal, and it is why several comparisons are run on fewer cells than
the raw capture suggests. This panel records what it removed, so the filter is
reported rather than assumed.

The part that needs watching is whether debris was removed evenly. A cell type
or a sample losing a much larger share than the rest means the filter is also
removing a biological population, and any count comparison involving it
inherits that. The debris rate is therefore plotted per cell type and per
sample, not as one overall number.

There is no robust/all split: this is the filtering record itself.

Usage:
  python plot_cell_counts_by_sample_tissue.py \
      --input results/rescue/label_qc/cell_counts_by_sample_tissue.tsv \
      --outdir figures/supp
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
    p.add_argument("--stem", default="supp_cell_counts_debris")
    p.add_argument("--flag_rate", type=float, default=0.5,
                   help="a group losing more than this share is flagged")
    return p.parse_args()


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    for c in ("n_labeled", "n_debris", "n_kept"):
        d[c] = pd.to_numeric(d[c], errors="coerce").fillna(0)
    if "n_kept" not in d.columns:
        d["n_kept"] = d["n_labeled"] - d["n_debris"]
    tot_lab, tot_deb = d["n_labeled"].sum(), d["n_debris"].sum()
    print(f"{int(tot_lab):,} labelled cells, {int(tot_deb):,} removed as debris "
          f"({tot_deb / max(tot_lab, 1):.1%}), {int(d['n_kept'].sum()):,} kept")

    by_ct = d.groupby("tissue")[["n_labeled", "n_debris", "n_kept"]].sum()
    by_ct["rate"] = by_ct["n_debris"] / by_ct["n_labeled"].replace(0, np.nan)
    by_ct = by_ct.sort_values("rate")
    by_s = d.groupby("sample")[["n_labeled", "n_debris", "n_kept"]].sum()
    by_s["rate"] = by_s["n_debris"] / by_s["n_labeled"].replace(0, np.nan)
    by_s = by_s.sort_values("rate")

    fig = plt.figure(figsize=(WIDE, max(2.8, 0.22 * len(by_ct) + 1.6)))
    gs = fig.add_gridspec(1, 3, wspace=0.18, width_ratios=[1.25, 1.0, 1.0])

    # ---- A: debris rate per cell type ----
    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(len(by_ct))
    ax.barh(y, by_ct["rate"], height=0.68, zorder=3,
            color=[COL["vermillion"] if r > args.flag_rate else COL["blue"]
                   for r in by_ct["rate"]])
    ax.axvline(tot_deb / max(tot_lab, 1), color=COL["black"], lw=0.9, ls="--")
    ax.set_yticks(y, ct_short(by_ct.index), fontsize=6)
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction removed as debris")
    ax.set_title("Was debris removed evenly?", pad=4)
    note(ax, "dashed line is the overall rate", y=-0.10)
    panel_letter(ax, "A", dx=-0.50)

    # ---- B: cells kept per cell type ----
    ax = fig.add_subplot(gs[0, 1])
    ax.barh(y, by_ct["n_kept"].clip(lower=0.5), height=0.68, color=COL["green"],
            zorder=3)
    ax.set_xscale("log")
    ax.set_yticks(y, [])
    ax.set_xlabel("cells kept")
    ax.set_title("What is left to analyse", pad=4)
    panel_letter(ax, "B", dx=-0.08)

    # ---- C: debris rate per sample ----
    ax = fig.add_subplot(gs[0, 2])
    ys = np.arange(len(by_s))
    ax.barh(ys, by_s["rate"], height=0.68, zorder=3,
            color=[COL["vermillion"] if r > args.flag_rate else COL["grey"]
                   for r in by_s["rate"]])
    ax.set_yticks(ys, [str(s).replace("_", " ") for s in by_s.index], fontsize=6)
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction removed")
    ax.set_title("By library", pad=4)
    panel_letter(ax, "C", dx=-0.52)

    save(fig, args.outdir, args.stem)

    print(by_ct.round(3).to_string())
    print("\nBy sample:")
    print(by_s.round(3).to_string())
    hot = by_ct[by_ct["rate"] > args.flag_rate]
    if len(hot):
        print(f"\n{len(hot)} cell types lose more than {args.flag_rate:.0%} to the "
              f"debris filter: {', '.join(ct_short(hot.index))}. The filter may "
              "be removing a real population there, so count comparisons "
              "involving them need that caveat.")
    spread = by_s["rate"].max() - by_s["rate"].min()
    if spread > 0.2:
        print(f"Debris rate varies by {spread:.0%} between libraries "
              f"({by_s.index[0]} to {by_s.index[-1]}). Any count compared across "
              "samples inherits that difference.")


if __name__ == "__main__":
    main()

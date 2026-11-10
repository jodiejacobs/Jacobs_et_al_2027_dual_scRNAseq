#!/usr/bin/env python
"""Supplementary: cell-type composition against the reference atlas.

Plots results/rescue/atlas_fidelity/atlas_proportions.tsv.

Cell-type labels here are transferred from a published atlas, so every
downstream per-cell-type claim inherits whatever that transfer got wrong. This
panel is the check: if the composition of the wild-type uninfected samples
departs sharply from the atlas, either the dissection sampled different tissue
or the labels are being assigned to the wrong cells, and either way a
per-cell-type result needs that caveat.

Departure is expected and is not by itself a problem. A dissection enriched for
germaria should over-represent early germline, and saying so is better than
pretending the composition should match. What matters is whether the departures
have an explanation, so the panel plots the full ladder rather than a single
agreement statistic.

There is no robust/all split: this is a composition comparison with no
per-cell-type statistics to filter on.

Usage:
  python plot_atlas_proportions.py \
      --input results/rescue/atlas_fidelity/atlas_proportions.tsv \
      --outdir figures/supp
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import COL, WIDE, apply_style, ct_short, note, panel_letter, save
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="supp_atlas_proportions")
    p.add_argument("--query_col", default=None,
                   help="which query column to compare; default is the reference "
                        "group column if present, else query_all")
    return p.parse_args()


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t", index_col=0)
    for c in d.columns:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    qcols = [c for c in d.columns if c.startswith("query_")]
    col = args.query_col or next((c for c in qcols if c != "query_all"), "query_all")
    if col not in d.columns:
        raise SystemExit(f"{col} not in {list(d.columns)}")
    print(f"Comparing {col} against atlas over {len(d)} cell types")

    shared = d[(d["atlas"] > 0) & (d[col] > 0)]
    rho = (stats.spearmanr(shared[col], shared["atlas"]) if len(shared) > 2
           else None)
    if rho is not None:
        print(f"Composition agreement over {len(shared)} shared types: "
              f"ρ = {rho.statistic:.3f} (p = {rho.pvalue:.3g})")

    fig = plt.figure(figsize=(WIDE, 2.9))
    gs = fig.add_gridspec(1, 2, wspace=0.34, width_ratios=[1.0, 1.25])

    # ---- A: query against atlas proportion ----
    ax = fig.add_subplot(gs[0, 0])
    lim = float(max(d[col].max(), d["atlas"].max())) * 1.2
    ax.plot([1e-5, lim], [1e-5, lim], color=COL["grey"], lw=0.6, ls="--", zorder=1)
    ax.scatter(d["atlas"], d[col], s=26, color=COL["blue"], zorder=3,
               edgecolors=COL["black"], linewidths=0.4)
    far = d.reindex(d["log2_ratio_vs_atlas"].abs().sort_values().tail(5).index)
    for ct, r in far.iterrows():
        ax.annotate(ct_short([ct])[0], (r["atlas"], r[col]),
                    textcoords="offset points", xytext=(4, 3), fontsize=5.5)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("atlas proportion")
    ax.set_ylabel("observed proportion")
    title = "Composition vs atlas"
    if rho is not None:
        title += f" (ρ = {rho.statistic:.2f})"
    ax.set_title(title, pad=4)
    panel_letter(ax, "A", dx=-0.28)

    # ---- B: the departures, ranked ----
    ax = fig.add_subplot(gs[0, 1])
    s = d.sort_values("log2_ratio_vs_atlas")
    y = np.arange(len(s))
    ax.barh(y, s["log2_ratio_vs_atlas"], height=0.68, zorder=3,
            color=[COL["vermillion"] if v > 0 else COL["blue"]
                   for v in s["log2_ratio_vs_atlas"]])
    ax.axvline(0, color=COL["black"], lw=0.7)
    for v in (-1, 1):
        ax.axvline(v, color=COL["grey"], lw=0.6, ls="--")
    ax.set_yticks(y, ct_short(s.index), fontsize=6)
    ax.set_xlabel("log$_2$ observed / atlas")
    ax.set_title("Where the dissection differs", pad=4)
    note(ax, "dashed lines at two-fold", y=-0.16)
    panel_letter(ax, "B", dx=-0.46)

    save(fig, args.outdir, args.stem)

    print(d[[col, "atlas", "log2_ratio_vs_atlas"]]
          .sort_values("log2_ratio_vs_atlas").round(4).to_string())
    big = s[s["log2_ratio_vs_atlas"].abs() > 1]
    if len(big):
        print(f"\n{len(big)} cell types differ from the atlas by more than "
              f"two-fold: {', '.join(ct_short(big.index))}. Give the dissection "
              "reason where there is one, and flag the rest as a label-transfer "
              "caveat on any per-cell-type claim in them.")


if __name__ == "__main__":
    main()

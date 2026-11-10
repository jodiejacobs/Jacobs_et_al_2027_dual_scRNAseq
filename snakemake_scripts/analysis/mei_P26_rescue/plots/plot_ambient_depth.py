#!/usr/bin/env python
"""Figure 2: whether the wMel signal in uninfected cells is ambient.

Plots results/rescue/uninfected_floor/ambient_depth.tsv.

The wMel reads in uninfected cells have to come from somewhere, and the two
candidates behave differently. Ambient RNA is captured in proportion to how
much a droplet captures of everything, so it predicts a positive correlation
between host UMIs and wMel UMIs within a sample. Index hopping and barcode
misassignment do not, because they depend on the library neighbourhood rather
than on droplet depth. The panel tests that prediction per sample, which is
the right unit: pooling samples would manufacture a correlation out of
differences in sequencing depth between them.

The result is also the justification for including host depth as a covariate
in the titer-program analysis, so this panel and Figure 7 argue together.

There is no robust/all split: each sample contributes one measured correlation
and dropping any of them would weaken the control.

Usage:
  python plot_ambient_depth.py \
      --input results/rescue/uninfected_floor/ambient_depth.tsv \
      --outdir figures/fig2
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import COL, WIDE, apply_style, fdr_stars, note, panel_letter, save
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig2_ambient_depth")
    p.add_argument("--alpha", type=float, default=0.05)
    return p.parse_args()


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    for c in ("n_cells", "median_host_umis", "median_wmel_umis", "pb_titer",
              "spearman_wmel_vs_host_umis", "pvalue"):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    d["label"] = d["sample"].astype(str).str.replace("_", " ")
    d = d.sort_values("spearman_wmel_vs_host_umis")
    sig = d["pvalue"] < args.alpha
    pos = d["spearman_wmel_vs_host_umis"] > 0
    print(f"{len(d)} uninfected samples; {int((sig & pos).sum())} show a positive "
          f"depth correlation at p < {args.alpha}")

    fig = plt.figure(figsize=(WIDE, 2.4))
    gs = fig.add_gridspec(1, 3, wspace=0.42)

    # ---- A: the depth correlation, per sample ----
    ax = fig.add_subplot(gs[0, 0])
    y = np.arange(len(d))
    ax.barh(y, d["spearman_wmel_vs_host_umis"], height=0.66, zorder=3,
            color=[COL["vermillion"] if s else COL["lightgrey"] for s in sig])
    for yi, (_, r) in zip(y, d.iterrows()):
        v = r["spearman_wmel_vs_host_umis"]
        ax.text(v + (0.01 if v >= 0 else -0.01), yi,
                fdr_stars(r["pvalue"]) or "n.s.", fontsize=5.5, va="center",
                ha="left" if v >= 0 else "right")
    ax.axvline(0, color=COL["black"], lw=0.7)
    ax.set_yticks(y, d["label"], fontsize=6)
    ax.set_xlabel("ρ, wMel UMIs vs host UMIs")
    ax.set_title("Does noise scale with depth?", pad=4)
    ax.margins(x=0.22)
    panel_letter(ax, "A", dx=-0.58)

    # ---- B: depth against the per-sample wMel level ----
    ax = fig.add_subplot(gs[0, 1])
    ax.scatter(d["median_host_umis"], d["pb_titer"], s=26, color=COL["blue"],
               edgecolors=COL["black"], linewidths=0.4, zorder=3)
    for _, r in d.iterrows():
        ax.annotate(r["label"], (r["median_host_umis"], r["pb_titer"]),
                    textcoords="offset points", xytext=(4, 3), fontsize=5.5)
    ax.set_xlabel("median host UMIs per cell")
    ax.set_ylabel("pooled wMel fraction")
    ax.set_title("Sample depth vs noise level", pad=4)
    ax.margins(0.18)
    panel_letter(ax, "B", dx=-0.40)

    # ---- C: how much wMel a typical uninfected cell carries ----
    ax = fig.add_subplot(gs[0, 2])
    ax.barh(y, d.set_index(d.index)["median_wmel_umis"].reindex(d.index),
            height=0.66, color=COL["grey"], zorder=3)
    ax.set_yticks(y, d["label"], fontsize=6)
    ax.set_xlabel("median wMel UMIs per cell")
    ax.set_title("Noise per cell", pad=4)
    ax.margins(x=0.18)
    panel_letter(ax, "C", dx=-0.58)

    save(fig, args.outdir, args.stem)

    cols = [c for c in ("sample", "n_cells", "median_host_umis", "median_wmel_umis",
                        "pb_titer", "spearman_wmel_vs_host_umis", "pvalue")
            if c in d.columns]
    print(d[cols].round(5).to_string(index=False))
    if (sig & pos).any():
        print("\nwMel signal in uninfected cells rises with droplet depth, which is "
              "the ambient signature. This is why host depth is a covariate in the "
              "titer-program analysis; without it a depth effect would be read as "
              "a titer effect.")
    else:
        print("\nNo sample shows wMel signal rising with droplet depth, so the "
              "uninfected signal does not look ambient. Index hopping or barcode "
              "misassignment is the better explanation, and depth adjustment "
              "will not remove it.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Figure 1: how many wMel loci DUOS-seq recovers, and at what cell support.

Plots results/rescue/detection_bias/wmel_breadth.tsv.

This is the method claim and it rests on one number, so the panel shows how
that number behaves as the evidence bar is raised rather than quoting it at a
single threshold. INVADEseq recovers one targeted 16S locus; the comparison
drawn in the text is against that, and the reference line makes it visible
instead of leaving it to the reader's memory.

There is no robust/all split here: the table has no per-locus statistics to
filter on, and the x-axis already is the stringency axis.

Usage:
  python plot_wmel_breadth.py \
      --input results/rescue/detection_bias/wmel_breadth.tsv \
      --outdir figures/fig1
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import COL, COLUMN, apply_style, panel_letter, save
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig1_wmel_breadth")
    p.add_argument("--reference_n", type=int, default=1,
                   help="loci recovered by the targeted method being compared against")
    p.add_argument("--reference_label", default="INVADEseq (1 targeted locus)")
    return p.parse_args()


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    for c in ("min_cells", "loci_detected", "frac_of_loci"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d.sort_values("min_cells")
    print(f"{len(d)} thresholds; {int(d['loci_detected'].iloc[0])} loci at the "
          f"loosest ({int(d['min_cells'].iloc[0])} cell), "
          f"{int(d['loci_detected'].iloc[-1])} at the strictest "
          f"({int(d['min_cells'].iloc[-1])} cells)")

    fig, ax = plt.subplots(figsize=(COLUMN, 2.3))
    x = np.arange(len(d))
    ax.bar(x, d["loci_detected"], color=COL["blue"], width=0.68, zorder=3)
    for xi, v, f in zip(x, d["loci_detected"], d["frac_of_loci"]):
        ax.text(xi, v, f"{int(v)}\n{f:.0%}", ha="center", va="bottom",
                fontsize=5.5, linespacing=1.1)
    if args.reference_n:
        ax.axhline(args.reference_n, color=COL["vermillion"], lw=0.9, ls="--",
                   zorder=4)
        ax.text(len(d) - 0.4, args.reference_n, f" {args.reference_label}",
                color=COL["vermillion"], fontsize=6, va="bottom", ha="right")
    ax.set_xticks(x, [f"≥{int(v)}" for v in d["min_cells"]])
    ax.set_xlabel("cells in which the locus is detected")
    ax.set_ylabel("wMel loci detected")
    ax.set_title("Breadth of wMel recovery", pad=4)
    ax.margins(y=0.20)
    panel_letter(ax, "A", dx=-0.19)
    save(fig, args.outdir, args.stem)

    print(d.to_string(index=False))
    print("\nPercentages are of loci that have coordinates in the annotation, "
          "not of all annotated features, so they do not describe loci the "
          "reference omits.")


if __name__ == "__main__":
    main()

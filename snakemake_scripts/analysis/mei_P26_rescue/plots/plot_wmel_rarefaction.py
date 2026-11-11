#!/usr/bin/env python
"""Figure 1: whether wMel locus recovery has saturated at the depth sequenced.

Plots results/rescue/detection_bias/wmel_rarefaction.tsv.

A breadth number means little on its own, because recovering many loci is
easy with enough cells. What matters for the method claim is whether the curve
has flattened: if it has, the locus count is a property of the assay, and if it
has not, it is a property of how many cells happened to be captured and a
larger run would report a larger number. The panel therefore plots the curve
and the slope over the final decile, and states which of the two conclusions
the data support rather than leaving the curve to speak for itself.

There is no robust/all split here: the table has no per-locus statistics to
filter on.

Usage:
  python plot_wmel_rarefaction.py \
      --input results/rescue/detection_bias/wmel_rarefaction.tsv \
      --outdir figures/fig1
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import COL, COLUMN, apply_style, note, panel_letter, save
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig1_wmel_rarefaction")
    p.add_argument("--total_loci", type=int, default=None,
                   help="wMel loci with coordinates in the annotation. Given it, "
                        "a curve reaching this value is reported as complete "
                        "recovery rather than as an unsaturated one")
    p.add_argument("--flat_threshold", type=float, default=0.02,
                   help="loci gained per 1%% more cells, below which the curve is "
                        "called saturated")
    return p.parse_args()


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t").sort_values("n_cells")
    d["n_cells"] = pd.to_numeric(d["n_cells"], errors="coerce")
    d["loci_detected"] = pd.to_numeric(d["loci_detected"], errors="coerce")
    total, final = d["loci_detected"].iloc[-1], d["n_cells"].iloc[-1]

    # Slope over the final interval, in loci gained per 1% more cells. The grid
    # is log spaced, so a fixed 10% window usually contains a single point and
    # would give no slope at all; the last two points always do.
    tail = d[d["n_cells"] >= 0.9 * final]
    if len(tail) < 2:
        tail = d.tail(2)
    gained = tail["loci_detected"].iloc[-1] - tail["loci_detected"].iloc[0]
    frac = (tail["n_cells"].iloc[-1] - tail["n_cells"].iloc[0]) / final * 100
    slope = gained / frac if frac else np.nan
    # Recovery can stop rising because the assay ran out of depth, or because
    # it has found every locus the annotation contains. Those support opposite
    # claims and the second is only knowable with the total.
    complete = args.total_loci is not None and total >= args.total_loci
    saturated = complete or (np.isfinite(slope) and slope < args.flat_threshold)
    span = f"{int(tail['n_cells'].iloc[0]):,} to {int(final):,} cells"
    print(f"{total} loci at {int(final):,} cells; over {span} the curve adds "
          f"{gained:.0f} loci ({slope:.3f} per 1% more cells)")
    if complete:
        print(f"-> COMPLETE: every one of the {args.total_loci:,} annotated loci "
              "is detected, so the curve stops at the annotation's ceiling and "
              "not at a limit of the assay. Do not call this unsaturated.")
    else:
        print(f"-> {'saturated' if saturated else 'NOT saturated'}")

    fig, ax = plt.subplots(figsize=(COLUMN, 2.3))
    ax.plot(d["n_cells"], d["loci_detected"], color=COL["blue"], lw=1.3, zorder=3)
    ax.scatter(d["n_cells"], d["loci_detected"], s=8, color=COL["blue"], zorder=4)
    ax.plot(tail["n_cells"], tail["loci_detected"], color=COL["vermillion"],
            lw=2.0, zorder=5, solid_capstyle="round")
    if args.total_loci:
        ax.axhline(args.total_loci, color=COL["vermillion"], lw=0.7, ls=":", zorder=1)
        ax.text(d["n_cells"].iloc[0], args.total_loci, " all annotated loci",
                fontsize=6, color=COL["vermillion"], va="bottom")
    ax.axhline(total, color=COL["grey"], lw=0.6, ls="--", zorder=1)
    ax.annotate(f"{int(total):,} loci", (final, total), textcoords="offset points",
                xytext=(-3, 4), ha="right", fontsize=6.5)
    ax.set_xlabel("cells sampled")
    ax.set_ylabel("wMel loci detected")
    ax.set_title("Complete recovery?" if complete else "Is recovery saturated?",
                 pad=4)
    ax.margins(y=0.14)
    tail_txt = ("all annotated loci recovered" if complete
                else f"{'flat' if saturated else 'still rising'}")
    note(ax, f"final interval (red): {slope:.2f} loci per 1% more cells, "
             f"{tail_txt}", y=-0.30)
    panel_letter(ax, "B", dx=-0.19)
    save(fig, args.outdir, args.stem)

    if not saturated:
        print("The curve has not flattened, so the locus count is a lower bound "
              "set by cell number and a deeper run would report more. Say that "
              "in the text rather than quoting the count as the assay's ceiling.")
    elif complete:
        print("The quotable claim is completeness, not saturation: DUOS-seq "
              "recovers the entire annotated wMel gene set at this depth.")


if __name__ == "__main__":
    main()

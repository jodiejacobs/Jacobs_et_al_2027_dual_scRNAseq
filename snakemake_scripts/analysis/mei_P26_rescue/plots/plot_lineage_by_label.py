#!/usr/bin/env python
"""Supplementary: do transferred labels agree with independent lineage calls?

Plots results/rescue/label_qc/lineage_by_label.tsv.

Atlas labels and lineage calls come from different evidence: the labels are
transferred by nearest-neighbour matching to the atlas, the lineage call from
germline and soma markers scored directly. Where they agree the label is
supported by something other than the transfer. Where they disagree, a cell
type's name is asserting more than the expression supports, and every result
in that cell type carries the doubt.

The germline labels matter most here, because the rescue claim rests on them.
A germline label that is only partly germline by marker score is a label whose
cells are mixed, and a per-cell-type effect in it can be a composition effect.

There is no robust/all split: this is the QC table itself, and filtering it
would hide the disagreements it exists to show.

Usage:
  python plot_lineage_by_label.py \
      --input results/rescue/label_qc/lineage_by_label.tsv \
      --outdir figures/supp
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import COL, PAL, WIDE, apply_style, ct_short, note, panel_letter, save
import matplotlib.pyplot as plt

GERMLINE_KEYS = ("germ", "cyst", "germline")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="supp_lineage_by_label")
    p.add_argument("--min_cells", type=int, default=50)
    p.add_argument("--purity", type=float, default=0.7,
                   help="a label below this fraction in its dominant lineage is "
                        "flagged as mixed")
    return p.parse_args()


def is_germline(name):
    return any(k in str(name).lower() for k in GERMLINE_KEYS)


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t", index_col=0)
    lin = [c for c in d.columns if c != "n_cells"]
    for c in d.columns:
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d[d["n_cells"].fillna(0) >= args.min_cells]
    d = d.loc[sorted(d.index, key=lambda r: (not is_germline(r), str(r)))]
    print(f"{len(d)} labels with ≥{args.min_cells} cells; lineages: {', '.join(lin)}")

    fig = plt.figure(figsize=(WIDE, max(2.6, 0.22 * len(d) + 1.4)))
    gs = fig.add_gridspec(1, 2, wspace=0.10, width_ratios=[1.45, 1.0])
    y = np.arange(len(d))

    # ---- A: lineage composition of each label ----
    ax = fig.add_subplot(gs[0, 0])
    left = np.zeros(len(d))
    for i, c in enumerate(lin):
        v = d[c].fillna(0).to_numpy()
        ax.barh(y, v, left=left, height=0.70, color=PAL[i % len(PAL)],
                label=str(c), zorder=3)
        left += v
    ticks = ct_short(d.index)
    ax.set_yticks(y, ticks, fontsize=6)
    for t, r in zip(ax.get_yticklabels(), d.index):
        t.set_color(COL["vermillion"] if is_germline(r) else COL["black"])
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of cells with this lineage call")
    ax.set_title("Label against independent lineage call", pad=4)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.02),
              ncol=min(len(lin), 5), fontsize=6)
    note(ax, "red label names are germline by atlas label", y=-0.10)
    panel_letter(ax, "A", dx=-0.50, dy=1.10)

    # ---- B: how many cells each label holds ----
    ax = fig.add_subplot(gs[0, 1])
    ax.barh(y, d["n_cells"], height=0.70, color=COL["grey"], zorder=3)
    ax.set_xscale("log")
    ax.set_yticks(y, [])
    ax.set_xlabel("cells (debris removed)")
    ax.set_title("Label size", pad=4)
    panel_letter(ax, "B", dx=-0.06)

    save(fig, args.outdir, args.stem)

    pur = d[lin].max(axis=1)
    dom = d[lin].idxmax(axis=1)
    out = pd.DataFrame({"n_cells": d["n_cells"].astype(int),
                        "dominant_lineage": dom, "purity": pur.round(3)})
    print(out.to_string())
    mixed = out[out["purity"] < args.purity]
    if len(mixed):
        print(f"\n{len(mixed)} labels are below {args.purity:.0%} pure in their "
              f"dominant lineage: {', '.join(ct_short(mixed.index))}. An effect "
              "in a mixed label can be a composition shift rather than a change "
              "within one cell type.")
    gl = out[[is_germline(i) for i in out.index]]
    bad_gl = gl[~gl["dominant_lineage"].astype(str).str.lower()
                .str.contains("germ")]
    if len(bad_gl):
        print(f"These germline labels are not dominantly germline by marker "
              f"score: {', '.join(ct_short(bad_gl.index))}. The rescue claim "
              "rests on germline labels, so this needs stating.")


if __name__ == "__main__":
    main()

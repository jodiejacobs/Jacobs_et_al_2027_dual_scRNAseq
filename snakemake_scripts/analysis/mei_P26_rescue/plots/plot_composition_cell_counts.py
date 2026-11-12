#!/usr/bin/env python
"""Figure 3: the raw cell counts behind every composition proportion.

Plots results/rescue/composition/composition_cell_counts.tsv.

Proportions hide two things that decide whether a composition result is real.
The first is total yield: if one group's libraries recovered half as many cells,
every proportion in them is estimated from half the data. The second is
absolute count, which can fall while a proportion rises, and a cell type
"enriched" in the mutant while fewer of its cells were recovered is a
compositional artifact rather than an expansion.

Both panels exist so a reader can check a proportion claim without going back
to another table. Counts below --min_cells are marked, because a proportion
from a handful of cells is noise whatever the test says.

There is no robust/all split: counts are the denominators, and filtering them
would defeat the purpose.

Usage:
  python plot_composition_cell_counts.py \
      --input results/rescue/composition/composition_cell_counts.tsv \
      --outdir figures/fig3
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (COL, LOW_CONFIDENCE, WIDE, apply_style, ct_short, note,
                      panel_letter, sample_colours, sample_label, save)
import matplotlib.pyplot as plt

META = ("sample", "genotype", "infection_status", "group", "replicate", "condition")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig3_composition_cell_counts")
    p.add_argument("--min_cells", type=int, default=20,
                   help="below this a per-sample proportion is noise")
    p.add_argument("--drop_labels", default=",".join(LOW_CONFIDENCE))
    return p.parse_args()


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t", index_col=0)
    d = d.drop(columns=[c for c in d.columns if c in META], errors="ignore")
    C = d.apply(pd.to_numeric, errors="coerce").fillna(0)
    C.index = C.index.astype(str)
    drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}
    C = C[[c for c in C.columns if c not in drop]]
    order = C.sum(axis=0).sort_values().index
    C = C[order]
    print(f"{C.shape[0]} samples x {C.shape[1]} cell types; "
          f"{int(C.to_numpy().sum()):,} cells total")

    fig = plt.figure(figsize=(WIDE, max(2.8, 0.24 * len(order) + 1.6)))
    gs = fig.add_gridspec(1, 2, wspace=0.12, width_ratios=[1.0, 1.35])

    # One colour per library, used in both panels, so a point in B can be traced
    # to its bar in A. Panel A is therefore the key and no legend is needed.
    scol, sorder = sample_colours(C.index)
    sorder = [s for s in sorder if s in set(C.index.astype(str))]

    # ---- A: total yield per sample ----
    ax = fig.add_subplot(gs[0, 0])
    tot = C.sum(axis=1)
    tot.index = tot.index.astype(str)
    tot = tot.reindex(sorder)
    ys = np.arange(len(tot))
    ax.barh(ys, tot, height=0.68, zorder=3,
            color=[scol[s] for s in tot.index],
            edgecolor=COL["black"], linewidth=0.4)
    ax.set_yticks(ys, [sample_label(s) for s in tot.index], fontsize=6)
    ax.set_xlabel("cells recovered")
    ax.set_title("Library yield", pad=4)
    lo, hi = tot.min(), tot.max()
    note(ax, f"{hi / max(lo, 1):.1f}x, largest to smallest library", y=-0.12)
    panel_letter(ax, "A", dx=-0.40)

    # ---- B: counts per cell type per sample ----
    ax = fig.add_subplot(gs[0, 1])
    y = np.arange(len(order))
    for i, s in enumerate(sorder):
        v = C.loc[s].to_numpy(dtype=float)
        thin = v < args.min_cells
        j = (np.random.default_rng(i).random(len(v)) - 0.5) * 0.3
        # Colour now carries the sample, so a thin count cannot also be grey.
        # It is marked by an open marker instead, which leaves both readable.
        for m, fc, ec in ((~thin, scol[s], "none"), (thin, "none", scol[s])):
            if m.any():
                ax.scatter(np.clip(v[m], 0.5, None), (y + j)[m], s=12, zorder=3,
                           facecolors=fc, edgecolors=ec, linewidths=0.6)
    ax.axvline(args.min_cells, color=COL["vermillion"], lw=0.9, ls="--", zorder=4)
    ax.set_xscale("log")
    ax.set_yticks(y, ct_short(order), fontsize=6)
    ax.set_xlabel("cells per sample")
    ax.set_title("Counts behind each proportion", pad=4)
    note(ax, f"colour is the sample, keyed by the bars in A; open points are "
             f"below {args.min_cells} cells", y=-0.12)
    panel_letter(ax, "B", dx=-0.10)

    save(fig, args.outdir, args.stem)

    thin_frac = (C < args.min_cells).mean(axis=0)
    bad = thin_frac[thin_frac > 0.5]
    print(C.T.astype(int).to_string())
    if len(bad):
        print(f"\n{len(bad)} cell types are below {args.min_cells} cells in more "
              f"than half the samples: {', '.join(ct_short(bad.index))}. "
              "Composition contrasts in these rest on counts too small to "
              "interpret, whatever their p-value.")


if __name__ == "__main__":
    main()

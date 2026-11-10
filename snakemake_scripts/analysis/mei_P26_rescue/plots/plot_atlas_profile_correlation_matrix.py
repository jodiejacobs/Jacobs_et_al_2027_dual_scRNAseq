#!/usr/bin/env python
"""Supplementary: the full query-by-atlas profile correlation matrix.

Plots results/rescue/atlas_fidelity/atlas_profile_correlation_matrix.tsv.

The per-label summary says whether each label's best match is itself. This
shows the whole matrix, which answers the question the summary cannot: when a
label is confusable, what is it confusable with? A bright diagonal means the
transfer is clean. An off-diagonal block means a group of labels are
interchangeable, and that is the situation where a per-cell-type result can be
an artifact of which of several similar labels a cell happened to receive.

Rows are ordered to put each query label next to its best atlas match, so a
clean transfer appears as a diagonal and the blocks are readable. The diagonal
is outlined rather than masked, because hiding it would make every transfer
look equally good.

There is no robust/all split: the matrix is the diagnostic, and filtering it
would remove the confusions it exists to show.

Usage:
  python plot_atlas_profile_correlation_matrix.py \
      --input results/rescue/atlas_fidelity/atlas_profile_correlation_matrix.tsv \
      --outdir figures/supp
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import COL, WIDE, apply_style, ct_short, note, panel_letter, save
import matplotlib.patches as mpatches
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="supp_atlas_correlation_matrix")
    p.add_argument("--margin", type=float, default=0.02,
                   help="a non-diagonal match within this of the diagonal counts "
                        "as a confusable pair")
    return p.parse_args()


def main():
    args = parse_args()
    apply_style()
    C = pd.read_csv(args.input, sep="\t", index_col=0)
    C = C.apply(pd.to_numeric, errors="coerce")
    print(f"{C.shape[0]} query labels x {C.shape[1]} atlas types")

    # Order query rows by their best atlas match, so a clean transfer is a diagonal.
    best = C.idxmax(axis=1)
    acols = list(C.columns)
    rows = sorted(C.index, key=lambda r: (acols.index(best[r]) if best[r] in acols
                                          else len(acols), str(r)))
    C = C.loc[rows]

    h = max(2.6, 0.17 * len(rows) + 1.5)
    fig, ax = plt.subplots(figsize=(WIDE, h))
    im = ax.imshow(C.to_numpy(), cmap="viridis", aspect="auto",
                   vmin=float(np.nanmin(C.to_numpy())), vmax=1.0)
    for i, r in enumerate(C.index):
        if r in C.columns:
            j = list(C.columns).index(r)
            ax.add_patch(mpatches.Rectangle((j - 0.5, i - 0.5), 1, 1, fill=False,
                                            edgecolor=COL["vermillion"], lw=1.1))
    ax.set_xticks(range(C.shape[1]), ct_short(C.columns), rotation=45,
                  ha="right", fontsize=5.5)
    ax.set_yticks(range(C.shape[0]), ct_short(C.index), fontsize=5.5)
    ax.set_xlabel("atlas cell type")
    ax.set_ylabel("transferred label")
    cb = fig.colorbar(im, ax=ax, fraction=0.022, pad=0.012)
    cb.set_label("profile correlation", fontsize=6)
    cb.ax.tick_params(labelsize=5.5)
    ax.set_title("Query profiles against every atlas type", pad=5)
    note(ax, "red boxes mark the diagonal: a label matching its own atlas type",
         y=-0.06)
    panel_letter(ax, "A", dx=-0.30, dy=1.02)
    save(fig, args.outdir, args.stem)

    conf = []
    for r in C.index:
        if r not in C.columns:
            continue
        self_r = C.at[r, r]
        others = C.loc[r].drop(labels=[r])
        if len(others) and others.max() >= self_r - args.margin:
            conf.append((r, others.idxmax(), self_r, others.max()))
    if conf:
        print(f"\n{len(conf)} confusable label pairs (best non-self match within "
              f"{args.margin} of the diagonal):")
        for a, b, sa, sb in conf:
            print(f"  {ct_short([a])[0]} ~ {ct_short([b])[0]} "
                  f"({sa:.3f} vs {sb:.3f})")
        print("Do not make claims that rest on distinguishing a confusable pair; "
              "pool them or say the labels cannot separate them.")
    else:
        print("\nEvery label is a clear best match for its own atlas type.")


if __name__ == "__main__":
    main()

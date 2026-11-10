#!/usr/bin/env python
"""Supplementary: do transferred labels match the atlas expression profile?

Plots results/rescue/atlas_fidelity/atlas_profile_correlations.tsv.

Composition agreement says the right proportions were recovered. It does not
say the right cells got the right names, and that is the assumption every
per-cell-type result depends on. The test here is stricter: for each
transferred label, does the mean expression profile of the cells carrying it
correlate best with its own atlas type, or with some other one? A label whose
profile matches a different atlas type better is a label that should not be
used, however plausible its name.

The panel reports both the self-correlation and whether the diagonal wins,
because they are different failures. A low self-correlation on a label that
still wins its diagonal is a noisy but usable label. A high self-correlation
that loses to a neighbour is a confusable pair, and claims that depend on
telling those two apart are not supported.

Two versions are written. "robust" keeps labels whose own atlas type is their
best match; "all" keeps every label, so the failures stay visible.

Usage:
  python plot_atlas_profile_correlations.py \
      --input results/rescue/atlas_fidelity/atlas_profile_correlations.tsv \
      --outdir figures/supp
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import empty_figure, COL, WIDE, apply_style, ct_short, note, panel_letter, save
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="supp_atlas_profile_correlations")
    p.add_argument("--min_r", type=float, default=0.5,
                   help="self-correlation below this is a weak label")
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def truthy(s):
    return s.astype(str).str.upper().isin(("TRUE", "T", "1", "YES"))


def build(d, args, suffix, note_txt):
    d = d.sort_values("self_correlation")
    y = np.arange(len(d))
    diag = truthy(d["diagonal_is_best"]) if "diagonal_is_best" in d.columns \
        else pd.Series(True, index=d.index)
    weak = d["self_correlation"] < args.min_r
    print(f"\n[{suffix}] {len(d)} labels, {int((~diag).sum())} match another atlas "
          f"type better, {int(weak.sum())} below r = {args.min_r}; {note_txt}")

    fig = plt.figure(figsize=(WIDE, 3.0))
    gs = fig.add_gridspec(1, 2, wspace=0.12, width_ratios=[1.3, 1.0])

    # ---- A: self vs best correlation ----
    ax = fig.add_subplot(gs[0, 0])
    if "best_correlation" in d.columns:
        for yi, (_, r) in zip(y, d.iterrows()):
            ax.plot([r["self_correlation"], r["best_correlation"]], [yi, yi],
                    color=COL["lightgrey"], lw=1.6, zorder=2)
        ax.scatter(d["best_correlation"], y, s=16, marker="s", zorder=3,
                   color=COL["grey"], label="best match")
    ax.scatter(d["self_correlation"], y, s=26, zorder=4,
               color=[COL["blue"] if ok else COL["vermillion"] for ok in diag],
               edgecolors=COL["black"], linewidths=0.4, label="own atlas type")
    ax.axvline(args.min_r, color=COL["vermillion"], lw=0.8, ls="--", zorder=1)
    ax.set_yticks(y, ct_short(d["cell_type"]), fontsize=6)
    ax.set_xlabel("profile correlation with atlas")
    ax.set_title("Does the label match its own atlas type?", pad=4)
    ax.legend(loc="lower left", fontsize=6)
    note(ax, "red points match a different atlas type better", y=-0.14)
    panel_letter(ax, "A", dx=-0.54)

    # ---- B: how many cells each label rests on ----
    ax = fig.add_subplot(gs[0, 1])
    if "n_cells_query" in d.columns:
        ax.barh(y, pd.to_numeric(d["n_cells_query"], errors="coerce"), height=0.68,
                zorder=3,
                color=[COL["blue"] if ok else COL["vermillion"] for ok in diag])
        ax.set_xscale("log")
    ax.set_yticks(y, [])
    ax.set_xlabel("cells with this label")
    ax.set_title("Label size", pad=4)
    panel_letter(ax, "B", dx=-0.06)

    fig.suptitle(note_txt, y=1.03, fontsize=6.5, color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")

    cols = [c for c in ("cell_type", "n_cells_query", "self_correlation",
                        "best_match", "best_correlation", "diagonal_is_best")
            if c in d.columns]
    print(d[cols].round(3).to_string(index=False))
    bad = d[~diag]
    for _, r in bad.iterrows():
        print(f"  {ct_short([r['cell_type']])[0]} matches "
              f"{ct_short([r.get('best_match')])[0]} better "
              f"({r.get('best_correlation', float('nan')):.3f} vs "
              f"{r['self_correlation']:.3f})")
    if len(bad):
        print("Treat these as confusable with their best match. A claim that "
              "depends on telling the pair apart is not supported by the labels.")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    for c in ("self_correlation", "best_correlation", "n_cells_query"):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d[d["self_correlation"].notna()]
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        if v == "robust" and "diagonal_is_best" in d.columns:
            sub = d[truthy(d["diagonal_is_best"])]
            txt = "labels whose own atlas type is their best expression match"
        else:
            sub = d
            txt = "every label with a self-correlation, best match or not"
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

#!/usr/bin/env python
"""Figure 2: whether wMel locus detection is biased by sequence features.

Plots results/rescue/detection_bias/wmel_feature_tests.tsv.

The figure exists to show the difference between the raw and the conditional
test, because they support opposite claims. Longer and more highly expressed
loci are detected in more cells, and expression and length are themselves
correlated, so a raw Spearman of detection against length mostly reports
expression. The partial coefficient holds expression rank fixed and asks
whether the feature adds anything. Plotting only the raw values would invent a
bias; plotting only the partial would hide why the question arose. Both are
drawn, paired per feature.

The two versions are the two tests: "all" is the raw Spearman, "robust" is the
coefficient conditional on expression, which is the one to quote.

Usage:
  python plot_wmel_feature_tests.py \
      --input results/rescue/detection_bias/wmel_feature_tests.tsv \
      --outdir figures/fig2
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import COL, COLUMN, apply_style, fdr_stars, note, panel_letter, save
import matplotlib.pyplot as plt

NICE = {"log10_length": "locus length\n(log$_{10}$)", "gc": "GC content",
        "rel_position": "position on\nthe chromosome", "plus_strand": "plus strand"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig2_wmel_feature_tests")
    p.add_argument("--alpha", type=float, default=0.05)
    return p.parse_args()


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    for c in ("spearman_raw", "pvalue_raw", "partial_spearman_given_cpm",
              "pvalue_partial", "n"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    d["label"] = [NICE.get(f, str(f).replace("_", " ")) for f in d["feature"]]
    d = d.reindex(d["partial_spearman_given_cpm"].abs().sort_values().index)
    print(f"{len(d)} features, n = {int(d['n'].max())} loci")

    fig, ax = plt.subplots(figsize=(COLUMN * 1.35, 2.4))
    y = np.arange(len(d))
    ax.barh(y - 0.19, d["spearman_raw"], height=0.34, color=COL["lightgrey"],
            label="raw", zorder=3)
    ax.barh(y + 0.19, d["partial_spearman_given_cpm"], height=0.34,
            color=[COL["vermillion"] if p < args.alpha else COL["blue"]
                   for p in d["pvalue_partial"]],
            label="conditional on expression", zorder=3)
    for yi, (_, r) in zip(y, d.iterrows()):
        v = r["partial_spearman_given_cpm"]
        s = fdr_stars(r["pvalue_partial"]) or "n.s."
        ax.text(v + (0.012 if v >= 0 else -0.012), yi + 0.19, s, fontsize=6,
                va="center", ha="left" if v >= 0 else "right")
    ax.axvline(0, color=COL["black"], lw=0.6)
    ax.set_yticks(y, d["label"], fontsize=6.5)
    ax.set_xlabel("Spearman with detection fraction")
    ax.set_title("Detection against sequence features", pad=4)
    ax.legend(loc="lower right", fontsize=6)
    ax.margins(x=0.22)
    note(ax, "stars on the conditional test only; raw values are shown to make "
             "the confound visible, not to be quoted", y=-0.26)
    panel_letter(ax, "A", dx=-0.30)
    save(fig, args.outdir, args.stem)

    print(d[["feature", "n", "spearman_raw", "pvalue_raw",
             "partial_spearman_given_cpm", "pvalue_partial"]]
          .round(4).to_string(index=False))
    bad = d.loc[d["pvalue_partial"] < args.alpha, "feature"].tolist()
    if bad:
        print(f"\n{', '.join(bad)} still predict detection after conditioning on "
              "expression, so recovery is not uniform across loci. Report this "
              "as a limitation and do not treat per-locus detection as "
              "interchangeable.")
    else:
        print("\nNo feature predicts detection once expression is held fixed, so "
              "the raw associations are expression in disguise.")


if __name__ == "__main__":
    main()

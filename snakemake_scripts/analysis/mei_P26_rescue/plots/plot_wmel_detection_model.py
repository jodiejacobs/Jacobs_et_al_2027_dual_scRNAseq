#!/usr/bin/env python
"""Figure 2: the joint model of wMel locus detection.

Plots results/rescue/detection_bias/wmel_detection_model.tsv, the OLS of
logit(detection fraction) on expression and sequence features.

The pairwise tests in the companion panel ask one feature at a time, which
leaves open whether two correlated features are reporting the same thing. This
model fits them together, so a coefficient here is the feature's contribution
with the others held fixed. Coefficients are plotted with their standard
errors rather than as bare point estimates, because the point of the panel is
whether a term is distinguishable from zero and a bar chart hides that.

The intercept is dropped: it is the logit detection of a locus at zero
expression and zero length, which is not a quantity anyone interprets.

There is no robust/all split here: a single model is fitted, and the error
bars are the uncertainty.

Usage:
  python plot_wmel_detection_model.py \
      --input results/rescue/detection_bias/wmel_detection_model.tsv \
      --outdir figures/fig2
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import COL, COLUMN, apply_style, note, panel_letter, save
import matplotlib.pyplot as plt

NICE = {"log10_cpm": "expression (log$_{10}$ CPM)", "log10_length": "length (log$_{10}$)",
        "gc": "GC content", "rel_position": "chromosomal position",
        "plus_strand": "plus strand", "const": "intercept"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig2_wmel_detection_model")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--keep_intercept", action="store_true")
    return p.parse_args()


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    for c in ("coef", "se", "t", "pvalue"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    if not args.keep_intercept:
        d = d[~d["term"].isin(("const", "Intercept"))]
    d["label"] = [NICE.get(t, str(t).replace("_", " ")) for t in d["term"]]
    d = d.reindex(d["coef"].abs().sort_values().index)
    sig = d["pvalue"] < args.alpha
    print(f"{len(d)} terms, {int(sig.sum())} with p < {args.alpha}")

    fig, ax = plt.subplots(figsize=(COLUMN * 1.35, 2.3))
    y = np.arange(len(d))
    ax.errorbar(d["coef"], y, xerr=1.96 * d["se"], fmt="none",
                ecolor=COL["grey"], elinewidth=0.9, capsize=1.6, zorder=2)
    ax.scatter(d["coef"], y, s=26, zorder=3,
               color=[COL["vermillion"] if s else COL["lightgrey"] for s in sig],
               edgecolors=COL["black"], linewidths=0.4)
    ax.axvline(0, color=COL["black"], lw=0.7)
    ax.set_yticks(y, d["label"], fontsize=6.5)
    ax.set_xlabel("coefficient on logit(detection fraction)")
    ax.set_title("Detection model, terms held jointly", pad=4)
    ax.margins(x=0.20, y=0.12)
    note(ax, f"bars are 95% CI; filled points p < {args.alpha}", y=-0.24)
    panel_letter(ax, "B", dx=-0.30)
    save(fig, args.outdir, args.stem)

    print(d[["term", "coef", "se", "t", "pvalue"]].round(4).to_string(index=False))
    keep = d.loc[sig & ~d["term"].isin(("log10_cpm",)), "term"].tolist()
    if keep:
        print(f"\n{', '.join(keep)} contribute beyond expression, so detection is "
              "structured by sequence and not only by abundance.")
    else:
        print("\nOnly expression predicts detection in the joint model, which is "
              "the cleanest version of this result: no sequence feature adds to it.")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""Figure 4: contrast tests on the germline metrics.

Plots either of the two tables written by 07_germline.py, which share a schema
(metric, contrast, estimate, se, t, pvalue, df):

  results/rescue/germline/germline_composition_tests.tsv
  results/rescue/germline/germline_state_tests.tsv

One script covers both because the tables are the same shape and the panel is
the same; --stem distinguishes the outputs. Splitting it would duplicate the
code without adding a figure.

These tests run on eight per-sample values with four groups, so the residual df
is 4. That is a real constraint and the panel states it rather than letting the
p-values imply more. At df = 4 a confidence interval is wide, and a contrast
that just clears 0.05 should be read as a direction worth following up, not an
established effect.

The interaction is the formal rescue test. The genotype contrasts describe the
defect and what is left of it, and showing them together is what makes the
interaction interpretable.

Two versions are written. "robust" keeps metrics with a significant
interaction; "all" keeps every metric.

Usage:
  python plot_germline_tests.py \
      --input results/rescue/germline/germline_composition_tests.tsv \
      --stem fig4_germline_composition_tests --outdir figures/fig4
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import empty_figure, COL, WIDE, apply_style, fdr_stars, note, panel_letter, save
import matplotlib.pyplot as plt

ORDER = ["genotype_in_uninf", "genotype_in_wMel", "interaction",
         "infection_in_OreR", "infection_in_MeiP26"]
CCOL = {"genotype_in_uninf": COL["vermillion"], "genotype_in_wMel": COL["blue"],
        "interaction": COL["green"], "infection_in_OreR": COL["skyblue"],
        "infection_in_MeiP26": COL["orange"]}
NICE = {"genotype_in_uninf": "the defect", "genotype_in_wMel": "what remains",
        "interaction": "rescue test", "infection_in_OreR": "wMel in OreR",
        "infection_in_MeiP26": "wMel in mei-P26"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig4_germline_tests")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def build(d, args, suffix, note_txt):
    cons = [c for c in ORDER if c in set(d["contrast"])] or \
        sorted(set(d["contrast"]))
    mets = list(dict.fromkeys(d["metric"]))
    df = int(pd.to_numeric(d["df"], errors="coerce").dropna().iloc[0]) \
        if "df" in d.columns and d["df"].notna().any() else None
    print(f"\n[{suffix}] {len(mets)} metrics x {len(cons)} contrasts"
          f"{f', residual df = {df}' if df is not None else ''}; {note_txt}")

    fig, ax = plt.subplots(figsize=(WIDE * 0.8,
                                    max(2.4, 0.42 * len(mets) + 1.2)))
    n = len(cons)
    y = np.arange(len(mets))
    for i, c in enumerate(cons):
        s = d[d["contrast"] == c].set_index("metric").reindex(mets)
        off = (i - (n - 1) / 2) * (0.74 / n)
        ax.errorbar(s["estimate"], y + off, xerr=1.96 * s["se"], fmt="none",
                    ecolor=COL["lightgrey"], elinewidth=0.8, capsize=1.2, zorder=2)
        sig = s["pvalue"] < args.alpha
        ax.scatter(s["estimate"], y + off, s=22, zorder=3,
                   color=[CCOL.get(c, COL["grey"]) if k else "white"
                          for k in sig.fillna(False)],
                   edgecolors=CCOL.get(c, COL["grey"]), linewidths=0.8,
                   label=NICE.get(c, c))
        for yi, (_, r) in zip(y, s.iterrows()):
            st = fdr_stars(r["pvalue"])
            if st:
                ax.text(r["estimate"], yi + off + 0.10, st, ha="center",
                        fontsize=5.5, color=CCOL.get(c, COL["grey"]))
    ax.axvline(0, color=COL["black"], lw=0.7)
    ax.set_yticks(y, [str(m).replace("_", " ") for m in mets], fontsize=6.5)
    ax.set_xlabel("contrast estimate (logit where a fraction)")
    ax.set_title("Germline contrasts", pad=4)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.02), ncol=min(n, 3),
              fontsize=6)
    tail = f"; residual df = {df}, so intervals are wide" if df is not None else ""
    note(ax, f"filled points p < {args.alpha}; bars 95% CI{tail}", y=-0.10)
    panel_letter(ax, "A", dx=-0.34, dy=1.12)
    save(fig, args.outdir, f"{args.stem}_{suffix}")

    cols = [c for c in ("metric", "contrast", "estimate", "se", "t", "pvalue", "df")
            if c in d.columns]
    print(d[cols].round(4).to_string(index=False))
    if df is not None and df <= 4:
        print(f"\nResidual df = {df}. A contrast near p = {args.alpha} here is a "
              "direction worth following up, not an established effect. Report "
              "the estimate and its interval rather than the p-value alone.")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    for c in ("estimate", "se", "t", "pvalue", "df"):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        if v == "robust":
            keep = set(d.loc[(d["contrast"] == "interaction")
                             & (d["pvalue"] < args.alpha), "metric"])
            sub = d[d["metric"].isin(keep)]
            txt = "metrics with a significant interaction, the formal rescue test"
        else:
            sub = d
            txt = "every metric tested"
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

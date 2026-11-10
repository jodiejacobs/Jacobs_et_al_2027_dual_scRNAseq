#!/usr/bin/env python
"""Figure 5: how many host genes respond, by contrast and cell type.

Plots results/rescue/de/host_de_summary.tsv.

This is the map of where anything happens at all, and it is read before any
individual gene. The mutant effect in uninfected animals sets the ceiling: a
cell type with no mutant DE genes has no defect to rescue, so a rescue claim
there is about noise whatever the downstream categories say. Plotting the
contrasts together keeps that comparison in view.

Counts are plotted on a log axis with the raw number annotated, because they
span orders of magnitude between germline and small somatic populations and a
linear axis would show only the largest cell type.

Two versions are written. "robust" keeps cell types with at least --min_de
genes in the mutant contrast, which is where downstream fractions mean
something; "all" keeps every cell type.

Usage:
  python plot_host_de_summary.py \
      --input results/rescue/de/host_de_summary.tsv \
      --outdir figures/fig5
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (empty_figure, COL, LOW_CONFIDENCE, WIDE, apply_style, ct_short, note,
                      panel_letter, save)
import matplotlib.pyplot as plt

ORDER = ["genotype_in_uninf", "genotype_in_wMel", "interaction",
         "infection_in_OreR", "infection_in_MeiP26", "infection_main"]
CCOL = {"genotype_in_uninf": COL["vermillion"], "genotype_in_wMel": COL["blue"],
        "interaction": COL["green"], "infection_in_OreR": COL["skyblue"],
        "infection_in_MeiP26": COL["orange"], "infection_main": COL["purple"]}
NICE = {"genotype_in_uninf": "the defect", "genotype_in_wMel": "what remains",
        "interaction": "rescue test", "infection_in_OreR": "wMel in OreR",
        "infection_in_MeiP26": "wMel in mei-P26", "infection_main": "wMel overall"}
REF = "genotype_in_uninf"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig5_host_de_summary")
    p.add_argument("--min_de", type=int, default=50)
    p.add_argument("--count_col", default=None,
                   help="column holding the gene count; detected if omitted")
    p.add_argument("--drop_labels", default=",".join(LOW_CONFIDENCE))
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def find_count_col(d):
    for c in ("n_sig", "n_de", "n_genes", "n", "count"):
        if c in d.columns:
            return c
    num = [c for c in d.columns if pd.api.types.is_numeric_dtype(d[c])]
    if not num:
        raise SystemExit(f"no numeric count column in {list(d.columns)}; "
                         "pass --count_col")
    return num[0]


def build(d, args, ncol, suffix, note_txt):
    cons = [c for c in ORDER if c in set(d["contrast"])] or \
        sorted(set(d["contrast"]))
    m = d.pivot_table(index="tissue", columns="contrast", values=ncol,
                      aggfunc="sum").reindex(columns=cons).fillna(0)
    key = REF if REF in m.columns else cons[0]
    m = m.sort_values(key)
    y = np.arange(len(m))
    print(f"\n[{suffix}] {len(m)} cell types x {len(cons)} contrasts; {note_txt}")

    fig, ax = plt.subplots(figsize=(WIDE, max(2.6, 0.26 * len(m) + 1.4)))
    n = len(cons)
    for i, c in enumerate(cons):
        off = (i - (n - 1) / 2) * (0.78 / n)
        ax.barh(y + off, m[c].clip(lower=0.4), height=0.74 / n,
                color=CCOL.get(c, COL["grey"]), label=NICE.get(c, c), zorder=3)
    ax.set_xscale("log")
    ax.axvline(args.min_de, color=COL["black"], lw=0.8, ls="--", zorder=4)
    ax.set_yticks(y, ct_short(m.index), fontsize=6)
    ax.set_xlabel("significant genes")
    ax.set_title("Where anything happens at all", pad=4)
    ax.legend(loc="lower left", bbox_to_anchor=(0, 1.02), ncol=min(n, 3),
              fontsize=6)
    note(ax, f"dashed line at {args.min_de} genes; bars at zero are drawn at 0.4 "
             "so the log axis can show them", y=-0.10)
    panel_letter(ax, "A", dx=-0.34, dy=1.10)
    save(fig, args.outdir, f"{args.stem}_{suffix}")

    print(m.astype(int).to_string())
    if key in m.columns:
        thin = m.index[m[key] < args.min_de]
        if len(thin):
            print(f"\n{', '.join(ct_short(thin))} have fewer than {args.min_de} "
                  "genes DE in the mutant, so there is little defect to rescue. "
                  "A rescue fraction in these is not interpretable.")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    if "tissue" not in d.columns or "contrast" not in d.columns:
        raise SystemExit(f"expected tissue and contrast columns, got {list(d.columns)}")
    ncol = args.count_col or find_count_col(d)
    d[ncol] = pd.to_numeric(d[ncol], errors="coerce")
    drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}
    d = d[~d["tissue"].isin(drop)]
    print(f"Counting column: {ncol}; labels excluded: {sorted(drop) or 'none'}")
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        if v == "robust":
            keep = set(d.loc[(d["contrast"] == REF) & (d[ncol] >= args.min_de),
                             "tissue"])
            sub = d[d["tissue"].isin(keep)]
            txt = f"cell types with at least {args.min_de} genes DE in the mutant"
        else:
            sub = d
            txt = "every cell type tested"
        if not len(sub):
            empty_figure(args.outdir, f"{args.stem}_{v}",
                         f"No rows survive the {v} filter.\n"
                         f"An empty set here is a result, not a failure:\n"
                         f"say so in the text rather than citing the\n"
                         f"unfiltered figure as though it passed.")
            continue
        build(sub, args, ncol, v, txt)


if __name__ == "__main__":
    main()

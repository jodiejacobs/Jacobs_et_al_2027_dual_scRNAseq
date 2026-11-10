#!/usr/bin/env python
"""Figure 6: how many wMel loci respond, by cell type and contrast.

Plots results/rescue/de/wolbachia_de_summary.tsv.

The gene-level panels show which loci respond. This one shows how much of the
genome responds at all, which is the context those panels need: a handful of
significant loci out of a thousand tested is a different claim from a broad
transcriptional shift, and only the second supports language about a
coordinated bacterial response.

The number tested is plotted beside the number significant, because they vary
together. A cell type where few loci pass the expression filter will show few
significant ones for that reason alone, and reading the significant count
without its denominator turns a coverage difference into a biological one.

Two versions are written. "robust" keeps cell types with at least --min_tested
loci tested; "all" keeps every cell type.

Usage:
  python plot_wolbachia_de_summary.py \
      --input results/rescue/de/wolbachia_de_summary.tsv \
      --outdir figures/fig6
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (empty_figure, COL, LOW_CONFIDENCE, NON_CELLTYPE, WIDE, apply_style,
                      ct_short, note, panel_letter, save)
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig6_wolbachia_de_summary")
    p.add_argument("--min_tested", type=int, default=15,
                   help="loci tested; the wMel DE expression filter leaves only "
                        "tens of loci per cell type, not the full genome")
    p.add_argument("--drop_labels", default=",".join(LOW_CONFIDENCE))
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def build(d, args, suffix, note_txt):
    d = d.groupby("tissue", as_index=False)[["n_tested", "n_sig", "n_up", "n_down"]] \
        .sum().sort_values("n_sig")
    y = np.arange(len(d))
    frac = d["n_sig"] / d["n_tested"].replace(0, np.nan)
    print(f"\n[{suffix}] {len(d)} cell types; {note_txt}")

    fig = plt.figure(figsize=(WIDE, 2.8))
    gs = fig.add_gridspec(1, 3, wspace=0.16, width_ratios=[1.3, 1.0, 1.0])

    # ---- A: up and down, stacked ----
    ax = fig.add_subplot(gs[0, 0])
    ax.barh(y, d["n_up"], height=0.68, color=COL["vermillion"],
            label="up in mei-P26", zorder=3)
    ax.barh(y, d["n_down"], left=d["n_up"], height=0.68, color=COL["blue"],
            label="down", zorder=3)
    ax.set_yticks(y, ct_short(d["tissue"]), fontsize=6)
    ax.set_xlabel("significant loci")
    ax.set_title("How many wMel loci respond", pad=4)
    ax.legend(loc="lower right", fontsize=6)
    panel_letter(ax, "A", dx=-0.50)

    # ---- B: the denominator ----
    ax = fig.add_subplot(gs[0, 1])
    ax.barh(y, d["n_tested"], height=0.68, color=COL["grey"], zorder=3)
    ax.axvline(args.min_tested, color=COL["vermillion"], lw=0.9, ls="--")
    ax.set_yticks(y, [])
    ax.set_xlabel("loci tested")
    ax.set_title("How many were testable", pad=4)
    panel_letter(ax, "B", dx=-0.08)

    # ---- C: the fraction, which is the claim ----
    ax = fig.add_subplot(gs[0, 2])
    ax.barh(y, frac, height=0.68, color=COL["green"], zorder=3)
    ax.set_yticks(y, [])
    ax.set_xlabel("significant / tested")
    ax.set_title("Share of the genome", pad=4)
    note(ax, "a broad response is a fraction, not a count", y=-0.18)
    panel_letter(ax, "C", dx=-0.08)

    fig.suptitle(note_txt, y=1.04, fontsize=6.5, color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")

    out = d.assign(frac_sig=frac.round(4))
    print(out.to_string(index=False))
    top = out.nlargest(1, "frac_sig")
    if len(top) and top["frac_sig"].iloc[0] < 0.05:
        print(f"\nNo cell type has more than 5% of tested loci responding "
              f"(max {top['frac_sig'].iloc[0]:.1%} in "
              f"{ct_short(top['tissue'])[0]}). Describe this as a targeted "
              "response in specific loci, not a coordinated genome-wide shift.")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    for c in ("n_tested", "n_sig", "n_up", "n_down"):
        if c not in d.columns:
            raise SystemExit(f"expected {c} in {list(d.columns)}")
        d[c] = pd.to_numeric(d[c], errors="coerce").fillna(0)
    drop = ({s.strip() for s in args.drop_labels.split(",") if s.strip()}
            | set(NON_CELLTYPE))
    d = d[~d["tissue"].isin(drop)]
    print(f"Excluded as not a cell type or low-confidence: {sorted(drop)}")
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        sub = d[d["n_tested"] >= args.min_tested] if v == "robust" else d
        txt = (f"cell types with at least {args.min_tested} loci tested"
               if v == "robust" else "every cell type, however few loci testable")
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

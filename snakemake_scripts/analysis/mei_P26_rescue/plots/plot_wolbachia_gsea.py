#!/usr/bin/env python
"""Supplementary: functional enrichment in the wMel response.

Plots results/rescue/wolbachia_followup/wolbachia_gsea.tsv.

The table mixes two kinds of test under a `family` column and they are not
comparable, so the panel separates them rather than ranking everything
together. For the gseapy families NES is a normalised enrichment score. For
the `category` family it is the median logFC difference between loci whose
product matches a keyword and the rest, which is a crude test on an annotation
that is mostly "hypothetical protein". Plotting them on one axis would imply
an equivalence that does not hold.

Set size is drawn alongside because a wMel category can reach significance on
three loci, and a three-gene enrichment is a description of three genes.

Two versions are written. "robust" keeps terms passing FDR with at least
--min_genes loci; "all" keeps every tested term.

Usage:
  python plot_wolbachia_gsea.py \
      --input results/rescue/wolbachia_followup/wolbachia_gsea.tsv \
      --outdir figures/supp
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import empty_figure, COL, WIDE, apply_style, note, panel_letter, save, shorten
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="supp_wolbachia_gsea")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--min_genes", type=int, default=5)
    p.add_argument("--top_n", type=int, default=10, help="terms shown per family")
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def build(d, args, suffix, note_txt):
    fams = [f for f in d["family"].dropna().unique()]
    print(f"\n[{suffix}] {len(d)} terms across {len(fams)} families; {note_txt}")
    if not fams:
        return

    fig = plt.figure(figsize=(WIDE, max(2.4, 0.22 * min(len(d), args.top_n *
                                                        len(fams)) + 1.4)))
    gs = fig.add_gridspec(1, len(fams), wspace=0.62)

    for k, fam in enumerate(fams):
        ax = fig.add_subplot(gs[0, k])
        s = d[d["family"] == fam]
        s = s.reindex(s["NES"].abs().sort_values().tail(args.top_n).index)
        y = np.arange(len(s))
        sig = s["FDR"] < args.fdr
        thin = s["n_genes"] < args.min_genes
        ax.barh(y, s["NES"], height=0.68, zorder=3,
                color=[COL["lightgrey"] if t else
                       (COL["vermillion"] if (v > 0 and g) else
                        (COL["blue"] if g else COL["lightgrey"]))
                       for v, g, t in zip(s["NES"], sig, thin)])
        for yi, (_, r) in zip(y, s.iterrows()):
            ax.text(r["NES"] + (0.02 if r["NES"] >= 0 else -0.02), yi,
                    f"{int(r['n_genes'])}", fontsize=5.5, va="center",
                    ha="left" if r["NES"] >= 0 else "right")
        ax.axvline(0, color=COL["black"], lw=0.7)
        ax.set_yticks(y, shorten(s["Term"], 30), fontsize=5.5)
        ax.set_xlabel("median Δ log$_2$FC" if fam == "category" else "NES",
                      fontsize=7)
        ax.set_title(str(fam), pad=4, fontsize=7.5)
        ax.margins(x=0.26)
        panel_letter(ax, "ABCDEF"[k], dx=-0.62)

    fig.suptitle(note_txt + "  ·  numbers are set sizes", y=1.03, fontsize=6.5,
                 color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")

    cols = [c for c in ("family", "Term", "n_genes", "NES", "p", "FDR")
            if c in d.columns]
    print(d[d["FDR"] < args.fdr][cols].round(4).to_string(index=False))
    small = d[(d["FDR"] < args.fdr) & (d["n_genes"] < args.min_genes)]
    if len(small):
        print(f"\n{len(small)} significant terms rest on fewer than "
              f"{args.min_genes} loci. Name the loci instead of citing the term.")
    if "category" in set(d["family"]):
        print("The category family's NES is a median logFC difference, not a "
              "normalised enrichment score, and it is computed on an annotation "
              "that is largely 'hypothetical protein'. Do not report it as GSEA.")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    for c in ("n_genes", "NES", "p", "FDR", "FDR_gseapy"):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d[d["NES"].notna()]
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        if v == "robust":
            sub = d[(d["FDR"] < args.fdr) & (d["n_genes"] >= args.min_genes)]
            txt = (f"terms at FDR < {args.fdr} resting on at least "
                   f"{args.min_genes} loci")
        else:
            sub = d
            txt = "every tested term, significant or not"
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

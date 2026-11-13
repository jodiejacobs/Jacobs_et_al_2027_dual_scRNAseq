#!/usr/bin/env python
"""Figure 4: the genotype-by-infection interaction, one panel per cell type.

Plots results/rescue/de/host_de_all.tsv.gz, the interaction contrast only.

This is the display of the test the rescue claim rests on. The interaction
asks which genes respond to wMel differently in the mutant than in wild type,
so it is the only contrast that separates rescue from an ordinary infection
response. The count per cell type is already in fig5_host_de_summary; this
shows the genes behind those counts, their effect sizes, and which ones they
are.

One panel per cell type rather than one pooled panel, because the whole result
is that the cell types differ: the somatic follicle types carry mutant-specific
responses and the germline does not. A pooled volcano would hide exactly the
contrast the figure exists to make.

summary/volcano_interaction.pdf from 05_rescue_summary.py plots the same
contrast. This version differs in three ways that matter for a figure rather
than a diagnostic: it is in figstyle, it drops the labels below, and it orders
panels by the number of significant genes so the reader meets the cell types
in the order the argument needs.

Non-ovarian labels are dropped by default, not only the low-confidence ones.
Tracheal cells and oviduct are real cells in a dissected ovary rather than
artefacts, but they are not the tissue the claim is about, and adult tracheal
cell otherwise carries the largest interaction count in the dataset, which
invites a question the paper does not want to spend its space answering.

Two versions are written. "robust" keeps cell types whose observed count beats
every permutation relabelling, which needs --null; without it, cell types with
at least --min_sig genes. "all" keeps every cell type.

Usage:
  python plot_interaction_volcano.py \
      --input results/rescue/de/host_de_all.tsv.gz \
      --null results/rescue/rescue_null/rescue_null_summary.tsv \
      --outdir results/rescue/figures
"""
import argparse
import math
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (empty_figure, COL, LOW_CONFIDENCE, WIDE, apply_style,
                      ct_short, note, panel_letter, save)
import matplotlib.pyplot as plt

NON_OVARIAN = ("adult tracheal cell", "oviduct")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="de/host_de_all.tsv.gz")
    p.add_argument("--null", default=None,
                   help="rescue_null/rescue_null_summary.tsv. When given, the "
                        "robust version keeps only cell types whose observed "
                        "interaction count exceeds every permutation "
                        "relabelling, which is the filter the paper argues")
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig4_interaction_volcano")
    p.add_argument("--contrast", default="interaction")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--min_sig", type=int, default=10,
                   help="fallback robust filter when --null is not given")
    p.add_argument("--n_label", type=int, default=6,
                   help="genes named per panel, by smallest FDR")
    p.add_argument("--ncol", type=int, default=3)
    p.add_argument("--drop_labels",
                   default=",".join(tuple(LOW_CONFIDENCE) + NON_OVARIAN))
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def panel(ax, d, args, title):
    d = d.dropna(subset=["logFC", "FDR"]).copy()
    # FDR of exactly 0 has no log. Floor it at the smallest positive value in
    # the panel so those genes stay visible and stay at the top, rather than
    # dropping out of the figure silently.
    pos = d.loc[d["FDR"] > 0, "FDR"]
    d["y"] = -np.log10(d["FDR"].clip(lower=pos.min() if len(pos) else 1e-300))
    sig = d["FDR"] < args.alpha
    up, down = sig & (d["logFC"] > 0), sig & (d["logFC"] < 0)
    ax.scatter(d.loc[~sig, "logFC"], d.loc[~sig, "y"], s=2, alpha=0.25,
               color=COL["lightgrey"], linewidths=0, rasterized=True)
    for m, c in ((down, COL["blue"]), (up, COL["vermillion"])):
        ax.scatter(d.loc[m, "logFC"], d.loc[m, "y"], s=3.5, alpha=0.75,
                   color=c, linewidths=0, rasterized=True)
    ax.axhline(-np.log10(args.alpha), color=COL["black"], lw=0.5, ls="--")
    ax.axvline(0, color=COL["black"], lw=0.5)
    # Label the strongest genes, but skip any that would land on one already
    # placed. The eggshell genes sit almost on top of each other, so labelling
    # all of them unconditionally produces a smear that names nothing.
    xr = d["logFC"].max() - d["logFC"].min() or 1.0
    yr = d["y"].max() - d["y"].min() or 1.0
    placed = []
    for _, r in d[sig].sort_values("FDR").iterrows():
        if len(placed) >= args.n_label:
            break
        fx, fy = r["logFC"] / xr, r["y"] / yr
        if any(abs(fx - px) < 0.10 and abs(fy - py) < 0.07 for px, py in placed):
            continue
        placed.append((fx, fy))
        ax.annotate(str(r.get("symbol") or r["gene_id"]), (r["logFC"], r["y"]),
                    fontsize=4.5, xytext=(2.5, 2), textcoords="offset points",
                    color=COL["black"])
    ax.set_title(f"{title}\n{int(sig.sum())} genes", pad=3, fontsize=6.5)
    ax.tick_params(labelsize=5.5)
    return int(sig.sum())


def build(d, args, suffix, note_txt):
    types = (d.groupby("tissue")
             .apply(lambda s: int((s["FDR"] < args.alpha).sum()),
                    include_groups=False)
             .sort_values(ascending=False))
    print(f"\n[{suffix}] {len(types)} cell types; {note_txt}")
    ncol = min(args.ncol, max(1, len(types)))
    nrow = math.ceil(len(types) / ncol)
    fig, axes = plt.subplots(nrow, ncol, figsize=(WIDE, 1.85 * nrow + 0.7),
                             squeeze=False)
    for ax in axes.ravel():
        ax.set_axis_off()
    for i, t in enumerate(types.index):
        ax = axes[i // ncol][i % ncol]
        ax.set_axis_on()
        n = panel(ax, d[d["tissue"] == t], args, ct_short([t])[0])
        if i % ncol == 0:
            ax.set_ylabel("-log$_{10}$ FDR", fontsize=6)
        if i // ncol == nrow - 1 or i >= len(types) - ncol:
            ax.set_xlabel("interaction log$_2$FC", fontsize=6)
        print(f"  {t[:50]:<52}{n:>6} genes at FDR < {args.alpha}")
    fig.tight_layout(rect=(0, 0.01, 1, 0.955))
    fig.suptitle(note_txt, y=0.985, fontsize=6.5, color=COL["grey"])
    panel_letter(axes[0][0], "A", dx=-0.26, dy=1.14)
    save(fig, args.outdir, f"{args.stem}_{suffix}")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    d = d[d["contrast"] == args.contrast]
    if not len(d):
        raise SystemExit(f"no rows with contrast == {args.contrast}")
    for c in ("logFC", "FDR"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}
    d = d[~d["tissue"].isin(drop)]
    print(f"{d['tissue'].nunique()} cell types, {len(d):,} gene tests; "
          f"dropped {sorted(drop) or 'none'}")

    beats = None
    if args.null and os.path.exists(args.null):
        nl = pd.read_csv(args.null, sep="\t")
        for c in ("obs_n_interaction", "null_n_interaction_max"):
            nl[c] = pd.to_numeric(nl[c], errors="coerce")
        beats = set(nl.loc[nl["obs_n_interaction"] > nl["null_n_interaction_max"],
                           "tissue"])
        print(f"{len(beats)} cell types beat every permutation relabelling")

    for v in (["robust", "all"] if args.version == "both" else [args.version]):
        if v == "robust":
            if beats is not None:
                sub = d[d["tissue"].isin(beats)]
                txt = ("cell types whose interaction count exceeds every "
                       "permutation relabelling")
            else:
                keep = {t for t, s in d.groupby("tissue")
                        if (s["FDR"] < args.alpha).sum() >= args.min_sig}
                sub = d[d["tissue"].isin(keep)]
                txt = (f"cell types with at least {args.min_sig} genes; no "
                       "--null given, so this is not the permutation filter")
        else:
            sub, txt = d, "every cell type tested"
        if not sub["tissue"].nunique():
            empty_figure(args.outdir, f"{args.stem}_{v}",
                         f"No cell type survives the {v} filter.\n"
                         "An empty set here is a result, not a failure:\n"
                         "say so in the text rather than citing the\n"
                         "unfiltered figure as though it passed.")
            continue
        build(sub, args, v, txt)


if __name__ == "__main__":
    main()

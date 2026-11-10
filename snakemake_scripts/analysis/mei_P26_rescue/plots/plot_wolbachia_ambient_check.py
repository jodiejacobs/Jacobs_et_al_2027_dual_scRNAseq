#!/usr/bin/env python
"""Figure 6: is there enough real wMel signal per cell type to run DE in it?

Plots results/rescue/de/wolbachia_ambient_check.tsv.

Wolbachia differential expression is only meaningful in cell types where the
wMel reads come from Wolbachia actually in those cells. The check is the ratio
of infected to uninfected pseudobulk titer in the same cell type: uninfected
animals have no Wolbachia, so their titer is that cell type's own noise, and a
ratio near 1 means the infected signal is indistinguishable from it.

This panel therefore decides which cell types any wMel DE result is allowed to
be reported in, and it should be read before Figure 6's gene panels rather
than as an afterthought. A cell type below the threshold can still produce
significant genes; those genes are differences between noise distributions.

Two versions are written. "robust" keeps cell types clearing --min_ratio;
"all" keeps every cell type, so what was excluded is visible.

Usage:
  python plot_wolbachia_ambient_check.py \
      --input results/rescue/de/wolbachia_ambient_check.tsv \
      --outdir figures/fig6
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


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig6_wolbachia_ambient_check")
    p.add_argument("--min_ratio", type=float, default=10.0,
                   help="infected titer must exceed the uninfected titer by this "
                        "factor before wMel DE is reportable in that cell type")
    p.add_argument("--drop_labels", default=",".join(LOW_CONFIDENCE))
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def pick(d, *names):
    for n in names:
        if n in d.columns:
            return n
    return None


def build(d, args, inf_c, un_c, suffix, note_txt):
    d = d.sort_values("infected_over_uninf")
    y = np.arange(len(d))
    ok = d["infected_over_uninf"] >= args.min_ratio
    print(f"\n[{suffix}] {len(d)} cell types, {int(ok.sum())} clear "
          f"{args.min_ratio:g}x; {note_txt}")

    fig = plt.figure(figsize=(WIDE, 2.9))
    gs = fig.add_gridspec(1, 2, wspace=0.12, width_ratios=[1.25, 1.0])

    # ---- A: the enrichment ratio ----
    # A cell type with no uninfected wMel at all gives an infinite ratio, which a
    # log axis cannot draw. Cap those at the largest finite value and mark them,
    # rather than dropping the cell types where the signal is cleanest.
    ax = fig.add_subplot(gs[0, 0])
    r = d["infected_over_uninf"].to_numpy(float)
    finite = np.isfinite(r)
    cap = (np.nanmax(r[finite]) * 1.6 if finite.any() else 1e3)
    shown = np.where(finite, r, cap)
    ax.barh(y, np.clip(shown, 0.1, None), height=0.68, zorder=3,
            color=[COL["green"] if k else COL["lightgrey"] for k in ok])
    for yi, f, v in zip(y, finite, shown):
        if not f:
            ax.text(v, yi, " no uninfected signal", fontsize=5, va="center",
                    color=COL["black"])
    ax.axvline(1, color=COL["black"], lw=0.8)
    ax.axvline(args.min_ratio, color=COL["vermillion"], lw=0.9, ls="--")
    ax.set_xscale("log")
    ax.set_yticks(y, ct_short(d["tissue"]), fontsize=6)
    ax.set_xlabel("infected / uninfected pseudobulk titer")
    ax.set_title("Is the wMel signal real here?", pad=4)
    note(ax, f"black line: no enrichment; red line: {args.min_ratio:g}x threshold",
         y=-0.16)
    panel_letter(ax, "A", dx=-0.52)

    # ---- B: the two levels behind the ratio ----
    ax = fig.add_subplot(gs[0, 1])
    ax.scatter(d[inf_c], y - 0.16, s=22, color=COL["green"], zorder=3,
               edgecolors=COL["black"], linewidths=0.3, label="infected")
    ax.scatter(d[un_c], y + 0.16, s=22, color=COL["grey"], zorder=3,
               edgecolors=COL["black"], linewidths=0.3, label="uninfected")
    for yi, (_, r) in zip(y, d.iterrows()):
        ax.plot([r[un_c], r[inf_c]], [yi + 0.16, yi - 0.16],
                color=COL["lightgrey"], lw=0.7, zorder=2)
    ax.set_xscale("log")
    ax.set_yticks(y, [])
    ax.set_xlabel("pseudobulk wMel fraction")
    ax.set_title("Signal and its floor", pad=4)
    ax.legend(loc="lower right", fontsize=6)
    panel_letter(ax, "B", dx=-0.06)

    fig.suptitle(note_txt, y=1.04, fontsize=6.5, color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")

    print(d[["tissue", inf_c, un_c, "infected_over_uninf"]].round(6)
          .to_string(index=False))
    bad = d.loc[~ok, "tissue"]
    if len(bad):
        print(f"\nwMel DE is not reportable in {', '.join(ct_short(bad))}: the "
              "infected titer does not clear the uninfected floor by "
              f"{args.min_ratio:g}x. Significant genes there are differences "
              "between noise distributions.")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    inf_c = pick(d, "titer.wMel", "titer_wMel", "infected_titer")
    un_c = pick(d, "titer.uninf", "titer_uninf", "uninf_titer")
    if not (inf_c and un_c):
        raise SystemExit(f"could not find infected/uninfected titer columns in "
                         f"{list(d.columns)}")
    for c in (inf_c, un_c, "infected_over_uninf"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    if "infected_over_uninf" not in d.columns:
        d["infected_over_uninf"] = d[inf_c] / d[un_c]
    drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}
    d = d[~d["tissue"].isin(drop) & d["infected_over_uninf"].notna()]
    print(f"Labels excluded as low-confidence: {sorted(drop) or 'none'}")
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        sub = d[d["infected_over_uninf"] >= args.min_ratio] if v == "robust" else d
        txt = (f"cell types where infected titer exceeds the uninfected floor by "
               f"{args.min_ratio:g}x" if v == "robust"
               else "every cell type, threshold cleared or not")
        if not len(sub):
            empty_figure(args.outdir, f"{args.stem}_{v}",
                         f"No rows survive the {v} filter.\n"
                         f"An empty set here is a result, not a failure:\n"
                         f"say so in the text rather than citing the\n"
                         f"unfiltered figure as though it passed.")
            continue
        build(sub, args, inf_c, un_c, v, txt)


if __name__ == "__main__":
    main()

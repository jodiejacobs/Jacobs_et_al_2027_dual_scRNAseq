#!/usr/bin/env python
"""Figure 7: per-cell wMel titer against host program usage, by cell type.

Plots results/rescue/titer_programs/titer_program_correlations.tsv.

This is the panel that distinguishes a biological relationship from ambient
contamination, and it does so by where the effect sits rather than by its size.
Ambient germline transcripts predict the strongest germline-program correlation
in somatic cells, where such transcripts can only be contamination, and none in
germ cells. The data show the reverse, so the layout puts germline and somatic
cell types side by side and lets the reader check the direction.

  A  cell type x program heatmap of the partial Spearman, infected samples,
     germline rows grouped first
  B  the decisive comparison: the programs with interpretable effects in
     16-cell cysts against the largest somatic population
  C  per-cell-type effect for one program, with each sample as a point, so
     replication is visible rather than asserted

Two versions are written. "robust" keeps only cell types where every infected
sample agrees in sign, which is the standard the paper leans on; "all" keeps
every group that passed --min_cells. Neither hides the other.

Usage:
  python plot_titer_program_correlations.py \
      --input results/rescue/titer_programs/titer_program_correlations.tsv \
      --outdir figures/fig7
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (empty_figure, COL, LOW_CONFIDENCE, WIDE, apply_style, ct_short,
                      panel_letter, save)
import matplotlib.pyplot as plt

GERMLINE = ("germ", "cyst", "germline")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig7_titer_programs_bycelltype")
    # Only programs that pass the uninfected control in titer_program_tests.tsv.
    # P2, P4, P8, P9, P10 and P13 show the same effect in uninfected samples, so
    # titer does not explain them and they do not belong in a rescue panel.
    p.add_argument("--programs", default="Program_7,Program_3,Program_1,Program_12,Program_5",
                   help="programs that pass the uninfected control, in panel order")
    p.add_argument("--focus", default="Program_7", help="program shown per sample in C")
    p.add_argument("--germline_label", default="16-cell germline cyst in germarium region 2a and 2b")
    p.add_argument("--somatic_label", default="central main body follicle cell ca. St. 6-8")
    p.add_argument("--drop_labels", default=",".join(LOW_CONFIDENCE))
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def is_germline(name):
    n = str(name).lower()
    return any(k in n for k in GERMLINE)


def build(d, args, suffix, note):
    progs = [p.strip() for p in args.programs.split(",") if p.strip()]
    inf = d[d["infection_status"] != "uninf"]
    print(f"\n[{suffix}] {len(inf)} group x program rows, "
          f"{inf['cell_type'].nunique()} cell types; {note}")

    fig = plt.figure(figsize=(WIDE, 4.4))
    gs = fig.add_gridspec(2, 2, hspace=0.46, wspace=0.34,
                          height_ratios=[1.0, 1.0], width_ratios=[1.45, 1.0])

    # ---- A: cell type x program heatmap, germline first ----
    ax = fig.add_subplot(gs[0, :])
    m = inf.pivot_table(index="cell_type", columns="program",
                        values="partial_spearman", aggfunc="mean")
    cols = sorted(m.columns, key=lambda c: int(str(c).split("_")[1]))
    rows = sorted(m.index, key=lambda r: (not is_germline(r), str(r)))
    m = m.reindex(index=rows, columns=cols)
    lim = float(np.nanmax(np.abs(m.to_numpy()))) or 0.1
    im = ax.imshow(m.to_numpy(), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    ax.set_xticks(range(len(cols)), [c.replace("Program_", "P") for c in cols])
    ax.set_yticks(range(len(rows)), ct_short(rows))
    n_germ = sum(is_germline(r) for r in rows)
    if 0 < n_germ < len(rows):
        ax.axhline(n_germ - 0.5, color=COL["black"], lw=1.0)
    # Colour the germline tick labels rather than bracketing them, so the
    # germline/somatic encoding is the same here as in panels B and C.
    for t, r in zip(ax.get_yticklabels(), rows):
        t.set_color(COL["vermillion"] if is_germline(r) else COL["blue"])
    cb = fig.colorbar(im, ax=ax, fraction=0.025, pad=0.015)
    cb.set_label("partial Spearman\n(titer, usage | host depth)", fontsize=6.5)
    cb.ax.tick_params(labelsize=6)
    ax.set_title("Titer against program usage, within cell type and sample", pad=5)
    panel_letter(ax, "A", dx=-0.13)

    # ---- B: germline against the largest somatic population ----
    ax = fig.add_subplot(gs[1, 0])
    pair = [args.germline_label, args.somatic_label]
    sub = (inf[inf["cell_type"].isin(pair) & inf["program"].isin(progs)]
           .pivot_table(index="program", columns="cell_type",
                        values="partial_spearman", aggfunc="mean")
           .reindex(index=progs, columns=pair))
    x = np.arange(len(progs))
    for i, (lab, c) in enumerate(zip(pair, (COL["vermillion"], COL["blue"]))):
        if lab in sub.columns:
            ax.bar(x + (i - 0.5) * 0.38, sub[lab], width=0.36, color=c,
                   label=ct_short([lab])[0])
    ax.axhline(0, color=COL["black"], lw=0.6)
    ax.set_xticks(x, [p.replace("Program_", "P") for p in progs])
    ax.set_ylabel("partial Spearman")
    ax.set_title("Germline against largest somatic", pad=4)
    ax.legend(loc="lower left", fontsize=6)
    panel_letter(ax, "B")

    # ---- C: one program, every cell type, each sample a point ----
    ax = fig.add_subplot(gs[1, 1])
    f = inf[inf["program"] == args.focus]
    if len(f):
        order = (f.groupby("cell_type")["partial_spearman"].mean()
                 .sort_values().index.tolist())
        for i, ct in enumerate(order):
            v = f.loc[f["cell_type"] == ct, "partial_spearman"].to_numpy()
            ax.scatter(v, np.full(len(v), i), s=9, alpha=0.9,
                       color=COL["vermillion"] if is_germline(ct) else COL["blue"],
                       edgecolors="none")
            ax.plot([v.mean()] * 2, [i - 0.3, i + 0.3], color=COL["black"], lw=1.0)
        ax.axvline(0, color=COL["grey"], lw=0.6, ls="--")
        ax.set_yticks(range(len(order)), ct_short(order), fontsize=5.5)
        ax.set_xlabel("partial Spearman")
        ax.set_title(f"{args.focus.replace('Program_', 'Program ')}, one point per sample",
                     pad=4)
        h = [plt.Line2D([], [], marker="o", ls="", color=c, markersize=4)
             for c in (COL["vermillion"], COL["blue"])]
        ax.legend(h, ["germline", "somatic"], loc="lower right", fontsize=6)
    panel_letter(ax, "C")

    fig.suptitle(note, y=1.003, fontsize=6.5, color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")

    if args.germline_label in m.index and args.somatic_label in m.index:
        for p in progs:
            if p in m.columns:
                print(f"  {p}: germline {m.loc[args.germline_label, p]:+.3f}  "
                      f"somatic {m.loc[args.somatic_label, p]:+.3f}")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}
    d = d[~d["cell_type"].isin(drop)].copy()
    print(f"Labels excluded as low-confidence: {sorted(drop) or 'none'}")

    # "robust" here means every infected sample in that group agrees in sign,
    # which is the standard the paper leans on at two replicates per genotype.
    inf = d[d["infection_status"] != "uninf"]
    sign_ok = (inf.groupby(["cell_type", "program"])["partial_spearman"]
               .apply(lambda s: len(set(np.sign(s))) == 1 and len(s) > 1))
    keep = set(sign_ok[sign_ok].index)
    versions = (["robust", "all"] if args.version == "both" else [args.version])
    for v in versions:
        if v == "robust":
            sub = d[[(ct, p) in keep for ct, p in zip(d["cell_type"], d["program"])]]
            note = "cell type x program groups where all infected samples agree in sign"
        else:
            sub = d
            note = "every group passing the minimum cell count"
        if not len(sub):
            empty_figure(args.outdir, f"{args.stem}_{v}",
                         f"No rows survive the {v} filter.\n"
                         f"An empty set here is a result, not a failure:\n"
                         f"say so in the text rather than citing the\n"
                         f"unfiltered figure as though it passed.")
            continue
        build(sub, args, v, note)


if __name__ == "__main__":
    main()

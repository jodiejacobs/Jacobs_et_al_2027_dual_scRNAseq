#!/usr/bin/env python
"""Eggshell figure: the eggshell program across the five rescue contrasts.

Reads de/host_de_all.tsv.gz. Gene set and stage classes from gene_sets.py.

  A  one row per follicle-cell type, one column per contrast. Colour is the
     median logFC over the eggshell genes tested in that cell type; the text
     is how many are significant out of how many tested. The argument is in
     reading across a row: down in the uninfected mutant, slightly down with
     wMel in wild type, strongly up with wMel in the mutant. The opposite
     signs of the two infection columns are what makes the interaction large,
     so they sit side by side.
  B  the interaction split by stage class. Vitelline membrane and yolk come
     on at stages 8-10, chorion at 11-14, so a class that recovers says the
     egg chambers reach that stage. If the early classes moved and chorion did
     not, the rescue would stop before choriogenesis.

The germline is left out. Germ cells do not make eggshell, and the eggshell
reads in germline droplets track follicle-cell levels across libraries
(18_eggshell.py, ambient check), so those rows would be carry-over.

Usage:
  python plot_eggshell_contrasts.py --input results/rescue/de/host_de_all.tsv.gz \
      --outdir results/rescue/figures
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (COL, LOW_CONFIDENCE, WIDE, apply_style, ct_short, note,
                      panel_letter, save)
from gene_sets import (EGGSHELL_ORDER, FOLLICLE_KEY, eggshell_class,
                       eggshell_symbols)
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

CONS = ["genotype_in_uninf", "infection_in_OreR", "infection_in_MeiP26",
        "interaction"]
NICE = {"genotype_in_uninf": "defect\n(mei-P26 vs OreR)",
        "infection_in_OreR": "wMel in\nOreR",
        "infection_in_MeiP26": "wMel in\nmei-P26",
        "interaction": "interaction\n(rescue test)"}
CMAP = LinearSegmentedColormap.from_list(
    "d", [COL["blue"], "#FFFFFF", COL["vermillion"]])


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="de/host_de_all.tsv.gz")
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig_eggshell_contrasts")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--vmax", type=float, default=3.0)
    p.add_argument("--min_genes", type=int, default=5,
                   help="a cell type testing fewer eggshell genes is left out")
    return p.parse_args()


def grid(ax, val, txt, rows, cols, vmax):
    im = ax.imshow(val, aspect="auto", cmap=CMAP, vmin=-vmax, vmax=vmax)
    for i in range(len(rows)):
        for j in range(len(cols)):
            v = val[i, j]
            if np.isfinite(v):
                ax.text(j, i, txt[i][j], ha="center", va="center", fontsize=5.5,
                        color="white" if abs(v) > vmax * 0.6 else COL["black"])
    ax.set_xticks(range(len(cols)))
    ax.set_yticks(range(len(rows)))
    ax.tick_params(length=0)
    for s in ax.spines.values():
        s.set_visible(False)
    return im


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    cls = eggshell_class()
    d = d[d["symbol"].isin(eggshell_symbols())
          & d["tissue"].str.lower().str.contains(FOLLICLE_KEY)
          & ~d["tissue"].isin(LOW_CONFIDENCE)].copy()
    d["class"] = d["symbol"].map(cls)
    for c in ("logFC", "FDR"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    n_t = d[d["contrast"] == "interaction"].groupby("tissue").size()
    tissues = [t for t in n_t.index if n_t[t] >= args.min_genes]
    # Order by the interaction, strongest at the top.
    med_it = (d[(d["contrast"] == "interaction") & d["tissue"].isin(tissues)]
              .groupby("tissue")["logFC"].median().sort_values(ascending=False))
    tissues = med_it.index.tolist()
    if not tissues:
        raise SystemExit("no follicle-cell type tests enough eggshell genes")

    def summary(sub, rows, cols, rkey, ckey):
        val = np.full((len(rows), len(cols)), np.nan)
        txt = [["" for _ in cols] for _ in rows]
        for i, r in enumerate(rows):
            for j, c in enumerate(cols):
                x = sub[(sub[rkey] == r) & (sub[ckey] == c)]
                if len(x):
                    val[i, j] = x["logFC"].median()
                    txt[i][j] = f"{int((x['FDR'] < args.fdr).sum())}/{len(x)}"
        return val, txt

    fig = plt.figure(figsize=(WIDE, 0.34 * len(tissues) + 1.9))
    gs = fig.add_gridspec(1, 3, wspace=0.08, width_ratios=[4.0, 3.0, 0.15])

    ax = fig.add_subplot(gs[0, 0])
    val, txt = summary(d[d["tissue"].isin(tissues)], tissues, CONS,
                       "tissue", "contrast")
    im = grid(ax, val, txt, tissues, CONS, args.vmax)
    ax.set_xticklabels([NICE[c] for c in CONS], fontsize=6)
    ax.set_yticklabels(ct_short(tissues), fontsize=6)
    ax.xaxis.tick_top()
    ax.set_title("Eggshell genes across the contrasts", pad=24, fontsize=7.5)
    panel_letter(ax, "A", dx=-0.38, dy=1.21)

    ax = fig.add_subplot(gs[0, 1])
    it = d[(d["contrast"] == "interaction") & d["tissue"].isin(tissues)]
    val, txt = summary(it, tissues, EGGSHELL_ORDER, "tissue", "class")
    grid(ax, val, txt, tissues, EGGSHELL_ORDER, args.vmax)
    ax.set_xticklabels(["vitelline memb.\n(St.8-10)", "yolk\n(St.8-10)",
                        "chorion\n(St.11-14)"], fontsize=6)
    ax.set_yticklabels([])
    ax.xaxis.tick_top()
    ax.set_title("Interaction, by stage class", pad=24, fontsize=7.5)
    panel_letter(ax, "B", dx=-0.04, dy=1.21)

    cax = fig.add_subplot(gs[0, 2])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("median log$_2$FC", fontsize=6)
    cb.ax.tick_params(labelsize=5.5)
    cb.outline.set_visible(False)
    note(fig.axes[0], f"text: genes at FDR < {args.fdr} / genes tested; "
                      "germline left out, its eggshell reads are ambient",
         y=-0.06)
    save(fig, args.outdir, args.stem)

    piv = (d[d["tissue"].isin(tissues)].groupby(["tissue", "contrast"])
           ["logFC"].median().unstack()[CONS].round(2))
    piv.index = ct_short(piv.index)
    print(piv.to_string())
    piv = it.groupby(["tissue", "class"])["logFC"].median().unstack()[EGGSHELL_ORDER]
    piv.index = ct_short(piv.index)
    print("\ninteraction by stage class:")
    print(piv.round(2).to_string())


if __name__ == "__main__":
    main()

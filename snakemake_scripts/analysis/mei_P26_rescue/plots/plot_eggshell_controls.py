#!/usr/bin/env python
"""Eggshell supplement: the ambient control and the titer link.

Plots eggshell_ambient.tsv and eggshell_titer.tsv from 18_eggshell.py.

  A  ambient. Per library, the mean eggshell score in marker-defined germline
     cells against the mean in follicle cells. Germ cells do not make
     eggshell, so if germline tracks follicle across libraries, the eggshell
     reads in germline droplets are carry-over from the same library. That is
     the justification for leaving germline out of the eggshell panels and for
     flagging eggshell genes in the germline panels of the interaction volcano.
  B  titer. Within each infected library, Spearman of per-cell eggshell score
     against per-cell wMel titer, per follicle-cell type. Exploratory only:
     titer and stage both change along oogenesis, so a correlation here would
     not separate wMel acting on the program from both tracking stage. In
     several follicle types the infected mutant sits slightly positive and
     infected wild type slightly negative; with two libraries per group that
     is a pattern to look at, not a result, and the title does not claim it.

Usage:
  python plot_eggshell_controls.py \
      --input results/rescue/eggshell/eggshell_ambient.tsv \
      --titer results/rescue/eggshell/eggshell_titer.tsv \
      --outdir results/rescue/figures
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (COL, GROUP_BASE, GROUP_ORDER, WIDE, apply_style,
                      ct_short, note, panel_letter, save)
from gene_sets import EGGSHELL_ORDER
import matplotlib.pyplot as plt

NICE = {"OreR_uninf": "OreR", "OreR_wMel": "OreR + wMel",
        "MeiP26_uninf": "mei-P26", "MeiP26_wMel": "mei-P26 + wMel"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="eggshell_ambient.tsv")
    p.add_argument("--titer", required=True, help="eggshell_titer.tsv")
    p.add_argument("--klass", default="all eggshell",
                   help="class drawn in both panels; every class's ambient "
                        "rho is printed in panel A regardless")
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="supp_eggshell_controls")
    return p.parse_args()


def main():
    args = parse_args()
    apply_style()
    a = pd.read_csv(args.input, sep="\t")
    t = pd.read_csv(args.titer, sep="\t")

    fig = plt.figure(figsize=(WIDE, 2.8))
    gs = fig.add_gridspec(1, 2, wspace=0.62, width_ratios=[1.0, 1.4])

    ax = fig.add_subplot(gs[0, 0])
    s = a[a["class"] == args.klass]
    for g in GROUP_ORDER:
        x = s[s["group"] == g]
        ax.scatter(x["follicle_mean"], x["germline_mean"], s=30, zorder=3,
                   color=GROUP_BASE[g], edgecolors=COL["black"], linewidths=0.4,
                   label=NICE[g])
    lines = []
    for c in EGGSHELL_ORDER + ["all eggshell"]:
        x = a[a["class"] == c]
        if len(x) >= 3:
            r = stats.spearmanr(x["follicle_mean"], x["germline_mean"]).statistic
            lines.append(f"{c}: ρ = {r:+.2f}")
    ax.text(0.03, 0.97, "\n".join(lines), transform=ax.transAxes, va="top",
            fontsize=5.5)
    ax.set_xlabel(f"follicle cells, mean {args.klass} score", fontsize=6.5)
    ax.set_ylabel(f"germline cells, mean {args.klass} score", fontsize=6.5)
    ax.set_title("Germline eggshell reads are ambient", pad=4, fontsize=7.5)
    ax.legend(loc="lower right", fontsize=5.5, frameon=False, handletextpad=0.2)
    note(ax, "one point per library; germline by vas/aub/AGO3/orb/osk", y=-0.24)
    panel_letter(ax, "A", dx=-0.30)

    ax = fig.add_subplot(gs[0, 1])
    s = t[t["class"] == args.klass]
    types = (s.groupby("tissue")["rho"].median().sort_values().index.tolist())
    rng = np.random.default_rng(0)
    for i, ty in enumerate(types):
        for g in ("OreR_wMel", "MeiP26_wMel"):
            v = s[(s["tissue"] == ty) & (s["group"] == g)]["rho"].to_numpy(float)
            off = -0.12 if g == "OreR_wMel" else 0.12
            ax.scatter(v, np.full(len(v), i + off) + rng.uniform(-0.04, 0.04, len(v)),
                       s=22, color=GROUP_BASE[g], edgecolors=COL["black"],
                       linewidths=0.4, zorder=3,
                       label=NICE[g] if i == 0 else None)
    ax.axvline(0, color=COL["black"], lw=0.7)
    ax.set_yticks(range(len(types)), ct_short(types), fontsize=6)
    ax.set_xlim(-1, 1)
    ax.set_xlabel(f"Spearman ρ, per-cell {args.klass} score vs titer",
                  fontsize=6.5)
    ax.set_title("Eggshell score against titer (exploratory)", pad=4,
                 fontsize=7.5)
    ax.legend(loc="lower right", fontsize=5.5, frameon=False, handletextpad=0.2)
    note(ax, "infected libraries only, one point each; titer and stage are "
             "confounded, so this cannot\nseparate wMel acting on the program "
             "from both tracking oogenesis", y=-0.24)
    panel_letter(ax, "B", dx=-0.40)
    save(fig, args.outdir, args.stem)

    print("\n".join(lines))
    print("\nmedian titer rho, by cell type and group:")
    print(s.groupby(["tissue", "group"])["rho"].median().unstack().round(2)
          .to_string())


if __name__ == "__main__":
    main()

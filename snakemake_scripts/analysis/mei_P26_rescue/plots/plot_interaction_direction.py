#!/usr/bin/env python
"""Figure 5: direction of the genotype x infection interaction, per cell type.

Plots the interaction rows of results/rescue/de/host_de_summary.tsv.

The interaction term asks whether wMel changes host expression differently in
mei-P26 than in OreR. Positive logFC means the infection effect is larger in
mei-P26, negative means smaller. A cell type whose significant genes split
evenly has no preferred direction; one that leans to a side has a consistent
genotype-specific response.

  A  significant interaction genes per cell type, down to the left and up to
     the right, ordered by how many there are
  B  the same as a percentage up, against the 50% a direction-blind split
     would give

Cell types with fewer than --min_sig significant genes are greyed in B and
marked, because a percentage over 12 genes is not worth reading. The artefact
and unannotated labels are shown, since they hold many of the hits, and are
marked as low-confidence labels.

Usage:
  python plot_interaction_direction.py \
      --input results/rescue/de/host_de_summary.tsv \
      --outdir results/rescue/figures
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (COL, DIRECTION, LOW_CONFIDENCE, WIDE, apply_style, ct_short,
                      note, panel_letter, save)
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig5_interaction_direction")
    p.add_argument("--contrast", default="interaction")
    p.add_argument("--min_sig", type=int, default=50,
                   help="significant genes below which the percentage is greyed")
    return p.parse_args()


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    d = d[(d["contrast"] == args.contrast) & (d["n_sig"] > 0)]
    d = d.sort_values("n_sig", ascending=False).reset_index(drop=True)
    d["pct_up"] = 100 * d["n_up"] / d["n_sig"]
    thin = d["n_sig"] < args.min_sig
    names = [f"{n}{'†' if t in LOW_CONFIDENCE else ''}"
             for n, t in zip(ct_short(d["tissue"]), d["tissue"])]
    print(d[["tissue", "n_sig", "n_up", "n_down", "pct_up"]].round(0).to_string(index=False))

    fig, (a, b) = plt.subplots(1, 2, figsize=(WIDE, 0.24 * len(d) + 0.9), sharey=True,
                               gridspec_kw={"width_ratios": [1.5, 1], "wspace": 0.08})
    y = np.arange(len(d))
    a.barh(y, -d["n_down"], height=0.72, color=DIRECTION["down"], zorder=3, label="down")
    a.barh(y, d["n_up"], height=0.72, color=DIRECTION["up"], zorder=3, label="up")
    a.axvline(0, color=COL["black"], lw=0.7, zorder=4)
    lim = d[["n_up", "n_down"]].to_numpy().max() * 1.1
    a.set_xlim(-lim, lim)
    a.set_xticks(a.get_xticks(), [f"{abs(int(t))}" for t in a.get_xticks()])
    a.set_xlabel("significant interaction genes")
    a.set_title("Direction of the interaction", pad=4)
    a.legend(loc="lower right", ncol=2, fontsize=6)
    a.set_yticks(y, [f"{n}  ({int(s)})" for n, s in zip(names, d["n_sig"])], fontsize=6)
    a.set_ylim(len(d) - 0.4, -0.6)
    panel_letter(a, "A", dx=-0.62, dy=1.10)

    b.barh(y, d["pct_up"], height=0.72, zorder=3,
           color=[COL["lightgrey"] if t else COL["grey"] for t in thin])
    b.axvline(50, color=COL["black"], lw=0.7, ls="--", zorder=4)
    for yi, p, t in zip(y, d["pct_up"], thin):
        b.text(104, yi, f"{p:.0f}%{' ✻' if t else ''}", fontsize=6, va="center",
               color=COL["grey"])
    b.set_xlim(0, 100)
    b.set_xticks([0, 25, 50, 75, 100])
    b.set_xlabel("% of genes up")
    b.set_title("Against an even split", pad=4)
    panel_letter(b, "B", dx=-0.06, dy=1.10)
    note(b, f"✻ fewer than {args.min_sig} genes   † low-confidence label", y=-0.22)

    save(fig, args.outdir, args.stem)


if __name__ == "__main__":
    main()

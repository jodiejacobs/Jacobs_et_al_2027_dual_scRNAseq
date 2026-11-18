#!/usr/bin/env python
"""Host pathway enrichment on a rescue contrast, as a term x cell-type dot plot.

Plots results/rescue/host_gsea/host_gsea.tsv from 17_host_gsea.py.

Colour is the normalised enrichment score, size is -log10 FDR, and a term that
does not reach --fdr in a cell type is drawn as a small open circle so its
direction stays visible without reading as a result. Two presets:

  top         the strongest terms: up to --per_type per cell type at FDR <
              --fdr, pooled and capped at --max_terms. For the supplement.
  secretory   the program the eggshell rescue needs: ribosome and
              translation, ER processing, protein export, ER to Golgi
              transport, and eggshell and vitelline membrane formation. For
              the eggshell figure, where the question is whether the
              machinery moves with the eggshell genes.

FDR is per library, cell type and contrast (see 17_host_gsea.py), and
ribosome and translation terms are the ones most exposed to depth and stage
mixture. The panel says so.

Usage:
  python plot_host_gsea.py --input results/rescue/host_gsea/host_gsea.tsv \
      --preset secretory --outdir results/rescue/figures --stem fig_eggshell_gsea
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (COL, WIDE, apply_style, ct_short, note, panel_letter,
                      save)
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

SECRETORY = [r"^Ribosome$", r"^cytoplasmic translation",
             r"^Protein processing in endoplasmic reticulum$",
             r"^Protein export$", r"^ER to Golgi vesicle-mediated transport",
             r"^chorion-containing eggshell formation",
             r"^vitelline membrane formation"]
CMAP = LinearSegmentedColormap.from_list(
    "nes", [COL["blue"], "#FFFFFF", COL["vermillion"]])


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="host_gsea.tsv")
    p.add_argument("--contrast", default="interaction")
    p.add_argument("--preset", choices=("top", "secretory"), default="top")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--per_type", type=int, default=3)
    p.add_argument("--max_terms", type=int, default=24)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="supp_host_gsea")
    return p.parse_args()


def label(term, lib):
    t = re.sub(r"\s*\(GO:\d+\)$", "", str(term))
    t = t if len(t) <= 52 else t[:50] + "…"
    return f"{t} (KEGG)" if str(lib).startswith("KEGG") else t


def main():
    args = parse_args()
    apply_style()
    g = pd.read_csv(args.input, sep="\t")
    g = g[g["contrast"] == args.contrast].copy()
    if not len(g):
        raise SystemExit(f"no rows for contrast {args.contrast}")
    g["key"] = g["library"] + "\t" + g["term"]

    if args.preset == "secretory":
        keys = []
        for pat in SECRETORY:
            hit = g.loc[g["term"].str.contains(pat, regex=True), "key"].unique()
            keys += [k for k in hit if k not in keys]
        title = "Secretory machinery moves with the eggshell genes"
    else:
        sig = g[g["fdr"] < args.fdr].sort_values("fdr")
        keys = []
        for _, s in sig.groupby("tissue"):
            for k in s["key"].head(args.per_type):
                if k not in keys:
                    keys.append(k)
        keys = keys[:args.max_terms]
        title = f"Strongest terms, {args.contrast}"
    if not keys:
        raise SystemExit("no terms selected; nothing passes for this preset")

    sub = g[g["key"].isin(keys)]
    types = (sub.groupby("tissue")["nes"].mean().sort_values(ascending=False)
             .index.tolist())
    order = (sub.groupby("key")["nes"].mean().sort_values().index.tolist())
    nes = sub.pivot_table(index="key", columns="tissue", values="nes").reindex(
        index=order, columns=types)
    fdr = sub.pivot_table(index="key", columns="tissue", values="fdr").reindex(
        index=order, columns=types)
    lib = sub.drop_duplicates("key").set_index("key")["library"]

    vmax = float(np.nanmax(np.abs(nes.to_numpy()))) if nes.notna().any().any() else 2.0
    fig, ax = plt.subplots(figsize=(WIDE * 0.86, 0.24 * len(order) + 1.9))
    for i, k in enumerate(order):
        for j, t in enumerate(types):
            v, q = nes.loc[k, t], fdr.loc[k, t]
            if not np.isfinite(v):
                continue
            if np.isfinite(q) and q < args.fdr:
                size = 14 + 22 * min(-np.log10(max(q, 1e-3)), 3)
                ax.scatter(j, i, s=size, color=CMAP((v + vmax) / (2 * vmax)),
                           edgecolors=COL["black"], linewidths=0.4, zorder=3)
            else:
                ax.scatter(j, i, s=10, facecolors="none",
                           edgecolors=CMAP((v + vmax) / (2 * vmax)),
                           linewidths=0.8, zorder=2)
    ax.set_xticks(range(len(types)), ct_short(types), rotation=35, ha="right",
                  fontsize=6)
    ax.set_yticks(range(len(order)),
                  [label(k.split("\t", 1)[1], lib[k]) for k in order], fontsize=6)
    ax.set_xlim(-0.6, len(types) - 0.4)
    ax.set_ylim(-0.6, len(order) - 0.4)
    ax.grid(color=COL["lightgrey"], lw=0.4, zorder=0)
    ax.set_title(title, pad=4, fontsize=7.5)
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(-vmax, vmax))
    cb = fig.colorbar(sm, ax=ax, fraction=0.035, pad=0.02)
    cb.set_label("NES", fontsize=6)
    cb.ax.tick_params(labelsize=5.5)
    cb.outline.set_visible(False)
    note(ax, f"filled: FDR < {args.fdr}, larger = smaller FDR; open: not "
             f"significant; contrast: {args.contrast}. Ribosome and "
             "translation terms\nmove with depth and stage mixture, both "
             "confounded with genotype here; read them as context",
         y=-0.30 if len(order) < 10 else -0.16)
    panel_letter(ax, "A", dx=-0.62 if args.preset == "top" else -0.50)
    save(fig, args.outdir, args.stem)

    show = nes.round(2)
    show.index = [label(k.split("\t", 1)[1], lib[k]) for k in show.index]
    show.columns = ct_short(show.columns)
    print(show.to_string())


if __name__ == "__main__":
    main()

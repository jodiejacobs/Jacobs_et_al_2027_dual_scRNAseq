#!/usr/bin/env python
"""Eggshell figure: eggshell genes across the eight libraries, per cell type.

Reads de/host_logcpm_<tissue>.tsv.gz for the values and de/host_de_all.tsv.gz
for the interaction FDR. Gene set and stage classes from gene_sets.py.

The panel the eye goes to first. Each block is one follicle-cell type; rows
are eggshell genes grouped by stage class, columns the eight libraries in
group order. Values are logCPM z-scored per gene within the block, so colour
says where a gene sits relative to its own range in that cell type, and the
pattern to look for is a cold mutant column pair turning hot with wMel while
the two wild-type pairs barely move. Replicates sit side by side so their
agreement is visible rather than averaged away.

Stars beside each row are that gene's interaction FDR in that cell type.

Usage:
  python plot_eggshell_heatmap.py --input results/rescue/de/host_de_all.tsv.gz \
      --logcpm_dir results/rescue/de --outdir results/rescue/figures
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (COL, GROUP_BASE, GROUP_ORDER, WIDE, apply_style,
                      ct_short, fdr_stars, note, panel_letter, sample_group,
                      save)
from gene_sets import EGGSHELL, EGGSHELL_ORDER, eggshell_class
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

DEFAULT_TISSUES = ("main body follicle cell ca. until St. 5",
                   "central main body follicle cell ca. St. 6-8",
                   "choriogenic main body follicle cell and corpus luteum")
# Blue through white to vermillion, the Okabe-Ito ends of DIRECTION.
CMAP = LinearSegmentedColormap.from_list(
    "z", [COL["blue"], "#FFFFFF", COL["vermillion"]])
# A gene not tested in a cell type is grey, so it cannot be read as "average".
CMAP.set_bad("#EBEBEB")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="de/host_de_all.tsv.gz")
    p.add_argument("--logcpm_dir", required=True)
    p.add_argument("--tissues", default="|".join(DEFAULT_TISSUES),
                   help="'|'-separated cell types, one block each")
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig_eggshell_heatmap")
    p.add_argument("--zmax", type=float, default=2.0)
    return p.parse_args()


def load_block(d, tissue, ids):
    """logCPM for one cell type, genes x samples, from whichever table has it."""
    for f in glob.glob(os.path.join(d, "host_logcpm_*.tsv.gz")):
        head = pd.read_csv(f, sep="\t", nrows=0).columns
        if any(str(c).endswith(f"|{tissue}") for c in head[1:]):
            t = pd.read_csv(f, sep="\t", index_col=0)
            t.columns = [str(c).split("|")[0] for c in t.columns]
            return t.reindex(ids)
    return None


def main():
    args = parse_args()
    apply_style()
    de = pd.read_csv(args.input, sep="\t")
    sym2id = (de.dropna(subset=["symbol"]).drop_duplicates("symbol")
              .set_index("symbol")["gene_id"])
    cls = eggshell_class()
    genes = [s for c in EGGSHELL_ORDER for s in EGGSHELL[c] if s in sym2id.index]
    ids = [sym2id[s] for s in genes]
    tissues = [t for t in args.tissues.split("|") if t]

    blocks = []
    for t in tissues:
        b = load_block(args.logcpm_dir, t, ids)
        if b is None:
            print(f"  no logCPM table for {t}; skipped")
            continue
        order = [s for g in GROUP_ORDER for s in sorted(b.columns)
                 if sample_group(s) == g]
        b = b[order]
        z = b.sub(b.mean(axis=1), axis=0).div(b.std(axis=1).replace(0, np.nan), axis=0)
        it = (de[(de["tissue"] == t) & (de["contrast"] == "interaction")]
              .set_index("gene_id")["FDR"].reindex(ids))
        blocks.append((t, z, it))
    if not blocks:
        raise SystemExit("no cell type had a logCPM table")

    nb = len(blocks)
    fig = plt.figure(figsize=(WIDE, 0.135 * len(genes) + 1.5))
    gs = fig.add_gridspec(1, nb + 1, wspace=0.24,
                          width_ratios=[1.0] * nb + [0.05])
    # Class boundaries, for separators and labels.
    bounds, start = [], 0
    for c in EGGSHELL_ORDER:
        n = sum(cls[s] == c for s in genes)
        if n:
            bounds.append((c, start, start + n))
            start += n

    for k, (t, z, it) in enumerate(blocks):
        ax = fig.add_subplot(gs[0, k])
        im = ax.imshow(z.to_numpy(float), aspect="auto", cmap=CMAP,
                       vmin=-args.zmax, vmax=args.zmax, interpolation="nearest")
        for _, s0, s1 in bounds[1:]:
            ax.axhline(s0 - 0.5, color=COL["black"], lw=0.7)
        # Group colour strip above the columns.
        for j, s in enumerate(z.columns):
            ax.add_patch(plt.Rectangle((j - 0.5, -1.45), 1, 0.75,
                                       color=GROUP_BASE[sample_group(s)],
                                       clip_on=False, lw=0))
        for j in (1.5, 3.5, 5.5):
            ax.axvline(j, color="white", lw=1.2)
        ax.set_xticks([])
        ax.set_ylim(len(genes) - 0.5, -1.6)
        if k == 0:
            ax.set_yticks(range(len(genes)), genes, fontsize=5.5,
                          fontstyle="italic")
            for c, s0, s1 in bounds:
                ax.text(-0.42, (s0 + s1 - 1) / 2, c, rotation=90, ha="center",
                        va="center", fontsize=6, transform=ax.get_yaxis_transform())
        else:
            ax.set_yticks([])
        for i, q in enumerate(it.to_numpy()):
            st = fdr_stars(q)
            if st:
                ax.text(len(z.columns) - 0.35, i, st, fontsize=5, va="center",
                        ha="left", color=COL["black"], clip_on=False)
        ax.set_title(ct_short([t])[0], pad=10, fontsize=7)
        for s in ax.spines.values():
            s.set_visible(False)
        ax.tick_params(length=0)
        panel_letter(ax, "ABCDEFG"[k], dx=-0.06 if k else -0.55, dy=1.06)

    cax = fig.add_subplot(gs[0, nb])
    cb = fig.colorbar(im, cax=cax)
    cb.set_label("logCPM, z-score within cell type", fontsize=6)
    cb.ax.tick_params(labelsize=5.5)
    cb.outline.set_visible(False)
    h = [plt.Rectangle((0, 0), 1, 1, color=GROUP_BASE[g]) for g in GROUP_ORDER]
    fig.legend(h, ["OreR", "OreR + wMel", "mei-P26", "mei-P26 + wMel"],
               loc="lower center", ncol=4, fontsize=6, frameon=False,
               bbox_to_anchor=(0.5, -0.02))
    note(fig.axes[0], "two libraries per group, side by side; stars: "
                      "interaction FDR in that cell type; grey: not tested "
                      "there", y=-0.06)
    save(fig, args.outdir, args.stem)

    print(f"{len(genes)} eggshell genes, {nb} cell types")
    for t, z, it in blocks:
        n = int((it < 0.05).sum())
        print(f"  {ct_short([t])[0]:<24} {n} of {it.notna().sum()} genes "
              "interaction FDR < 0.05")


if __name__ == "__main__":
    main()

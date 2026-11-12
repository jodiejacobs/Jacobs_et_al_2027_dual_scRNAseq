#!/usr/bin/env python
"""Figure 2: the wMel signal is specific, and detection is unbiased by sequence.

Composes the actual figure from several tables, rather than one table per
figure. The per-table scripts in this directory stay as they are and feed the
supplement; this is the version that carries the argument.

Both claims here are sameness claims: titer separates infection in every cell
type, and sequence predicts nothing about detection once expression is held
fixed. Bar charts show sameness badly, because the reader's eye goes to the
small differences between bars rather than to the fact that they are all alike.
Every panel below is therefore a distribution or a point cloud.

  A  per-cell titer, infected against uninfected. Two parts, because the
     result has two: a bar for the cells with no wMel UMIs at all, which is
     most of the uninfected set and a fifth of the infected one, and a
     histogram for the cells that carry signal. The uninfected libraries
     contain no Wolbachia, so their distribution is the assay's noise and the
     gap between the two is the whole specificity claim.
  B  AUROC per cell type as one point each, against chance. The claim is that
     they cluster high, not that any one of them is highest, so they are not
     ranked bars.
  C  detection against expression, loci split into length terciles. If length
     drove detection the terciles would separate vertically. They lie on top of
     each other, which shows the null as data rather than as a coefficient near
     zero, and replaces three separate coefficient panels.
Panel D, wMel against host depth in uninfected cells, is available behind
--with_ambient but off by default. Measured on this data the relationship is
weak and non-monotonic (rho 0.124, and the top depth decile falls), so a panel
titled "what is left is ambient" would claim more than it shows. The
per-sample version in fig2_ambient_depth reports it without that framing.

Usage:
  python plot_figure2.py \
      --counts results/rescue/pseudobulk/counts_labeled.h5ad \
      --separation results/rescue/uninfected_floor/titer_separation.tsv \
      --gene_detection results/rescue/detection_bias/wmel_gene_detection.tsv.gz \
      --outdir results/rescue/figures
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (COL, LOW_CONFIDENCE, WIDE, apply_style, ct_short, note,
                      panel_letter, save)
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts", required=True,
                   help="counts_labeled.h5ad; only .obs is read")
    p.add_argument("--separation", required=True)
    p.add_argument("--gene_detection", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig2_signal_specificity")
    p.add_argument("--floor_quantile", type=float, default=0.99)
    p.add_argument("--drop_labels", default=",".join(LOW_CONFIDENCE))
    p.add_argument("--max_points", type=int, default=40000)
    p.add_argument("--with_ambient", action="store_true",
                   help="add a fourth panel of wMel against host depth in "
                        "uninfected cells. Off by default: the relationship is "
                        "weak and non-monotonic (rho 0.124, binned means wobble "
                        "and the top decile falls), so it does not carry the "
                        "claim its title would make. The per-sample version in "
                        "fig2_ambient_depth is the honest place for it")
    return p.parse_args()


def load_obs(path):
    import anndata as ad
    o = ad.read_h5ad(path, backed="r").obs.copy()
    for c in ("infection_status", "tissue", "sample"):
        o[c] = o[c].astype(str)
    for c in ("titer", "host_umis", "wmel_umis"):
        if c in o.columns:
            o[c] = pd.to_numeric(o[c], errors="coerce")
    return o


def main():
    args = parse_args()
    apply_style()
    obs = load_obs(args.counts)
    drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}
    obs = obs[~obs["tissue"].isin(drop)]
    un = obs["infection_status"] == "uninf"
    floor = float(np.quantile(obs.loc[un, "titer"], args.floor_quantile))
    print(f"{len(obs):,} cells, {int(un.sum()):,} uninfected; "
          f"floor at q{args.floor_quantile} = {floor:.3e}")

    four = args.with_ambient
    fig = plt.figure(figsize=(WIDE, 5.0 if four else 2.7))
    gs = (fig.add_gridspec(2, 2, hspace=0.70, wspace=0.34) if four
          else fig.add_gridspec(1, 3, wspace=0.58,
                                width_ratios=[1.62, 1.0, 1.2]))
    cell = (lambda i: gs[i // 2, i % 2]) if four else (lambda i: gs[0, i])

    # ---- A: the separation, as distributions, with the zero class drawn ----
    # Most uninfected cells carry zero wMel UMIs. Binning those zeros with the
    # rest puts a spike at the axis edge taller than everything else and hides
    # the comparison, but leaving them out of the drawing hides the larger half
    # of the result, and a caption is not a panel. They therefore get their own
    # narrow axis beside the histogram: every cell is drawn somewhere, and
    # neither part crushes the other.
    sub = cell(0).subgridspec(1, 2, width_ratios=[0.36, 1.0], wspace=0.60)
    axz = fig.add_subplot(sub[0, 0])
    ax = fig.add_subplot(sub[0, 1])
    zu = float((obs.loc[un, "titer"] <= 0).mean())
    zi = float((obs.loc[~un, "titer"] <= 0).mean())
    axz.bar([0, 1], [zu * 100, zi * 100], width=0.66,
            color=[COL["grey"], COL["green"]], alpha=0.65,
            edgecolor=COL["black"], linewidth=0.5)
    for x, v in ((0, zu), (1, zi)):
        axz.text(x, v * 100 + 2.0, f"{v:.0%}", ha="center", fontsize=5.5)
    axz.set_xticks([0, 1], ["uninf", "inf"], fontsize=5.5)
    axz.set_ylim(0, 110)
    axz.set_yticks([0, 50, 100])
    axz.set_ylabel("% of cells with no\nwMel UMIs at all", fontsize=6)
    axz.set_title("No signal", pad=4, fontsize=6.5)
    panel_letter(axz, "A", dx=-0.58)
    pos = obs[obs["titer"] > 0]
    pu = pos["infection_status"] == "uninf"
    bins = np.linspace(np.log10(pos["titer"].min()), np.log10(pos["titer"].max()), 55)
    for m, c, lab in ((pu, COL["grey"], "uninfected"), (~pu, COL["green"], "infected")):
        ax.hist(np.log10(pos.loc[m, "titer"]), bins=bins, density=True, color=c,
                alpha=0.65, edgecolor="none", label=lab)
    ax.axvline(np.log10(floor), color=COL["vermillion"], lw=0.9, ls="--")
    ax.set_xlabel("per-cell wMel titer (log$_{10}$)")
    ax.set_ylabel("density of cells")
    ax.set_title("Infection, not background", pad=4)
    ax.legend(loc="upper left", fontsize=6)
    note(ax, "the cells that carry any signal; the rest are in the bar at left",
         y=-0.30)

    # ---- B: AUROC per cell type, one point each ----
    ax = fig.add_subplot(cell(1))
    sep = pd.read_csv(args.separation, sep="\t")
    sep = sep[~sep["cell_type"].isin(drop)]
    sep["auroc"] = pd.to_numeric(sep["auroc"], errors="coerce")
    sep = sep.sort_values("auroc")
    y = np.arange(len(sep))
    ax.scatter(sep["auroc"], y, s=16, color=COL["blue"], zorder=3,
               edgecolors=COL["black"], linewidths=0.3)
    ax.axvline(0.5, color=COL["black"], lw=0.8, zorder=2)
    ax.text(0.5, len(sep) - 0.3, " chance", fontsize=5.5, color=COL["grey"], va="top")
    ax.set_yticks(y, ct_short(sep["cell_type"]), fontsize=5)
    ax.set_xlim(0.45, 1.03)
    ax.set_xlabel("AUROC, titer vs infection")
    ax.set_title("In every cell type", pad=4)
    panel_letter(ax, "B", dx=-0.42)

    # ---- C: detection vs expression, by length tercile ----
    ax = fig.add_subplot(cell(2))
    g = pd.read_csv(args.gene_detection, sep="\t", index_col=0)
    for c in ("detection_frac", "log10_cpm", "log10_length"):
        g[c] = pd.to_numeric(g[c], errors="coerce")
    g = g.dropna(subset=["detection_frac", "log10_cpm", "log10_length"])
    g = g[g["detection_frac"] > 0]
    g["tercile"] = pd.qcut(g["log10_length"], 3, labels=["short", "mid", "long"])
    cols = {"short": COL["skyblue"], "mid": COL["orange"], "long": COL["purple"]}
    for t, c in cols.items():
        s_ = g[g["tercile"] == t]
        ax.scatter(s_["log10_cpm"], s_["detection_frac"], s=2.5, alpha=0.28,
                   color=c, linewidths=0, rasterized=True)
    edges = np.quantile(g["log10_cpm"], np.linspace(0, 1, 9))
    for t, c in cols.items():
        s_ = g[g["tercile"] == t]
        mids, med = [], []
        for lo, hi in zip(edges[:-1], edges[1:]):
            m = (s_["log10_cpm"] >= lo) & (s_["log10_cpm"] < hi)
            if m.sum() >= 5:
                mids.append((lo + hi) / 2)
                med.append(s_.loc[m, "detection_frac"].median())
        ax.plot(mids, med, color=c, lw=1.5, zorder=4, label=t)
    # Detection spans four orders of magnitude; on a linear axis every locus but
    # the two rRNA genes sits on the floor.
    ax.set_yscale("log")
    ax.set_xlabel("expression (log$_{10}$ CPM)")
    ax.set_ylabel("fraction of cells detected")
    ax.set_title("Expression, not sequence", pad=4)
    ax.legend(title="locus length", loc="upper left", fontsize=5.5,
              title_fontsize=5.5, ncol=3, columnspacing=0.8)
    note(ax, "terciles coincide: length adds nothing", y=-0.30)
    panel_letter(ax, "C", dx=-0.17)

    # ---- D: ambient rises with depth ----
    if four:
        ax = fig.add_subplot(cell(3))
        u = obs[un]
        # Most uninfected cells carry zero wMel UMIs, so a scatter is a solid
        # line along the bottom with a thin haze above it, which reads as no
        # relationship or the wrong one. Binned means show the trend the
        # correlation is reporting.
        q = np.quantile(u["host_umis"], np.linspace(0, 1, 11))
        mids, mean, lo, hi = [], [], [], []
        for a, b in zip(q[:-1], q[1:]):
            m = (u["host_umis"] >= a) & (u["host_umis"] < b)
            if m.sum() < 20:
                continue
            v = u.loc[m, "wmel_umis"].to_numpy(float)
            se = v.std(ddof=1) / np.sqrt(len(v))
            mids.append(np.sqrt(max(a, 1) * max(b, 1)))
            mean.append(v.mean()); lo.append(v.mean() - 1.96 * se)
            hi.append(v.mean() + 1.96 * se)
        ax.fill_between(mids, lo, hi, color=COL["lightgrey"], zorder=2)
        ax.plot(mids, mean, "o-", color=COL["grey"], markersize=3, lw=1.2, zorder=3)
        ax.set_xscale("log")
        ax.set_xlabel("host UMIs per cell (decile midpoint)")
        ax.set_ylabel("mean wMel UMIs per cell")
        ax.set_title("What is left is ambient", pad=4)
        rs = stats.spearmanr(u["host_umis"], u["wmel_umis"])
        note(ax, f"uninfected cells only; \u03c1 = {rs.statistic:.3f}, band is 95% CI",
             y=-0.30)
        panel_letter(ax, "D", dx=-0.42)
        print(f"ambient: rho(host UMIs, wMel UMIs) in uninfected = "
              f"{rs.statistic:.4f} (p = {rs.pvalue:.3g})")

    save(fig, args.outdir, args.stem)
    print(f"\nA: zero wMel UMIs in {zu:.1%} of uninfected cells against "
          f"{zi:.1%} of infected; median titer among signal-carrying cells "
          f"{pos.loc[pu, 'titer'].median():.3e} against "
          f"{pos.loc[~pu, 'titer'].median():.3e}")
    print(f"B: AUROC {sep['auroc'].min():.3f} to {sep['auroc'].max():.3f} "
          f"over {len(sep)} cell types, all above chance")


if __name__ == "__main__":
    main()

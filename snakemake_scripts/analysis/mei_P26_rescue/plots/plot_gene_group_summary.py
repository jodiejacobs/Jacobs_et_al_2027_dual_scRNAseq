#!/usr/bin/env python
"""Figure 1: wMel capture across all samples, as one summary panel.

Combines the per-sample outputs of rule gene_group_stats, which otherwise
exist only as eight separate four-panel QC figures under
results/qc/gene_group_stats/<sample_id>/. Those are diagnostics, one library at
a time; this is the version for the paper, where the comparison between
samples is the point.

  A  per-cell wMel share, one violin per sample. The infected and uninfected
     groups are the whole argument for the method: uninfected libraries have
     no Wolbachia, so their distribution is the measured noise floor and
     anything the infected libraries show above it is signal.
  B  share of each library that is Wolbachia against Drosophila, from
     group_summary.csv. The number quoted in the methods comes from here.
  C  per-gene detection, genes pooled across samples and coloured by group,
     so the wMel genes can be seen against the host background rather than
     in isolation.
  D  per-cell wMel share against library size. This is the panel that shows
     whether the signal is a depth artefact: if wMel share rises with total
     counts in uninfected cells, that is ambient capture.

The per-sample figures are not replaced. They are the right thing to look at
when one library misbehaves, and this panel cannot show that.

Usage:
  python plot_gene_group_summary.py \
      --stats_dir results/qc/gene_group_stats \
      --samples OreR_wMel-1_10x,OreR_uninf-1_10x,... \
      --exclude_cells results/rescue/label_qc/debris_cells.tsv \
      --outdir results/rescue/figures
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import COL, WIDE, apply_style, note, panel_letter, save
import matplotlib.pyplot as plt

GCOL = {"OreR_uninf": COL["skyblue"], "OreR_wMel": COL["blue"],
        "MeiP26_uninf": COL["orange"], "MeiP26_wMel": COL["vermillion"]}
ORDER = ["OreR_uninf", "OreR_wMel", "MeiP26_uninf", "MeiP26_wMel"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--stats_dir", required=True,
                   help="results/qc/gene_group_stats")
    p.add_argument("--samples", required=True,
                   help="comma-separated sample ids, one subdirectory each")
    p.add_argument("--group", default="Wolbachia",
                   help="the gene group plotted as titer; matches the --prefix "
                        "label given to gene_group_stats.py")
    p.add_argument("--host_group", default="Drosophila")
    p.add_argument("--exclude_cells", default=None,
                   help="TSV with an obs_name column (label_qc/debris_cells.tsv). "
                        "cell_stats.csv is indexed by bare barcode, so the sample "
                        "id is prepended before matching")
    p.add_argument("--points_per_sample", type=int, default=1500,
                   help="cells drawn over each box in panel A; the box itself "
                        "is computed on every cell, so quartiles are exact")
    p.add_argument("--max_points", type=int, default=60000,
                   help="cells drawn in panel D; all are used for the statistics")
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig1_wmel_capture_summary")
    return p.parse_args()


def group_of(sample):
    g = "MeiP26" if "Mei" in sample else "OreR"
    i = "uninf" if "uninf" in sample else "wMel"
    return f"{g}_{i}"


def load(args):
    samples = [s.strip() for s in args.samples.split(",") if s.strip()]
    drop = set()
    if args.exclude_cells and os.path.exists(args.exclude_cells):
        d = pd.read_csv(args.exclude_cells, sep="\t")
        col = "obs_name" if "obs_name" in d.columns else d.columns[0]
        drop = set(d[col].astype(str))

    cells, summ, genes, dropped = [], [], [], 0
    for s in samples:
        base = os.path.join(args.stats_dir, s)
        c = pd.read_csv(os.path.join(base, "cell_stats.csv"), index_col=0)
        c["sample"] = s
        if drop:
            keyed = pd.Index([f"{s}__{b}" for b in c.index])
            keep = ~keyed.isin(drop)
            dropped += int((~keep).sum())
            c = c[keep]
        cells.append(c)
        gs = pd.read_csv(os.path.join(base, "group_summary.csv"))
        gs["sample"] = s
        summ.append(gs)
        gp = os.path.join(base, "gene_stats.csv")
        if os.path.exists(gp):
            g = pd.read_csv(gp)
            g["sample"] = s
            genes.append(g)

    cells = pd.concat(cells)
    cells["group"] = cells["sample"].map(group_of)
    if drop and dropped == 0:
        raise SystemExit(
            f"none of the {len(drop):,} excluded barcodes matched any "
            "cell_stats.csv after prepending the sample id. Either these "
            "stats are already debris-free, in which case drop "
            "--exclude_cells, or the barcode format differs. Compare a few "
            "names on each side before rerunning.")
    print(f"{len(samples)} samples, {len(cells):,} cells"
          + (f", {dropped:,} debris removed" if drop else ""))
    return cells, pd.concat(summ), (pd.concat(genes) if genes else None)


def main():
    args = parse_args()
    apply_style()
    cells, summ, genes = load(args)
    pct, host = f"pct_{args.group}", f"pct_{args.host_group}"
    if pct not in cells.columns:
        raise SystemExit(f"{pct} not in cell_stats.csv ({list(cells.columns)}); "
                         "check --group matches the label given to "
                         "gene_group_stats.py")

    order = [s for g in ORDER for s in sorted(
        {x for x in cells["sample"] if group_of(x) == g})]
    fig = plt.figure(figsize=(WIDE, 5.0))
    gs = fig.add_gridspec(2, 2, hspace=0.62, wspace=0.30)

    # ---- A: per-cell wMel share, by sample ----
    ax = fig.add_subplot(gs[0, 0])
    data = [cells.loc[cells["sample"] == s, pct].dropna().to_numpy() for s in order]
    floor = 1e-3
    logd = [np.log10(np.clip(v, floor, None)) for v in data]
    # Points first, box over them. Every cell cannot be drawn (20,000+ a sample
    # renders as a solid bar and hides the distribution), so a fixed number is
    # sampled per sample and the count is stated under the panel. The box is
    # computed on all cells regardless, so the quartiles are not a subsample.
    rng = np.random.default_rng(0)
    for i, (s, v) in enumerate(zip(order, logd)):
        d = v if len(v) <= args.points_per_sample else rng.choice(
            v, args.points_per_sample, replace=False)
        ax.scatter(i + rng.uniform(-0.22, 0.22, len(d)), d, s=0.8, alpha=0.25,
                   linewidths=0, color=GCOL.get(group_of(s), COL["grey"]),
                   rasterized=True, zorder=2)
    bp = ax.boxplot(logd, positions=np.arange(len(order)), widths=0.52,
                    showfliers=False, patch_artist=True, zorder=3,
                    medianprops=dict(color=COL["black"], linewidth=1.1),
                    whiskerprops=dict(color=COL["black"], linewidth=0.7),
                    capprops=dict(color=COL["black"], linewidth=0.7))
    for patch, s in zip(bp["boxes"], order):
        patch.set_facecolor(GCOL.get(group_of(s), COL["grey"]))
        patch.set_alpha(0.45)
        patch.set_edgecolor(COL["black"])
        patch.set_linewidth(0.7)
    ax.set_xticks(range(len(order)),
                  [s.replace("_10x", "").replace("_", " ") for s in order],
                  rotation=40, ha="right", fontsize=5.5)
    ax.set_ylabel(f"% {args.group} per cell (log$_{{10}}$)")
    ax.set_title("Per-cell wMel share, by sample", pad=4)
    note(ax, f"box on all cells; {args.points_per_sample:,} drawn per sample; "
             f"zeros shown at {floor}%", y=-0.36)
    panel_letter(ax, "A", dx=-0.22)

    # ---- B: share of the library, per sample ----
    ax = fig.add_subplot(gs[0, 1])
    piv = summ.pivot_table(index="sample", columns="group",
                           values="pct_of_library", aggfunc="sum").reindex(order)
    x = np.arange(len(order))
    if args.group in piv.columns:
        ax.bar(x, piv[args.group], width=0.7,
               color=[GCOL.get(group_of(s), COL["grey"]) for s in order])
    ax.set_xticks(x, [s.replace("_10x", "").replace("_", " ") for s in order],
                  rotation=40, ha="right", fontsize=5.5)
    ax.set_ylabel(f"% of library that is {args.group}")
    ax.set_title("Library composition", pad=4)
    panel_letter(ax, "B", dx=-0.22)

    # ---- C: per-gene detection, pooled ----
    ax = fig.add_subplot(gs[1, 0])
    if genes is not None:
        agg = (genes.groupby(["gene", "group"], observed=True)
               .agg(total_counts=("total_counts", "sum"),
                    pct_cells=("pct_cells_expressing", "mean")).reset_index())
        for lab, col in ((args.host_group, COL["lightgrey"]),
                         (args.group, COL["vermillion"])):
            m = agg["group"] == lab
            if m.any():
                ax.scatter(agg.loc[m, "total_counts"] + 1, agg.loc[m, "pct_cells"],
                           s=3, alpha=0.45, color=col, linewidths=0,
                           label=f"{lab} ({int(m.sum()):,})", rasterized=True)
        ax.set_xscale("log")
        ax.set_xlabel("total counts + 1, summed across samples")
        ax.set_ylabel("% of cells expressing")
        ax.legend(loc="upper left", fontsize=6, markerscale=2.5)
    else:
        ax.text(0.5, 0.5, "gene_stats.csv not found", ha="center", va="center",
                transform=ax.transAxes, color=COL["grey"])
        ax.set_axis_off()
    ax.set_title("Per-gene detection", pad=4)
    panel_letter(ax, "C", dx=-0.22)

    # ---- D: wMel share against library size ----
    ax = fig.add_subplot(gs[1, 1])
    sub = cells
    if len(sub) > args.max_points:
        sub = sub.sample(args.max_points, random_state=0)
    for g in ORDER:
        m = sub["group"] == g
        if m.any():
            ax.scatter(sub.loc[m, "total_counts"] + 1,
                       np.clip(sub.loc[m, pct], floor, None),
                       s=1.5, alpha=0.35, linewidths=0, color=GCOL[g],
                       label=g.replace("_", " "), rasterized=True)
    ax.set_xscale("log")
    ax.set_yscale("log")
    ax.set_xlabel("total counts per cell + 1")
    ax.set_ylabel(f"% {args.group} per cell")
    ax.set_title("Is the signal a depth artefact?", pad=4)
    # The lower left is the densest part of this panel, so the legend needs a
    # background or it sits unreadable on top of the points.
    leg = ax.legend(loc="lower left", fontsize=5.5, markerscale=5, ncol=2,
                    frameon=True, framealpha=0.92, borderpad=0.4)
    leg.get_frame().set_edgecolor("none")
    leg.get_frame().set_facecolor("white")
    if len(sub) < len(cells):
        note(ax, f"{len(sub):,} of {len(cells):,} cells drawn", y=-0.30)
    panel_letter(ax, "D", dx=-0.22)

    save(fig, args.outdir, args.stem)

    # The numbers a methods section quotes, printed so the log records them.
    out = (cells.groupby("sample")[pct]
           .agg(n="size", median="median", mean="mean",
                q99=lambda v: v.quantile(0.99)).reindex(order))
    out["group"] = [group_of(s) for s in out.index]
    print(out.round(4).to_string())
    inf = cells.loc[~cells["group"].str.contains("uninf"), pct]
    un = cells.loc[cells["group"].str.contains("uninf"), pct]
    if len(un) and len(inf):
        print(f"\nuninfected median {un.median():.4f}%, 99th {un.quantile(.99):.4f}%; "
              f"infected median {inf.median():.4f}%, 99th {inf.quantile(.99):.4f}%")
        print("Quote the uninfected 99th percentile as the floor, not zero: "
              "uninfected libraries carry no Wolbachia, so whatever they show "
              "is the assay's own noise.")


if __name__ == "__main__":
    main()

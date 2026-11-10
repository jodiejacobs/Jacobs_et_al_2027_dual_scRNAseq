#!/usr/bin/env python
"""Figure 2: per-locus wMel detection against expression and sequence features.

Plots results/rescue/detection_bias/wmel_gene_detection.tsv.gz, the per-locus
table the feature tests and the detection model are computed from.

The two companion panels report coefficients. This one shows the data they
came from, which is the only way a reader can see that the expression-detection
relationship is a saturating curve rather than a line, and that the residual
spread at fixed expression is wide. Both matter: the curve is why a raw test
against any expression-correlated feature returns a significant result, and
the spread is why per-locus detection should not be treated as interchangeable
between loci.

  A  detection against expression, the dominant relationship, rRNA loci marked
     because they carry most of the read budget and sit at the ceiling
  B  detection against length, coloured by expression, so the confound is
     visible as a colour gradient rather than asserted
  C  residual detection after expression, against each remaining feature

Two versions are written. "robust" keeps loci detected in at least --min_cells
cells, where the detection fraction is estimated from enough observations to
mean something; "all" keeps every locus with counts and coordinates.

Usage:
  python plot_wmel_gene_detection.py \
      --input results/rescue/detection_bias/wmel_gene_detection.tsv.gz \
      --outdir figures/fig2
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import empty_figure, COL, WIDE, apply_style, panel_letter, save
import matplotlib.pyplot as plt

NICE = {"log10_length": "length (log$_{10}$ bp)", "gc": "GC content",
        "rel_position": "position on chromosome"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig2_wmel_gene_detection")
    p.add_argument("--min_cells", type=int, default=10,
                   help="cells a locus must be seen in for the robust version")
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def residualise(y, x):
    """Rank residuals of y on x, matching how the feature tests condition."""
    ry, rx = stats.rankdata(y).astype(float), stats.rankdata(x).astype(float)
    Z = np.column_stack([np.ones_like(rx), rx])
    return ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]


def build(d, args, suffix, note_txt):
    feats = [f for f in ("log10_length", "gc", "rel_position") if f in d.columns
             and d[f].notna().any()]
    print(f"\n[{suffix}] {len(d)} loci, features: {', '.join(feats) or 'none'}; "
          f"{note_txt}")

    fig = plt.figure(figsize=(WIDE, 2.5))
    gs = fig.add_gridspec(1, 3, wspace=0.40, width_ratios=[1.0, 1.0, 1.25])

    # ---- A: detection against expression ----
    ax = fig.add_subplot(gs[0, 0])
    is_rrna = (d["product"].astype(str).str.contains("ribosomal RNA|rRNA|23S|16S|5S",
                                                     case=False, na=False)
               if "product" in d.columns else pd.Series(False, index=d.index))
    ax.scatter(d.loc[~is_rrna, "log10_cpm"], d.loc[~is_rrna, "detection_frac"],
               s=5, alpha=0.5, color=COL["blue"], edgecolors="none", zorder=2)
    if is_rrna.any():
        ax.scatter(d.loc[is_rrna, "log10_cpm"], d.loc[is_rrna, "detection_frac"],
                   s=22, color=COL["vermillion"], edgecolors=COL["black"],
                   linewidths=0.4, zorder=3, label=f"rRNA ({int(is_rrna.sum())})")
        ax.legend(loc="lower right", fontsize=6)
    r = stats.spearmanr(d["log10_cpm"], d["detection_frac"])
    ax.set_xlabel("expression (log$_{10}$ CPM)")
    ax.set_ylabel("fraction of cells detected")
    ax.set_title(f"Detection tracks expression\n(ρ = {r.statistic:.2f})",
                 pad=4, fontsize=7.5)
    panel_letter(ax, "A", dx=-0.32)

    # ---- B: detection against length, coloured by expression ----
    ax = fig.add_subplot(gs[0, 1])
    if "log10_length" in d.columns:
        sc = ax.scatter(d["log10_length"], d["detection_frac"], s=7,
                        c=d["log10_cpm"], cmap="viridis", edgecolors="none")
        cb = fig.colorbar(sc, ax=ax, fraction=0.045, pad=0.02)
        cb.set_label("log$_{10}$ CPM", fontsize=6)
        cb.ax.tick_params(labelsize=5.5)
        rl = stats.spearmanr(d["log10_length"], d["detection_frac"])
        ax.set_xlabel("length (log$_{10}$ bp)")
        ax.set_ylabel("fraction of cells detected")
        ax.set_title(f"Length is expression in disguise\n(raw ρ = {rl.statistic:.2f})",
                     pad=4, fontsize=7.5)
    panel_letter(ax, "B", dx=-0.32)

    # ---- C: residual detection against each feature ----
    ax = fig.add_subplot(gs[0, 2])
    res = residualise(d["detection_frac"], d["log10_cpm"])
    rows = []
    for f in feats:
        v = d[f].to_numpy(float)
        ok = np.isfinite(v)
        if ok.sum() < 10:
            continue
        rf = residualise(v[ok], d["log10_cpm"].to_numpy()[ok])
        rr, pp = stats.pearsonr(rf, res[ok])
        rows.append({"feature": f, "partial": rr, "p": pp, "n": int(ok.sum())})
    t = pd.DataFrame(rows)
    if len(t):
        t = t.reindex(t["partial"].abs().sort_values().index)
        y = np.arange(len(t))
        ax.barh(y, t["partial"], height=0.6,
                color=[COL["vermillion"] if p < 0.05 else COL["lightgrey"]
                       for p in t["p"]])
        ax.axvline(0, color=COL["black"], lw=0.6)
        ax.set_yticks(y, [NICE.get(f, f) for f in t["feature"]], fontsize=6.5)
        ax.set_xlabel("partial ρ with detection,\nexpression held fixed", fontsize=7)
        ax.set_title("What survives conditioning", pad=4, fontsize=7.5)
        ax.margins(x=0.26)
    panel_letter(ax, "C", dx=-0.42)

    fig.suptitle(note_txt, y=1.05, fontsize=6.5, color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")
    if len(t):
        print(t.round(4).to_string(index=False))


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t", index_col=0)
    for c in ("detection_frac", "log10_cpm", "log10_length", "gc", "rel_position",
              "n_cells", "cpm"):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    d = d[d["detection_frac"].notna() & d["log10_cpm"].notna()]
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        if v == "robust":
            sub = d[d["n_cells"] >= args.min_cells]
            txt = (f"loci detected in ≥{args.min_cells} cells, where the detection "
                   "fraction rests on enough observations to interpret")
        else:
            sub = d
            txt = "every locus with counts and coordinates"
        if len(sub) < 10:
            empty_figure(args.outdir, f"{args.stem}_{v}",
                         f"No rows survive the {v} filter.\n"
                         f"An empty set here is a result, not a failure:\n"
                         f"say so in the text rather than citing the\n"
                         f"unfiltered figure as though it passed.")
            continue
        build(sub, args, v, txt)


if __name__ == "__main__":
    main()

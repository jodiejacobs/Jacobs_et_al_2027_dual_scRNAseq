#!/usr/bin/env python
"""Figure 1: atlas-projected UMAP and per-sample composition, real cells only.

Unlike every other script in this directory, this one reads an AnnData object
rather than a results table, because the UMAP coordinates exist nowhere else:
they live in results/integrated/integrated.h5ad as obsm["X_umap_atlas"].

Why it exists. The equivalent panels in figures/ (umap_annotation.pdf,
composition_by_sample_annotation.pdf) are written by plot_diagnostics() in
integrate_via_atlas_projection.py, which runs BEFORE 06_label_qc.py has called
debris, so they show all barcodes. Everything downstream of rescue_pseudobulk
is already debris-free, because that rule passes --exclude_cells, which leaves
those two panels as the only figures in the paper that still include debris.
This script regenerates them from the same object with the debris list applied,
so the figure and the analysis describe the same set of cells.

The all-cells versions are not replaced. Keep both: the difference between them
is the debris filter's effect, which is worth being able to show.

Usage:
  python plot_umap_annotation.py \
      --integrated results/integrated/integrated.h5ad \
      --exclude_cells results/rescue/label_qc/debris_cells.tsv \
      --label_col cell_type_annotation \
      --outdir results/rescue/figures
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import COL, WIDE, apply_style, ct_short, note, panel_letter, save
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--integrated", required=True,
                   help="results/integrated/integrated.h5ad")
    p.add_argument("--exclude_cells", default=None,
                   help="TSV with an obs_name column (label_qc/debris_cells.tsv)")
    p.add_argument("--label_col", default="cell_type_annotation",
                   help="matches rescue_label_col in config.yaml; the object's "
                        "column is atlas_<label_col>")
    p.add_argument("--sample_col", default="source_file")
    p.add_argument("--umap_key", default="X_umap_atlas")
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig1_umap_annotation")
    p.add_argument("--comp_stem", default="fig1_composition_by_sample")
    p.add_argument("--min_frac", type=float, default=0.0,
                   help="pool labels below this share of cells into 'other'. "
                        "Default 0, which pools nothing: on the fine annotation "
                        "even 0.005 would hide 11 of 21 cell types, escort cells "
                        "and germ cell stage 4+ among them. Raise it only for the "
                        "broad annotation, or when a legend is genuinely unreadable")
    return p.parse_args()


def load(args):
    import anndata as ad
    a = ad.read_h5ad(args.integrated)
    n0 = a.n_obs
    dropped = 0
    if args.exclude_cells and os.path.exists(args.exclude_cells):
        d = pd.read_csv(args.exclude_cells, sep="\t")
        col = "obs_name" if "obs_name" in d.columns else d.columns[0]
        names = set(d[col].astype(str))
        keep = ~a.obs_names.isin(names)
        dropped = int((~keep).sum())
        if dropped == 0:
            # The object may already be debris-free, or the names may not match.
            # Those look identical here and mean very different things, so say
            # which one it is rather than reporting a silent no-op.
            overlap = len(set(a.obs_names) & names)
            raise SystemExit(
                f"none of the {len(names):,} excluded barcodes are in "
                f"{os.path.basename(args.integrated)} ({n0:,} cells, "
                f"{overlap} name overlap). Either this object is already "
                "debris-free, in which case drop --exclude_cells, or the "
                "obs_name formats differ between the two files. Check "
                f"a few names on each side before rerunning.")
        a = a[keep].copy()
    print(f"{n0:,} cells, {dropped:,} debris removed, {a.n_obs:,} real cells kept")
    return a


def main():
    args = parse_args()
    apply_style()
    a = load(args)

    lab = f"atlas_{args.label_col}"
    if lab not in a.obs.columns:
        lab = args.label_col if args.label_col in a.obs.columns else None
    if lab is None:
        raise SystemExit(f"neither atlas_{args.label_col} nor {args.label_col} "
                         f"is in obs: {list(a.obs.columns)[:20]}")
    if args.umap_key not in a.obsm:
        raise SystemExit(f"{args.umap_key} not in obsm ({list(a.obsm.keys())}); "
                         "this is not the integrated object")

    xy = np.asarray(a.obsm[args.umap_key])[:, :2]
    lv = a.obs[lab].astype(str)
    frac = lv.value_counts(normalize=True)
    big = [k for k in frac.index if frac[k] >= args.min_frac]
    shown = lv.where(lv.isin(big), "other")
    order = [k for k in frac.index if k in big] + (["other"] if (shown == "other").any() else [])
    # tab20 repeats after 20 colours and the fine annotation has 21 cell types,
    # so the 21st would be indistinguishable from the 1st. Two maps give 40.
    pal = [c for cm in ("tab20", "tab20b") for c in plt.get_cmap(cm).colors]
    colour = {k: (COL["lightgrey"] if k == "other" else pal[i % len(pal)])
              for i, k in enumerate(order)}

    # ---- UMAP ----
    fig, ax = plt.subplots(figsize=(WIDE * 0.62, 4.2))
    # Draw most abundant first so the rarest land on top. Shuffling would bury
    # them: several fine labels hold under 100 cells against 35,000 in the
    # largest, and a population drawn under that is simply not visible.
    sv = shown.to_numpy()
    for k in order:
        m = sv == k
        if not m.any():
            continue
        ax.scatter(xy[m, 0], xy[m, 1], s=1.2, linewidths=0, alpha=0.6,
                   color=colour[k], rasterized=True)
    ax.set_xticks([])
    ax.set_yticks([])
    for s in ("left", "bottom"):
        ax.spines[s].set_visible(False)
    ax.set_xlabel("UMAP 1")
    ax.set_ylabel("UMAP 2")
    ax.set_title(f"Atlas-projected cell types, {a.n_obs:,} real cells", pad=4)
    h = [plt.Line2D([], [], marker="o", ls="", markersize=4, color=colour[k])
         for k in order]
    ax.legend(h, ct_short(order), loc="center left", bbox_to_anchor=(1.01, 0.5),
              fontsize=5.5, handletextpad=0.3, labelspacing=0.28, frameon=False)
    note(ax, "debris removed; the all-cells version is figures/umap_annotation.pdf",
         y=-0.07)
    panel_letter(ax, "A", dx=-0.10)
    save(fig, args.outdir, args.stem)

    # ---- composition per sample ----
    samp = (a.obs[args.sample_col].astype(str) if args.sample_col in a.obs.columns
            else pd.Series("all", index=a.obs_names))
    counts = pd.crosstab(samp, shown)
    comp = counts.div(counts.sum(axis=1), axis=0)
    comp = comp.reindex(columns=[k for k in order if k in comp.columns])
    os.makedirs(args.outdir, exist_ok=True)
    csv = os.path.join(args.outdir, f"{args.comp_stem}.csv")
    comp.to_csv(csv)
    print(f"  wrote {csv}")

    fig, ax = plt.subplots(figsize=(WIDE * 0.7, 3.4))
    left = np.zeros(len(comp))
    y = np.arange(len(comp))
    for k in comp.columns:
        v = comp[k].to_numpy()
        ax.barh(y, v, left=left, height=0.72, color=colour[k], label=str(k))
        left += v
    ax.set_yticks(y, [s.replace("_", " ") for s in comp.index], fontsize=6)
    ax.set_xlim(0, 1)
    ax.set_xlabel("fraction of real cells")
    ax.set_title("Cell-type composition by sample", pad=4)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=5.5,
              handletextpad=0.3, labelspacing=0.28, frameon=False)
    panel_letter(ax, "B", dx=-0.34)
    save(fig, args.outdir, args.comp_stem)

    print(f"\nLabels drawn: {len(big)} at or above {args.min_frac:.1%} of cells"
          + (f", {int((shown == 'other').sum()):,} cells pooled into 'other'"
             if (shown == "other").any() else ""))
    print(comp.round(4).to_string())


if __name__ == "__main__":
    main()

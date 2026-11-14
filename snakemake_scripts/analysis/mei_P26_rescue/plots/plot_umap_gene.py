#!/usr/bin/env python
"""Figure 4: a gene's expression on the atlas UMAP, one panel per group.

Reads results/integrated/integrated.h5ad, because the UMAP coordinates live
only in obsm["X_umap_atlas"]. Like plot_umap_annotation.py it reads the stored
projection and never recomputes one, so the cells sit where the atlas put them
and the four panels are directly comparable.

Default gene is mei-P26. The point of splitting by genotype and infection
rather than overlaying them is that the claim is about where expression is
lost: the hypomorph depletes mei-P26 across the whole tissue and not only in
the germline, and wMel does not bring it back in either genotype. All four
panels share one colour scale and one pair of axis limits, because the whole
comparison is between panels and any per-panel rescaling would destroy it.

Every cell is drawn in light grey behind each panel, so a panel with little
expression still shows the tissue and cannot be mistaken for missing cells.
Expressing cells are drawn last, brightest on top, since at this density a
few positive cells underneath 30,000 negative ones are invisible.

Expression scale. integrated.h5ad may hold raw counts, normalised values or
scaled ones depending on how the object was written, and the figure means
different things in each case. --normalize auto inspects the matrix: if it is
non-negative integers it is treated as counts and normalised per cell then
log1p transformed, otherwise it is used as it stands. Whatever it decides is
printed, and the colour bar is labelled accordingly. Pass --normalize
explicitly if the object is something the check cannot recognise.

Usage:
  python plot_umap_gene.py \
      --integrated results/integrated/integrated.h5ad \
      --exclude_cells results/rescue/label_qc/debris_cells.tsv \
      --gene FBgn0026206 --gene_label mei-P26 \
      --outdir results/rescue/figures
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (COL, WIDE, apply_style, note, panel_letter, sample_group,
                      save)
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

GROUP_ORDER = ["OreR_uninf", "OreR_wMel", "MeiP26_uninf", "MeiP26_wMel"]
NICE = {"OreR_uninf": "OreR", "OreR_wMel": "OreR + wMel",
        "MeiP26_uninf": "mei-P26", "MeiP26_wMel": "mei-P26 + wMel"}
# Light grey through to a strong blue. A sequential ramp that starts near the
# background keeps zero-expressing cells from reading as a category of their
# own while still letting the high end carry.
CMAP = LinearSegmentedColormap.from_list(
    "expr", ["#E8E8E8", "#BFD9EC", "#56B4E9", "#0072B2", "#08306B"])


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--integrated", required=True)
    p.add_argument("--exclude_cells", default=None,
                   help="TSV with an obs_name column "
                        "(label_qc/debris_cells.tsv). Required: the integrated "
                        "object still contains debris, because it is written "
                        "before 06_label_qc.py calls it")
    p.add_argument("--allow_debris", action="store_true",
                   help="run without --exclude_cells. Only for deliberately "
                        "showing what the debris filter removes; the figure "
                        "then does not describe the same cells as every other "
                        "panel in the paper")
    p.add_argument("--gene", default="FBgn0026206",
                   help="FlyBase id, or a symbol if var carries one")
    p.add_argument("--gene_label", default="mei-P26")
    p.add_argument("--umap_key", default="X_umap_atlas")
    p.add_argument("--sample_col", default="source_file",
                   help="used to derive the genotype by infection group when "
                        "genotype/infection_status are not both in obs")
    p.add_argument("--expression", default=None,
                   help="h5ad to read expression from, when the integrated "
                        "object does not carry the gene. integrated.h5ad is "
                        "subset to the highly variable genes during "
                        "integration, so most genes are absent from it; "
                        "results/rescue/pseudobulk/counts_labeled.h5ad keeps "
                        "all of them. Cells become the intersection of the "
                        "two objects, matched on obs_names")
    p.add_argument("--layer", default=None,
                   help="layer to read expression from; default .X")
    p.add_argument("--normalize", choices=("auto", "none", "lognorm"),
                   default="auto")
    p.add_argument("--vmax_quantile", type=float, default=0.995,
                   help="colour scale top, as a quantile of the expressing "
                        "cells. A handful of very high cells otherwise set the "
                        "scale and flatten everything else to the floor")
    p.add_argument("--point_size", type=float, default=1.1)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig4_umap_meip26")
    return p.parse_args()


def load(args):
    import anndata as ad
    import scipy.sparse as sp
    if not args.exclude_cells and not args.allow_debris:
        raise SystemExit(
            "--exclude_cells is required. integrated.h5ad is written by "
            "integrate_via_atlas_projection.py, which runs before "
            "06_label_qc.py calls debris, so it still holds every barcode. "
            "Pass results/rescue/label_qc/debris_cells.tsv, or --allow_debris "
            "if you deliberately want the unfiltered version.")
    if args.exclude_cells and not os.path.exists(args.exclude_cells):
        raise SystemExit(f"--exclude_cells {args.exclude_cells} does not exist")
    a = ad.read_h5ad(args.integrated, backed="r")
    n0 = a.n_obs
    if args.exclude_cells:
        d = pd.read_csv(args.exclude_cells, sep="\t")
        col = "obs_name" if "obs_name" in d.columns else d.columns[0]
        names = set(d[col].astype(str))
        keep = ~a.obs_names.isin(names)
        if int((~keep).sum()) == 0:
            raise SystemExit(
                f"none of the {len(names):,} excluded barcodes are in this "
                f"object ({n0:,} cells). Either it is already debris-free, in "
                "which case drop --exclude_cells, or the obs_name formats "
                "differ. Compare a few names on each side before rerunning.")
        a = a[keep]
    print(f"{n0:,} cells, {n0 - a.n_obs:,} debris removed, {a.n_obs:,} kept")

    if args.umap_key not in a.obsm:
        raise SystemExit(f"{args.umap_key} not in obsm ({list(a.obsm.keys())}); "
                         "this is not the integrated object")
    xy = np.asarray(a.obsm[args.umap_key])[:, :2]
    obs = a.obs.copy()
    obs_names = pd.Index(a.obs_names.astype(str))

    # Where the expression comes from. integrated.h5ad is cut down to the
    # highly variable genes when the integration runs, so most genes are
    # simply not in it and a lookup there fails for reasons that have nothing
    # to do with the gene being absent from the data. --expression points at
    # an object that kept every gene, usually
    # rescue/pseudobulk/counts_labeled.h5ad, and the cells become the
    # intersection of the two.
    e = a
    if args.expression:
        if not os.path.exists(args.expression):
            raise SystemExit(f"--expression {args.expression} does not exist")
        e = ad.read_h5ad(args.expression, backed="r")
        en = pd.Index(e.obs_names.astype(str))
        common = obs_names.intersection(en)
        if len(common) < 0.05 * min(len(obs_names), len(en)):
            raise SystemExit(
                f"only {len(common):,} barcodes are shared between the two "
                f"objects ({len(obs_names):,} with coordinates, {len(en):,} "
                "with expression). The obs_name formats probably differ; "
                f"compare {list(obs_names[:2])} against {list(en[:2])}.")
        print(f"expression from {os.path.basename(args.expression)}: "
              f"{len(en):,} cells x {e.n_vars:,} genes; "
              f"{len(common):,} shared with the coordinates")
        ia = obs_names.get_indexer(common)
        xy = xy[ia]
        obs = obs.iloc[ia]
        e = e[en.get_indexer(common)]
        obs_names = common

    # The gene, by id or by whatever symbol column var carries.
    name = args.gene
    if name not in e.var_names:
        for c in ("symbol", "gene_symbol", "gene_name", "Symbol"):
            if c in e.var.columns:
                hit = e.var.index[e.var[c].astype(str) == name]
                if len(hit):
                    name = hit[0]
                    break
    if name not in e.var_names:
        where = args.expression or args.integrated
        extra = ("" if args.expression else
                 "\nintegrated.h5ad is usually subset to the highly variable "
                 "genes, so most genes are missing from it. Pass --expression "
                 "results/rescue/pseudobulk/counts_labeled.h5ad, which keeps "
                 "every gene.")
        raise SystemExit(f"{args.gene} is not in var_names of "
                         f"{os.path.basename(where)} and no symbol column "
                         f"matched it. var has {e.n_vars:,} entries, columns "
                         f"{list(e.var.columns)[:8]}.{extra}")

    def dense(m):
        return (np.asarray(m.todense()) if sp.issparse(m)
                else np.asarray(m))

    sub = e[:, [e.var_names.get_loc(name)]]
    v = dense(sub.layers[args.layer] if args.layer else sub.X).ravel()

    mode = args.normalize
    if mode == "auto":
        # Judge the matrix, not this one gene: a lowly expressed gene is all
        # small integers under any encoding, and a magnitude threshold would
        # misread it. Counts are non-negative integers; normalised values are
        # not integers and scaled values go negative. A random block is enough,
        # and going through AnnData indexing keeps this working on a
        # backed object, where X cannot be sliced twice directly.
        rng = np.random.default_rng(0)
        rows = np.sort(rng.choice(e.n_obs, min(2000, e.n_obs), replace=False))
        cols = np.sort(rng.choice(e.n_vars, min(200, e.n_vars), replace=False))
        blk_ad = e[rows][:, cols]
        blk = dense(blk_ad.layers[args.layer] if args.layer else blk_ad.X).ravel()
        blk = blk[np.isfinite(blk)]
        looks_counts = bool(len(blk) and blk.min() >= 0
                            and np.array_equal(blk, np.round(blk)))
        mode = "lognorm" if looks_counts else "none"
        what = "raw counts" if looks_counts else "already transformed values"
        print(f"--normalize auto: sampled {len(blk):,} entries, matrix looks "
              f"like {what}, so using {mode}")

    eobs = e.obs
    if mode == "lognorm":
        # Per-cell totals. Summing the whole matrix on a backed object reads
        # every byte off disk, so use a stored depth column when the object
        # carries one. counts_labeled.h5ad has host_umis, which is exactly
        # this quantity for the host side.
        depth = None
        for c in ("host_umis", "total_counts", "n_counts"):
            if c in eobs.columns:
                depth = pd.to_numeric(eobs[c], errors="coerce").to_numpy(float)
                print(f"per-cell depth from obs['{c}']")
                break
        if depth is None:
            print("no depth column in obs, summing the matrix (slow on a "
                  "backed object)")
            depth = np.asarray(dense(e.X).sum(axis=1)).ravel()
        depth = np.where(np.isfinite(depth) & (depth > 0), depth, 1.0)
        v = np.log1p(v / depth * 1e4)
        unit = "log1p(CP10K)"
    else:
        unit = "value in the object"

    # Group: prefer the explicit columns, fall back to the sample id. Look in
    # the expression object first when there is one, since it is the object
    # the per-cell values came from.
    o = eobs if (args.expression and
                 ({"genotype", "infection_status"} <= set(eobs.columns)
                  or args.sample_col in eobs.columns)) else obs
    if {"genotype", "infection_status"} <= set(o.columns):
        g = (o["genotype"].astype(str).str.replace("-", "", regex=False) + "_"
             + o["infection_status"].astype(str))
    elif args.sample_col in o.columns:
        g = o[args.sample_col].astype(str).map(sample_group)
    else:
        raise SystemExit(
            "cannot work out the group: obs has neither genotype and "
            f"infection_status nor {args.sample_col}. obs columns: "
            f"{list(o.columns)[:15]}")
    return xy, v, np.asarray(g), unit


def main():
    args = parse_args()
    apply_style()
    xy, v, grp, unit = load(args)

    groups = [g for g in GROUP_ORDER if (grp == g).any()]
    groups += [g for g in pd.unique(grp) if g not in groups]
    if not groups:
        raise SystemExit("no cells in any recognised group")
    pos = v[v > 0]
    vmax = float(np.quantile(pos, args.vmax_quantile)) if len(pos) else 1.0
    if vmax <= 0:
        vmax = float(v.max()) or 1.0
    print(f"{args.gene_label}: {len(pos):,} of {len(v):,} cells nonzero; "
          f"colour scale 0 to {vmax:.2f} ({unit}, "
          f"q{args.vmax_quantile} of expressing cells)")

    nr, nc = (2, 2) if len(groups) == 4 else (1, len(groups))
    fig, axes = plt.subplots(nr, nc, figsize=(WIDE * 0.78, 3.0 * nr),
                             squeeze=False)
    xpad = (xy[:, 0].max() - xy[:, 0].min()) * 0.03
    ypad = (xy[:, 1].max() - xy[:, 1].min()) * 0.03
    for i, g in enumerate(groups):
        ax = axes[i // nc][i % nc]
        m = grp == g
        # All cells behind, so a quiet panel still shows the tissue.
        ax.scatter(xy[:, 0], xy[:, 1], s=args.point_size * 0.7, linewidths=0,
                   color="#F0F0F0", rasterized=True, zorder=1)
        o = np.argsort(v[m])          # brightest last, or they vanish
        sc = ax.scatter(xy[m][o, 0], xy[m][o, 1], s=args.point_size,
                        c=np.clip(v[m][o], 0, vmax), cmap=CMAP, vmin=0,
                        vmax=vmax, linewidths=0, rasterized=True, zorder=2)
        ax.set_xticks([]); ax.set_yticks([])
        for s in ax.spines.values():
            s.set_visible(False)
        ax.set_xlim(xy[:, 0].min() - xpad, xy[:, 0].max() + xpad)
        ax.set_ylim(xy[:, 1].min() - ypad, xy[:, 1].max() + ypad)
        frac = float((v[m] > 0).mean())
        ax.set_title(f"{NICE.get(g, g)}\n{int(m.sum()):,} cells, "
                     f"{frac:.0%} expressing", pad=3, fontsize=7)
        panel_letter(ax, "ABCD"[i] if i < 4 else str(i + 1), dx=-0.04, dy=1.02)
    for j in range(len(groups), nr * nc):
        axes[j // nc][j % nc].set_axis_off()

    cb = fig.colorbar(sc, ax=axes, fraction=0.030, pad=0.02, aspect=28)
    cb.set_label(f"{args.gene_label} expression, {unit}", fontsize=6)
    cb.ax.tick_params(labelsize=5.5)
    cb.outline.set_visible(False)
    note(axes[nr - 1][0],
         "one colour scale and one pair of axis limits across all panels; "
         "grey is every cell, for shape", y=-0.06)
    save(fig, args.outdir, args.stem)

    out = pd.DataFrame({
        "group": groups,
        "n_cells": [int((grp == g).sum()) for g in groups],
        "frac_expressing": [round(float((v[grp == g] > 0).mean()), 3)
                            for g in groups],
        "mean_expr": [round(float(v[grp == g].mean()), 3) for g in groups],
        "mean_expr_in_expressing": [
            round(float(v[(grp == g) & (v > 0)].mean()), 3)
            if ((grp == g) & (v > 0)).any() else float("nan") for g in groups]})
    print()
    print(out.to_string(index=False))
    print("\nfrac_expressing moves with sequencing depth as well as with "
          "biology, so read it beside mean_expr_in_expressing rather than on "
          "its own. The per-cell-type test is fig4_meip26_expression; this "
          "figure is for showing where the expression sits, not for testing.")


if __name__ == "__main__":
    main()

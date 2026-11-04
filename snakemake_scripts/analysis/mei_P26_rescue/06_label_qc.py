#!/usr/bin/env python
"""What are the 'artefact', 'unannotated' and 'adult tracheal cell' labels?

Atlas projection gave three labels that are hard to read biologically and
dominate the host DE counts: 'artefact' (17% -> 42% of Mei-P26 cells with
wMel), 'unannotated' (40% of OreR) and 'adult tracheal cell' (~20% of all
cells). This script checks every label against canonical ovary markers and
re-clusters the three problem labels on their own.

1. Lineage scores (scanpy score_genes on log-CP10k host genes) for germline,
   early germline, somatic follicle, late follicle (vitelline/chorion),
   muscle, trachea and hemocyte gene sets. Each cell gets the top-scoring
   lineage, or 'ambiguous' if no score exceeds --min_score.
2. Marker heatmap: mean log-CP10k per label x condition.
3. Leiden subclusters within each problem label (all cells, not subsampled),
   with per-subcluster QC (host UMIs, genes, % mito), lineage mix, top markers
   and cell counts per sample. Per-sample counts matter: with 2 replicates a
   subcluster driven by one library is a library effect, not biology.

Outputs (in --outdir):
  lineage_by_label.tsv        label x lineage cell fractions
  lineage_by_label_sample.tsv label x sample x lineage cell counts
  marker_heatmap.pdf/.png     mean marker expression per label x condition
  subclusters_<label>.tsv     per-subcluster QC, lineage mix, top markers
  subcluster_samples_<label>.tsv  cells per subcluster x sample
  subclusters_<label>.pdf     UMAP coloured by subcluster, lineage, sample, titer
  cell_lineage.tsv.gz         per-cell lineage call + scores
  debris_cells.tsv            obs_names of cells in debris subclusters (median
                              genes < --debris_max_genes AND median % mito >
                              --debris_min_mito); 01_build_pseudobulk.py
                              --exclude_cells drops them

Usage:
  python 06_label_qc.py --counts results/rescue/pseudobulk/counts_labeled.h5ad \
      --annotation reference/fbgn_annotation_ID_fb_2025_04.tsv.gz \
      --outdir results/rescue/label_qc
"""
import argparse
import gzip
import os

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns

# Canonical ovary markers. Genes absent from the matrix are dropped silently
# (printed once); scores need >= 2 genes present.
LINEAGES = {
    "germline":       ["vas", "aub", "AGO3", "piwi", "zuc", "bgcn", "orb", "osk", "BicD", "egl", "otu", "stau", "grk"],
    "early_germline": ["bam", "Dad", "bgcn", "mei-P26", "pum"],
    "follicle":       ["tj", "Fas3", "eya", "cas", "bab1", "bab2", "upd1", "ImpL2"],
    "late_follicle":  ["Vm26Aa", "Vm26Ab", "Vm26Ac", "Cp36", "Cp38", "Cp15", "Cp16", "Cp19", "br"],
    "muscle":         ["Mhc", "Act57B", "Mef2", "up", "Prm"],
    "trachea":        ["btl", "Gasp"],
    "hemocyte":       ["Hml", "srp", "Pxn"],
}
HEATMAP_GENES = ["vas", "aub", "bam", "Dad", "bgcn", "mei-P26", "orb", "osk", "grk",
                 "tj", "Fas3", "upd1", "eya", "cas", "ptc", "hh", "Vm26Aa", "Cp36", "br",
                 "Mhc", "Act57B", "btl", "Gasp", "Hml", "Pxn", "abd-A", "Abd-B"]
PROBLEM_LABELS = ["artefact", "unannotated", "adult tracheal cell"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts", required=True, help="counts_labeled.h5ad from 01_build_pseudobulk.py")
    p.add_argument("--annotation", required=True,
                   help="FlyBase fbgn_annotation_ID tsv(.gz): col1 symbol, col3 primary FBgn")
    p.add_argument("--outdir", required=True)
    p.add_argument("--labels", nargs="+", default=PROBLEM_LABELS)
    p.add_argument("--min_score", type=float, default=0.1,
                   help="min top lineage score to call a lineage (else 'ambiguous')")
    p.add_argument("--resolution", type=float, default=0.4)
    p.add_argument("--n_markers", type=int, default=12)
    p.add_argument("--debris_max_genes", type=float, default=120,
                   help="subcluster is debris if its median genes/cell is below this ...")
    p.add_argument("--debris_min_mito", type=float, default=4,
                   help="... and its median %% mito is above this")
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def slug(x):
    return "".join(c if c.isalnum() else "_" for c in x).strip("_")


def load_symbols(path):
    opener = gzip.open if path.endswith(".gz") else open
    sym = {}
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                sym[parts[2]] = parts[0]
    return sym


def load_host(path, symbols):
    a = ad.read_h5ad(path)
    a = a[:, a.var_names.str.startswith("FB")].copy()
    a.X = a.X.astype(np.float32)
    a.var["symbol"] = [symbols.get(g, g) for g in a.var_names]
    a.var["mt"] = a.var["symbol"].str.startswith("mt:")
    a.obs["condition"] = (a.obs["genotype"].astype(str).str.replace("-", "") + "_" +
                          a.obs["infection_status"].astype(str))
    a.obs["n_genes"] = np.asarray((a.X > 0).sum(axis=1)).ravel()
    a.obs["pct_mito"] = 100 * np.asarray(a.X[:, a.var["mt"].values].sum(axis=1)).ravel() / \
        np.maximum(np.asarray(a.X.sum(axis=1)).ravel(), 1)
    a.layers["counts"] = a.X.copy()
    sc.pp.normalize_total(a, target_sum=1e4)
    sc.pp.log1p(a)
    # score_genes / heatmaps look genes up by symbol; keep FBgn in var
    a.var["fbgn"] = a.var_names
    a.var_names = a.var["symbol"].astype(str)
    a.var_names_make_unique()
    return a


def score_lineages(a, min_score, seed):
    present = {k: [g for g in v if g in a.var_names] for k, v in LINEAGES.items()}
    missing = sorted({g for v in LINEAGES.values() for g in v} - set(a.var_names))
    if missing:
        print("Markers not in matrix:", ", ".join(missing))
    cols = []
    for k, genes in present.items():
        if len(genes) < 2:
            print(f"  skip lineage {k}: {len(genes)} marker(s) present")
            continue
        sc.tl.score_genes(a, genes, score_name=f"score_{k}", random_state=seed)
        cols.append(f"score_{k}")
    # early_germline is a sub-state of germline: call lineage on the others
    call_cols = [c for c in cols if c != "score_early_germline"]
    s = a.obs[call_cols]
    lin = s.idxmax(axis=1).str.replace("score_", "")
    lin[s.max(axis=1) < min_score] = "ambiguous"
    a.obs["lineage"] = lin.values
    return cols


def marker_heatmap(a, outdir):
    genes = [g for g in HEATMAP_GENES if g in a.var_names]
    grp = a.obs["tissue"].astype(str) + " | " + a.obs["condition"]
    X = a[:, genes].X
    X = X.toarray() if hasattr(X, "toarray") else X
    df = pd.DataFrame(X, columns=genes).groupby(grp.values).mean()
    n = grp.value_counts()
    df = df.loc[[i for i in df.index if n[i] >= 50]].sort_index()
    df.index = [f"{i} (n={n[i]:,})" for i in df.index]
    fig, ax = plt.subplots(figsize=(0.35 * len(genes) + 5, 0.22 * len(df) + 2))
    sns.heatmap(df, cmap="magma", ax=ax, cbar_kws={"label": "mean log1p CP10k"},
                xticklabels=True, yticklabels=True)
    ax.tick_params(labelsize=6)
    ax.set_title("Marker expression by atlas label x condition (groups with >= 50 cells)", fontsize=8)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(outdir, f"marker_heatmap.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)


def subcluster(a, label, args):
    sub = a[a.obs["tissue"] == label].copy()
    if sub.n_obs < 200:
        print(f"  {label}: {sub.n_obs} cells, skipped")
        return None
    sc.pp.highly_variable_genes(sub, n_top_genes=2000, flavor="seurat", batch_key="sample")
    sc.pp.pca(sub, n_comps=30, mask_var="highly_variable", random_state=args.seed)
    sc.pp.neighbors(sub, n_neighbors=15, random_state=args.seed)
    sc.tl.leiden(sub, resolution=args.resolution, key_added="subcluster",
                 flavor="igraph", n_iterations=2, random_state=args.seed)
    sc.tl.umap(sub, random_state=args.seed)
    sc.tl.rank_genes_groups(sub, "subcluster", method="wilcoxon", n_genes=args.n_markers)

    rg = sub.uns["rank_genes_groups"]["names"]
    obs = sub.obs
    rows = []
    for cl in sorted(obs["subcluster"].unique(), key=int):
        o = obs[obs["subcluster"] == cl]
        lin = o["lineage"].value_counts(normalize=True)
        rows.append({
            "subcluster": cl, "n_cells": len(o),
            "median_host_umis": o["host_umis"].median(),
            "median_genes": o["n_genes"].median(),
            "median_pct_mito": round(o["pct_mito"].median(), 2),
            "median_titer_wMel_cells": o.loc[o["infection_status"] == "wMel", "titer"].median(),
            "top_lineage": lin.index[0], "top_lineage_frac": round(lin.iloc[0], 3),
            "lineage_mix": "; ".join(f"{k}={v:.2f}" for k, v in lin.head(3).items()),
            "max_sample_frac": round(o["sample"].value_counts(normalize=True).iloc[0], 3),
            "top_markers": ", ".join(rg[cl][: args.n_markers]),
        })
    tab = pd.DataFrame(rows)
    tab["debris"] = (tab["median_genes"] < args.debris_max_genes) & \
        (tab["median_pct_mito"] > args.debris_min_mito)
    tag = slug(label)
    tab.to_csv(os.path.join(args.outdir, f"subclusters_{tag}.tsv"), sep="\t", index=False)
    pd.crosstab(obs["subcluster"], obs["sample"]).to_csv(
        os.path.join(args.outdir, f"subcluster_samples_{tag}.tsv"), sep="\t")
    print(f"\n== {label}: {sub.n_obs:,} cells, {len(tab)} subclusters")
    print(tab.drop(columns="top_markers").to_string(index=False))

    fig, axes = plt.subplots(1, 4, figsize=(22, 5))
    sc.pl.umap(sub, color="subcluster", ax=axes[0], show=False, legend_loc="on data",
               legend_fontsize=7, title=f"{label}: subcluster")
    sc.pl.umap(sub, color="lineage", ax=axes[1], show=False, title="lineage call")
    sc.pl.umap(sub, color="condition", ax=axes[2], show=False, title="condition")
    sc.pl.umap(sub, color="titer", ax=axes[3], show=False, vmax="p99", title="wMel titer")
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(args.outdir, f"subclusters_{tag}.{ext}"), dpi=150, bbox_inches="tight")
    plt.close(fig)
    debris_cl = set(tab.loc[tab["debris"], "subcluster"].astype(str))
    sub_ids = obs["subcluster"].astype(str)
    return sub_ids, sub_ids.isin(debris_cl)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    sc.settings.verbosity = 1

    a = load_host(args.counts, load_symbols(args.annotation))
    print(f"{a.n_obs:,} cells x {a.n_vars:,} host genes")
    score_cols = score_lineages(a, args.min_score, args.seed)
    marker_heatmap(a, args.outdir)

    frac = pd.crosstab(a.obs["tissue"], a.obs["lineage"], normalize="index").round(3)
    frac.insert(0, "n_cells", a.obs["tissue"].value_counts())
    frac.sort_values("n_cells", ascending=False).to_csv(
        os.path.join(args.outdir, "lineage_by_label.tsv"), sep="\t")
    print("\nLineage call by atlas label (row fractions):")
    print(frac.sort_values("n_cells", ascending=False).to_string())
    a.obs.groupby(["tissue", "sample", "lineage"], observed=True).size().rename("n_cells") \
        .reset_index().to_csv(os.path.join(args.outdir, "lineage_by_label_sample.tsv"),
                              sep="\t", index=False)

    out = a.obs[["sample", "condition", "tissue", "lineage"] + score_cols].copy()
    out["subcluster"] = ""
    out["debris"] = False
    for label in args.labels:
        res = subcluster(a, label, args)
        if res is not None:
            sub, deb = res
            out.loc[sub.index, "subcluster"] = sub.values
            out.loc[deb.index, "debris"] = deb.values
    out.to_csv(os.path.join(args.outdir, "cell_lineage.tsv.gz"), sep="\t")

    deb = out[out["debris"]]
    pd.Series(deb.index, name="obs_name").to_csv(
        os.path.join(args.outdir, "debris_cells.tsv"), sep="\t", index=False)
    print(f"\nDebris cells flagged: {len(deb):,}")
    print(pd.crosstab(deb["tissue"] + ":" + deb["subcluster"], deb["sample"]).to_string())


if __name__ == "__main__":
    main()

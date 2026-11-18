#!/usr/bin/env python
"""Expression of the marker genes per cell type and group, as tables.

Companion to the marker UMAPs (rule rescue_plot_umap_markers): the UMAP shows
where a gene is expressed, this gives the numbers. Genes are MARKERS in
plots/gene_sets.py: mei-P26, Vm26Aa, bam, orb, Sxl and the wMel loci TomO and
CifA.

Reads results/rescue/pseudobulk/counts_labeled.h5ad (debris-free, every gene,
obs host_umis and wmel_umis). For each cell type, genotype by infection group
and gene:

  n_cells          cells of that type in that group
  n_with_value     cells the value is defined for. For host genes, every cell;
                   for wMel genes, cells carrying any wMel read
  frac_expressing  fraction of those with a nonzero count
  cpm              pseudobulk counts per million: summed counts over summed
                   depth across all cells of the group (host_umis for host
                   genes, wmel_umis for wMel genes), times 1e6
  cpm_rep1/2       the same per library, so the replicate spread is visible
  counts           summed counts behind cpm. Read this before any wMel number:
                   TomO has about a hundred counts in the whole dataset, so a
                   large TomO cpm in one cell type can rest on a few reads
  mean_log1p_cp10k mean per-cell log1p(CP10K), the scale of the UMAP colour

wMel genes are normalised to wMel UMIs, so their cpm is per unit of wMel
captured, the same normalisation as the wMel DE. wMel UMIs are mostly 23S and
16S rRNA, so this is in practice per unit of Wolbachia rRNA. In uninfected
groups the wMel numbers come from the few cells with ambient wMel reads and
mean nothing; they are written, flagged by n_with_value, and not interpreted.

These are descriptive. The tests for host genes are in de/host_de_all.tsv.gz
and for wMel genes in wolbachia_followup/wolbachia_de_checked.tsv.gz.

Outputs, in --outdir:
  marker_expression_long.tsv   one row per cell type x group x gene
  marker_expression_cpm.tsv    wide: cell type x gene, one column per group,
                               pseudobulk CPM; the table to read first
  marker_expression_frac.tsv   wide, fraction of cells expressing
  marker_expression_counts.tsv wide, summed counts behind each cpm

Usage:
  python 19_marker_expression.py \
      --counts results/rescue/pseudobulk/counts_labeled.h5ad \
      --outdir results/rescue/marker_expression
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "plots"))
from gene_sets import MARKERS

GROUP_ORDER = ["OreR_uninf", "OreR_wMel", "MeiP26_uninf", "MeiP26_wMel"]
DROP = ("artefact", "unannotated", "adult tracheal cell", "oviduct")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts", required=True, help="counts_labeled.h5ad")
    p.add_argument("--outdir", required=True)
    p.add_argument("--drop_labels", default=",".join(DROP),
                   help="labels left out, comma-separated; '' keeps all")
    p.add_argument("--min_cells", type=int, default=20,
                   help="cell type x group cells needed to report a row")
    return p.parse_args()


def main():
    args = parse_args()
    import anndata as ad
    import scipy.sparse as sp
    os.makedirs(args.outdir, exist_ok=True)
    a = ad.read_h5ad(args.counts, backed="r")
    missing = [k for k, v in MARKERS.items() if v not in a.var_names]
    if missing:
        print(f"not in {os.path.basename(args.counts)}, skipped: {missing}")
    genes = {k: v for k, v in MARKERS.items() if v in a.var_names}
    gi = np.sort(a.var_names.get_indexer(list(genes.values())))
    sub = a[:, gi]                      # one index operation on the backed file
    X = sub.X
    X = X.toarray() if sp.issparse(X) else np.asarray(X)
    col = {g: i for i, g in enumerate(a.var_names[gi])}

    obs = a.obs.copy()
    drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}
    keep = ~obs["tissue"].astype(str).isin(drop).to_numpy()
    obs = obs.loc[keep]
    X = X[keep]
    obs["group"] = (obs["genotype"].astype(str).str.replace("-", "", regex=False)
                    + "_" + obs["infection_status"].astype(str))
    host = pd.to_numeric(obs["host_umis"], errors="coerce").fillna(0).to_numpy()
    wmel = pd.to_numeric(obs["wmel_umis"], errors="coerce").fillna(0).to_numpy()
    print(f"{len(obs):,} cells after dropping {sorted(drop) or 'nothing'}")

    tissues = obs["tissue"].astype(str).to_numpy()
    groups = obs["group"].to_numpy()
    samples = obs["sample"].astype(str).to_numpy()
    rows = []
    for t in list(pd.unique(tissues)) + ["all cells"]:
        mt = np.ones(len(obs), bool) if t == "all cells" else tissues == t
        for g in GROUP_ORDER:
            m = mt & (groups == g)
            if m.sum() < args.min_cells:
                continue
            reps = sorted(pd.unique(samples[m]))
            for lab, gid in genes.items():
                is_w = gid.startswith("WD_")
                depth = wmel if is_w else host
                x = X[:, col[gid]]
                has = m & (depth > 0)
                v = np.log1p(x[has] / depth[has] * 1e4)
                dsum = depth[m].sum()
                rep_cpm = []
                for r in reps[:2]:
                    mr = m & (samples == r)
                    d = depth[mr].sum()
                    rep_cpm.append(x[mr].sum() / d * 1e6 if d > 0 else np.nan)
                rep_cpm += [np.nan] * (2 - len(rep_cpm))
                rows.append({
                    "tissue": t, "group": g, "gene": lab, "gene_id": gid,
                    "source": "wMel" if is_w else "host",
                    "n_cells": int(m.sum()), "n_with_value": int(has.sum()),
                    "frac_expressing": float((x[has] > 0).mean()) if has.any()
                    else np.nan,
                    "cpm": float(x[m].sum() / dsum * 1e6) if dsum > 0 else np.nan,
                    "cpm_rep1": rep_cpm[0], "cpm_rep2": rep_cpm[1],
                    "counts": float(x[m].sum()),
                    "mean_log1p_cp10k": float(v.mean()) if has.any() else np.nan,
                })
    long = pd.DataFrame(rows)
    long.to_csv(os.path.join(args.outdir, "marker_expression_long.tsv"),
                sep="\t", index=False, float_format="%.4g")

    order_t = (long[long["gene"] == "mei-P26"].groupby("tissue")["n_cells"]
               .sum().sort_values(ascending=False).index.tolist()
               if "mei-P26" in genes else list(pd.unique(long["tissue"])))
    order_t = [t for t in order_t if t != "all cells"] + ["all cells"]
    for val, name in (("cpm", "marker_expression_cpm.tsv"),
                      ("frac_expressing", "marker_expression_frac.tsv"),
                      ("counts", "marker_expression_counts.tsv")):
        w = long.pivot_table(index=["tissue", "gene"], columns="group",
                             values=val)
        w = w.reindex(columns=[g for g in GROUP_ORDER if g in w.columns])
        w = w.reindex(pd.MultiIndex.from_product(
            [order_t, list(genes)], names=["tissue", "gene"])).dropna(how="all")
        w.to_csv(os.path.join(args.outdir, name), sep="\t",
                 float_format="%.4g")

    pd.set_option("display.width", 200)
    w = pd.read_csv(os.path.join(args.outdir, "marker_expression_cpm.tsv"),
                    sep="\t")
    w["tissue"] = w["tissue"].str[:40]
    print("\nPseudobulk CPM (host genes per host UMI, wMel genes per wMel UMI):")
    print(w.round(1).to_string(index=False))
    c = long[(long["source"] == "wMel") & long["group"].str.endswith("wMel")]
    thin = c[(c["counts"] < 10) & (c["counts"] > 0)]
    if len(thin):
        print(f"\n{len(thin)} infected wMel rows rest on fewer than 10 counts "
              "(see marker_expression_counts.tsv); do not read their cpm.")
    print("\nwMel values in uninfected groups come from ambient reads in a "
          "handful of cells; see n_with_value in the long table.")


if __name__ == "__main__":
    main()

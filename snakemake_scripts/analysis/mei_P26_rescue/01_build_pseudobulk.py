#!/usr/bin/env python
"""Build a labeled raw-count object and sample x tissue pseudobulk matrices.

Why this step exists: results/integrated/integrated.h5ad carries the atlas
cell-type labels, but its .X is the atlas HVG panel and its .raw is
log1p-normalised, so neither can feed a count-based DE model. This script
re-reads raw counts from each results/filtered_h5ad/<sample>.h5ad (.raw.X),
keeps only cells that survived integration, and attaches the labels by
obs_name ("<sample>__<barcode>", the convention integrate_via_atlas_projection
uses).

Outputs (in --outdir):
  counts_labeled.h5ad              raw counts, all genes (FB* host + GQ* wMel),
                                   one row per labeled cell
  pseudobulk_host_counts.tsv.gz    host genes x pseudobulk (sample|tissue)
  pseudobulk_wolbachia_counts.tsv.gz  wMel genes x pseudobulk
  pseudobulk_meta.tsv              one row per pseudobulk: design factors,
                                   n_cells, host/wMel UMIs, titer summaries

Every sample x tissue combination is written, however few cells it has;
downstream scripts apply --min_cells themselves so the threshold lives in one
place per analysis.
"""
import argparse
import os

import anndata as ad
import numpy as np
import pandas as pd
import scipy.sparse as sp

HOST_PREFIX = "FB"
WMEL_PREFIX = "GQ"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--integrated", required=True,
                   help="results/integrated/integrated.h5ad (labels source)")
    p.add_argument("--filtered", required=True, nargs="+",
                   help="results/filtered_h5ad/<sample>.h5ad files (raw counts in .raw)")
    p.add_argument("--outdir", required=True)
    p.add_argument("--label_col", default="cell_type_annotation_broad",
                   help="obs column in integrated.h5ad used as 'tissue'")
    p.add_argument("--min_confidence", type=float, default=0.5,
                   help="drop cells whose atlas label KNN confidence is below this "
                        "(0 keeps everything)")
    p.add_argument("--drop_qc_outliers", action="store_true",
                   help="also drop cells flagged qc_outlier during integration")
    return p.parse_args()


def load_labels(path, label_col, min_confidence, drop_qc_outliers):
    obs = ad.read_h5ad(path, backed="r").obs.copy()
    needed = [label_col, "source_file", "genotype", "infection_status", "replicate"]
    missing = [c for c in needed if c not in obs.columns]
    if missing:
        raise ValueError(f"{path} obs is missing {missing}; available: {list(obs.columns)}")

    n0 = len(obs)
    conf_col = label_col.replace("cell_type_", "atlas_") + "_confidence"
    if min_confidence > 0:
        if conf_col not in obs.columns:
            raise ValueError(f"--min_confidence set but {conf_col} not in obs")
        obs = obs[obs[conf_col].astype(float) >= min_confidence]
    print(f"Label confidence >= {min_confidence}: kept {len(obs):,}/{n0:,} cells")

    if drop_qc_outliers and "qc_outlier" in obs.columns:
        n1 = len(obs)
        obs = obs[~obs["qc_outlier"].astype(bool)]
        print(f"Dropped qc_outlier cells: kept {len(obs):,}/{n1:,}")

    out = pd.DataFrame({
        "sample": obs["source_file"].astype(str),
        "genotype": obs["genotype"].astype(str),
        "infection_status": obs["infection_status"].astype(str),
        "replicate": obs["replicate"].astype(str),
        "tissue": obs[label_col].astype(str),
        "label_confidence": obs[conf_col].astype(float) if conf_col in obs.columns else np.nan,
    }, index=obs.index)
    return out


def load_raw_counts(paths, keep_obs):
    adatas = []
    for path in paths:
        a = ad.read_h5ad(path)
        if a.raw is None:
            raise ValueError(f"{path} has no .raw (expected raw counts from filter_h5ad)")
        base = os.path.splitext(os.path.basename(path))[0]
        X = a.raw.X
        X = X.tocsr() if sp.issparse(X) else sp.csr_matrix(X)
        b = ad.AnnData(X=X, var=pd.DataFrame(index=a.raw.var_names.astype(str)))
        b.obs_names = [f"{base}__{bc}" for bc in a.obs_names]
        b = b[b.obs_names.isin(keep_obs)].copy()
        print(f"  {base}: {b.n_obs:,} labeled cells")
        adatas.append(b)
    # outer join: samples can differ in which genes survived filtering
    combined = ad.concat(adatas, join="outer", fill_value=0)
    combined.X = combined.X.tocsr()
    return combined


def sum_by_group(X, codes, n_groups):
    """Sum rows of sparse X into n_groups via an indicator matrix."""
    ind = sp.csr_matrix((np.ones(len(codes)), (codes, np.arange(len(codes)))),
                        shape=(n_groups, len(codes)))
    return np.asarray((ind @ X).todense())


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    labels = load_labels(args.integrated, args.label_col,
                         args.min_confidence, args.drop_qc_outliers)
    counts = load_raw_counts(args.filtered, set(labels.index))
    missing = labels.index.difference(counts.obs_names)
    if len(missing):
        print(f"WARNING: {len(missing):,} labeled cells not found in filtered files "
              f"(e.g. {list(missing[:3])}); check obs_name convention")
    counts.obs = labels.loc[counts.obs_names]

    genes = counts.var_names.astype(str)
    host = np.asarray(genes.str.startswith(HOST_PREFIX))
    wmel = np.asarray(genes.str.startswith(WMEL_PREFIX))
    print(f"Genes: {host.sum():,} host ({HOST_PREFIX}*), {wmel.sum():,} wMel ({WMEL_PREFIX}*), "
          f"{(~host & ~wmel).sum():,} other (ignored)")

    host_umis = np.asarray(counts.X[:, host].sum(axis=1)).ravel()
    wmel_umis = np.asarray(counts.X[:, wmel].sum(axis=1)).ravel()
    denom = host_umis + wmel_umis
    counts.obs["host_umis"] = host_umis
    counts.obs["wmel_umis"] = wmel_umis
    counts.obs["titer"] = np.where(denom > 0, wmel_umis / np.where(denom == 0, 1, denom), 0.0)

    # plain-string copy: write_h5ad converts string columns to categoricals in place
    obs = counts.obs.copy()
    counts.write_h5ad(os.path.join(args.outdir, "counts_labeled.h5ad"), compression="gzip")

    # ---- pseudobulk: one column per sample x tissue ----
    pb_id = obs["sample"] + "|" + obs["tissue"]
    pb_cat = pd.Categorical(pb_id)
    codes, levels = pb_cat.codes, list(pb_cat.categories)

    for mask, name in [(host, "host"), (wmel, "wolbachia")]:
        mat = sum_by_group(counts.X[:, mask], codes, len(levels)).T
        df = pd.DataFrame(mat.astype(np.int64), index=genes[mask], columns=levels)
        df.index.name = "gene_id"
        df.to_csv(os.path.join(args.outdir, f"pseudobulk_{name}_counts.tsv.gz"), sep="\t")

    obs = obs.assign(pb_id=pb_id.values)
    meta = obs.groupby("pb_id", observed=True).agg(
        sample=("sample", "first"),
        genotype=("genotype", "first"),
        infection_status=("infection_status", "first"),
        replicate=("replicate", "first"),
        tissue=("tissue", "first"),
        n_cells=("titer", "size"),
        host_umis=("host_umis", "sum"),
        wmel_umis=("wmel_umis", "sum"),
        mean_cell_titer=("titer", "mean"),
        median_cell_titer=("titer", "median"),
    ).reindex(levels)
    meta["pb_titer"] = meta["wmel_umis"] / (meta["wmel_umis"] + meta["host_umis"]).replace(0, np.nan)
    meta.to_csv(os.path.join(args.outdir, "pseudobulk_meta.tsv"), sep="\t")

    print("\nCells per tissue x condition:")
    print(pd.crosstab(obs["tissue"], obs["genotype"] + "_" + obs["infection_status"]).to_string())
    print(f"\nWrote {len(levels)} pseudobulk columns to {args.outdir}")


if __name__ == "__main__":
    main()

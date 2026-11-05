#!/usr/bin/env python
"""Does the host side of DUOS-seq look like published ovary scRNA-seq?

With no matched standard-chemistry run, the host-fidelity argument has to come
from the reference atlas the labels were transferred from. Two comparisons,
and the second is the one that matters.

  proportions   cell-type composition in the uninfected OreR samples against
                the atlas. Wild-type uninfected ovaries should compose roughly
                like a published wild-type ovary; large departures point at
                dissociation bias, not at biology. Mutant and infected samples
                are reported for context but are expected to differ.

  profiles      per-cell-type mean expression in the query against the atlas,
                over shared genes, as a full query x atlas correlation matrix.
                The fidelity test is whether each query cell type correlates
                best with its own atlas cell type, i.e. whether the matrix
                diagonal is the row maximum. A high mean correlation with an
                off-diagonal argmax means labels were transferred onto the
                wrong neighbours, which no proportion table would reveal.

Because labels came from this atlas by projection, the comparison is not
independent: the test is whether transcriptomes agree, not whether labels were
assigned correctly from scratch. State that in the methods. The diagonal check
is still informative, because projection constrains neighbourhoods in PCA
space and not per-gene expression.

Usage:
  python 15_atlas_fidelity.py \
      --query results/integrated/integrated_clean.h5ad \
      --atlas ovary_atlas/atlas/h5ad_unprocessed/s_fca_biohub_ovary_10x.h5ad \
      --flybase_annotation reference/fbgn_annotation_ID_fb_2025_04.tsv.gz \
      --outdir results/rescue/atlas_fidelity
"""
import argparse
import os
import sys

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import seaborn as sns
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
try:
    from integrate_via_atlas_projection import harmonise_atlas_gene_ids
except ImportError as e:  # pragma: no cover
    raise SystemExit("cannot import harmonise_atlas_gene_ids from "
                     f"integrate_via_atlas_projection.py ({e}); this script expects "
                     "to live in snakemake_scripts/analysis/mei_P26_rescue/")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--query", required=True, help="integrated_clean.h5ad")
    p.add_argument("--atlas", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--flybase_annotation", default=None,
                   help="fbgn_annotation_ID tsv.gz, for symbol->FBgn on the atlas")
    p.add_argument("--label_col", default="cell_type_annotation")
    p.add_argument("--atlas_label_col", default="annotation")
    p.add_argument("--reference_group", default="OreR_uninf",
                   help="query condition whose composition is compared to the atlas")
    p.add_argument("--min_cells", type=int, default=20,
                   help="minimum cells for a cell type to get an expression profile")
    p.add_argument("--n_genes", type=int, default=3000,
                   help="shared highly variable genes used for the profile correlation")
    p.add_argument("--max_atlas_cells", type=int, default=0,
                   help="subsample the atlas to this many cells (0 = all)")
    p.add_argument("--min_shared_genes", type=int, default=500,
                   help="fail if fewer genes are shared after harmonisation")
    return p.parse_args()


def normalise(a):
    a = a.copy()
    if a.raw is not None and a.raw.shape[1] == a.shape[1]:
        a.X = a.raw.X.copy()
    mx = float(a.X[:50].max()) if a.n_obs else 0.0
    if mx > 50:  # looks like counts, not already logged
        sc.pp.normalize_total(a, target_sum=1e4)
        sc.pp.log1p(a)
    return a


def profiles(a, label_col, genes, min_cells):
    """mean expression per cell type, over `genes`."""
    a = a[:, genes]
    rows, keep = {}, []
    lab = a.obs[label_col].astype(str)
    for ct, idx in lab.groupby(lab).groups.items():
        if len(idx) < min_cells:
            continue
        X = a[list(idx)].X
        rows[ct] = np.asarray(X.mean(axis=0)).ravel()
        keep.append((ct, len(idx)))
    return pd.DataFrame(rows, index=genes).T, dict(keep)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    print(f"Reading {args.query}")
    q = ad.read_h5ad(args.query)
    print(f"  {q.n_obs:,} cells x {q.n_vars:,} genes")
    print(f"Reading {args.atlas}")
    at = ad.read_h5ad(args.atlas)
    at = harmonise_atlas_gene_ids(at, args.flybase_annotation)
    at.var_names_make_unique()
    if args.max_atlas_cells and at.n_obs > args.max_atlas_cells:
        sc.pp.subsample(at, n_obs=args.max_atlas_cells, random_state=0)
        print(f"  atlas subsampled to {at.n_obs:,} cells")
    if args.atlas_label_col not in at.obs.columns:
        raise SystemExit(f"--atlas_label_col {args.atlas_label_col!r} not in the atlas; "
                         f"present: {sorted(at.obs.columns)[:40]}")

    # ---- composition ----
    cond = "condition" if "condition" in q.obs.columns else None
    ql = q.obs[args.label_col].astype(str)
    comp = pd.DataFrame({"query_all": ql.value_counts(normalize=True)})
    if cond and (q.obs[cond].astype(str) == args.reference_group).any():
        sub = ql[q.obs[cond].astype(str) == args.reference_group]
        comp[f"query_{args.reference_group}"] = sub.value_counts(normalize=True)
        print(f"\n{len(sub):,} cells in the reference group {args.reference_group!r}")
    else:
        print(f"\nreference group {args.reference_group!r} not found in obs[{cond!r}]; "
              "comparing all query cells")
    comp["atlas"] = at.obs[args.atlas_label_col].astype(str).value_counts(normalize=True)
    comp = comp.fillna(0.0)
    comp["log2_ratio_vs_atlas"] = np.log2((comp.iloc[:, -2] + 1e-4)
                                          / (comp["atlas"] + 1e-4))
    comp.index.name = "cell_type"
    comp.to_csv(os.path.join(args.outdir, "atlas_proportions.tsv"), sep="\t")

    shared_ct = comp.index[(comp["atlas"] > 0) & (comp["query_all"] > 0)]
    rho = None
    if len(shared_ct) > 2:
        col = comp.columns[-2]
        if comp.loc[shared_ct, col].nunique() > 1 and comp.loc[shared_ct, "atlas"].nunique() > 1:
            rho = stats.spearmanr(comp.loc[shared_ct, col], comp.loc[shared_ct, "atlas"])
            print(f"Composition agreement over {len(shared_ct)} shared cell types: "
                  f"Spearman rho = {rho.statistic:.3f} (p = {rho.pvalue:.3g})")
    print("\nLargest composition departures from the atlas:")
    print(comp.reindex(comp["log2_ratio_vs_atlas"].abs().sort_values(ascending=False).index)
          .head(8).round(4).to_string())

    # ---- expression profiles ----
    shared_genes = q.var_names.intersection(at.var_names)
    print(f"\n{len(shared_genes):,} shared genes after harmonisation")
    if len(shared_genes) < args.min_shared_genes:
        raise SystemExit("too few shared genes; check --flybase_annotation and that the "
                         "atlas var_names harmonised to FBgn IDs")
    qn, an = normalise(q[:, shared_genes]), normalise(at[:, shared_genes])
    sc.pp.highly_variable_genes(an, n_top_genes=min(args.n_genes, len(shared_genes) - 1))
    genes = an.var_names[an.var["highly_variable"]].tolist()
    print(f"using {len(genes)} atlas-variable shared genes for the profiles")

    qp, qn_cells = profiles(qn, args.label_col, genes, args.min_cells)
    ap, an_cells = profiles(an, args.atlas_label_col, genes, args.min_cells)
    print(f"{len(qp)} query cell types and {len(ap)} atlas cell types with "
          f">= {args.min_cells} cells")

    C = pd.DataFrame(
        np.corrcoef(np.vstack([qp.to_numpy(), ap.to_numpy()]))[:len(qp), len(qp):],
        index=qp.index, columns=ap.index)
    C.to_csv(os.path.join(args.outdir, "atlas_profile_correlation_matrix.tsv"), sep="\t")

    rows = []
    for ct in C.index:
        best = C.loc[ct].idxmax()
        rows.append({"cell_type": ct, "n_cells_query": qn_cells.get(ct),
                     "self_correlation": float(C.loc[ct, ct]) if ct in C.columns else np.nan,
                     "best_match": best, "best_correlation": float(C.loc[ct, best]),
                     "diagonal_is_best": ct == best})
    per_ct = pd.DataFrame(rows).sort_values("self_correlation", ascending=False)
    per_ct.to_csv(os.path.join(args.outdir, "atlas_profile_correlations.tsv"),
                  sep="\t", index=False)

    matched = per_ct[per_ct["self_correlation"].notna()]
    n_diag = int(matched["diagonal_is_best"].sum())
    print(f"\nSelf-correlation available for {len(matched)} cell types; "
          f"{n_diag} of {len(matched)} match their own atlas type best "
          f"(median self r = {matched['self_correlation'].median():.3f})")
    print(matched.round(3).to_string(index=False))
    off = matched[~matched["diagonal_is_best"]]
    if len(off):
        print("\nCell types whose best match is another atlas type (check these labels):")
        print(off[["cell_type", "self_correlation", "best_match",
                   "best_correlation"]].round(3).to_string(index=False))

    with open(os.path.join(args.outdir, "atlas_fidelity_summary.txt"), "w") as fh:
        fh.write(f"query\t{args.query}\natlas\t{args.atlas}\n")
        fh.write(f"shared_genes\t{len(shared_genes)}\nprofile_genes\t{len(genes)}\n")
        fh.write(f"cell_types_compared\t{len(matched)}\n")
        fh.write(f"diagonal_is_best\t{n_diag}\n")
        fh.write(f"median_self_correlation\t{matched['self_correlation'].median():.4f}\n")
        if rho is not None:
            fh.write(f"composition_spearman\t{rho.statistic:.4f}\n")
            fh.write(f"composition_spearman_p\t{rho.pvalue:.4g}\n")

    # ---- figures ----
    col = comp.columns[-2]
    fig, axes = plt.subplots(1, 2, figsize=(13, 5.5))
    ax = axes[0]
    ax.scatter(comp.loc[shared_ct, "atlas"], comp.loc[shared_ct, col], s=30, c="#4C72B0")
    lim = [1e-4, max(comp[["atlas", col]].to_numpy().max() * 1.6, 1e-3)]
    ax.plot(lim, lim, ls="--", c="grey", lw=1)
    for ct in shared_ct:
        ax.annotate(str(ct)[:22], (comp.loc[ct, "atlas"], comp.loc[ct, col]),
                    fontsize=6, alpha=0.8)
    ax.set(xscale="log", yscale="log", xlim=lim, ylim=lim,
           xlabel="atlas proportion", ylabel=f"{col} proportion",
           title="Cell-type composition against the atlas")

    ax = axes[1]
    order = [c for c in C.index if c in C.columns]
    sns.heatmap(C.loc[order, order], cmap="viridis", ax=ax,
                cbar_kws={"label": "Pearson r of mean expression"})
    ax.set(title="Query (rows) against atlas (columns)", xlabel="", ylabel="")
    ax.tick_params(labelsize=6)
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "atlas_fidelity.pdf"))
    fig.savefig(os.path.join(args.outdir, "atlas_fidelity.png"), dpi=200)
    plt.close(fig)
    print(f"\nWrote {args.outdir}")


if __name__ == "__main__":
    main()

"""
integrate_via_atlas_projection.py
==================================
Replaces integrate.py's joint BBKNN re-clustering with a single projection
of every filtered sample onto a FROZEN, pre-annotated reference atlas --
same mechanism as Lum et al. 2027's integrate_via_atlas_projection.py
(which itself replaced a Harmony/BBKNN step with scanpy.tl.ingest reference
projection against the Flysta3D-v2 atlas). Here the reference is the Fly
Cell Atlas (FCA) BioHub ovary object instead, since this project's cells
are ovary/ovary-cell-culture, not whole embryo.

Why project instead of jointly re-cluster (ported rationale from Lum's
script -- same mechanism, same reasoning)
-------------------------------------------------------------------------
Reference projection means every query cell's coordinates come from ONE
fixed, external transform (the atlas's own PCA/UMAP) -- no query file ever
sees or influences another query file's placement, so there is no "batch"
axis between your own samples (genotype x infection status x replicate)
for anything to leak into. Every sample lands in the same coordinate
system because they're all measured against the same fixed yardstick, not
because they were corrected to agree with each other. This also makes
adding a new sample later an O(1) operation: project it onto the same
frozen reference, no need to recompute anything for samples already
projected.

What this object is NOT for: the atlas was fit on its OWN biology (adult
ovary cell types from the Fly Cell Atlas), not on your experimental axis
(Wolbachia titer, infection status, genotype). It has no reason to be
sensitive to the subtle within-cell-type expression shifts a titer/
infection-status analysis cares about -- that signal was never part of
what the atlas's PCs were fit to capture. Use THIS object (atlas_<label> +
atlas UMAP) to answer "what ovary cell type is this cell, and how does
cell-type composition / titer differ by genotype and infection status."
For "does expression change with titer within a cell type," a Harmony
embedding fit on your own cells alone (the old integrate.py path) is still
the right tool -- two different questions, two different embeddings.

Differences from Lum et al. 2027's version
-------------------------------------------------------------------------
  - No Dsim->Dmel ortholog remapping: every sample here (OreR, Mei-P26
    genotypes) is D. melanogaster, so that step (and --ortholog_map) is
    dropped entirely rather than carried over unused.
  - No embryo/primary_cells/cell_culture 3-way split: this project's
    experimental axis is genotype x infection status
    (OreR/Mei-P26 x wMel/uninf), parsed from each sample's condition
    string instead (see add_sample_metadata below). No is_embryo/
    sample_type columns are written.
  - The query files are dual host+Wolbachia transcriptomes (FBgn-prefixed
    Dmel genes + GQ-prefixed wMel genes in the same object, see
    kallisto_bustools_qcfilter_adata_no_pybiomart.py). The ovary atlas is
    host-only, so this needs no special handling: project_query_onto_
    reference restricts to genes shared with the atlas's HVG panel, which
    naturally excludes GQ* Wolbachia genes (the atlas has none). Wolbachia
    genes -- and the per-cell wolbachia_titer column already computed
    upstream -- survive untouched in the full-gene .raw snapshot and in
    .obs, same as every other obs column load_all_query_files keeps.
  - plot_diagnostics adds a titer-by-atlas-cell-type panel (boxplot of
    wolbachia_titer per atlas_<label>, split by infection status) in
    addition to Lum's UMAP/composition plots, since titer-vs-cell-type is
    this project's actual question for this object.

Everything else -- frozen-atlas embedding fit (build_reference_embedding),
ingest projection + KNN confidence (project_query_onto_reference /
_add_knn_confidence), and atlas gene-ID harmonisation (symbol->FBgn via
--flybase_annotation) -- is the same mechanism as Lum's script, ported
directly rather than re-derived.

Run with:
    mamba activate scanpy
    python snakemake_scripts/analysis/integrate_via_atlas_projection.py \\
        --atlas ovary_atlas/atlas/h5ad_unprocessed/s_fca_biohub_ovary_10x.h5ad \\
        --query results/filtered_h5ad/*.h5ad \\
        --out_path results/integrated/integrated.h5ad \\
        --fig_dir results/integrated/figures \\
        --flybase_annotation reference/fbgn_annotation_ID_fb_2025_04.tsv.gz

First run without --label_cols to print the atlas's candidate obs columns
and exit before doing any expensive work -- same discovery behaviour as
Lum's script, needed here since this atlas's obs schema hasn't been
inspected yet.
"""

import os
import re
import glob
import argparse

import numpy as np
import pandas as pd
import matplotlib
matplotlib.use("Agg")  # headless -- no display on SLURM compute nodes
import matplotlib.pyplot as plt
import scipy.sparse
import anndata as ad
import scanpy as sc


FBGN_RE = re.compile(r"^FBgn\d{7,8}$")


# -----------------------------------------------------------------------------
# Small helpers
# -----------------------------------------------------------------------------

def _savefig(fig, path):
    fig.savefig(path, bbox_inches="tight", dpi=150)
    plt.close(fig)
    print(f"   Saved: {path}")


# -----------------------------------------------------------------------------
# Gene ID harmonisation: the FCA ovary atlas may be symbol-indexed; your own
# kallisto|bustools output is FBgn-indexed (dual host+Wolbachia object, see
# kallisto_bustools_qcfilter_adata_no_pybiomart.py). Ported from Lum et al.
# 2027's annotate_with_flysta3d.py unchanged -- same file format, same fix.
# -----------------------------------------------------------------------------

def load_symbol_to_fbgn(flybase_path):
    """Build a gene_symbol -> primary FBgn# map from a FlyBase
    fbgn_annotation_ID_fb_*.tsv.gz file (column 0 = gene_symbol, column 2 =
    primary_FBgn#). Ambiguous symbols (mapping to >1 distinct FBgn) are
    dropped rather than picking one arbitrarily.
    """
    import gzip
    from io import StringIO

    with gzip.open(flybase_path, "rt") as f:
        lines = [line for line in f if not line.startswith("#")]
    df = pd.read_csv(StringIO("".join(lines)), sep="\t", header=None)
    df = df[[0, 2]].rename(columns={0: "symbol", 2: "fbgn"}).dropna()

    n_raw = len(df)
    df = df.drop_duplicates(subset="symbol", keep=False)
    n_kept = len(df)
    if n_kept < n_raw:
        print(f"  FlyBase symbol->FBgn map: dropped {n_raw - n_kept}/{n_raw} "
              f"ambiguous (non-unique) symbol rows from {flybase_path}")

    mapping = dict(zip(df["symbol"], df["fbgn"]))
    print(f"  Loaded {len(mapping)} unique symbol->FBgn mappings")
    return mapping


def _recover_var_index(var, label="var"):
    """Some h5ad exports end up with a purely positional integer var index
    (0, 1, 2, ...) while the real gene identifiers survive only as an
    ordinary column. Detects a numeric-looking index and, if found,
    promotes whichever column looks like real (non-numeric) gene IDs back
    to the index."""
    idx_is_numeric = pd.Series(var.index.astype(str)).str.match(r"^\d+$").mean() > 0.9
    if not idx_is_numeric:
        return var

    for col in var.columns:
        vals = var[col].astype(str)
        if vals.str.match(r"^\d+$").mean() < 0.5:
            print(f"   {label}: index is positional (0,1,2,...) -- "
                  f"promoting column '{col}' (e.g. {vals.iloc[0]!r}) to be "
                  "the real var_names")
            var = var.set_index(col, drop=True)
            var.index.name = None
            return var

    print(f"   WARNING: {label} index looks positional (0,1,2,...) and no "
          f"column looks like real gene IDs either (columns: "
          f"{list(var.columns)}) -- leaving var_names as row numbers, gene "
          "ID matching below will fail")
    return var


def harmonise_atlas_gene_ids(atlas, flybase_annotation=None):
    """Make atlas.var_names FBgn IDs so they line up with your query files.

    Priority:
      1. var_names already look like FBgn IDs -> leave as-is.
      2. atlas.var has a column that looks like it holds FBgn/FlyBase IDs
         -> use that column directly.
      3. Fall back to symbol->FBgn remap via --flybase_annotation.
    """
    frac_fbgn = atlas.var_names.str.match(FBGN_RE).mean()
    print(f"  Atlas var_names matching FBgn pattern: {frac_fbgn*100:.1f}%")
    print(f"  Atlas var_names sample (first 20): {atlas.var_names[:20].tolist()}")
    print(f"  Atlas var.columns: {atlas.var.columns.tolist()}")
    if frac_fbgn > 0.5:
        print("  Atlas already FBgn-indexed -- no remap needed")
        return atlas

    for col in atlas.var.columns:
        if re.search(r"fbgn|flybase|gene_?id", col, re.IGNORECASE):
            vals = atlas.var[col].astype(str)
            frac = vals.str.match(FBGN_RE).mean()
            if frac > 0.5:
                print(f"  Using atlas.var['{col}'] as FBgn ID "
                      f"({frac*100:.1f}% match FBgn pattern)")
                atlas = atlas.copy()
                atlas.var_names = vals.values
                atlas.var_names_make_unique()
                return atlas

    if flybase_annotation is None:
        raise ValueError(
            "Atlas var_names don't look like FBgn IDs and no --flybase_annotation "
            "was given to remap them. Either pass --flybase_annotation "
            "reference/fbgn_annotation_ID_fb_2025_04.tsv.gz, or inspect "
            "atlas.var.columns yourself and tell the script which column holds "
            "FlyBase IDs."
        )

    print("  Atlas appears symbol-indexed -- remapping to FBgn via "
          f"{flybase_annotation}")
    symbol_to_fbgn = load_symbol_to_fbgn(flybase_annotation)
    mapped = atlas.var_names.map(symbol_to_fbgn)
    keep = mapped.notna()
    n_total, n_kept = atlas.n_vars, int(keep.sum())
    print(f"  Symbol->FBgn remap: {n_kept}/{n_total} atlas genes matched "
          f"({n_total - n_kept} dropped -- no unique FBgn for that symbol)")
    # BUG FIX vs. the Lum et al. 2027 version this was ported from: that
    # script used an absolute floor (n_kept < 1000) to decide the symbol
    # match "went wrong". That's only a safe proxy when the reference atlas
    # is known to have >>1000 genes (true for the whole-embryo Flysta3D-v2
    # atlas it was written against); a smaller reference panel can match
    # 100% of its genes and still trip an absolute floor. Use a fraction of
    # n_total instead, which is what the check is actually trying to catch.
    if n_total > 0 and (n_kept / n_total) < 0.1:
        sample_symbols = list(symbol_to_fbgn.keys())[:20]
        print(f"  Atlas var_names sample (unmatched): {atlas.var_names[:20].tolist()}")
        print(f"  FlyBase table symbol sample (what we're matching against): "
              f"{sample_symbols}")
        raise ValueError(
            f"Only {n_kept}/{n_total} ({100*n_kept/n_total:.1f}%) atlas genes "
            "remapped to FBgn -- something is wrong with the symbol matching "
            "(check for a species prefix or case mismatches in "
            "atlas.var_names, or whether atlas.var_names are actually "
            "symbols at all -- compare the two samples printed above)."
        )
    atlas = atlas[:, keep].copy()
    atlas.var_names = mapped[keep].astype(str).values
    atlas.var_names_make_unique()
    return atlas


# -----------------------------------------------------------------------------
# Step 1 -- Load atlas reference
# -----------------------------------------------------------------------------

def load_atlas_reference(atlas_path, label_cols, flybase_annotation=None,
                          subsample_ref=None, random_state=42):
    print(f"\n-- Loading reference atlas: {atlas_path} --")
    atlas = sc.read_h5ad(atlas_path)
    print(f"   {atlas.n_obs:,} cells x {atlas.n_vars:,} genes")

    if not label_cols:
        candidates = [c for c in atlas.obs.columns
                      if atlas.obs[c].dtype == object
                      or str(atlas.obs[c].dtype).startswith("category")]
        print("\n   No --label_cols given. Candidate annotation columns in "
              "the atlas obs (dtype object/category):")
        for c in candidates:
            n_uniq = atlas.obs[c].nunique()
            print(f"     {c!r}  ({n_uniq} unique values)")
        print("\n   All obs columns:")
        print(f"     {list(atlas.obs.columns)}")
        raise SystemExit(
            "\nPick one or more of the columns above and re-run with "
            "--label_cols <col1> [<col2> ...]."
        )

    missing = [c for c in label_cols if c not in atlas.obs.columns]
    if missing:
        raise ValueError(
            f"--label_cols {missing} not found in atlas obs. "
            f"Available: {list(atlas.obs.columns)}"
        )

    for c in label_cols:
        print(f"\n   Atlas '{c}' distribution (top 15):")
        print(atlas.obs[c].astype(str).value_counts().head(15).to_string())

    if subsample_ref and atlas.n_obs > subsample_ref:
        print(f"\n   Subsampling atlas to {subsample_ref:,} cells "
              f"(--subsample_ref) for a faster first pass")
        sc.pp.subsample(atlas, n_obs=subsample_ref, random_state=random_state)

    # Prefer raw counts for renormalisation consistency with the query data.
    if atlas.raw is not None:
        X = atlas.raw.X
        var = atlas.raw.var.copy()
        print("   Using atlas.raw.X as counts source")
    elif "counts" in atlas.layers:
        X = atlas.layers["counts"]
        var = atlas.var.copy()
        print("   Using atlas.layers['counts'] as counts source")
    else:
        X = atlas.X
        var = atlas.var.copy()
        print("   WARNING: no .raw or 'counts' layer found in the atlas -- "
              "using .X as-is. If this atlas object is already "
              "normalised/log1p'd, renormalising it here is not strictly "
              "correct (values won't be integer counts), but the relative "
              "structure used for PCA/neighbors/UMAP is only mildly "
              "affected. Check atlas.X.max() if this matters for your "
              "analysis.")

    # A few public atlas exports carry a purely positional var index with
    # the real gene IDs stranded in an ordinary column instead -- fix that
    # up before var_names is used for anything downstream.
    var = _recover_var_index(var, label="atlas var")

    if scipy.sparse.issparse(X):
        X = X.tocsr()
    else:
        X = scipy.sparse.csr_matrix(X)
    X.data = X.data.astype(np.float32)

    # JUDGMENT CALL: drop any all-zero cells (no RNA signal to cluster on --
    # left in, sc.pp.scale's per-gene z-scoring would pull these into a
    # tight, artificial cluster with no biological meaning). This atlas
    # isn't known to co-embed a second assay the way Flysta3D-v2 does, but
    # the check itself is generic and cheap, so it stays as a safety net.
    total_counts = np.asarray(X.sum(axis=1)).ravel()
    keep_cell = total_counts > 0
    n_dropped = int((~keep_cell).sum())
    if n_dropped:
        print(f"   Dropping {n_dropped:,} / {X.shape[0]:,} atlas cells with "
              f"zero total counts (no RNA signal)")
        X = X[keep_cell]
        atlas = atlas[keep_cell].copy()
    else:
        print("   No zero-count atlas cells found (nothing dropped)")

    ref = ad.AnnData(X=X, obs=atlas.obs[label_cols].copy(), var=var)
    ref.obs_names = atlas.obs_names
    ref.obs["dataset"] = "ovary_atlas"
    ref = harmonise_atlas_gene_ids(ref, flybase_annotation)

    print(f"\n   Reference ready: {ref.n_obs:,} cells x {ref.n_vars:,} genes")
    return ref


# -----------------------------------------------------------------------------
# Step 2 -- Fit the reference space ONCE, on the atlas alone. Frozen from
# here on: nothing downstream is allowed to move these coordinates. Ported
# unchanged from Lum et al. 2027's annotate_with_flysta3d_ingest.py.
# -----------------------------------------------------------------------------

def build_reference_embedding(ref, n_pcs=30, n_top_genes=3000, n_neighbors=30):
    print("\n-- Building frozen reference embedding (atlas only, no query) --")
    ref = ref.copy()

    if scipy.sparse.issparse(ref.X):
        ref.X = ref.X.tocsr()
    else:
        ref.X = scipy.sparse.csr_matrix(ref.X)

    print(f"   Highly variable genes (seurat_v3, atlas-only, top {n_top_genes}) ...")
    sc.pp.highly_variable_genes(ref, flavor="seurat_v3", n_top_genes=n_top_genes)

    print("   Normalising (1e4 per cell) + log1p ...")
    sc.pp.normalize_total(ref, target_sum=1e4)
    sc.pp.log1p(ref)

    ref = ref[:, ref.var["highly_variable"]].copy()
    print(f"   Reference HVG panel: {ref.n_vars:,} genes")

    # Deliberately NOT scaling (no sc.pp.scale) here -- ingest re-applies
    # the reference's own PCA loadings directly to the query's log-
    # normalised values; scaling the reference with its own per-gene mean/
    # std and then not re-applying those exact stats to the query is a
    # known footgun (mirrors scanpy's own "Integrating data using ingest"
    # tutorial recipe).
    print(f"   PCA ({n_pcs} components) ...")
    sc.pp.pca(ref, n_comps=n_pcs, svd_solver="arpack")

    print(f"   Neighbors (k={n_neighbors}) + UMAP ...")
    sc.pp.neighbors(ref, n_neighbors=n_neighbors, n_pcs=n_pcs)
    sc.tl.umap(ref)

    print(f"   Reference embedding fit: {ref.n_obs:,} cells x {ref.n_vars:,} HVGs, "
          f"{n_pcs} PCs -- this UMAP is now frozen.")
    return ref


# -----------------------------------------------------------------------------
# Step 3 -- Project query cells into the frozen reference via scanpy.tl.
# ingest, then a KNN confidence score in that same space. Ported unchanged
# from Lum et al. 2027's annotate_with_flysta3d_ingest.py. Note this is
# also what naturally keeps GQ*-prefixed Wolbachia genes out of the
# projection: `common` below is genes shared with ref_genes (host-only
# atlas HVGs), so Wolbachia genes are simply never part of `q`.
# -----------------------------------------------------------------------------

def project_query_onto_reference(ref, query, label_cols, k=30):
    """Project ALL query cells into ref's frozen PCA/UMAP space.

    Nothing about `ref` is modified by this call. Returns a NEW AnnData
    (the query cells) carrying obsm['X_pca']/obsm['X_umap'] in the
    reference's own coordinate system, plus atlas_<col>/
    atlas_<col>_confidence obs columns for every label_col.
    """
    ref_genes = ref.var_names
    n_present = int(query.var_names.isin(ref_genes).sum())
    print(f"\n-- Projecting {query.n_obs:,} query cells onto the frozen "
          f"reference embedding --")
    print(f"   Query genes overlapping reference HVG panel: "
          f"{n_present}/{len(ref_genes)}")
    if n_present < 0.5 * len(ref_genes):
        print("   WARNING: fewer than half the reference HVGs are present in "
              "the query -- projection quality may be poor. Check gene ID "
              "harmonisation (FBgn vs. symbol) above.")

    common = ref_genes.intersection(query.var_names)
    q = query[:, common].copy()
    missing = ref_genes.difference(common)
    if len(missing):
        print(f"   Zero-filling {len(missing)} reference HVGs absent from "
              "this query set")
        pad = ad.AnnData(
            X=scipy.sparse.csr_matrix((q.n_obs, len(missing)), dtype=np.float32),
            obs=pd.DataFrame(index=q.obs_names.copy()),
            var=pd.DataFrame(index=missing),
        )
        q = ad.concat([q, pad], axis=1, join="outer")
    q = q[:, ref_genes].copy()
    q.obs = query.obs.copy()

    if scipy.sparse.issparse(q.X):
        q.X = q.X.tocsr()
    else:
        q.X = scipy.sparse.csr_matrix(q.X)

    print("   Normalising (1e4 per cell) + log1p (same recipe as reference) ...")
    sc.pp.normalize_total(q, target_sum=1e4)
    sc.pp.log1p(q)

    print(f"   Running sc.tl.ingest (label_cols={label_cols}) ...")
    sc.tl.ingest(q, ref, obs=label_cols, embedding_method=["umap", "pca"])

    # ingest overwrites obs[col] in place with the transferred value; rename
    # to atlas_<col> so it's unambiguous once this sits next to other obs
    # columns on the same object.
    q.obs = q.obs.rename(columns={col: f"atlas_{col}" for col in label_cols})

    _add_knn_confidence(q, ref, label_cols, k=k)

    return q


def _add_knn_confidence(query, ref, label_cols, k=30):
    """ingest only writes a hard label transfer -- compute a majority-vote
    confidence score (winning fraction among k nearest reference neighbours
    in the same frozen PCA space ingest used)."""
    from sklearn.neighbors import NearestNeighbors

    print(f"\n-- KNN confidence scoring (k={k}, in reference PCA space) --")
    nbrs = NearestNeighbors(n_neighbors=k, metric="euclidean", n_jobs=-1)
    nbrs.fit(ref.obsm["X_pca"])
    _, indices = nbrs.kneighbors(query.obsm["X_pca"])

    for col in label_cols:
        labels_arr = ref.obs[col].astype(str).values
        confidence = []
        for row_idx in indices:
            neigh = labels_arr[row_idx]
            counts = pd.Series(neigh).value_counts()
            confidence.append(counts.iloc[0] / k)
        query.obs[f"atlas_{col}_confidence"] = confidence

        print(f"\n   Transferred 'atlas_{col}' distribution:")
        print(query.obs[f"atlas_{col}"].value_counts().head(20).to_string())
        mean_conf = float(np.mean(confidence))
        low_conf = float(np.mean(np.array(confidence) < 0.5))
        print(f"   Mean confidence: {mean_conf:.3f}")
        if low_conf > 0.2:
            print(f"   WARNING: {low_conf*100:.1f}% of cells have confidence "
                  f"< 0.5 for '{col}' -- treat this column cautiously, or "
                  "increase --k, or check gene-ID harmonisation upstream")


# -----------------------------------------------------------------------------
# Step 4 -- Load all query files with EVERY obs column kept (wolbachia_titer,
# condition, method, replicate, ...), not just a cherry-picked subset --
# this object is meant to be the final analysis-ready integrated object.
# Ported from Lum et al. 2027's integrate_via_atlas_projection.py, minus
# the Dsim->Dmel remap branch (not needed -- every sample here is Dmel).
# -----------------------------------------------------------------------------

def load_all_query_files(query_paths):
    adatas = []
    for path in query_paths:
        print(f"\n   Loading query: {path}")
        adata = sc.read_h5ad(path)
        print(f"   {adata.n_obs} cells x {adata.n_vars} genes")

        if adata.raw is None:
            raise ValueError(
                f"{path} has no .raw -- expected filtered h5ad output with "
                "adata.raw set to pre-normalisation counts."
            )

        raw_X = adata.raw.X
        raw_X = raw_X.tocsr() if scipy.sparse.issparse(raw_X) else scipy.sparse.csr_matrix(raw_X)
        raw_X.data = raw_X.data.astype(np.float32)

        basename = os.path.splitext(os.path.basename(path))[0]
        a = ad.AnnData(X=raw_X, obs=adata.obs.copy(), var=adata.raw.var.copy())
        a.obs["dataset"]     = "query"
        a.obs["source_file"] = basename

        a.obs_names = [f"{basename}__{bc}" for bc in a.obs_names]
        adatas.append(a)

    print(f"\n-- Concatenating {len(adatas)} query files --")
    query = ad.concat(adatas, join="outer", index_unique=None)
    query.obs_names_make_unique()
    if scipy.sparse.issparse(query.X):
        query.X = query.X.tocsr()
    print(f"   Total query cells: {query.n_obs:,}")
    return query


# -----------------------------------------------------------------------------
# Step 5 -- condition/replicate/method/genotype/infection_status, parsed
# from source_file (the per-sample h5ad basename, i.e. samples.csv's
# "{condition}-{replicate}_{seq_platform}" sample_id -- see Snakefile).
# condition/replicate/method regex ported unchanged from Lum et al. 2027's
# add_sample_metadata (same "<condition>-<replicate>_<10x|pipseq>" sample_id
# convention). genotype/infection_status is new here, splitting condition
# on its trailing "_wMel"/"_uninf" suffix -- this project's experimental
# axis (OreR/Mei-P26 genotype x infection status), in place of Lum's embryo
# vs. primary_cells vs. cell_culture sample_type split, which doesn't apply
# to this dataset.
# -----------------------------------------------------------------------------

_INFECTION_SUFFIXES = ("wMel", "uninf")


def add_sample_metadata(adata, batch_key="source_file"):
    parsed = adata.obs[batch_key].astype(str).str.extract(
        r"^(?P<condition>.+)-(?P<replicate>\d+)_(?P<method>10x|pipseq)$"
    )
    unparsed = parsed["condition"].isna()
    if unparsed.any():
        bad = sorted(adata.obs.loc[unparsed, batch_key].unique())
        print(f"  WARNING: {int(unparsed.sum())} cells have a {batch_key!r} "
              "value that doesn't match '<condition>-<replicate>_<10x|pipseq>': "
              f"{bad}. condition/replicate/method left as 'unknown' for these.")

    adata.obs["condition"] = parsed["condition"].fillna(adata.obs[batch_key]).astype(str)
    adata.obs["replicate"] = parsed["replicate"].fillna("unknown").astype(str)
    adata.obs["method"]    = parsed["method"].fillna("unknown").astype(str)

    # JUDGMENT CALL: genotype/infection_status split on condition's trailing
    # "_wMel" / "_uninf" suffix (matches samples.csv today: OreR_wMel,
    # OreR_uninf, Mei-P26_wMel, Mei-P26_uninf). A condition that doesn't end
    # in one of these is left with infection_status='unknown' and genotype
    # equal to the full condition string, rather than guessed -- check the
    # WARNING below if a future sample's condition doesn't follow this
    # pattern (e.g. a third infection label besides wMel/uninf).
    cond = adata.obs["condition"].astype(str)
    suffix_pattern = r"^(?P<genotype>.+)_(?P<infection_status>" + "|".join(_INFECTION_SUFFIXES) + r")$"
    split = cond.str.extract(suffix_pattern)
    no_match = split["genotype"].isna()
    if no_match.any():
        bad = sorted(cond[no_match].unique())
        print(f"  WARNING: {int(no_match.sum())} cells have a condition that "
              f"doesn't end in {_INFECTION_SUFFIXES} -- genotype/"
              f"infection_status left as 'unknown' for: {bad}")
    adata.obs["genotype"] = split["genotype"].fillna(cond).astype(str)
    adata.obs["infection_status"] = split["infection_status"].fillna("unknown").astype(str)

    # Kept for parity with the old integrate.py's output schema.
    adata.obs["bio_condition"] = adata.obs["condition"]
    return adata


# -----------------------------------------------------------------------------
# Diagnostics
# -----------------------------------------------------------------------------

def plot_diagnostics(query, label_cols, fig_dir):
    os.makedirs(fig_dir, exist_ok=True)
    sc.settings.figdir = fig_dir
    query_plot = query.copy()
    query_plot.obsm["X_umap"] = query_plot.obsm["X_umap_atlas"]

    print(f"\n-- Diagnostic plots -- writing to {fig_dir}/ --")
    sc.pl.umap(query_plot, color="genotype", save="_genotype.pdf",
               title="Genotype, all projected onto the atlas")
    sc.pl.umap(query_plot, color="infection_status", save="_infection_status.pdf",
               title="Infection status, all projected onto the atlas")
    sc.pl.umap(query_plot, color="source_file", save="_source_file.pdf",
               title="Every sample, atlas-projected")

    for col in label_cols:
        atlas_col = f"atlas_{col}"
        if atlas_col in query_plot.obs.columns:
            sc.pl.umap(query_plot, color=atlas_col, save=f"_{col}.pdf",
                       title=f"'{col}' transferred onto every cell")

    if "wolbachia_titer" in query_plot.obs.columns:
        sc.pl.umap(query_plot, color="wolbachia_titer", save="_wolbachia_titer.pdf",
                    title="Wolbachia titer, atlas-projected", cmap="viridis")

    # Composition comparison: cell-type proportions per sample.
    for col in label_cols:
        atlas_col = f"atlas_{col}"
        if atlas_col not in query.obs.columns:
            continue
        counts = query.obs.groupby(["source_file", atlas_col], observed=True).size().unstack(fill_value=0)
        comp = counts.div(counts.sum(axis=1), axis=0)
        comp.to_csv(os.path.join(fig_dir, f"composition_by_sample_{col}.csv"))

        fig, ax = plt.subplots(figsize=(max(8, len(comp) * 0.6), 6))
        comp.plot(kind="bar", stacked=True, ax=ax, colormap="tab20", legend=False)
        ax.set_ylabel(f"Fraction of cells ({atlas_col})")
        ax.set_title(f"Cell-type composition by sample -- {col}")
        ax.legend(bbox_to_anchor=(1.02, 1), loc="upper left", fontsize=7)
        plt.xticks(rotation=45, ha="right")
        _savefig(fig, os.path.join(fig_dir, f"composition_by_sample_{col}.pdf"))

    # Titer by atlas cell type, split by infection status -- this project's
    # actual question for this object (see module docstring).
    if "wolbachia_titer" in query.obs.columns:
        for col in label_cols:
            atlas_col = f"atlas_{col}"
            if atlas_col not in query.obs.columns:
                continue
            df = query.obs[[atlas_col, "infection_status", "wolbachia_titer"]].copy()
            order = df.groupby(atlas_col, observed=True)["wolbachia_titer"].median().sort_values(ascending=False).index
            fig, ax = plt.subplots(figsize=(max(8, len(order) * 0.5), 5))
            data = [df.loc[df[atlas_col] == ct, "wolbachia_titer"].dropna().values for ct in order]
            ax.boxplot(data, labels=order, showfliers=False)
            ax.set_ylabel("Wolbachia titer")
            ax.set_title(f"Wolbachia titer by atlas cell type -- {col}")
            plt.xticks(rotation=45, ha="right")
            _savefig(fig, os.path.join(fig_dir, f"titer_by_{col}.pdf"))
            df.groupby(atlas_col, observed=True)["wolbachia_titer"].agg(
                ["mean", "median", "std", "count"]
            ).to_csv(os.path.join(fig_dir, f"titer_by_{col}_summary.csv"))

    print(f"   Diagnostic plots complete -- see {fig_dir}/")


# -----------------------------------------------------------------------------
# Main
# -----------------------------------------------------------------------------

def main():
    parser = argparse.ArgumentParser(
        description="Project every filtered per-sample h5ad file onto a "
                     "FROZEN reference atlas PCA/UMAP in one shot (scanpy."
                     "tl.ingest), replacing integrate.py's joint BBKNN "
                     "re-clustering with pure reference projection."
    )
    parser.add_argument("--atlas", required=True,
                         help="Path to the reference atlas h5ad, e.g. "
                              "ovary_atlas/atlas/h5ad_unprocessed/"
                              "s_fca_biohub_ovary_10x.h5ad")
    parser.add_argument("--query", required=True, nargs="+",
                         help="ALL your filtered per-sample h5ad files, e.g. "
                              "results/filtered_h5ad/*.h5ad")
    parser.add_argument("--out_path", required=True)
    parser.add_argument("--label_cols", nargs="+", default=None,
                         help="obs column(s) in the atlas to transfer. Omit "
                              "to print candidate columns from the atlas "
                              "and exit.")
    parser.add_argument("--flybase_annotation", default=None,
                         help="Used to remap atlas gene symbols to FBgn IDs "
                              "if the atlas isn't already FBgn-indexed.")
    parser.add_argument("--k", type=int, default=30)
    parser.add_argument("--n_pcs", type=int, default=30)
    parser.add_argument("--n_top_genes", type=int, default=3000)
    parser.add_argument("--subsample_ref", type=int, default=None)
    parser.add_argument("--fig_dir", default=None)
    args = parser.parse_args()

    query_paths = []
    for pattern in args.query:
        matches = glob.glob(pattern)
        query_paths.extend(matches if matches else [pattern])
    query_paths = sorted(set(query_paths))
    print(f"Query files ({len(query_paths)}):")
    for p in query_paths:
        print(f"  {p}")

    ref = load_atlas_reference(
        args.atlas, args.label_cols,
        flybase_annotation=args.flybase_annotation,
        subsample_ref=args.subsample_ref,
    )

    query = load_all_query_files(query_paths)

    print("\n-- Extracting sample metadata (condition/replicate/method/"
          "genotype/infection_status) from source_file --")
    query = add_sample_metadata(query, batch_key="source_file")
    print(query.obs["condition"].value_counts().to_string())

    # Full-gene, log1p-normalised copy BEFORE the projection step restricts
    # query down to the atlas's HVG panel -- stashed as .raw afterwards,
    # same convention the old integrate.py's preprocess() uses. This keeps
    # Wolbachia (GQ*) genes available downstream even though they play no
    # role in the atlas projection itself.
    print("\n-- Building full-gene log1p-normalised .raw layer --")
    raw_full = query.copy()
    sc.pp.normalize_total(raw_full, target_sum=1e4)
    sc.pp.log1p(raw_full)

    ref = build_reference_embedding(
        ref, n_pcs=args.n_pcs, n_top_genes=args.n_top_genes, n_neighbors=args.k,
    )

    projected = project_query_onto_reference(ref, query, args.label_cols, k=args.k)

    # Reattach by obs_names (not position) -- robust regardless of any
    # internal reordering project_query_onto_reference's gene-padding/concat
    # steps may have done.
    projected.raw = raw_full[projected.obs_names].copy()

    # Named explicitly (not the scanpy-default 'X_umap'/'X_pca' keys) so
    # it's unambiguous if this ever sits next to another embedding on the
    # same object -- but ALSO keep the conventional 'X_umap'/'X_pca' keys
    # so generic scanpy tooling works against this object unmodified.
    projected.obsm["X_umap_atlas"] = projected.obsm["X_umap"]
    projected.obsm["X_pca_atlas"]  = projected.obsm["X_pca"]

    if args.label_cols:
        for col in args.label_cols:
            atlas_col = f"atlas_{col}"
            if atlas_col in projected.obs.columns:
                projected.obs[f"cell_type_{col}"] = projected.obs[atlas_col]

    out_dir = os.path.dirname(args.out_path)
    if out_dir:
        os.makedirs(out_dir, exist_ok=True)
    projected.write(args.out_path)
    print(f"\nWrote combined atlas-projected object -> {args.out_path}  "
          f"({projected.n_obs:,} cells x {projected.n_vars:,} genes)")

    if args.fig_dir:
        plot_diagnostics(projected, args.label_cols, args.fig_dir)

    print("\n" + "=" * 60)
    print("COMPLETE (atlas-projected integration)")
    print("=" * 60)
    print(f"-> {args.out_path}")
    print("Every cell now carries atlas_<label> / cell_type_<label> / "
          "atlas_<label>_confidence, genotype, infection_status, a "
          "full-gene .raw layer (host + Wolbachia genes), and "
          "obsm['X_umap']/['X_umap_atlas']. No 'leiden' column is written "
          "(no clustering step here -- use atlas_<label> for cell identity "
          "instead) and no per-cell-type Harmony re-embedding is done for "
          "your own titer/infection-status axis -- keep using the old "
          "integrate.py's Harmony output for that question.")


if __name__ == "__main__":
    main()

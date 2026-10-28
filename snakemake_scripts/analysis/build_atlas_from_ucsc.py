"""
build_atlas_from_ucsc.py
Build an h5ad reference atlas from a UCSC Cell Browser download
(exprMatrix.tsv.gz = genes x cells, meta.tsv = per-cell annotations),
in the shape integrate_via_atlas_projection.py expects.

    python snakemake_scripts/analysis/build_atlas_from_ucsc.py \
        --expr ovary_atlas/ucsc/exprMatrix.tsv.gz \
        --meta ovary_atlas/ucsc/meta.tsv \
        --out  ovary_atlas/atlas/fca_ovary_ucsc.h5ad
"""
import argparse
import gzip

import numpy as np
import pandas as pd
import scipy.sparse as sp
import anndata as ad

p = argparse.ArgumentParser()
p.add_argument("--expr", required=True)
p.add_argument("--meta", required=True)
p.add_argument("--out", required=True)
p.add_argument("--chunksize", type=int, default=1000, help="genes per chunk")
a = p.parse_args()

# Header: first field is the gene column name, the rest are cell barcodes
with gzip.open(a.expr, "rt") as f:
    header = f.readline().rstrip("\n").split("\t")
gene_col, cells = header[0], header[1:]
print(f"{len(cells):,} cells in expression matrix")

# Read the dense genes x cells TSV in chunks, keeping each chunk sparse
dtypes = {gene_col: str, **{c: np.float32 for c in cells}}
genes, blocks = [], []
for i, chunk in enumerate(pd.read_csv(a.expr, sep="\t", index_col=0,
                                      dtype=dtypes, chunksize=a.chunksize)):
    genes.extend(chunk.index.astype(str))
    blocks.append(sp.csr_matrix(chunk.to_numpy()))
    print(f"  read {len(genes):,} genes", flush=True)
X = sp.vstack(blocks).T.tocsr()  # -> cells x genes
del blocks

# Cell Browser gene IDs are sometimes "symbol|geneID"; keep the symbol
genes = pd.Index([g.split("|")[0] for g in genes])

meta = pd.read_csv(a.meta, sep="\t", index_col=0)
n_missing = (~pd.Index(cells).isin(meta.index)).sum()
if n_missing:
    print(f"WARNING: {n_missing:,} cells have no row in meta.tsv")
meta = meta.reindex(cells)
for c in meta.columns:
    if meta[c].dtype == object:
        meta[c] = meta[c].astype("category")

adata = ad.AnnData(X=X, obs=meta, var=pd.DataFrame(index=genes))
adata.var_names_make_unique()

# Raw counts or already normalised? Decides where the matrix is stored.
sample = adata.X.data[:1_000_000]
is_counts = sample.size > 0 and np.allclose(sample, np.round(sample))
if is_counts:
    adata.layers["counts"] = adata.X.copy()
    print("Matrix looks like integer counts -> also stored in layers['counts']")
else:
    print(f"NOTE: matrix is NOT integer counts (max {adata.X.max():.2f}); "
          "it is likely normalised/log1p. Stored as .X only.")

print(adata)
print("\nobs columns (pick label columns for config atlas_label_cols):")
for c in adata.obs.columns:
    print(f"  {c!r}: {adata.obs[c].nunique()} unique")
print(f"\nvar_names sample: {adata.var_names[:10].tolist()}")

adata.write_h5ad(a.out, compression="gzip")
print(f"Written -> {a.out}")

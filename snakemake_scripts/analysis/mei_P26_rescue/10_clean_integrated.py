#!/usr/bin/env python
"""Write a copy of the integrated object without the debris cells.

Debris is only identifiable after integration (06_label_qc.py works on the
atlas labels), so integrated.h5ad still contains it. Any analysis that
compares conditions from the integrated object (NMF, GSEA, cluster_stats,
summary_statistics, UMAP plots) should read integrated_clean.h5ad instead.

The embedding, labels, .raw and everything else are kept as is; only the cells
listed in debris_cells.tsv (obs_name column) are dropped.

Usage:
  python 10_clean_integrated.py --integrated results/integrated/integrated.h5ad \
      --exclude_cells results/rescue/label_qc/debris_cells.tsv \
      --out results/integrated/integrated_clean.h5ad
"""
import argparse

import anndata as ad
import pandas as pd


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--integrated", required=True)
    p.add_argument("--exclude_cells", required=True)
    p.add_argument("--out", required=True)
    args = p.parse_args()

    drop = set(pd.read_csv(args.exclude_cells, sep="\t")["obs_name"].astype(str))
    a = ad.read_h5ad(args.integrated)
    hit = a.obs_names.isin(drop)
    if len(drop) and not hit.any():
        raise ValueError("no debris obs_names found in the integrated object; check the "
                         "obs_name convention (<sample>__<barcode>)")
    print(f"Dropping {hit.sum():,} of {a.n_obs:,} cells "
          f"({len(drop) - hit.sum():,} listed cells not in the integrated object)")
    sample_col = "source_file" if "source_file" in a.obs.columns else None
    if sample_col:
        print(pd.crosstab(a.obs[sample_col].astype(str), hit).rename(columns={False: "kept", True: "dropped"})
              .to_string())
    a = a[~hit].copy()
    a.write_h5ad(args.out, compression="gzip")
    print(f"Wrote {a.n_obs:,} cells to {args.out}")


if __name__ == "__main__":
    main()

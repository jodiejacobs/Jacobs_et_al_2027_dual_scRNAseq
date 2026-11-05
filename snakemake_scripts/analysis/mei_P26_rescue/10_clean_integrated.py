#!/usr/bin/env python
"""Write the analysis cell set from the integrated object.

integrated.h5ad holds every barcode that survived knee filtering, most of which
the atlas could not label confidently: in the Nov 2026 data 534,168 of 758,471
cells have label confidence below 0.5, with a median of 47 genes. Those are
ambient droplets, not cells, and they swamp anything run on the whole object.

Two filters are applied, giving exactly the cell set the rescue analysis uses:

  --min_confidence  drop cells whose atlas label confidence is below this
                    (same cutoff as 01_build_pseudobulk.py, config
                    rescue_min_confidence)
  --exclude_cells   drop the debris clusters flagged by 06_label_qc.py

Anything that compares conditions from the integrated object (NMF, GSEA,
cluster_stats, summary_statistics, UMAP plots) should read the output,
integrated_clean.h5ad, rather than integrated.h5ad. The embedding, labels,
.raw and everything else are kept as they are.

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
    p.add_argument("--label_col", default="cell_type_annotation",
                   help="label column whose confidence gates the cells")
    p.add_argument("--min_confidence", type=float, default=0.5,
                   help="drop cells below this atlas label confidence (0 keeps everything)")
    args = p.parse_args()

    drop = set(pd.read_csv(args.exclude_cells, sep="\t")["obs_name"].astype(str))
    a = ad.read_h5ad(args.integrated)
    hit = a.obs_names.isin(drop)
    if len(drop) and not hit.any():
        raise ValueError("no debris obs_names found in the integrated object; check the "
                         "obs_name convention (<sample>__<barcode>)")
    print(f"Dropping {hit.sum():,} of {a.n_obs:,} cells "
          f"({len(drop) - hit.sum():,} listed cells not in the integrated object)")

    conf_col = args.label_col.replace("cell_type_", "atlas_") + "_confidence"
    if args.min_confidence > 0:
        if conf_col not in a.obs.columns:
            raise ValueError(f"--min_confidence set but {conf_col} is not in obs; "
                             f"available: {list(a.obs.columns)}")
        low = a.obs[conf_col].astype(float) < args.min_confidence
        print(f"Dropping {int((low & ~hit).sum()):,} further cells below "
              f"{conf_col} {args.min_confidence} "
              f"(median {int(a.obs.loc[low, 'n_genes'].median()) if 'n_genes' in a.obs else '?'} genes)")
        hit = hit | low.to_numpy()
    sample_col = "source_file" if "source_file" in a.obs.columns else None
    if sample_col:
        print(pd.crosstab(a.obs[sample_col].astype(str), hit).rename(columns={False: "kept", True: "dropped"})
              .to_string())
    a = a[~hit].copy()
    a.write_h5ad(args.out, compression="gzip")
    print(f"Wrote {a.n_obs:,} cells to {args.out}")


if __name__ == "__main__":
    main()

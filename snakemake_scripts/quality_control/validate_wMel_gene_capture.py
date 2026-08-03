import scanpy as sc
import numpy as np
import pandas as pd
from scipy.sparse import issparse
from argparse import ArgumentParser

parser = ArgumentParser()
parser.add_argument("--input", help="Input h5ad file")
parser.add_argument("--output", help="Output file for results (optional)", default=None)
args = parser.parse_args()

if args.output:
    with open(args.output, "w") as f:
        f.write("Validation results for wMel gene capture\n")

adata = sc.read_h5ad(args.input)

var_names = adata.var_names
wMel_mask = var_names.str.startswith("G")
Dmel_mask = var_names.str.startswith("F")

# Total genes in the reference, by prefix
print(f"G-prefixed genes in var: {wMel_mask.sum()}")
print(f"F-prefixed genes in var: {Dmel_mask.sum()}")

# Genes actually detected (nonzero count in >=1 cell)
counts_per_gene = adata.X.sum(axis=0)
counts_per_gene = np.asarray(counts_per_gene).flatten() if issparse(adata.X) else np.asarray(counts_per_gene).flatten()
detected = counts_per_gene > 0

print(f"wMel genes detected (>0 counts): {(wMel_mask & detected).sum()}")
print(f"Dmel genes detected (>0 counts): {(Dmel_mask & detected).sum()}")

# Total UMI counts assigned to each group
print((f'Percent of total UMIs assigned to wMel genes: {counts_per_gene[wMel_mask].sum() / counts_per_gene.sum() * 100:.2f}%'))
print(f"Total UMIs -> wMel genes: {counts_per_gene[wMel_mask].sum()}")
print(f"Total UMIs -> Dmel genes: {counts_per_gene[Dmel_mask].sum()}")

# --- List of detected wMel genes, with counts, sorted descending ---
wMel_detected_mask = wMel_mask & detected
wMel_gene_names = var_names[wMel_detected_mask]
wMel_gene_counts = counts_per_gene[wMel_detected_mask]

wMel_df = pd.DataFrame({
    "gene": wMel_gene_names,
    "counts": wMel_gene_counts
}).sort_values("counts", ascending=False).reset_index(drop=True)

print("\nDetected wMel genes (sorted by counts):")
print(wMel_df.to_string(index=False))

# Write results to output file if specified
if args.output:
    with open(args.output, "a") as f:
        f.write(f"G-prefixed genes in var: {wMel_mask.sum()}\n")
        f.write(f"F-prefixed genes in var: {Dmel_mask.sum()}\n")
        f.write(f"wMel genes detected (>0 counts): {(wMel_mask & detected).sum()}\n")
        f.write(f"Dmel genes detected (>0 counts): {(Dmel_mask & detected).sum()}\n")
        f.write(f'Percent of total UMIs assigned to wMel genes: {counts_per_gene[wMel_mask].sum() / counts_per_gene.sum() * 100:.2f}%\n')
        f.write(f"Total UMIs -> wMel genes: {counts_per_gene[wMel_mask].sum()}\n")
        f.write(f"Total UMIs -> Dmel genes: {counts_per_gene[Dmel_mask].sum()}\n")
        f.write("\nDetected wMel genes (sorted by counts):\n")
        f.write(wMel_df.to_string(index=False))
        f.write("\n")

    # Sidecar CSV for easy reloading (e.g. for plotting top wMel genes)
    csv_path = args.output.rsplit(".", 1)[0] + "_wMel_genes.csv"
    wMel_df.to_csv(csv_path, index=False)
    print(f"\nWrote detected wMel gene list to: {csv_path}")
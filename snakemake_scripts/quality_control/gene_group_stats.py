#!/usr/bin/env python
"""
Compute per-gene and per-cell count statistics for gene groups defined by
var_name prefixes in an AnnData object -- e.g. Wolbachia (GQ...) vs
Drosophila (FBgn...) gene IDs in a combined host+symbiont alignment.

Outputs three CSVs into --outdir:
    gene_stats.csv     per-gene total counts, n_cells_expressing, group
    group_summary.csv  per-group totals, % of library, and per-cell pct stats
    cell_stats.csv      per-cell total counts and per-group counts/pct

Also writes a 2x2 QC figure (PNG + PDF, supplement-ready) into --outdir:
    {sample}_gene_group_qc.png / .pdf
        A. per-cell titer distribution for --titer-group
        B. group composition (n genes, % of library)
        C. per-gene detection (total counts vs. % cells expressing)
        D. per-cell total counts vs. titer

Usage:
    python gene_group_stats.py --input adata.h5ad --outdir stats/ \
        --prefix Wolbachia:GQ --prefix Drosophila:FB

    # if raw counts are in a layer instead of adata.raw:
    python gene_group_stats.py --input adata.h5ad --outdir stats/ \
        --source layer --layer-name counts
"""

import argparse
import sys
from pathlib import Path

import matplotlib
matplotlib.use("Agg")  # headless -- no display on SLURM compute nodes
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
from scipy.sparse import issparse

# Colorblind-safe categorical palette (Okabe-Ito), assigned in --prefix order
# so colors stay consistent across all four panels. Extra groups beyond this
# cycle fall back to gray, same as the 'other' bucket.
_PALETTE_CYCLE = ['#D55E00', '#0072B2', '#009E73', '#E69F00', '#CC79A7', '#56B4E9']
_OTHER_COLOR = '#999999'


def get_group_colors(labels):
    """Assign a stable colorblind-safe color to each group label, plus 'other'."""
    colors = {label: _PALETTE_CYCLE[i % len(_PALETTE_CYCLE)] for i, label in enumerate(labels)}
    colors['other'] = _OTHER_COLOR
    return colors


def parse_prefixes(prefix_args):
    groups = {}
    for item in prefix_args:
        if ':' not in item:
            raise ValueError(f"--prefix must be LABEL:PREFIX, got '{item}'")
        label, prefix = item.split(':', 1)
        groups[label] = prefix
    return groups


def get_counts_matrix(adata, source, layer_name):
    """Return (X, var_names) for the requested counts source."""
    if source == 'raw':
        if adata.raw is None:
            sys.exit("Error: --source raw requested but adata.raw is None.")
        return adata.raw.X, adata.raw.var_names
    if source == 'layer':
        if layer_name not in adata.layers:
            sys.exit(
                f"Error: layer '{layer_name}' not found. "
                f"Available layers: {list(adata.layers.keys())}"
            )
        return adata.layers[layer_name], adata.var_names
    return adata.X, adata.var_names


def to_flat_array(mat):
    """Flatten a sparse or dense matrix-sum result to a 1D numpy array."""
    return np.asarray(mat).flatten()


def make_qc_plots(gene_stats, cell_stats, group_summary, groups, titer_group, outdir, sample_name):
    """Build a 2x2 QC figure and save as PNG (quick look) + PDF (supplement-ready).

    Panels: A. per-cell titer distribution for `titer_group`; B. group
    composition (n genes, % of library); C. per-gene detection (total counts
    vs. % cells expressing); D. per-cell total counts vs. titer.
    """
    colors = get_group_colors(list(groups.keys()))
    pct_col = f'pct_{titer_group}'

    fig, axes = plt.subplots(2, 2, figsize=(8, 6.5))
    fig.suptitle(sample_name, fontsize=10, fontweight='bold')

    # A. Per-cell titer distribution
    ax = axes[0, 0]
    vals = cell_stats[pct_col].dropna()
    ax.hist(vals, bins=50, color=colors.get(titer_group, _OTHER_COLOR), edgecolor='none')
    median_val = vals.median()
    ax.axvline(median_val, color='black', linestyle='--', linewidth=0.8)
    ax.text(0.98, 0.95, f'median = {median_val:.1f}%', transform=ax.transAxes,
            ha='right', va='top', fontsize=7)
    ax.set_xlabel(f'% {titer_group} of transcripts per cell')
    ax.set_ylabel('Number of cells')
    ax.set_title('A. Per-cell titer distribution', loc='left', fontsize=8)

    # B. Group composition -- n genes (left axis) and % of library (right axis)
    ax = axes[0, 1]
    ax2 = ax.twinx()
    labels = group_summary['group'].tolist()
    x = np.arange(len(labels))
    width = 0.35
    ax.bar(x - width / 2, group_summary['n_genes'], width, color=_OTHER_COLOR)
    ax2.bar(x + width / 2, group_summary['pct_of_library'], width,
            color=[colors.get(g, _OTHER_COLOR) for g in labels])
    ax.set_xticks(x)
    ax.set_xticklabels(labels)
    ax.set_ylabel('Number of genes')
    ax2.set_ylabel('% of library (counts)')
    ax.set_title('B. Group composition', loc='left', fontsize=8)
    ax.bar([np.nan], [np.nan], color=_OTHER_COLOR, label='n genes')
    ax.bar([np.nan], [np.nan], color='black', label='% of library')
    ax.legend(fontsize=6, loc='upper right')

    # C. Per-gene detection: total counts vs. % cells expressing, by group
    ax = axes[1, 0]
    for label in list(groups.keys()) + ['other']:
        sub = gene_stats[gene_stats['group'] == label]
        if sub.empty:
            continue
        ax.scatter(sub['total_counts'] + 1, sub['pct_cells_expressing'],
                   s=4, alpha=0.5, color=colors.get(label, _OTHER_COLOR),
                   label=label, linewidths=0, rasterized=True)
    ax.set_xscale('log')
    ax.set_xlabel('Total counts + 1 (log scale)')
    ax.set_ylabel('% cells expressing')
    ax.set_title('C. Per-gene detection', loc='left', fontsize=8)
    ax.legend(fontsize=6, loc='lower right', markerscale=2)

    # D. Per-cell total counts vs. titer
    ax = axes[1, 1]
    ax.scatter(cell_stats['total_counts'] + 1, cell_stats[pct_col],
               s=3, alpha=0.3, color=colors.get(titer_group, _OTHER_COLOR),
               linewidths=0, rasterized=True)
    ax.set_xscale('log')
    ax.set_xlabel('Total counts per cell + 1 (log scale)')
    ax.set_ylabel(f'% {titer_group} of transcripts')
    ax.set_title('D. Titer vs. library size', loc='left', fontsize=8)

    fig.tight_layout(rect=[0, 0, 1, 0.96])

    png_path = outdir / f'{sample_name}_gene_group_qc.png'
    pdf_path = outdir / f'{sample_name}_gene_group_qc.pdf'
    fig.savefig(png_path, dpi=300)
    fig.savefig(pdf_path, dpi=300)  # dpi controls the rasterized scatter layers only
    plt.close(fig)
    return png_path, pdf_path


def main():
    parser = argparse.ArgumentParser(
        description=__doc__, formatter_class=argparse.RawDescriptionHelpFormatter
    )
    parser.add_argument('--input', required=True, help='Path to .h5ad AnnData file')
    parser.add_argument('--outdir', required=True, help='Directory to write stats files into')
    parser.add_argument(
        '--prefix', action='append', default=[],
        help='LABEL:PREFIX pair for a gene group, e.g. Wolbachia:GQ. '
             'Repeat for multiple groups. Default: Wolbachia:GQ Drosophila:FB'
    )
    parser.add_argument(
        '--source', choices=['raw', 'layer', 'X'], default='raw',
        help="Where to pull counts from: adata.raw.X ('raw', default), "
             "adata.layers[--layer-name] ('layer'), or adata.X ('X')."
    )
    parser.add_argument(
        '--layer-name', default='counts',
        help="Layer name to use when --source layer (default: 'counts')"
    )
    parser.add_argument(
        '--titer-group', default=None,
        help="Group label to use for the titer-distribution and titer-vs-library-size "
             "plots (panels A and D). Defaults to whichever --prefix label contains "
             "'wolbachia' (case-insensitive); falls back to the first group given."
    )
    args = parser.parse_args()

    prefix_args = args.prefix or ['Wolbachia:GQ', 'Drosophila:FB']
    groups = parse_prefixes(prefix_args)

    # JUDGMENT CALL: which group counts as "titer" for panels A/D. Auto-detects
    # a label containing "wolbachia"; pass --titer-group explicitly to override.
    if args.titer_group:
        if args.titer_group not in groups:
            sys.exit(f"Error: --titer-group '{args.titer_group}' not among group labels: {list(groups)}")
        titer_group = args.titer_group
    else:
        titer_group = next((g for g in groups if 'wolbachia' in g.lower()), next(iter(groups)))

    outdir = Path(args.outdir)
    outdir.mkdir(parents=True, exist_ok=True)

    print(f"Loading {args.input} ...")
    adata = sc.read_h5ad(args.input)

    X, var_names = get_counts_matrix(adata, args.source, args.layer_name)
    var_names = pd.Index(var_names)

    # Sanity check: warn if this looks like scaled/normalized data (negative
    # values), which would make count-based stats meaningless.
    sample = X[:1000] if X.shape[0] > 1000 else X
    sample_arr = sample.toarray() if issparse(sample) else np.asarray(sample)
    if (sample_arr < 0).any():
        print(
            f"WARNING: negative values detected in the '{args.source}' matrix -- "
            "this looks like scaled/transformed data, not raw counts. "
            "Stats below may not be meaningful as counts. "
            "Consider --source raw or --source layer with --layer-name counts.",
            file=sys.stderr,
        )

    # Assign each gene to a group by prefix (first match wins; else 'other')
    gene_group = pd.Series('other', index=var_names)
    for label, prefix in groups.items():
        gene_group[var_names.str.startswith(prefix)] = label

    # --- per-gene stats ---
    total_counts = to_flat_array(X.sum(axis=0))
    n_cells_expr = to_flat_array((X > 0).sum(axis=0))
    pct_cells_expr = n_cells_expr / X.shape[0] * 100

    gene_stats = pd.DataFrame({
        'gene': var_names,
        'group': gene_group.values,
        'total_counts': total_counts,
        'n_cells_expressing': n_cells_expr,
        'pct_cells_expressing': pct_cells_expr,
    }).sort_values(['group', 'total_counts'], ascending=[True, False])

    gene_stats_path = outdir / 'gene_stats.csv'
    gene_stats.to_csv(gene_stats_path, index=False)
    print(f"Wrote per-gene stats: {gene_stats_path}")

    # --- per-cell stats ---
    cell_stats = pd.DataFrame(index=adata.obs_names)
    cell_stats['total_counts'] = to_flat_array(X.sum(axis=1))
    for label, prefix in groups.items():
        mask = var_names.str.startswith(prefix)
        group_counts = to_flat_array(X[:, mask].sum(axis=1))
        cell_stats[f'counts_{label}'] = group_counts
        cell_stats[f'pct_{label}'] = np.where(
            cell_stats['total_counts'] > 0,
            group_counts / cell_stats['total_counts'] * 100,
            np.nan,
        )

    cell_stats_path = outdir / 'cell_stats.csv'
    cell_stats.to_csv(cell_stats_path)
    print(f"Wrote per-cell stats: {cell_stats_path}")

    # --- group-level summary ---
    # Gene-level totals (sum of raw counts across all cells, per group)
    group_summary = gene_stats.groupby('group').agg(
        n_genes=('gene', 'count'),
        total_counts=('total_counts', 'sum'),
        mean_counts_per_gene=('total_counts', 'mean'),
    ).reset_index()
    library_total = total_counts.sum()
    group_summary['pct_of_library'] = group_summary['total_counts'] / library_total * 100

    # Per-cell pct stats (e.g. average of pct_GQ across cells), skipping the
    # NaNs assigned above for cells with zero total counts.
    pct_rows = []
    for label in groups:
        pct_col = cell_stats[f'pct_{label}']
        pct_rows.append({
            'group': label,
            'mean_pct_per_cell': pct_col.mean(),
            'median_pct_per_cell': pct_col.median(),
            'std_pct_per_cell': pct_col.std(),
        })
    pct_summary = pd.DataFrame(pct_rows)

    group_summary = group_summary.merge(pct_summary, on='group', how='left')

    summary_path = outdir / 'group_summary.csv'
    group_summary.to_csv(summary_path, index=False)
    print(f"Wrote group summary: {summary_path}")

    # --- QC figure (PNG + PDF) ---
    sample_name = outdir.name
    png_path, pdf_path = make_qc_plots(
        gene_stats, cell_stats, group_summary, groups, titer_group, outdir, sample_name
    )
    print(f"Wrote QC figure: {png_path}")
    print(f"Wrote QC figure: {pdf_path}")

    print("Done.")


if __name__ == '__main__':
    main()

# mei-P26 rescue analysis

Pseudobulk analysis by atlas cell type (`cell_type_annotation_broad` by default) for the four-condition design: OreR / Mei-P26 × wMel / uninfected, 2 replicates each.

| Question | Script | Main output (`results/rescue/`) |
|---|---|---|
| Wolbachia titer across tissues, Mei-P26 vs OreR | `02_wolbachia_titer.py` | `titer/titer_genotype_by_tissue.tsv`, `titer_by_tissue.pdf`, `titer_ratio_forest.pdf` |
| Wolbachia gene expression by tissue, Mei-P26 vs OreR | `03_edger_de.R` | `de/wolbachia_de_all.tsv.gz` (per tissue + across-tissue F tests) |
| Host: infected vs uninfected, per tissue | `03_edger_de.R` | `de/host_de_all.tsv.gz` contrasts `infection_in_OreR`, `infection_in_MeiP26`, `infection_main` |
| Host: Mei-P26 vs OreR, uninfected | `03_edger_de.R` | contrast `genotype_in_uninf` |
| ~infection + genotype + infection:genotype by tissue | `03_edger_de.R`, `05_rescue_summary.py` | contrast `interaction`; `summary/rescue_genes.tsv`, `rescue_scatter.pdf` |
| Cell-type composition | `04_cell_composition.py` | `composition/composition_tests.tsv`, `composition.pdf` |

## Run

```bash
# one-time: R env for edgeR
mamba create -n edger -c conda-forge -c bioconda bioconductor-edger r-base
snakemake rescue_all --executor slurm --default-resources slurm_partition=medium runtime=120 mem_mb=8000 -j 16
```

## Settings (config.yaml, all optional)

| key | default | effect |
|---|---|---|
| `edger_env` | `edger` | conda env with R + edgeR |
| `rescue_label_col` | `cell_type_annotation_broad` | tissue resolution; `cell_type_annotation` is finer but splits cells thinner |
| `rescue_min_confidence` | 0.5 | drop cells with atlas label confidence below this |
| `rescue_min_cells` | 20 | min cells for a sample × tissue pseudobulk to be used |
| `rescue_min_wmel_umis` | 2000 | min wMel UMIs for a pseudobulk in Wolbachia DE |
| `rescue_block_replicate` | false | add replicate to the host design (costs 1 of 4 residual df) |
| `rescue_drop_qc_outliers` | false | also drop `qc_outlier` cells from integration |

## Interpretation notes

- **Replication.** Each tissue model has 8 pseudobulks and 4 residual df. edgeR QL with robust dispersion is the right tool at this n, but expect modest power, especially for the interaction.
- **Wolbachia DE is composition, not titer.** wMel genes are TMM-normalized among themselves, so a logFC means a gene makes up a larger share of the Wolbachia transcriptome. Titer is in `02`.
- **Ambient wMel reads.** All libraries share a lane, so uninfected samples carry wMel UMIs from ambient RNA and index hopping. `02` uses the uninfected 99th percentile as a per-cell floor, and `de/wolbachia_ambient_check.tsv` gives infected/uninfected titer per tissue. Treat Wolbachia DE in tissues with a low ratio with caution.
- **Titer p-values** come from a quasi-binomial GLM. Pseudobulks from the same sample are correlated across tissues, and a sample random effect isn't estimable with 2 replicates, so check replicate agreement in `titer_by_tissue.pdf`.
- **Genotype contrasts** compare Mei-P26 to OreR, so background differences beyond the mei-P26 allele also land in `genotype_*`.
- **Rescue categories** (`05`) use descriptive cutoffs (rescue fraction 0.25 / 0.75). `FDR_interaction` is the formal test.

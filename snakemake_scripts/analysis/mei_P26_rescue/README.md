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
| What are the artefact / unannotated / tracheal labels? Debris flag | `06_label_qc.py` | `label_qc/marker_heatmap.pdf`, `subclusters_<label>.tsv`, `debris_cells.tsv` |
| Marker-defined germline: abundance, subclusters, germline pseudobulks | `07_germline.py` (+ `03`, `08` on its pseudobulks) | `germline/germline_composition_tests.tsv`, `germline/de/`, `germline/rescue_null/` |
| Is rescue larger than selection bias? Interaction-supported rescue genes | `08_rescue_null.R` | `rescue_null/rescue_null_summary.tsv`, `rescue_supported_genes.tsv` |
| Wolbachia DE robustness (replicates, normalisation, pooled 2 vs 2) | `09_wolbachia_followup.R` | `wolbachia_followup/wolbachia_candidates.tsv` |

## Run

```bash
# one-time: R env for edgeR
mamba create -n edger -c conda-forge -c bioconda bioconductor-edger r-base
snakemake rescue_all rescue_followup_all --executor slurm --default-resources slurm_partition=medium runtime=120 mem_mb=8000 -j 16
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
| `rescue_drop_debris` | true | build pseudobulks without the debris cells flagged by `06` (see below) |
| `wmel_gff` | (none) | RefSeq wMel GFF3; adds gene products to `09` output |

## Debris filter

`rescue_pseudobulk_all` runs `01` on every labeled cell; `06` subclusters the
artefact, unannotated and adult tracheal labels and flags a subcluster as debris
when its median cell has < 120 genes and > 4% mito (`--debris_max_genes`,
`--debris_min_mito`); `rescue_pseudobulk` then reruns `01 --exclude_cells
label_qc/debris_cells.tsv`, and every downstream rule uses that output.
On the Nov 2026 data this removes artefact c1 and tracheal c1/c2 (46,400 cells;
30,460 from Mei-P26_wMel-1). The filter works on subclusters, not single cells:
a per-cell cutoff (e.g. < 200 genes and >= 3% mito) catches the same debris but
also ~50% of every Mei-P26 sample and ~5% of OreR, which would confound genotype.

## Interpretation notes

- **Replication.** Each tissue model has 8 pseudobulks and 4 residual df. edgeR QL with robust dispersion is the right tool at this n, but expect modest power, especially for the interaction.
- **Wolbachia DE is composition, not titer.** wMel genes are TMM-normalized among themselves, so a logFC means a gene makes up a larger share of the Wolbachia transcriptome. Titer is in `02`.
- **Ambient wMel reads.** All libraries share a lane, so uninfected samples carry wMel UMIs from ambient RNA and index hopping. `02` uses the uninfected 99th percentile as a per-cell floor, and `de/wolbachia_ambient_check.tsv` gives infected/uninfected titer per tissue. Treat Wolbachia DE in tissues with a low ratio with caution.
- **Titer p-values** come from a quasi-binomial GLM. Pseudobulks from the same sample are correlated across tissues, and a sample random effect isn't estimable with 2 replicates, so check replicate agreement in `titer_by_tissue.pdf`.
- **Genotype contrasts** compare Mei-P26 to OreR, so background differences beyond the mei-P26 allele also land in `genotype_*`.
- **Rescue categories** (`05`) use descriptive cutoffs (rescue fraction 0.25 / 0.75). `FDR_interaction` is the formal test.

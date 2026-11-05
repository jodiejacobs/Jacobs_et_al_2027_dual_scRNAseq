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
| Integrated object without debris, for NMF / GSEA / UMAP / summary scripts | `10_clean_integrated.py` | `results/integrated/integrated_clean.h5ad` |
| Which wMel genes could act in host cells (secretion class + evidence + DE) | `11_wolbachia_secretion.py` | `wolbachia_followup/wolbachia_secretion_all.tsv`, `_shortlist.tsv` |
| Is the wMel response a functional shift? GO / KEGG / product-category enrichment | `12_wolbachia_gsea.py` | `wolbachia_followup/wolbachia_gsea.tsv` |

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
when its median cell has < 150 genes and > 4.5% mito (`--debris_max_genes`,
`--debris_min_mito`); `rescue_pseudobulk` then reruns `01 --exclude_cells
label_qc/debris_cells.tsv`, and every downstream rule uses that output.
Thresholds were set from two independent clusterings of the Nov 2026 data:
the junk clusters have median 77-124 genes and 4.9-7.2% mito (mito, rRNA and
ambient yolk/chorion markers, >90% Mei-P26), while the nearest real clusters
are ovarian sheath muscle (124 genes, 4.1% mito) and tracheal c8 (176 genes).
Leiden output differs slightly between environments, so check the `debris`
column of `label_qc/subclusters_*.tsv` after each run. The filter works on subclusters, not single cells:
a per-cell cutoff (e.g. < 200 genes and >= 3% mito) catches the same debris but
also ~50% of every Mei-P26 sample and ~5% of OreR, which would confound genotype.

Every analysis that compares conditions runs without the debris: `02`-`05`
and `07`-`09` read the filtered pseudobulks, `06` builds its marker heatmap and
lineage tables from non-debris cells (its subcluster tables keep the debris,
since that is where it is flagged), and `10` writes
`results/integrated/integrated_clean.h5ad` for scripts that read the integrated
object. `integrated.h5ad` itself and the `integrate` rule's diagnostic figures
still contain the debris; use the rescue outputs or `integrated_clean.h5ad`
instead. Per-sample counts are in `label_qc/cell_counts_by_sample.tsv`.

## Interpretation notes

- **Replication.** Each tissue model has 8 pseudobulks and 4 residual df. edgeR QL with robust dispersion is the right tool at this n, but expect modest power, especially for the interaction.
- **Wolbachia DE is composition, not titer.** wMel genes are TMM-normalized among themselves, so a logFC means a gene makes up a larger share of the Wolbachia transcriptome. Titer is in `02`.
- **Ambient wMel reads.** All libraries share a lane, so uninfected samples carry wMel UMIs from ambient RNA and index hopping. `02` uses the uninfected 99th percentile as a per-cell floor, and `de/wolbachia_ambient_check.tsv` gives infected/uninfected titer per tissue. Treat Wolbachia DE in tissues with a low ratio with caution.
- **Titer p-values** come from a quasi-binomial GLM. Pseudobulks from the same sample are correlated across tissues, and a sample random effect isn't estimable with 2 replicates, so check replicate agreement in `titer_by_tissue.pdf`.
- **Genotype contrasts** compare Mei-P26 to OreR, so background differences beyond the mei-P26 allele also land in `genotype_*`.
- **Rescue categories** (`05`) use descriptive cutoffs (rescue fraction 0.25 / 0.75). `FDR_interaction` is the formal test.

## Wolbachia secretion annotation (11)

Run as `rule rescue_secretion` and `rule rescue_gsea`, both part of
`rescue_followup_all`:

```bash
snakemake rescue_followup_all --executor slurm -j 16
```

The protein FASTA and the signal-peptide calls are params, not inputs. When
`wmel_proteins` / `wmel_signalp` are unset in config.yaml, or point at files
that do not exist, the rule still runs: `secretion_class` reads "unknown" for
every locus, the "uncharacterised" tier drops out, and ranking falls back to
published evidence, eukaryote-like domains and expression. `wmel_gff` IS
required; without it the rule stops with a message.

| config key | default | effect |
|---|---|---|
| `wmel_gff` | (none) | annotation with Wolbachia locus tags; required |
| `wmel_proteins` | `reference/wMel_proteins/wmel_proteins.faa` | protein lengths + the transmembrane screen |
| `wmel_signalp` | `reference/wMel_proteins/wmel_deepsig.gff` | signal peptides (DeepSig GFF, or TSV of protein_id + Y/N) |
| `rescue_secretion_min_aa` | 400 | length floor for the "uncharacterised" tier |
| `rescue_secretion_min_cpm` | 100 | germline CPM floor for the shortlist |

`wmel_proteins.faa` is the RefSeq wMel proteome (GCF_000008025.1, same assembly
as the GTF) and `wmel_deepsig.gff` holds DeepSig signal-peptide calls for it.
Both ship in `reference/wMel_proteins/`. To redo the predictions, run the script
directly with `--write_faa out.faa` to dump the proteins named in the GTF for a
web predictor (DeepTMHMM, SignalP), or use `run_deepsig.py`, which runs DeepSig
locally under Keras 3.

A Sec signal peptide means the bacterial envelope, not the host cytoplasm;
host-cytoplasmic delivery goes through the T4SS, whose signal is C-terminal and
not predictable from sequence. Rank by the evidence tier, not the class, and
never drop a candidate because the Kyte-Doolittle screen calls it membrane: it
calls TomO and WalE1 membrane, and both act in host cells.

## Wolbachia gene-set enrichment (12)

`rule rescue_gsea` scores three set families against the pooled Wolbachia DE
ranking: GO terms parsed from `wmel_gff`'s go_* attributes, KEGG pathways for
`wmel_kegg_org`, and product-name categories defined in the script. KEGG comes
from the cache in `wmel_kegg_dir` (shipped in `reference/wMel_proteins/`), so
the rule runs with no network; clear that key to fetch from the KEGG REST API.
Needs `gseapy` in the scanpy env.

On the Nov 2026 data every family is null (best p = 0.014, FDR 0.14, for mobile
elements). That is a statement about the annotation more than the biology: only
~210 of 1,154 wMel genes survive filterByExpr, and GO and KEGG annotate
housekeeping rather than host interaction, so 271 of 273 hypothetical proteins
and all of CifA, TomO and WalE1 carry no GO term at all. Treat the gene-level
results as the stronger read, and rerun this if the DE gains power.

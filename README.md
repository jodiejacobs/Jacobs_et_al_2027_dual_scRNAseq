# Jacobs et al. 2027 — Dual scRNA-seq of *Wolbachia*-Infected *Drosophila melanogaster* Ovaries

Snakemake pipeline for dual scRNA-seq analysis of *w*Mel-infected *D. melanogaster* ovaries, covering pseudoalignment, QC/filtering, cell cycle annotation, host-symbiont integration, and *Wolbachia* 16S rRNA read validation.

## Pipeline overview

Raw FASTQs (10x or PIPseq) are pseudoaligned with `kb count` (kallisto|bustools), QC'd and filtered with scanpy, annotated for cell cycle phase with cyclum, then integrated across samples. In parallel, reads are aligned to a combined *Dmel*/*w*Mel rRNA reference to validate 16S capture and quantify bacterial genus composition by BLAST.

```
FASTQ (10x)
  └─ map_10x (kb count)                     → h5ad, unfiltered BUS
       ├─ inspect_10x_corrected / summarize_inspect_10x  → barcode/UMI QC table
       ├─ validate_wMel_gene_capture                     → wMel vs Dmel gene capture report
       └─ filter_h5ad (scanpy)                            → filtered h5ad
            └─ annotate_cell_cycle (cyclum)                → annotated h5ad
                 └─ integrate                              → integrated.h5ad

FASTQ (R2)
  └─ align_gene_reads (bwa mem)             → target-gene BAM
       ├─ calculate_coverage (samtools depth)     → coverage.tsv → plot_coverage_by_group
       └─ extract_16s_sequences → blast_16s → summarize_blast   → plot_blast_by_group
                                             └─ extract_abundant_16s (optional)
```

Gene program / NMF and GSEA rules (`nmf_programs`, `nmf_continuous_var`, `nmf_categorical_var`, `nmf_annotate_programs`) are implemented in `snakemake_scripts/analysis/` but currently commented out in the `Snakefile`.

## Repo layout

```
Snakefile                   # main pipeline
config/
  config.yaml                # paths, resources, SLURM settings
  samples.csv                # sample sheet (condition, platform, replicate, FASTQ paths)
  scanpy_env.yml              # scanpy conda/mamba env
  env.yml                     # legacy scanpy_ipynb env
  envs/cyclum_env.yml         # cyclum (cell cycle) conda/mamba env
snakemake_scripts/
  filtering/                  # QC filtering of raw h5ad
  analysis/                   # cell cycle, integration, NMF/GSEA
    mei_P26_wMel_scripts/     # exploratory NMF/GSEA/titer analyses
  quality_control/            # cluster QC, wMel gene capture validation
  plotting/                   # UMAPs, QC histograms, summary stats
  rRNA_analysis/               # 16S read extraction, BLAST, coverage
```

## Requirements

- Snakemake ≥ 9.0
- mamba/conda environments: `scanpy` (config/scanpy_env.yml), `cyclum` (config/envs/cyclum_env.yml), plus `kallisto_bustools` and `sra-tools` (bwa/samtools/blast) environments referenced by shell rules but not bundled here
- SLURM cluster (rules are written with `slurm_partition`/`mem_mb`/`slurm_time` resources)

Set up environments:

```bash
mamba env create -f config/scanpy_env.yml
mamba env create -f config/envs/cyclum_env.yml
```

Then update `config/config.yaml` with the absolute paths to those envs (`scanpy_env`, `cyclum_env`), your FASTQ directory, and kallisto index.

## Configuration

Edit `config/config.yaml`:
- `samples_file` — path to the sample sheet
- `fastq_dir`, `kallisto_index`, `transcripts_to_genes` — alignment inputs
- `target_gene`, `rRNA_regions`, `ref_fasta`, `ref_bed`, `blast_db` — rRNA analysis
- per-rule thread/memory/time/partition settings for SLURM

`config/samples.csv` has no header; columns are:

```
condition, seq_platform (10x|pipseq), replicate, L005_R1, L005_R2, L006_R1, L006_R2
```

## Running

```bash
mamba activate snakemake

# dry run
snakemake -n

# submit to SLURM
snakemake --executor slurm \
  --default-resources slurm_partition=medium slurm_time="2:00:00" runtime=120 mem_mb=8000 \
  -j 16
```

Outputs land under `results/`, logs under `logs/`.

## Citation

Jacobs et al. 2027 (in preparation).

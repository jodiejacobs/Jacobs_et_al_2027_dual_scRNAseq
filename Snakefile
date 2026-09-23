# Snakefile for processing scRNA-seq data with Snakemake and kallisto bustools
# mamba activate snakemake #Needs snakemake>=9.0
# snakemake --executor slurm --default-resources slurm_partition=medium runtime=120 mem_mb=8000 -j 16 -n
# (rules set their own "runtime" resource in minutes; slurm_time is not read by the SLURM executor plugin)

import pandas as pd
import os

# Configuration
configfile: "config/config.yaml"
SCANPY_ENV = config["scanpy_env"]
CYCLUM_ENV = config["cyclum_env"]

# Run small, quick summary rules locally instead of submitting them as slurm jobs
localrules: summarize_inspect_10x


# Load samples information
samples_df = pd.read_csv(config["samples_file"], header=None, sep=',')
samples_df = samples_df.apply(lambda x: x.str.strip() if x.dtype == "object" else x)

# Change index format to use hyphen between condition and replicate
samples_df.index = [(row[0] + "-" + str(row[2]) + "_" + row[1]) for _, row in samples_df.iterrows()]

# Get actual sample IDs that exist
SAMPLE_IDS = samples_df.index.tolist()

CONDITIONS = samples_df[0].unique().tolist()
REPLICATES = samples_df[2].unique().tolist()
SEQUENCING_PLATFORMS = samples_df[1].unique().tolist()

# 10x sample IDs (used for the bustools inspect QC rules below)
SAMPLE_IDS_10X = [s for s in SAMPLE_IDS if s.endswith("_10x")]

# Create a dictionary to track which condition-seq_platform combinations exist
# and what replicates are available for each
CONDITION_PLATFORM_REPS = {}
for sample_id in SAMPLE_IDS:
    condition = samples_df.loc[sample_id][0]
    seq_platform = samples_df.loc[sample_id][1]
    replicate = samples_df.loc[sample_id][2]

    key = (condition, seq_platform)
    if key not in CONDITION_PLATFORM_REPS:
        CONDITION_PLATFORM_REPS[key] = []
    CONDITION_PLATFORM_REPS[key].append(replicate)

# Get list of condition-seq_platform combinations that actually exist
CONDITION_PLATFORM_COMBOS = list(CONDITION_PLATFORM_REPS.keys())

# Find conditions that have BOTH 10x and pipseq data for method comparison
CONDITIONS_WITH_BOTH_METHODS = []
for condition in CONDITIONS:
    has_10x = any(p == '10x' for c, p in CONDITION_PLATFORM_COMBOS if c == condition)
    has_pipseq = any(p == 'pipseq' for c, p in CONDITION_PLATFORM_COMBOS if c == condition)
    if has_10x and has_pipseq:
        CONDITIONS_WITH_BOTH_METHODS.append(condition)

print(f"Conditions with both methods for comparison: {CONDITIONS_WITH_BOTH_METHODS}")

# Sample lookup by sample name; samples_df.loc['sample_name'][3] gives the R1 path
print(f"Conditions: {CONDITIONS}")
print(f"Replicates: {REPLICATES}")
print(f"Sequencing Platforms: {SEQUENCING_PLATFORMS}")
print(f"Sample IDs: {SAMPLE_IDS}")
print(f"Condition-Platform combinations: {CONDITION_PLATFORM_COMBOS}")

# Helper function to get fastq files for a sample
def get_fastq_files(sample_id):
    """Get R1/R2 fastq paths for a sample.

    Each of the R1/R2 columns in samples.csv normally holds a single path.
    To add a second (or third, ...) set of reads for a sample - e.g. a
    top-up sequencing run or an extra lane - list the paths separated by
    ';' in both columns, in matching order (R1[i] pairs with R2[i]):
        cond,10x,1,run1_R1.fastq.gz;run2_R1.fastq.gz,run1_R2.fastq.gz;run2_R2.fastq.gz
    Returns two lists (r1_files, r2_files), each of length 1 for a
    single-lane sample.
    """
    sample_info = samples_df.loc[sample_id]
    r1_files = [p.strip() for p in str(sample_info[3]).split(";") if p.strip()]
    r2_files = [p.strip() for p in str(sample_info[4]).split(";") if p.strip()]
    if len(r1_files) != len(r2_files):
        raise ValueError(
            f"{sample_id}: R1 column lists {len(r1_files)} file(s) but R2 column "
            f"lists {len(r2_files)} - samples.csv must list the same number of "
            f"';'-separated R1 and R2 files, in matching order, for each sample."
        )
    return r1_files, r2_files

# Interleave a sample's R1/R2 files the way kb count expects multiple
# lanes/runs to be passed: R1_a R2_a R1_b R2_b ...
def get_kb_reads(sample_id):
    r1_files, r2_files = get_fastq_files(sample_id)
    interleaved = []
    for r1, r2 in zip(r1_files, r2_files):
        interleaved += [r1, r2]
    return interleaved

# snakemake-executor-plugin-slurm reads wall-time from the "runtime"
# resource (integer minutes) - it does NOT recognize "slurm_time", so a
# rule that only sets slurm_time submits with no --time at all and falls
# back to the cluster/partition default. This converts the "H:MM:SS" (or
# "M:SS") strings used throughout this Snakefile/config.yaml into the
# minutes the plugin actually reads, so rules keep their human-readable
# time strings but the SLURM submission gets the wall time we intend.
def hms_to_minutes(time_str):
    parts = [int(p) for p in str(time_str).split(":")]
    while len(parts) < 3:
        parts.insert(0, 0)
    hours, minutes, seconds = parts
    total_minutes = hours * 60 + minutes + (1 if seconds else 0)
    return max(total_minutes, 1)

# Helper function to get replicates for a condition-seq_platform combo
def get_replicates_for_combo(condition, seq_platform):
    """Get the list of replicates that exist for a given condition-seq_platform combo."""
    key = (condition, seq_platform)
    return CONDITION_PLATFORM_REPS.get(key, [])

# Main rule that defines the final output
rule all:
    input:
        # Use actual sample IDs for rRNA analysis.
        # Note: {gene}_aligned.bam is intentionally NOT requested here - it's
        # a temp() intermediate (see align_gene_reads) that Snakemake deletes
        # once calculate_coverage and extract_16s_sequences have both used
        # it. coverage.tsv and blast.summary below still pull it through the
        # DAG transitively.
        expand("results/rRNA_analysis/coverage/{sample_id}/{gene}_coverage.tsv",
               sample_id=SAMPLE_IDS,
               gene=config.get("target_gene", ["GQX67_05945"])),
        expand("results/rRNA_analysis/blast/{sample_id}/{gene}.blast.summary",
               sample_id=SAMPLE_IDS,
               gene=config.get("target_gene", ["GQX67_05945"])),
        expand("results/rRNA_analysis/plots/coverage_{condition}_{seq_platform}",
               zip,
               condition=[c for c, p in CONDITION_PLATFORM_COMBOS],
               seq_platform=[p for c, p in CONDITION_PLATFORM_COMBOS]),
        expand("results/rRNA_analysis/plots/blast_{condition}_{seq_platform}",
               zip,
               condition=[c for c, p in CONDITION_PLATFORM_COMBOS],
               seq_platform=[p for c, p in CONDITION_PLATFORM_COMBOS]),
        # Integration
        "results/integrated/integrated.h5ad",

        # Barcode/UMI QC across 10x samples (post barcode-correction)
        "results/10x/inspect_summary.tsv",

        # wMel gene capture validation, run on the raw (unfiltered) h5ad
        expand("results/qc/wMel_gene_capture/{sample_id}_wMel_gene_capture.txt",
               sample_id=SAMPLE_IDS_10X),

        # Wolbachia (GQ) vs Dmel (FB) gene-group count stats, run on the
        # filtered h5ad (needs adata.raw, which filter_h5ad populates)
        expand("results/qc/gene_group_stats/{sample_id}/group_summary.csv",
               sample_id=SAMPLE_IDS_10X),

        # # Gene program and pathway analysis
        # "results/nmf_programs/.done",
        # "results/nmf_continuous_var/.done",
        # "results/nmf_categorical_var/.done",
        # "results/nmf_annotate_programs/program_cellcycle_overlap.csv",


# Process 10X samples with kallisto bustools
rule map_10x:
    input:
        reads1 = lambda wildcards: get_fastq_files(wildcards.sample_id)[0],
        reads2 = lambda wildcards: get_fastq_files(wildcards.sample_id)[1],
    output:
        h5ad = "results/h5ad_results/{sample_id}.h5ad",
        bus = "results/10x/{sample_id}/output.unfiltered.bus",
        ec = "results/10x/{sample_id}/matrix.ec",
        transcripts = "results/10x/{sample_id}/transcripts.txt"
    params:
        sample_id = "{sample_id}",
        outdir = "results/10x/{sample_id}",
        kallisto_index = config["kallisto_index"],
        transcripts_to_genes = config["transcripts_to_genes"],
        kb_reads = lambda wildcards: " ".join(get_kb_reads(wildcards.sample_id))
    wildcard_constraints:
        sample_id = ".*_10x"  # Only match samples ending with _10x
    log:
        "logs/10x/{sample_id}.log"
    threads:
        config["pseudoalign_threads"]
    resources:
        slurm_partition = config["pseudoalign_partition"],
        mem_mb = config["pseudoalign_mem"],
        runtime = hms_to_minutes(config["pseudoalign_time"])
    shell:
        """
        exec > {log} 2>&1
        echo "Starting 10x processing for {params.sample_id}"
        echo "Input files (R1 R2 pairs, in order): {params.kb_reads}"
        echo "Output directory: {params.outdir}"

        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate kallisto_bustools

        kb count \
            --kallisto /private/home/jomojaco/kallisto/build/src/kallisto \
            -i {params.kallisto_index} \
            -g {params.transcripts_to_genes} \
            --mm \
            --keep-tmp \
            -x 0,0,16:0,16,28:1,0,0 \
            -o {params.outdir} \
            -t {threads} \
            --h5ad \
            {params.kb_reads}

        echo "Moving h5ad file to final location"
        mv {params.outdir}/counts_unfiltered/adata.h5ad {output.h5ad}
        echo "10x processing complete for {params.sample_id}"
        """

##################################################################
# Barcode/UMI QC rules
##################################################################
# Rule: Inspect the corrected + sorted BUS file for each 10x sample.
# kb count writes the final corrected+resorted bus directly to
# {outdir}/output.unfiltered.bus (NOT into tmp/ - the tmp dir only holds
# the intermediate output.s.bus / output.s.c.bus files along the way).
# This gives post-correction barcode/UMI/on-list stats, as opposed to
# inspect.json (which kb count generates automatically on the
# pre-correction sorted bus file).
rule inspect_10x_corrected:
    input:
        bus = "results/10x/{sample_id}/output.unfiltered.bus"
    output:
        json = "results/10x/{sample_id}/inspect_corrected.json"
    params:
        whitelist = "results/10x/{sample_id}/whitelist.txt"
    wildcard_constraints:
        sample_id = ".*_10x"
    log:
        "logs/inspect_10x/{sample_id}.log"
    threads: 1
    resources:
        slurm_partition = "medium",
        mem_mb = 8000,
        runtime = hms_to_minutes("30:00")
    shell:
        """
        exec > {log} 2>&1
        echo "Inspecting corrected BUS file for {wildcards.sample_id}"

        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate kallisto_bustools

        bustools inspect -w {params.whitelist} -o {output.json} {input.bus}

        echo "Inspect complete for {wildcards.sample_id}"
        """

# Rule: Summarize inspect_corrected.json across all 10x samples into one table
rule summarize_inspect_10x:
    input:
        jsons = expand("results/10x/{sample_id}/inspect_corrected.json", sample_id=SAMPLE_IDS_10X)
    output:
        summary = "results/10x/inspect_summary.tsv"
    log:
        "logs/inspect_10x/summarize.log"
    run:
        import json

        rows = []
        for sample_id, jf in zip(SAMPLE_IDS_10X, input.jsons):
            with open(jf) as fh:
                d = json.load(fh)
            rows.append({
                "sample_id": sample_id,
                "numBarcodes": d.get("numBarcodes"),
                "numReads": d.get("numReads"),
                "meanReadsPerBarcode": d.get("meanReadsPerBarcode"),
                "medianReadsPerBarcode": d.get("medianReadsPerBarcode"),
                "numBarcodeUMIs": d.get("numBarcodeUMIs"),
                "meanUMIsPerBarcode": d.get("meanUMIsPerBarcode"),
                "medianUMIsPerBarcode": d.get("medianUMIsPerBarcode"),
                "percentageBarcodesOnOnlist": d.get("percentageBarcodesOnOnlist"),
                "percentageReadsOnOnlist": d.get("percentageReadsOnOnlist"),
            })

        df = pd.DataFrame(rows)
        df.to_csv(output.summary, sep="\t", index=False)
        print(df.to_string(index=False))

# Rule: Validate wMel (Wolbachia) vs Dmel gene capture on the RAW
# (unfiltered) h5ad, before filter_h5ad has a chance to gzip it.
rule validate_wMel_gene_capture:
    input:
        h5ad = "results/h5ad_results/{sample_id}.h5ad"
    output:
        report = "results/qc/wMel_gene_capture/{sample_id}_wMel_gene_capture.txt",
        csv = "results/qc/wMel_gene_capture/{sample_id}_wMel_gene_capture_wMel_genes.csv"
    params:
        script = config.get(
            "validate_wMel_script",
            "/private/groups/russelllab/jodie/scRNAseq/Jacobs_et_al_2027_dual_scRNAseq/snakemake_scripts/quality_control/validate_wMel_gene_capture.py"
        )
    wildcard_constraints:
        sample_id = ".*_10x"
    log:
        "logs/validate_wMel/{sample_id}.log"
    threads: 1
    resources:
        slurm_partition = "medium",
        mem_mb = 8000,
        runtime = hms_to_minutes("30:00")
    shell:
        """
        exec > {log} 2>&1
        echo "Validating wMel gene capture for {wildcards.sample_id}"

        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}

        python {params.script} \
            --input {input.h5ad} \
            --output {output.report}

        echo "wMel gene capture validation complete for {wildcards.sample_id}"
        """

# Rule: Gene-group (Wolbachia GQ vs Dmel FB) count stats, run on the
# filtered h5ad (post doublet-removal/QC-filter, pre-normalization) since
# that's the object with adata.raw populated with raw counts.
rule gene_group_stats:
    input:
        h5ad = "results/filtered_h5ad/{sample_id}.h5ad"
    output:
        gene_stats = "results/qc/gene_group_stats/{sample_id}/gene_stats.csv",
        group_summary = "results/qc/gene_group_stats/{sample_id}/group_summary.csv",
        cell_stats = "results/qc/gene_group_stats/{sample_id}/cell_stats.csv",
        qc_plot_png = "results/qc/gene_group_stats/{sample_id}/{sample_id}_gene_group_qc.png",
        qc_plot_pdf = "results/qc/gene_group_stats/{sample_id}/{sample_id}_gene_group_qc.pdf"
    params:
        script = config.get(
            "gene_group_stats_script",
            "/private/groups/russelllab/jodie/scRNAseq/Jacobs_et_al_2027_dual_scRNAseq/snakemake_scripts/quality_control/gene_group_stats.py"
        ),
        outdir = "results/qc/gene_group_stats/{sample_id}"
    wildcard_constraints:
        sample_id = ".*_10x"
    log:
        "logs/gene_group_stats/{sample_id}.log"
    threads: 1
    resources:
        slurm_partition = "medium",
        mem_mb = lambda wildcards, attempt: 16000 * attempt,  # was 8000; OOM'd on Mei-P26_wMel-1_10x, scales up on retry
        runtime = hms_to_minutes("30:00")
    shell:
        """
        exec > {log} 2>&1
        echo "Computing Wolbachia/Dmel gene-group stats for {wildcards.sample_id}"

        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}

        python {params.script} \
            --input {input.h5ad} \
            --outdir {params.outdir} \
            --prefix Wolbachia:GQ \
            --prefix Drosophila:FB

        echo "Gene-group stats complete for {wildcards.sample_id}"
        """

# Filter h5ad output and output qc:
rule filter_h5ad:
    input: "results/h5ad_results/{sample_id}.h5ad"
    output:
        filtered_h5ad = "results/filtered_h5ad/{sample_id}.h5ad"
    params:
        script = config["filter_script"]
    log:
        "logs/filter/{sample_id}.log"
    threads:
        config["filter_threads"]
    resources:
        slurm_partition = config["filter_partition"],
        mem_mb = config["filter_mem"],
        runtime = hms_to_minutes(config["filter_time"])
    shell:
        """
        exec > {log} 2>&1
        echo "Starting filtering for {wildcards.sample_id}"
        echo "Input file: {input}"

        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}

        python {params.script} \
            --input {input} \
            --output {output.filtered_h5ad}

        echo "Compressing original h5ad file"
        gzip {input}
        echo "Filtering complete for {wildcards.sample_id}"
        """

# Cell cycle Annotation
rule annotate_cell_cycle: # This needs the cyclum conda environment
    input:
        h5ad = "results/filtered_h5ad/{sample_id}.h5ad" # Output of filtered script
    output:
        annotated_h5ad = "results/annotated_h5ad/{sample_id}.h5ad"
    params:
        script = config["cell_cycle_script"]
    log:
        "logs/annotate/{sample_id}.log"
    threads:
        config["cell_cycle_threads"]
    resources:
        slurm_partition = config["cell_cycle_partition"],
        mem_mb = config["cell_cycle_mem"],
        runtime = hms_to_minutes(config["cell_cycle_time"])
    shell:
        """
        exec > {log} 2>&1
        echo "Starting cell cycle annotation for {wildcards.sample_id}"
        echo "Input file: {input.h5ad}"

        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {CYCLUM_ENV}

        python {params.script} \
            --input {input.h5ad} \
            --output {output.annotated_h5ad}

        echo "Compressing filtered h5ad file"
        gzip {input.h5ad}
        echo "Cell cycle annotation complete for {wildcards.sample_id}"
        """

##################################################################
# rRNA analysis rules
##################################################################
# Rule 1: Align extracted gene reads to rRNA reference with BWA
rule align_gene_reads:
    input:
        # r1 = "results/gene_extracted/{sample_id}/{gene}_R1.fastq.gz",
        r2 = lambda wildcards: get_fastq_files(wildcards.sample_id)[1],
        ref = config["ref_fasta"],
        regions = config["rRNA_regions"]
    output:
        bam = temp("results/rRNA_analysis/alignment/{sample_id}/{gene}_aligned.bam"),
        bai = temp("results/rRNA_analysis/alignment/{sample_id}/{gene}_aligned.bam.bai")
    log:
        "logs/align_gene/{sample_id}_{gene}.log"
    threads: 8
    resources:
        slurm_partition = "medium",
        mem_mb = 1000000,
        runtime = hms_to_minutes("12:00:00")
    shell:
        """
        # exec > {log} 2>&1
        set -euo pipefail
        echo "Starting BWA alignment for {wildcards.sample_id} - {wildcards.gene}"

        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate sra-tools

        # Make the directory
        mkdir -p results/rRNA_analysis/alignment/{wildcards.sample_id}

        # {input.r2} is one file per lane/run for this sample (see
        # get_fastq_files in the Snakefile); concatenating gzipped fastqs
        # is valid (multi-member gzip streams decompress fine), so this
        # merges however many R2 files the sample has into one before
        # alignment - unchanged behavior for single-lane samples. Streamed
        # through process substitution instead of written to a temp file
        # first, so large samples don't need a second full-size copy of R2
        # sitting on disk during alignment.
        #
        # {input.r2} covers the whole library, but {input.ref} is only the
        # small rRNA/wMel reference, so the large majority of reads here
        # won't map. The previous version kept every read (mapped and
        # unmapped) through sort, which meant runtime/memory/disk scaled
        # with the full library size instead of the much smaller mapped
        # fraction - this is almost certainly why large samples were
        # struggling. `-F 4` drops unmapped reads right after alignment,
        # before the expensive sort step.
        bwa mem -t {threads} {input.ref} <(cat {input.r2}) | \
            samtools view -Sb -F 4 | \
            samtools sort -@ {threads} -o results/rRNA_analysis/alignment/{wildcards.sample_id}/all_aligned.bam

        samtools index results/rRNA_analysis/alignment/{wildcards.sample_id}/all_aligned.bam

        # Extract reads for this specific gene using the full chromosome name
        CHROM=$(samtools idxstats results/rRNA_analysis/alignment/{wildcards.sample_id}/all_aligned.bam | grep "^{wildcards.gene}::" | cut -f1)

        if [ -z "$CHROM" ]; then
            echo "ERROR: Gene {wildcards.gene} not found"
            samtools idxstats results/rRNA_analysis/alignment/{wildcards.sample_id}/all_aligned.bam
            exit 1
        fi

        samtools view -b results/rRNA_analysis/alignment/{wildcards.sample_id}/all_aligned.bam "$CHROM" -o {output.bam}
        samtools index {output.bam}

        # Clean up the whole-sample intermediate now that the gene-specific
        # slice has been pulled out of it - even filtered to mapped reads
        # only, it's not a declared output and nothing else needs it.
        # NOTE: this file is shared by sample_id only, not by gene, so if
        # target_gene is ever a list of more than one value, deleting it
        # here will break any other in-flight/queued gene job for this same
        # sample that still expects to read it - fine for the current
        # single-gene config, not safe if that changes.
        rm -f results/rRNA_analysis/alignment/{wildcards.sample_id}/all_aligned.bam
        rm -f results/rRNA_analysis/alignment/{wildcards.sample_id}/all_aligned.bam.bai

        echo "Alignment complete for {wildcards.sample_id} - {wildcards.gene}"
        echo "Filtered to $(samtools view -c {output.bam}) reads in target regions"
        """

# Rule 2: Calculate coverage depth
rule calculate_coverage:
    input:
        bam = "results/rRNA_analysis/alignment/{sample_id}/{gene}_aligned.bam",
        bai = "results/rRNA_analysis/alignment/{sample_id}/{gene}_aligned.bam.bai"
    output:
        cov = "results/rRNA_analysis/coverage/{sample_id}/{gene}_coverage.tsv"
    log:
        "logs/coverage/{sample_id}_{gene}.log"
    threads: 1
    resources:
        slurm_partition = "medium",
        mem_mb = 4000,
        runtime = hms_to_minutes("30:00")
    shell:
        """
        exec > {log} 2>&1
        echo "Calculating coverage for {wildcards.sample_id} - {wildcards.gene}"

        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate sra-tools

        samtools depth {input.bam} > {output.cov}

        echo "Coverage calculation complete"
        """

# Rule 3: Extract 16S sequences from aligned BAM
rule extract_16s_sequences:
    input:
        bam = "results/rRNA_analysis/alignment/{sample_id}/{gene}_aligned.bam",
        bai = "results/rRNA_analysis/alignment/{sample_id}/{gene}_aligned.bam.bai"
    output:
        fasta = "results/rRNA_analysis/extracted_16S/{sample_id}/{gene}_16S.fasta"
    log:
        "logs/extract_16s/{sample_id}_{gene}.log"
    threads: 1
    resources:
        slurm_partition = "medium",
        mem_mb = 4000,
        runtime = hms_to_minutes("30:00")
    shell:
        """
        exec > {log} 2>&1
        echo "Extracting 16S sequences for {wildcards.sample_id} - {wildcards.gene}"

        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate sra-tools

        # Extract reads mapping to 16S region
        samtools fasta {input.bam} > {output.fasta}

        echo "16S extraction complete"
        """

# Rule 4: BLAST 16S sequences against database
rule blast_16s:
    input:
        fasta = "results/rRNA_analysis/extracted_16S/{sample_id}/{gene}_16S.fasta"
    output:
        blast = "results/rRNA_analysis/blast/{sample_id}/{gene}.blast"
    params:
        db = config["blast_db"]
    log:
        "logs/blast/{sample_id}_{gene}.log"
    threads: 16
    resources:
        slurm_partition = "medium",
        mem_mb = 32000,
        runtime = hms_to_minutes("4:00:00")
    shell:
        """
        exec > {log} 2>&1
        echo "Running BLAST for {wildcards.sample_id} - {wildcards.gene}"

        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate sra-tools

        blastn -db {params.db} \
            -query {input.fasta} \
            -out {output.blast} \
            -num_threads {threads} \
            -outfmt 6 \
            -max_target_seqs 5 \
            -evalue 1e-5

        echo "BLAST complete"
        """

# Rule 5: Summarize BLAST results
rule summarize_blast:
    input:
        blast = "results/rRNA_analysis/blast/{sample_id}/{gene}.blast"
    output:
        summary = "results/rRNA_analysis/blast/{sample_id}/{gene}.blast.summary"
    log:
        "logs/summarize_blast/{sample_id}_{gene}.log"
    threads: 1
    resources:
        slurm_partition = "medium",
        mem_mb = 2000,
        runtime = hms_to_minutes("15:00")
    shell:
        """
        exec > {log} 2>&1
        echo "Summarizing BLAST results for {wildcards.sample_id} - {wildcards.gene}"

        # Summarize: count hits per subject, track best identity
        awk '{{count[$2]++; if($3 > best[$2]) best[$2]=$3}}
             END {{for(s in count) printf "%6d %s (%.3f%% identity)\\n", count[s], s, best[s]}}' \
            {input.blast} | \
            sort -k1,1nr -k3,3nr > {output.summary}

        echo "BLAST summarization complete"
        """

# Rule 6: Plot coverage for samples grouped by condition and seq_platform
rule plot_coverage_by_group:
    input:
        coverage_files = lambda wildcards: expand(
            "results/rRNA_analysis/coverage/{condition}-{replicate}_{seq_platform}/{gene}_coverage.tsv",
            condition=wildcards.condition,
            replicate=get_replicates_for_combo(wildcards.condition, wildcards.seq_platform),
            seq_platform=wildcards.seq_platform,
            gene=config.get("target_gene", ["GQX67_05945"])
        )
    output:
        plot_dir = directory("results/rRNA_analysis/plots/coverage_{condition}_{seq_platform}")
    params:
        script = config["plot_coverage_script"]
    log:
        "logs/plot_coverage/{condition}_{seq_platform}.log"
    threads: 1
    resources:
        slurm_partition = "medium",
        mem_mb = 8000,
        runtime = hms_to_minutes("1:00:00")
    shell:
        """
        exec > {log} 2>&1
        echo "Plotting coverage for {wildcards.condition}_{wildcards.seq_platform}"

        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate sra-tools

        mkdir -p {output.plot_dir}

        # Create file list
        echo {input.coverage_files} | tr ' ' '\\n' > {output.plot_dir}/coverage_files.txt

        python {params.script} {output.plot_dir}/coverage_files.txt \
            --output-dir {output.plot_dir}

        echo "Coverage plotting complete"
        """

# Rule 7: Plot BLAST pie charts for samples grouped by condition and seq_platform
rule plot_blast_by_group:
    input:
        blast_summaries = lambda wildcards: expand(
            "results/rRNA_analysis/blast/{condition}-{replicate}_{seq_platform}/{gene}.blast.summary",
            condition=wildcards.condition,
            replicate=get_replicates_for_combo(wildcards.condition, wildcards.seq_platform),
            seq_platform=wildcards.seq_platform,
            gene=config.get("target_gene", ["GQX67_05945"])
        )
    output:
        plot_dir = directory("results/rRNA_analysis/plots/blast_{condition}_{seq_platform}")
    params:
        script = config["plot_blast_script"]
    log:
        "logs/plot_blast/{condition}_{seq_platform}.log"
    threads: 1
    resources:
        slurm_partition = "medium",
        mem_mb = 8000,
        runtime = hms_to_minutes("1:00:00")
    shell:
        """
        exec > {log} 2>&1
        echo "Plotting BLAST results for {wildcards.condition}_{wildcards.seq_platform}"

        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate sra-tools

        mkdir -p {output.plot_dir}

        # Plot each summary file
        for summary in {input.blast_summaries}; do
            echo "Processing $(basename $summary)..."
            python {params.script} "$summary" --output-dir {output.plot_dir}
        done

        echo "BLAST plotting complete"
        """

# Optional: Extract highly represented 16S sequences
rule extract_abundant_16s:
    input:
        blast = "results/rRNA_analysis/blast/{sample_id}/{gene}.blast",
        summary = "results/rRNA_analysis/blast/{sample_id}/{gene}.blast.summary"
    output:
        ids = "results/rRNA_analysis/abundant_16S/{sample_id}/{gene}_abundant_ids.txt",
        fasta = "results/rRNA_analysis/abundant_16S/{sample_id}/{gene}_abundant.fasta"
    params:
        db = config["blast_db"],
        min_reads = config.get("min_blast_reads", 100)
    log:
        "logs/extract_abundant/{sample_id}_{gene}.log"
    threads: 1
    resources:
        slurm_partition = "medium",
        mem_mb = 4000,
        runtime = hms_to_minutes("30:00")
    shell:
        """
        exec > {log} 2>&1
        echo "Extracting abundant 16S sequences for {wildcards.sample_id} - {wildcards.gene}"

        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate sra-tools

        # Extract IDs with at least min_reads reads
        awk -v min={params.min_reads} '$1 >= min {{print $2}}' {input.summary} > {output.ids}

        # Extract sequences from BLAST database
        if [ -s {output.ids} ]; then
            blastdbcmd -db {params.db} \
                -entry_batch {output.ids} \
                -out {output.fasta}
        else
            echo "No abundant sequences found" > {output.fasta}
        fi

        echo "Abundant sequence extraction complete"
        """

##################################################################
# Integration via frozen-atlas projection (now the default path)
##################################################################
# Produces results/integrated/integrated.h5ad by projecting every filtered
# sample onto the SAME frozen ovary reference atlas embedding (Fly Cell
# Atlas / BioHub ovary object), via integrate_via_atlas_projection.py --
# same mechanism as Lum et al. 2027's integrate_via_atlas_projection.py
# (scanpy.tl.ingest reference projection), replacing this rule's old
# BBKNN joint-clustering path (integrate.py, kept in
# snakemake_scripts/analysis/ if you want to compare or revert -- point
# integrate_atlas_script at it and add back its own CLI flags to switch).
# See integrate_via_atlas_projection.py's module docstring for what this
# object is/isn't good for (atlas_<label> cell-type identity and
# wolbachia_titer-by-cell-type composition, not fine within-cell-type
# expression shifts along the titer/infection-status axis -- keep using
# integrate.py's own Harmony output, run on your cells alone, for that).
rule integrate:
    input:
        files              = expand("results/filtered_h5ad/{sample_id}.h5ad", sample_id=SAMPLE_IDS),
        atlas              = config.get("ovary_atlas",
                                  "ovary_atlas/atlas/h5ad_unprocessed/s_fca_biohub_ovary_10x.h5ad"),
        flybase_annotation = config["flybase_annotation"],
    output:
        integrated = "results/integrated/integrated.h5ad"
    params:
        script   = config.get("integrate_atlas_script",
                       "snakemake_scripts/analysis/integrate_via_atlas_projection.py"),
        fig_dir  = "results/integrated/figures",
        k        = config.get("atlas_k", 30),
        n_pcs    = config.get("atlas_n_pcs", 30),
        label_cols_flag = (
            "--label_cols " + " ".join(config["atlas_label_cols"])
            if config.get("atlas_label_cols") else ""
        ),
        subsample_flag = (
            f"--subsample_ref {config['atlas_subsample_ref']}"
            if config.get("atlas_subsample_ref") else ""
        ),
    log:
        "logs/integrate/integrate.log"
    threads:
        config.get("integrate_threads", 16)
    resources:
        slurm_partition = config.get("integrate_partition", "medium"),
        mem_mb          = config.get("integrate_mem", 128000),
        runtime      = hms_to_minutes(config.get("integrate_time", "8:00:00"))
    shell:
        """
        exec > {log} 2>&1
        echo "Starting atlas-projected integration"

        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}

        python {params.script} \
            --atlas {input.atlas} \
            --query {input.files} \
            --out_path {output.integrated} \
            --fig_dir {params.fig_dir} \
            --flybase_annotation {input.flybase_annotation} \
            --k {params.k} \
            --n_pcs {params.n_pcs} \
            {params.label_cols_flag} \
            {params.subsample_flag}

        echo "Atlas-projected integration complete"
        """

# ##################################################################
# # Gene program and pathway analysis rules
# ##################################################################
# rule nmf_programs:
#     input:
#         h5ad = "results/integrated/integrated_with_cellcycle.h5ad"
#     output:
#         adata_with_programs = "results/nmf_programs/adata_with_programs.h5ad",
#         summary = "results/nmf_programs/SUMMARY.txt",
#         flag = touch("results/nmf_programs/.done")
#     params:
#         script = config.get("nmf_script"),
#         output_dir = "results/nmf_programs",
#         n_programs = config.get("n_programs", 15),
#         n_top_genes = config.get("n_top_genes", 2000),
#         organism = config.get("organism", "Fly"),
#         gene_id_type = config.get("gene_id_type", "flybase"),
#         flybase_annotation = config.get("flybase_annotation",
#             "reference/fbgn_annotation_ID_fb_2025_04.tsv.gz")
#     log:
#         "logs/nmf/nmf_programs.log"
#     threads:
#         config.get("nmf_threads", 8)
#     resources:
#         slurm_partition = config.get("nmf_partition", "medium"),
#         mem_mb          = config.get("nmf_mem", 64000),
#         runtime      = hms_to_minutes(config.get("nmf_time", "4:00:00"))
#     shell:
#         """
#         exec > {log} 2>&1
#         echo "Starting NMF program discovery"

#         source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
#         conda activate {SCANPY_ENV}

#         python {params.script} \
#             --input {input.h5ad} \
#             --output_dir {params.output_dir} \
#             --n_programs {params.n_programs} \
#             --n_top_genes {params.n_top_genes} \
#             --organism {params.organism} \
#             --gene_id_type {params.gene_id_type} \
#             --flybase_annotation {params.flybase_annotation}

#         echo "NMF program discovery complete"
#         """

# rule nmf_continuous_var:
#     input:
#         adata_with_programs = "results/nmf_programs/adata_with_programs.h5ad"
#     output:
#         correlations = "results/nmf_continuous_var/program_correlations.csv",
#         summary = "results/nmf_continuous_var/SUMMARY.txt",
#         flag = touch("results/nmf_continuous_var/.done")
#     params:
#         script = config.get("nmf_continuous_script"),
#         output_dir = "results/nmf_continuous_var",
#         continuous_var = config.get("continuous_var", None),  # Auto-detect if None
#         flybase_annotation = config.get("flybase_annotation",
#             "reference/fbgn_annotation_ID_fb_2025_04.tsv.gz")
#     log:
#         "logs/nmf/nmf_continuous_var.log"
#     threads:
#         config.get("nmf_continuous_threads", 4)
#     resources:
#         slurm_partition = config.get("nmf_continuous_partition", "medium"),
#         mem_mb = config.get("nmf_continuous_mem", 32000),
#         runtime = hms_to_minutes(config.get("nmf_continuous_time", "2:00:00"))
#     shell:
#         """
#         exec > {log} 2>&1
#         echo "Starting NMF continuous variable analysis"

#         source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
#         conda activate {SCANPY_ENV}

#         python {params.script} \
#             --input {input.adata_with_programs} \
#             --output_dir {params.output_dir} \
#             --continuous_var {params.continuous_var} \
#             --flybase_annotation {params.flybase_annotation}

#         echo "NMF continuous variable analysis complete"
#         """

# rule nmf_categorical_var:
#     input:
#         adata_with_programs = "results/nmf_programs/adata_with_programs.h5ad"
#     output:
#         comparison = "results/nmf_categorical_var/program_comparison.csv",
#         flag = touch("results/nmf_categorical_var/.done")
#     params:
#         script = config.get("nmf_categorical_script"),
#         output_dir = "results/nmf_categorical_var",
#         categorical_var = config.get("categorical_var", None)  # Auto-detect if None
#     log:
#         "logs/nmf/nmf_categorical_var.log"
#     threads:
#         config.get("nmf_categorical_threads", 4)
#     resources:
#         slurm_partition = config.get("nmf_categorical_partition", "medium"),
#         mem_mb = config.get("nmf_categorical_mem", 32000),
#         runtime = hms_to_minutes(config.get("nmf_categorical_time", "2:00:00"))
#     shell:
#         """
#         exec > {log} 2>&1
#         echo "Starting NMF categorical variable analysis"

#         source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
#         conda activate {SCANPY_ENV}

#         python {params.script} \
#             --input {input.adata_with_programs} \
#             --output_dir {params.output_dir} \
#             --categorical_var {params.categorical_var}

#         echo "NMF categorical variable analysis complete"
#         """

# rule nmf_annotate_programs:
#     input:
#         adata_with_programs = "results/nmf_programs/adata_with_programs.h5ad",
#         nmf_done            = "results/nmf_programs/.done",
#         mapping             = config["transcripts_to_genes"]
#     output:
#         output  = "results/nmf_annotate_programs/program_cellcycle_overlap.csv",
#     params:
#         script       = config.get("nmf_annotate_script"),
#         output_dir   = "results/nmf_annotate_programs",
#         program_dir  = "results/nmf_programs",
#         sample_name  = config.get("nmf_annotate_sample", "wolbachia_infection"),
#         titer_var    = config.get("continuous_var", "wolbachia_titer"),
#         cc_s_var     = config.get("cc_s_var",      "S_score"),
#         cc_g2m_var   = config.get("cc_g2m_var",    "G2M_score"),
#         cc_phase_var = config.get("cc_phase_var",  "phase"),
#         top_genes    = config.get("nmf_annotate_top_genes", 200),
#         skip_gsea    = "--skip_gsea" if config.get("nmf_annotate_skip_gsea", False) else "",
#         skip_fly     = "--skip_flyenrichr" if config.get("nmf_annotate_skip_flyenrichr", False) else ""
#     log:
#         "logs/nmf/nmf_annotate_programs.log"
#     threads:
#         config.get("nmf_annotate_threads", 8)
#     resources:
#         slurm_partition = config.get("nmf_annotate_partition", "medium"),
#         mem_mb          = config.get("nmf_annotate_mem",       32000),
#         runtime      = hms_to_minutes(config.get("nmf_annotate_time",      "4:00:00"))
#     shell:
#         """
#         exec > {log} 2>&1
#         echo "Starting NMF program annotation"

#         source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
#         conda activate {SCANPY_ENV}

#         python {params.script} \
#             --input        {input.adata_with_programs} \
#             --program_dir  {params.program_dir} \
#             --output_dir   {params.output_dir} \
#             --mapping      {input.mapping} \
#             --titer_var    {params.titer_var} \
#             --cc_s_var     {params.cc_s_var} \
#             --cc_g2m_var   {params.cc_g2m_var} \
#             --cc_phase_var {params.cc_phase_var} \
#             --top_genes    {params.top_genes} \
#             {params.skip_gsea} \
#             {params.skip_fly}

#         echo "NMF program annotation complete"
#         """
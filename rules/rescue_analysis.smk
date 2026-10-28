# mei-P26 fertility rescue analysis (pseudobulk by atlas cell type)
#
# Run everything:   snakemake rescue_all --executor slurm ...
# Inputs come from rule integrate (labels) and filter_h5ad (raw counts).
# Scripts: snakemake_scripts/analysis/mei_P26_rescue/
#   01_build_pseudobulk.py   labeled raw counts + sample x tissue pseudobulk
#   02_wolbachia_titer.py    titer by tissue, Mei-P26 vs OreR
#   03_edger_de.R            host genotype x infection DE; Wolbachia DE
#   04_cell_composition.py   cell-type proportions ~ genotype * infection
#   05_rescue_summary.py     rescue classification + plots

RESCUE_DIR     = "results/rescue"
RESCUE_SCRIPTS = "snakemake_scripts/analysis/mei_P26_rescue"
EDGER_ENV      = config.get("edger_env", "/private/groups/russelllab/jodie/miniforge3/envs/edger")
RESCUE_LABEL   = config.get("rescue_label_col", "cell_type_annotation_broad")
RESCUE_MINCONF = config.get("rescue_min_confidence", 0.5)
RESCUE_MINCELL = config.get("rescue_min_cells", 20)

rule rescue_all:
    input:
        f"{RESCUE_DIR}/titer/titer_genotype_by_tissue.tsv",
        f"{RESCUE_DIR}/composition/composition_tests.tsv",
        f"{RESCUE_DIR}/summary/rescue_genes.tsv",

rule rescue_pseudobulk:
    input:
        integrated = "results/integrated/integrated.h5ad",
        filtered   = expand("results/filtered_h5ad/{sample_id}.h5ad", sample_id=SAMPLE_IDS),
    output:
        counts = f"{RESCUE_DIR}/pseudobulk/counts_labeled.h5ad",
        host   = f"{RESCUE_DIR}/pseudobulk/pseudobulk_host_counts.tsv.gz",
        wmel   = f"{RESCUE_DIR}/pseudobulk/pseudobulk_wolbachia_counts.tsv.gz",
        meta   = f"{RESCUE_DIR}/pseudobulk/pseudobulk_meta.tsv",
    params:
        outdir = f"{RESCUE_DIR}/pseudobulk",
        qc_flag = "--drop_qc_outliers" if config.get("rescue_drop_qc_outliers", False) else "",
    log: "logs/rescue/pseudobulk.log"
    resources:
        slurm_partition = "medium", mem_mb = 64000, runtime = hms_to_minutes("2:00:00")
    shell:
        """
        exec > {log} 2>&1
        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}
        {SCANPY_ENV}/bin/python {RESCUE_SCRIPTS}/01_build_pseudobulk.py \
            --integrated {input.integrated} --filtered {input.filtered} \
            --outdir {params.outdir} --label_col {RESCUE_LABEL} \
            --min_confidence {RESCUE_MINCONF} {params.qc_flag}
        """

rule rescue_titer:
    input:
        counts = rules.rescue_pseudobulk.output.counts,
        meta   = rules.rescue_pseudobulk.output.meta,
    output:
        f"{RESCUE_DIR}/titer/titer_genotype_by_tissue.tsv",
        f"{RESCUE_DIR}/titer/titer_by_sample_tissue.tsv",
    params:
        indir = f"{RESCUE_DIR}/pseudobulk", outdir = f"{RESCUE_DIR}/titer",
    log: "logs/rescue/titer.log"
    resources:
        slurm_partition = "medium", mem_mb = 32000, runtime = hms_to_minutes("1:00:00")
    shell:
        """
        exec > {log} 2>&1
        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}
        {SCANPY_ENV}/bin/python {RESCUE_SCRIPTS}/02_wolbachia_titer.py \
            --pseudobulk_dir {params.indir} --outdir {params.outdir} --min_cells {RESCUE_MINCELL}
        """

rule rescue_de:
    input:
        host = rules.rescue_pseudobulk.output.host,
        wmel = rules.rescue_pseudobulk.output.wmel,
        meta = rules.rescue_pseudobulk.output.meta,
        annotation = config["flybase_annotation"],
    output:
        host = f"{RESCUE_DIR}/de/host_de_all.tsv.gz",
        host_summary = f"{RESCUE_DIR}/de/host_de_summary.tsv",
        wmel = f"{RESCUE_DIR}/de/wolbachia_de_all.tsv.gz",
    params:
        indir = f"{RESCUE_DIR}/pseudobulk", outdir = f"{RESCUE_DIR}/de",
        block = "--block_replicate" if config.get("rescue_block_replicate", False) else "",
        min_wmel = config.get("rescue_min_wmel_umis", 2000),
    log: "logs/rescue/de.log"
    resources:
        slurm_partition = "medium", mem_mb = 16000, runtime = hms_to_minutes("1:00:00")
    shell:
        """
        exec > {log} 2>&1
        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {EDGER_ENV}
        {EDGER_ENV}/bin/Rscript {RESCUE_SCRIPTS}/03_edger_de.R \
            --pseudobulk_dir {params.indir} --outdir {params.outdir} \
            --min_cells {RESCUE_MINCELL} --min_wmel_umis {params.min_wmel} \
            --annotation {input.annotation} {params.block}
        """

rule rescue_composition:
    input:
        meta = rules.rescue_pseudobulk.output.meta,
    output:
        f"{RESCUE_DIR}/composition/composition_tests.tsv",
    params:
        indir = f"{RESCUE_DIR}/pseudobulk", outdir = f"{RESCUE_DIR}/composition",
    log: "logs/rescue/composition.log"
    resources:
        slurm_partition = "medium", mem_mb = 8000, runtime = hms_to_minutes("0:30:00")
    shell:
        """
        exec > {log} 2>&1
        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}
        {SCANPY_ENV}/bin/python {RESCUE_SCRIPTS}/04_cell_composition.py \
            --pseudobulk_dir {params.indir} --outdir {params.outdir}
        """

rule rescue_summary:
    input:
        host = rules.rescue_de.output.host,
        host_summary = rules.rescue_de.output.host_summary,
        wmel = rules.rescue_de.output.wmel,
    output:
        f"{RESCUE_DIR}/summary/rescue_genes.tsv",
    params:
        indir = f"{RESCUE_DIR}/de", outdir = f"{RESCUE_DIR}/summary",
    log: "logs/rescue/summary.log"
    resources:
        slurm_partition = "medium", mem_mb = 8000, runtime = hms_to_minutes("0:30:00")
    shell:
        """
        exec > {log} 2>&1
        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}
        {SCANPY_ENV}/bin/python {RESCUE_SCRIPTS}/05_rescue_summary.py --de_dir {params.indir} --outdir {params.outdir}
        """

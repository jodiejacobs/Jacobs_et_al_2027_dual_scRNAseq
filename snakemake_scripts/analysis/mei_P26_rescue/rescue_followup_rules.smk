
# ---------------------------------------------------------------------------
# Rescue follow-up (06-09), added 2026-11-03
#   06_label_qc.py            marker check of atlas labels + subclusters of
#                             artefact / unannotated / adult tracheal cell
#   07_germline.py            marker-defined germline: composition, pseudotime,
#                             germline pseudobulks (01 format)
#   03_edger_de.R (reused)    host + Wolbachia DE on germline pseudobulks
#   08_rescue_null.R          permutation null for rescue fraction +
#                             interaction-supported rescue genes
#   09_wolbachia_followup.R   replicate / normalisation / pooled checks of
#                             Wolbachia DE; --gff adds gene products
# ---------------------------------------------------------------------------
WMEL_GFF = config.get("wmel_gff", "")

rule rescue_followup_all:
    input:
        f"{RESCUE_DIR}/label_qc/lineage_by_label.tsv",
        f"{RESCUE_DIR}/germline/germline_composition_tests.tsv",
        f"{RESCUE_DIR}/germline/de/host_de_all.tsv.gz",
        f"{RESCUE_DIR}/germline/rescue_null/rescue_null_summary.tsv",
        f"{RESCUE_DIR}/rescue_null/rescue_null_summary.tsv",
        f"{RESCUE_DIR}/wolbachia_followup/wolbachia_followup_summary.tsv",

rule rescue_label_qc:
    input:
        counts = rules.rescue_pseudobulk.output.counts,
        annotation = config["flybase_annotation"],
    output:
        f"{RESCUE_DIR}/label_qc/lineage_by_label.tsv",
        f"{RESCUE_DIR}/label_qc/cell_lineage.tsv.gz",
    params:
        outdir = f"{RESCUE_DIR}/label_qc",
    log: "logs/rescue/label_qc.log"
    resources:
        slurm_partition = "medium", mem_mb = 128000, runtime = hms_to_minutes("3:00:00")
    shell:
        """
        exec > {log} 2>&1
        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}
        {SCANPY_ENV}/bin/python {RESCUE_SCRIPTS}/06_label_qc.py \
            --counts {input.counts} --annotation {input.annotation} --outdir {params.outdir}
        """

rule rescue_germline:
    input:
        counts = rules.rescue_pseudobulk.output.counts,
        annotation = config["flybase_annotation"],
    output:
        tests = f"{RESCUE_DIR}/germline/germline_composition_tests.tsv",
        host  = f"{RESCUE_DIR}/germline/pseudobulk/pseudobulk_host_counts.tsv.gz",
        wmel  = f"{RESCUE_DIR}/germline/pseudobulk/pseudobulk_wolbachia_counts.tsv.gz",
        meta  = f"{RESCUE_DIR}/germline/pseudobulk/pseudobulk_meta.tsv",
    params:
        outdir = f"{RESCUE_DIR}/germline",
    log: "logs/rescue/germline.log"
    resources:
        slurm_partition = "medium", mem_mb = 128000, runtime = hms_to_minutes("2:00:00")
    shell:
        """
        exec > {log} 2>&1
        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}
        {SCANPY_ENV}/bin/python {RESCUE_SCRIPTS}/07_germline.py \
            --counts {input.counts} --annotation {input.annotation} --outdir {params.outdir}
        """

rule rescue_germline_de:
    input:
        host = rules.rescue_germline.output.host,
        wmel = rules.rescue_germline.output.wmel,
        meta = rules.rescue_germline.output.meta,
        annotation = config["flybase_annotation"],
    output:
        host = f"{RESCUE_DIR}/germline/de/host_de_all.tsv.gz",
        wmel = f"{RESCUE_DIR}/germline/de/wolbachia_de_all.tsv.gz",
    params:
        indir = f"{RESCUE_DIR}/germline/pseudobulk", outdir = f"{RESCUE_DIR}/germline/de",
        min_wmel = config.get("rescue_min_wmel_umis", 2000),
    log: "logs/rescue/germline_de.log"
    resources:
        slurm_partition = "medium", mem_mb = 32000, runtime = hms_to_minutes("1:00:00")
    shell:
        """
        exec > {log} 2>&1
        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {EDGER_ENV}
        {EDGER_ENV}/bin/Rscript {RESCUE_SCRIPTS}/03_edger_de.R \
            --pseudobulk_dir {params.indir} --outdir {params.outdir} \
            --min_cells {RESCUE_MINCELL} --min_wmel_umis {params.min_wmel} \
            --annotation {input.annotation}
        """

rule rescue_null:
    input:
        host = rules.rescue_pseudobulk.output.host,
        meta = rules.rescue_pseudobulk.output.meta,
        de   = rules.rescue_de.output.host,
    output:
        f"{RESCUE_DIR}/rescue_null/rescue_null_summary.tsv",
        f"{RESCUE_DIR}/rescue_null/rescue_supported_genes.tsv",
    params:
        indir = f"{RESCUE_DIR}/pseudobulk", dedir = f"{RESCUE_DIR}/de", outdir = f"{RESCUE_DIR}/rescue_null",
    log: "logs/rescue/rescue_null.log"
    resources:
        slurm_partition = "medium", mem_mb = 32000, runtime = hms_to_minutes("2:00:00")
    shell:
        """
        exec > {log} 2>&1
        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {EDGER_ENV}
        {EDGER_ENV}/bin/Rscript {RESCUE_SCRIPTS}/08_rescue_null.R \
            --pseudobulk_dir {params.indir} --de_dir {params.dedir} --outdir {params.outdir} \
            --min_cells {RESCUE_MINCELL}
        """

rule rescue_germline_null:
    input:
        meta = rules.rescue_germline.output.meta,
        de   = rules.rescue_germline_de.output.host,
    output:
        f"{RESCUE_DIR}/germline/rescue_null/rescue_null_summary.tsv",
    params:
        indir = f"{RESCUE_DIR}/germline/pseudobulk", dedir = f"{RESCUE_DIR}/germline/de",
        outdir = f"{RESCUE_DIR}/germline/rescue_null",
    log: "logs/rescue/germline_null.log"
    resources:
        slurm_partition = "medium", mem_mb = 32000, runtime = hms_to_minutes("1:00:00")
    shell:
        """
        exec > {log} 2>&1
        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {EDGER_ENV}
        {EDGER_ENV}/bin/Rscript {RESCUE_SCRIPTS}/08_rescue_null.R \
            --pseudobulk_dir {params.indir} --de_dir {params.dedir} --outdir {params.outdir} \
            --min_cells {RESCUE_MINCELL}
        """

rule rescue_wolbachia_followup:
    input:
        wmel = rules.rescue_pseudobulk.output.wmel,
        meta = rules.rescue_pseudobulk.output.meta,
    output:
        f"{RESCUE_DIR}/wolbachia_followup/wolbachia_followup_summary.tsv",
        f"{RESCUE_DIR}/wolbachia_followup/wolbachia_candidates.tsv",
    params:
        indir = f"{RESCUE_DIR}/pseudobulk", outdir = f"{RESCUE_DIR}/wolbachia_followup",
        min_wmel = config.get("rescue_min_wmel_umis", 2000),
        gff = f"--gff {WMEL_GFF}" if WMEL_GFF else "",
    log: "logs/rescue/wolbachia_followup.log"
    resources:
        slurm_partition = "medium", mem_mb = 16000, runtime = hms_to_minutes("0:30:00")
    shell:
        """
        exec > {log} 2>&1
        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {EDGER_ENV}
        {EDGER_ENV}/bin/Rscript {RESCUE_SCRIPTS}/09_wolbachia_followup.R \
            --pseudobulk_dir {params.indir} --outdir {params.outdir} \
            --min_cells {RESCUE_MINCELL} --min_wmel_umis {params.min_wmel} {params.gff}
        """

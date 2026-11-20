# ---------------------------------------------------------------------------
# Hypothetical wMel proteins that are up in MeiP26 ovaries, added 2026-11-20
#   20_hypo_select.py   significant-up hypothetical loci (+ robustness flags) -> hits.tsv/.faa
#   InterProScan, eggNOG-mapper   sequence annotation, whole wMel proteome (cheap, reusable)
#   ColabFold + Foldseek          structure and structural homology, hits only
#   21_hypo_context.py  neighbouring genes, mobile-element / co-up neighbours
#   22_hypo_merge.py    one evidence table (also pulls rescue_secretion's class and euk_domain)
#
# snakemake --executor slurm --default-resources slurm_partition=medium -j 4 hypo_all
#
# Tools that need data or GPUs are set in config.yaml (hypo_* keys); a rule with
# an unset path stops with a message saying which key to fill.
# ---------------------------------------------------------------------------
HYPO_DIR = f"{RESCUE_DIR}/hypothetical"
HYPO_ENV = config.get("hypo_env", "/private/groups/russelllab/jodie/miniforge3/envs/hypo")
FOLDSEEK_DBS = config.get("hypo_foldseek_dbs") or {}
CONDA_SH = "source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh"

rule hypo_all:
    input:
        f"{HYPO_DIR}/hypothetical_evidence.tsv",

# wMel proteome keyed by locus tag (headers cut to the locus tag so the annotation
# tools return it unchanged).
rule hypo_proteome_faa:
    output: f"{HYPO_DIR}/wmel_proteome.faa"
    params: gtf = WMEL_GFF, faa = WMEL_FAA, de_dir = RESCUE_DIR, outdir = HYPO_DIR
    log: "logs/rescue/hypo_proteome_faa.log"
    resources: slurm_partition = "medium", mem_mb = 4000, runtime = hms_to_minutes("0:10:00")
    shell:
        """
        exec > {log} 2>&1
        {CONDA_SH}
        conda activate {SCANPY_ENV}
        {SCANPY_ENV}/bin/python {RESCUE_SCRIPTS}/11_wolbachia_secretion.py \\
            --gtf {params.gtf} --de_dir {params.de_dir} --outdir {params.outdir} \\
            --proteins {params.faa} --write_faa {output}.tmp
        sed '/^>/ s/|.*//' {output}.tmp > {output} && rm {output}.tmp
        """

rule hypo_select:
    input:
        de = rules.rescue_de.output.wmel,
        checked = rules.rescue_wolbachia_followup.output[2],
        faa = rules.hypo_proteome_faa.output,
    output:
        hits = f"{HYPO_DIR}/hits.tsv",
        by_tissue = f"{HYPO_DIR}/hits_by_tissue.tsv",
        up = f"{HYPO_DIR}/up_loci.txt",
        faa = f"{HYPO_DIR}/hits.faa",
    params:
        outdir = HYPO_DIR, contrast = config.get("hypo_contrast", "MeiP26_vs_OreR"),
        fdr = config.get("hypo_fdr", 0.05),
        regex = config.get("hypo_product_regex", "hypothetical|DUF|uncharacteri[sz]ed"),
    log: "logs/rescue/hypo_select.log"
    resources: slurm_partition = "medium", mem_mb = 4000, runtime = hms_to_minutes("0:10:00")
    shell:
        """
        exec > {log} 2>&1
        {CONDA_SH}
        conda activate {SCANPY_ENV}
        {SCANPY_ENV}/bin/python {RESCUE_SCRIPTS}/20_hypo_select.py \\
            --de {input.de} --checked {input.checked} --proteome_faa {input.faa} \\
            --outdir {params.outdir} --contrast {params.contrast} --fdr {params.fdr} \\
            --regex "{params.regex}"
        """

rule hypo_interproscan:
    input: rules.hypo_proteome_faa.output
    output: f"{HYPO_DIR}/interproscan/wmel_proteome.tsv"
    params:
        ipr = config.get("interproscan_dir", ""),
        apps = f"--applications {config['interproscan_apps']}" if config.get("interproscan_apps") else "",
        tmp = f"{HYPO_DIR}/interproscan/tmp",
    threads: 16
    log: "logs/rescue/hypo_interproscan.log"
    resources: slurm_partition = "medium", mem_mb = 32000, runtime = hms_to_minutes("4:00:00")
    shell:
        """
        exec > {log} 2>&1
        if [ -z "{params.ipr}" ]; then echo "ERROR: set interproscan_dir in config.yaml"; exit 1; fi
        mkdir -p {params.tmp}
        {params.ipr}/interproscan.sh -i {input} -f tsv -o {output} -dp \\
            -cpu {threads} -T {params.tmp} {params.apps}
        """

rule hypo_eggnog:
    input: rules.hypo_proteome_faa.output
    output: f"{HYPO_DIR}/eggnog/wmel.emapper.annotations"
    params: data = config.get("eggnog_data_dir", ""), outdir = f"{HYPO_DIR}/eggnog"
    threads: 16
    log: "logs/rescue/hypo_eggnog.log"
    resources: slurm_partition = "medium", mem_mb = 32000, runtime = hms_to_minutes("2:00:00")
    shell:
        """
        exec > {log} 2>&1
        if [ -z "{params.data}" ]; then echo "ERROR: set eggnog_data_dir in config.yaml"; exit 1; fi
        {CONDA_SH}
        conda activate {HYPO_ENV}
        emapper.py -i {input} --itype proteins -o wmel --output_dir {params.outdir} \\
            --data_dir {params.data} --cpu {threads} --override
        """

# Hits only (a handful of proteins). Rank-1 model per protein is copied to
# <id>.pdb and its scores to <id>.json so Foldseek and the merge see clean names.
rule hypo_colabfold:
    input: rules.hypo_select.output.faa
    output:
        pdb = directory(f"{HYPO_DIR}/colabfold/rank1_pdb"),
        scores = directory(f"{HYPO_DIR}/colabfold/rank1_scores"),
    params:
        bin = config.get("colabfold_bin", "colabfold_batch"),
        extra = config.get("colabfold_extra", ""),
        raw = f"{HYPO_DIR}/colabfold/raw",
    log: "logs/rescue/hypo_colabfold.log"
    resources:
        slurm_extra = "'--gres=gpu:1'",
        mem_mb = 64000, runtime = hms_to_minutes("4:00:00")
    shell:
        """
        exec > {log} 2>&1
        mkdir -p {output.pdb} {output.scores} {params.raw}
        if [ -s {input} ]; then
            {params.bin} {params.extra} {input} {params.raw}
            for f in {params.raw}/*_unrelaxed_rank_001_*.pdb; do
                id=$(basename $f | sed 's/_unrelaxed_rank_001.*//')
                cp $f {output.pdb}/$id.pdb
                cp {params.raw}/${{id}}_scores_rank_001_*.json {output.scores}/$id.json
            done
        else
            echo "no hit proteins, nothing to fold"
        fi
        """

# Pre-build each database once: foldseek databases PDB / Alphafold/Swiss-Prot / Alphafold/Proteome.
rule hypo_foldseek:
    input: rules.hypo_colabfold.output.pdb
    output: f"{HYPO_DIR}/foldseek/{{db}}.m8"
    params:
        db = lambda wc: FOLDSEEK_DBS[wc.db],
        tmp = f"{HYPO_DIR}/foldseek/tmp_{{db}}",
    threads: 16
    log: "logs/rescue/hypo_foldseek_{db}.log"
    resources: slurm_partition = "medium", mem_mb = 64000, runtime = hms_to_minutes("2:00:00")
    shell:
        """
        exec > {log} 2>&1
        {CONDA_SH}
        conda activate {HYPO_ENV}
        if ls {input}/*.pdb >/dev/null 2>&1; then
            foldseek easy-search {input} {params.db} {output} {params.tmp} -e 1e-3 \\
                --threads {threads} \\
                --format-output "query,target,fident,alnlen,evalue,bits,alntmscore,prob,theader"
        else
            : > {output}
        fi
        """

rule hypo_context:
    input:
        hits = rules.hypo_select.output.hits,
        up = rules.hypo_select.output.up,
    output: f"{HYPO_DIR}/hits_context.tsv"
    params: gtf = WMEL_GFF, window = config.get("hypo_context_window", 5)
    log: "logs/rescue/hypo_context.log"
    resources: slurm_partition = "medium", mem_mb = 4000, runtime = hms_to_minutes("0:10:00")
    shell:
        """
        exec > {log} 2>&1
        {CONDA_SH}
        conda activate {SCANPY_ENV}
        {SCANPY_ENV}/bin/python {RESCUE_SCRIPTS}/21_hypo_context.py --gtf {params.gtf} \\
            --hits {input.hits} --up_loci {input.up} --window {params.window} --out {output}
        """

rule hypo_merge:
    input:
        hits = rules.hypo_select.output.hits,
        secretion = rules.rescue_secretion.output.all_loci,
        interpro = rules.hypo_interproscan.output,
        eggnog = rules.hypo_eggnog.output,
        scores = rules.hypo_colabfold.output.scores,
        context = rules.hypo_context.output,
        foldseek = expand(f"{HYPO_DIR}/foldseek/{{db}}.m8", db=list(FOLDSEEK_DBS)),
    output: f"{HYPO_DIR}/hypothetical_evidence.tsv"
    params: fs = " ".join(f"{db}={HYPO_DIR}/foldseek/{db}.m8" for db in FOLDSEEK_DBS)
    log: "logs/rescue/hypo_merge.log"
    resources: slurm_partition = "medium", mem_mb = 4000, runtime = hms_to_minutes("0:10:00")
    shell:
        """
        exec > {log} 2>&1
        {CONDA_SH}
        conda activate {SCANPY_ENV}
        {SCANPY_ENV}/bin/python {RESCUE_SCRIPTS}/22_hypo_merge.py --hits {input.hits} \\
            --secretion {input.secretion} --interpro {input.interpro} --eggnog {input.eggnog} \\
            --scores_dir {input.scores} --context {input.context} \\
            --foldseek {params.fs} --out {output}
        """

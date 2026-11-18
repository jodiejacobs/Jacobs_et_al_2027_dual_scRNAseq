#!/usr/bin/env python
"""Patch the Snakefile for the eggshell figure and the host GSEA.

  - removes rule cluster_gsea (one-vs-rest GSEA; cluster_gsea.py stays in the
    repo, unwired) and its line in the section header
  - adds rule rescue_host_gsea  (17_host_gsea.py, rescue contrasts per cell type)
  - adds rule rescue_eggshell   (18_eggshell.py, per-cell eggshell scores)
  - adds both to rescue_followup_all
  - adds the eggshell and GSEA panels to the _S plot registry
  - flags eggshell genes as ambient in fig4_interaction_volcano
  - adds rule rescue_plot_umap_eggshell and puts it in rescue_plots_all

Writes Snakefile.bak_<date> first. Safe to run twice: each step checks
whether it has already been applied.

Usage (repo root):  python add_eggshell_rules.py
"""
import datetime
import re
import shutil
import sys

SF = "Snakefile"
s = open(SF).read()
orig = s


def once(marker, fn):
    global s
    if marker in s:
        print(f"  already applied: {marker}")
        return
    s = fn(s)


# 1. remove cluster_gsea -----------------------------------------------------
if "rule cluster_gsea:" in s:
    m = re.search(r"\nrule cluster_gsea:\n.*?\n        \"\"\"\n\n", s, re.S)
    if not m:
        sys.exit("could not find the end of rule cluster_gsea; nothing written")
    s = s[:m.start()] + "\n" + s[m.end():]
    s = s.replace("#   snakemake cluster_gsea   one-vs-rest GSEA per cell type\n",
                  "#   (cluster_gsea, one-vs-rest GSEA, was replaced by "
                  "rescue_host_gsea)\n")
    print("  removed rule cluster_gsea")
else:
    print("  already applied: cluster_gsea removed")

# 2. analysis rules, inserted before rescue_followup_all ----------------------
RULES = '''
# Host pathway enrichment on the rescue contrasts, per cell type: genes ranked
# by signed -log10 p from rescue_de, gseapy prerank against KEGG_2019 and
# GO_Biological_Process_2018. Replaces the one-vs-rest cluster_gsea, which
# pooled the four groups and so could not speak to rescue. The gmt cache is the
# one cluster_gsea downloaded; 17_host_gsea.py downloads only what is missing.
rule rescue_host_gsea:
    input:
        de = rules.rescue_de.output.host,
        script = f"{RESCUE_SCRIPTS}/17_host_gsea.py",
        gene_sets = f"{RESCUE_SCRIPTS}/plots/gene_sets.py",
    output:
        gsea = f"{RESCUE_DIR}/host_gsea/host_gsea.tsv",
        summary = f"{RESCUE_DIR}/host_gsea/host_gsea_summary.tsv",
    params:
        outdir = f"{RESCUE_DIR}/host_gsea",
        gmt_dir = config.get("host_gsea_gmt_dir", "results/cluster_gsea/gene_sets"),
    log: "logs/rescue/host_gsea.log"
    threads: 4
    resources:
        slurm_partition = "medium", mem_mb = 16000, runtime = hms_to_minutes("1:30:00")
    shell:
        """
        exec > {log} 2>&1
        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}
        {SCANPY_ENV}/bin/python {input.script} \\
            --de {input.de} --gmt_dir {params.gmt_dir} \\
            --outdir {params.outdir} --threads {threads}
        """

# Eggshell program per cell: vitelline membrane, yolk and chorion scores,
# fraction of follicle cells above each library's germline ambient floor and
# score among those cells, tested on the same 4-df per-sample model as the
# composition panels. Gene sets live in plots/gene_sets.py.
rule rescue_eggshell:
    input:
        counts = rules.rescue_pseudobulk.output.counts,
        annotation = config["flybase_annotation"],
        script = f"{RESCUE_SCRIPTS}/18_eggshell.py",
        gene_sets = f"{RESCUE_SCRIPTS}/plots/gene_sets.py",
    output:
        genes = f"{RESCUE_DIR}/eggshell/eggshell_genes.tsv",
        by_sample = f"{RESCUE_DIR}/eggshell/eggshell_by_sample.tsv",
        tests = f"{RESCUE_DIR}/eggshell/eggshell_tests.tsv",
        floor = f"{RESCUE_DIR}/eggshell/eggshell_floor.tsv",
        ambient = f"{RESCUE_DIR}/eggshell/eggshell_ambient.tsv",
        titer = f"{RESCUE_DIR}/eggshell/eggshell_titer.tsv",
    params:
        outdir = f"{RESCUE_DIR}/eggshell",
    log: "logs/rescue/eggshell.log"
    resources:
        slurm_partition = "medium", mem_mb = 64000, runtime = hms_to_minutes("1:00:00")
    shell:
        """
        exec > {log} 2>&1
        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}
        {SCANPY_ENV}/bin/python {input.script} \\
            --counts {input.counts} --annotation {input.annotation} \\
            --outdir {params.outdir}
        """

'''
once("rule rescue_host_gsea:",
     lambda t: t.replace("\nrule rescue_followup_all:\n",
                         RULES + "rule rescue_followup_all:\n", 1))

# 3. followup_all targets -----------------------------------------------------
once('f"{RESCUE_DIR}/host_gsea/host_gsea.tsv",\n        f"{RESCUE_DIR}/eggshell',
     lambda t: t.replace(
         'rule rescue_followup_all:\n    input:\n',
         'rule rescue_followup_all:\n    input:\n'
         '        f"{RESCUE_DIR}/host_gsea/host_gsea.tsv",\n'
         '        f"{RESCUE_DIR}/eggshell/eggshell_tests.tsv",\n', 1))

# 4. volcano: eggshell genes drawn as ambient in germline panels -------------
VOLC_OLD = '''        f"--null {RESCUE_DIR}/rescue_null/rescue_null_summary.tsv"),'''
VOLC_NEW = '''        f"--null {RESCUE_DIR}/rescue_null/rescue_null_summary.tsv"
        " --ambient_genes eggshell"),'''
if "--ambient_genes eggshell" in s:
    print("  already applied: volcano ambient flag")
elif VOLC_OLD in s:
    s = s.replace(VOLC_OLD, VOLC_NEW, 1)
else:
    sys.exit("fig4_interaction_volcano entry not in the expected form; "
             "nothing written")

# 5. _S registry --------------------------------------------------------------
REG = '''    # Eggshell figure. The heatmap and contrast summary read the edgeR table;
    # the score, controls and UMAP read rescue_eggshell; the GSEA panels read
    # rescue_host_gsea. Second tables passed in extra come from the same rule
    # as the declared one, so the dependency is still tracked.
    "fig_eggshell_heatmap": (
        "plot_eggshell_heatmap.py",
        f"{RESCUE_DIR}/de/host_de_all.tsv.gz",
        f"--logcpm_dir {RESCUE_DIR}/de"),
    "fig_eggshell_contrasts": (
        "plot_eggshell_contrasts.py",
        f"{RESCUE_DIR}/de/host_de_all.tsv.gz", ""),
    "fig_eggshell_score": (
        "plot_eggshell_score.py",
        f"{RESCUE_DIR}/eggshell/eggshell_by_sample.tsv",
        f"--tests {RESCUE_DIR}/eggshell/eggshell_tests.tsv"),
    "supp_eggshell_controls": (
        "plot_eggshell_controls.py",
        f"{RESCUE_DIR}/eggshell/eggshell_ambient.tsv",
        f"--titer {RESCUE_DIR}/eggshell/eggshell_titer.tsv"),
    "fig_eggshell_gsea": (
        "plot_host_gsea.py",
        f"{RESCUE_DIR}/host_gsea/host_gsea.tsv", "--preset secretory"),
    "supp_host_gsea": (
        "plot_host_gsea.py",
        f"{RESCUE_DIR}/host_gsea/host_gsea.tsv", "--preset top"),
}

PLOTS = {**_V, **_S}'''
once('"fig_eggshell_heatmap": (',
     lambda t: t.replace("}\n\nPLOTS = {**_V, **_S}", REG, 1))

# 6. eggshell UMAP rule -------------------------------------------------------
UMAP = '''
rule rescue_plot_umap_eggshell:
    """Summed eggshell expression on the atlas UMAP, one panel per group.

    Same script and the same debris rule as rescue_plot_umap_gene, with the
    22 eggshell genes from rescue_eggshell summed per cell instead of one gene.
    """
    input:
        integrated = "results/integrated/integrated.h5ad",
        expression = f"{RESCUE_DIR}/pseudobulk/counts_labeled.h5ad",
        gene_set = f"{RESCUE_DIR}/eggshell/eggshell_genes.tsv",
        exclude = RESCUE_DEBRIS if RESCUE_DEBRIS else [],
        script = f"{PLOT_SCRIPTS}/plot_umap_gene.py",
        style = f"{PLOT_SCRIPTS}/figstyle.py",
    output:
        f"{FIG_DIR}/fig_eggshell_umap.pdf",
    params:
        outdir = FIG_DIR,
        exclude_flag = f"--exclude_cells {RESCUE_DEBRIS}" if RESCUE_DEBRIS else "",
    log: "logs/rescue/plots/fig_eggshell_umap.log"
    resources:
        slurm_partition = "medium", mem_mb = 64000, runtime = hms_to_minutes("1:00:00")
    shell:
        """
        exec > {log} 2>&1
        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}
        {SCANPY_ENV}/bin/python {input.script} \\
            --integrated {input.integrated} \\
            --expression {input.expression} --outdir {params.outdir} \\
            --gene_set {input.gene_set} --gene_label eggshell \\
            --stem fig_eggshell_umap {params.exclude_flag}
        """
'''
once("rule rescue_plot_umap_eggshell:",
     lambda t: t.replace("\nrule rescue_plot_gene_group_summary:",
                         UMAP + "\nrule rescue_plot_gene_group_summary:", 1))
ALL_MARK = ('f"{FIG_DIR}/fig1_umap_annotation.pdf",\n'
            '        f"{FIG_DIR}/fig_eggshell_umap.pdf",\n')
once(ALL_MARK,
     lambda t: t.replace(
         '        f"{FIG_DIR}/fig1_umap_annotation.pdf",\n',
         '        f"{FIG_DIR}/fig1_umap_annotation.pdf",\n'
         '        f"{FIG_DIR}/fig_eggshell_umap.pdf",\n', 1))

if s == orig:
    print("nothing to change")
    sys.exit(0)
for need in ("rule rescue_host_gsea:", "rule rescue_eggshell:",
             '"fig_eggshell_heatmap": (', "rule rescue_plot_umap_eggshell:",
             ALL_MARK, "--ambient_genes eggshell"):
    if need not in s:
        sys.exit(f"patch incomplete ({need!r} missing); nothing written")
bak = f"{SF}.bak_{datetime.date.today():%Y%m%d}_eggshell"
shutil.copy(SF, bak)
open(SF, "w").write(s)
print(f"wrote {SF} (backup {bak})")

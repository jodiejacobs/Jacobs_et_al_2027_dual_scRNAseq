#!/usr/bin/env python
"""Patch the Snakefile for the marker UMAPs and the marker expression table.

  - rule rescue_marker_expression (19_marker_expression.py), added to
    rescue_followup_all
  - rule rescue_plot_umap_markers, one UMAP figure per gene in MARKERS
    (plots/gene_sets.py) except mei-P26, which keeps rule
    rescue_plot_umap_gene and fig4_umap_meip26; added to rescue_plots_all

Writes Snakefile.bak_<date>_markers first; safe to run twice.

Usage (repo root):  python add_marker_rules.py
"""
import datetime
import shutil
import sys

SF = "Snakefile"
s = open(SF).read()
orig = s

RULES = '''
rule rescue_marker_expression:
    input:
        counts = rules.rescue_pseudobulk.output.counts,
        script = f"{RESCUE_SCRIPTS}/19_marker_expression.py",
        gene_sets = f"{PLOT_SCRIPTS}/gene_sets.py",
    output:
        long = f"{RESCUE_DIR}/marker_expression/marker_expression_long.tsv",
        cpm = f"{RESCUE_DIR}/marker_expression/marker_expression_cpm.tsv",
        frac = f"{RESCUE_DIR}/marker_expression/marker_expression_frac.tsv",
        counts_out = f"{RESCUE_DIR}/marker_expression/marker_expression_counts.tsv",
    params:
        outdir = f"{RESCUE_DIR}/marker_expression",
    log: "logs/rescue/marker_expression.log"
    resources:
        slurm_partition = "medium", mem_mb = 32000, runtime = hms_to_minutes("0:30:00")
    shell:
        """
        exec > {log} 2>&1
        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}
        {SCANPY_ENV}/bin/python {input.script} \\
            --counts {input.counts} --outdir {params.outdir}
        """

rule rescue_plot_umap_markers:
    """One marker gene on the atlas UMAP, one panel per group.

    Same script and debris rule as rescue_plot_umap_gene. wMel genes (WD_ ids)
    are normalised to the cell's wMel UMIs; cells with no wMel reads are drawn
    in the grey background only.
    """
    input:
        integrated = "results/integrated/integrated.h5ad",
        expression = f"{RESCUE_DIR}/pseudobulk/counts_labeled.h5ad",
        exclude = RESCUE_DEBRIS if RESCUE_DEBRIS else [],
        script = f"{PLOT_SCRIPTS}/plot_umap_gene.py",
        style = f"{PLOT_SCRIPTS}/figstyle.py",
    output:
        f"{FIG_DIR}/fig_umap_{{marker}}.pdf",
    wildcard_constraints:
        marker = "|".join(UMAP_MARKER_PLOTS),
    params:
        outdir = FIG_DIR,
        gene = lambda w: UMAP_MARKER_PLOTS[w.marker],
        exclude_flag = f"--exclude_cells {RESCUE_DEBRIS}" if RESCUE_DEBRIS else "",
    log: "logs/rescue/plots/fig_umap_{marker}.log"
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
            --gene {params.gene} --gene_label {wildcards.marker} \\
            --stem fig_umap_{wildcards.marker} {params.exclude_flag}
        """
'''

DEFS = '''
# Marker genes on the UMAP and per cell type. The list is MARKERS in
# plots/gene_sets.py, so adding a gene there adds a figure and a table row.
# mei-P26 is drawn by rescue_plot_umap_gene (fig4_umap_meip26) and only
# tabulated here.
import sys as _sys
_sys.path.insert(0, PLOT_SCRIPTS)
from gene_sets import MARKERS as UMAP_MARKERS
UMAP_MARKER_PLOTS = {k: v for k, v in UMAP_MARKERS.items() if k != "mei-P26"}

'''
REG = "\nPLOTS = {**_V, **_S}\n"
if "UMAP_MARKER_PLOTS = " in s:
    print("  already applied: marker definitions")
elif REG in s:
    # Before rescue_plots_all, which lists the marker figures.
    s = s.replace(REG, REG + DEFS, 1)
else:
    sys.exit("PLOTS registry line not found; nothing written")

ANCHOR = "\nrule rescue_plot_gene_group_summary:"
if "rule rescue_plot_umap_markers:" in s:
    print("  already applied: marker rules")
elif ANCHOR in s:
    s = s.replace(ANCHOR, RULES + ANCHOR, 1)
else:
    sys.exit("rule rescue_plot_gene_group_summary not found; nothing written")

FU = 'rule rescue_followup_all:\n    input:\n'
FU_ADD = '        f"{RESCUE_DIR}/marker_expression/marker_expression_cpm.tsv",\n'
if FU_ADD in s:
    print("  already applied: followup target")
elif FU in s:
    s = s.replace(FU, FU + FU_ADD, 1)
else:
    sys.exit("rule rescue_followup_all not found; nothing written")

PA = '        f"{FIG_DIR}/fig1_umap_annotation.pdf",\n'
PA_ADD = ('        expand(f"{FIG_DIR}/fig_umap_{{marker}}.pdf", '
          'marker=list(UMAP_MARKER_PLOTS)),\n')
if PA_ADD in s:
    print("  already applied: plots_all target")
elif PA in s:
    s = s.replace(PA, PA + PA_ADD, 1)
else:
    sys.exit("fig1_umap_annotation target not found in rescue_plots_all; "
             "nothing written")

if s == orig:
    print("nothing to change")
    sys.exit(0)
bak = f"{SF}.bak_{datetime.date.today():%Y%m%d}_markers"
shutil.copy(SF, bak)
open(SF, "w").write(s)
print(f"wrote {SF} (backup {bak})")

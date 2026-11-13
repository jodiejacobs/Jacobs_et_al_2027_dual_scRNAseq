#!/usr/bin/env python3
"""Add rule rescue_plot_umap_gene after rescue_plot_umap_annotation. Idempotent."""
import re, shutil, sys
P = sys.argv[1] if len(sys.argv) > 1 else "Snakefile"

RULE = '''
rule rescue_plot_umap_gene:
    """mei-P26 on the atlas UMAP, one panel per genotype by infection group.

    Reads the integrated AnnData for the same reason rescue_plot_umap_annotation
    does: obsm["X_umap_atlas"] exists nowhere else. The script refuses to run
    without --exclude_cells, because that object is written before
    06_label_qc.py calls debris and still holds every barcode, so an
    unfiltered version would not describe the same cells as the rest of the
    paper. If RESCUE_DEBRIS is empty the script fails with that message rather
    than quietly including debris.
    """
    input:
        integrated = "results/integrated/integrated.h5ad",
        exclude = RESCUE_DEBRIS if RESCUE_DEBRIS else [],
        script = f"{PLOT_SCRIPTS}/plot_umap_gene.py",
        style = f"{PLOT_SCRIPTS}/figstyle.py",
    output:
        f"{FIG_DIR}/fig4_umap_meip26.pdf",
    params:
        outdir = FIG_DIR,
        exclude_flag = f"--exclude_cells {RESCUE_DEBRIS}" if RESCUE_DEBRIS else "",
        gene = config.get("rescue_umap_gene", "FBgn0026206"),
        gene_label = config.get("rescue_umap_gene_label", "mei-P26"),
    log: "logs/rescue/plots/fig4_umap_meip26.log"
    resources:
        slurm_partition = "medium", mem_mb = 64000, runtime = hms_to_minutes("1:00:00")
    shell:
        """
        exec > {log} 2>&1
        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}
        {SCANPY_ENV}/bin/python {input.script} \\
            --integrated {input.integrated} --outdir {params.outdir} \\
            --gene {params.gene} --gene_label {params.gene_label} \\
            {params.exclude_flag}
        """

'''

s = open(P).read()
if "rule rescue_plot_umap_gene:" in s:
    print("already present"); raise SystemExit
m = re.search(r"\nrule rescue_plot_umap_annotation:\n", s)
if not m: raise SystemExit("anchor rule not found")
nxt = re.search(r"\nrule \w+:\n", s[m.end():])
if not nxt: raise SystemExit("could not find the rule after it")
at = m.end() + nxt.start() + 1
shutil.copyfile(P, P + ".bak5")
open(P, "w").write(s[:at] + RULE.lstrip("\n") + s[at:])
print("added rule rescue_plot_umap_gene")

#!/usr/bin/env python3
"""Point rescue_plot_umap_gene at counts_labeled.h5ad for expression.

integrated.h5ad is subset to the highly variable genes during integration, so
mei-P26 is not in it. Coordinates still come from there; expression comes from
the pseudobulk input, which keeps every gene. Idempotent.
"""
import sys, shutil
P = sys.argv[1] if len(sys.argv) > 1 else "Snakefile"
s = open(P).read()
if "--expression {input.expression}" in s:
    print("already patched"); raise SystemExit

OLD_IN = '''        integrated = "results/integrated/integrated.h5ad",
        exclude = RESCUE_DEBRIS if RESCUE_DEBRIS else [],
        script = f"{PLOT_SCRIPTS}/plot_umap_gene.py",'''
NEW_IN = '''        integrated = "results/integrated/integrated.h5ad",
        # integrated.h5ad is cut to the highly variable genes, so most genes
        # are not in it. Coordinates come from there, expression from here.
        expression = f"{RESCUE_DIR}/pseudobulk/counts_labeled.h5ad",
        exclude = RESCUE_DEBRIS if RESCUE_DEBRIS else [],
        script = f"{PLOT_SCRIPTS}/plot_umap_gene.py",'''
OLD_SH = '''            --integrated {input.integrated} --outdir {params.outdir} \\
            --gene {params.gene} --gene_label {params.gene_label} \\
            {params.exclude_flag}'''
NEW_SH = '''            --integrated {input.integrated} \\
            --expression {input.expression} --outdir {params.outdir} \\
            --gene {params.gene} --gene_label {params.gene_label} \\
            {params.exclude_flag}'''

k = s.find("rule rescue_plot_umap_gene:")
if k == -1: raise SystemExit("rule rescue_plot_umap_gene not found")
end = s.find("\nrule ", k + 10)
end = len(s) if end == -1 else end
blk = s[k:end]
for old, new, what in ((OLD_IN, NEW_IN, "input"), (OLD_SH, NEW_SH, "shell")):
    if old not in blk: raise SystemExit(f"could not find the {what} block to patch")
    blk = blk.replace(old, new, 1)
shutil.copyfile(P, P + ".bak6")
open(P, "w").write(s[:k] + blk + s[end:])
print("patched rescue_plot_umap_gene to pass --expression")

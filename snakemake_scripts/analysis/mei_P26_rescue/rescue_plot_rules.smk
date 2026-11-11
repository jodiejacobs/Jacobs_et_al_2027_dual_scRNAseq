# ---------------------------------------------------------------------------
# Publication figure panels, added 2026-11-09
#
# One plotting script per claim-bearing table under results/rescue. Each script
# reads exactly one table, writes PDF (vector, for Illustrator) and PNG at
# 400 dpi, and prints the numbers behind its panels to the log so a run record
# exists. Shared style lives in plots/figstyle.py.
#
# Most scripts write two versions of their figure:
#   <stem>_robust.*   the filtered view the paper's claims rest on
#   <stem>_all.*      everything, so what the filter removed stays visible
# The filter differs by table (the composite `robust` DE flag, sign agreement
# across infected samples, a minimum gene or cell count) and each script's
# docstring says which it applies. Scripts whose table has nothing to filter on
# write a single figure.
#
# Plotting is cheap: seconds, a few hundred MB. Everything runs on short.
# ---------------------------------------------------------------------------
PLOT_SCRIPTS = f"{RESCUE_SCRIPTS}/plots"
FIG_DIR = config.get("rescue_figure_dir", f"{RESCUE_DIR}/figures")

# stem -> (script, table, extra CLI args)
_V = {  # two versions: <stem>_robust and <stem>_all
    "fig2_wmel_gene_detection": (
        "plot_wmel_gene_detection.py",
        f"{RESCUE_DIR}/detection_bias/wmel_gene_detection.tsv.gz", ""),
    "fig2_titer_separation": (
        "plot_titer_separation.py",
        f"{RESCUE_DIR}/uninfected_floor/titer_separation.tsv", ""),
    "fig2_ambient_profile": (
        "plot_ambient_profile.py",
        f"{RESCUE_DIR}/uninfected_floor/ambient_profile.tsv.gz", ""),
    "fig3_composition_proportions": (
        "plot_composition_proportions.py",
        f"{RESCUE_DIR}/composition/composition_proportions.tsv", ""),
    "fig3_composition_tests": (
        "plot_composition_tests.py",
        f"{RESCUE_DIR}/composition/composition_tests.tsv", ""),
    "fig3_titer_by_sample_tissue": (
        "plot_titer_by_sample_tissue.py",
        f"{RESCUE_DIR}/titer/titer_by_sample_tissue.tsv", ""),
    "fig3_titer_genotype": (
        "plot_titer_genotype_by_tissue.py",
        f"{RESCUE_DIR}/titer/titer_genotype_by_tissue.tsv", ""),
    "fig4_germline_composition": (
        "plot_germline_composition_by_sample.py",
        f"{RESCUE_DIR}/germline/germline_composition_by_sample.tsv", ""),
    "fig4_germline_state": (
        "plot_germline_state_by_sample.py",
        f"{RESCUE_DIR}/germline/germline_state_by_sample.tsv", ""),
    "fig4_germline_subclusters": (
        "plot_germline_subclusters.py",
        f"{RESCUE_DIR}/germline/germline_subclusters.tsv", ""),
    "fig4_germline_composition_tests": (
        "plot_germline_tests.py",
        f"{RESCUE_DIR}/germline/germline_composition_tests.tsv", ""),
    "fig4_germline_state_tests": (
        "plot_germline_tests.py",
        f"{RESCUE_DIR}/germline/germline_state_tests.tsv", ""),
    "fig5_rescue_categories": (
        "plot_rescue_category_counts.py",
        f"{RESCUE_DIR}/summary/rescue_category_counts.tsv", ""),
    "fig5_rescue_null": (
        "plot_rescue_null_summary.py",
        f"{RESCUE_DIR}/rescue_null/rescue_null_summary.tsv", ""),
    "fig5_host_de_summary": (
        "plot_host_de_summary.py",
        f"{RESCUE_DIR}/de/host_de_summary.tsv", ""),
    "fig6_wolbachia_ambient_check": (
        "plot_wolbachia_ambient_check.py",
        f"{RESCUE_DIR}/de/wolbachia_ambient_check.tsv", ""),
    "fig6_wolbachia_de_summary": (
        "plot_wolbachia_de_summary.py",
        f"{RESCUE_DIR}/de/wolbachia_de_summary.tsv", ""),
    "fig6_wolbachia_de": (
        "plot_wolbachia_de_sig_checked.py",
        f"{RESCUE_DIR}/wolbachia_followup/wolbachia_de_sig_checked.tsv", ""),
    # Recounts contrasts over real cell types from the DE table, because the
    # counts shipped in the secretion table include the pooled contrast.
    "fig6_secretion_candidates": (
        "plot_wolbachia_secretion_all.py",
        f"{RESCUE_DIR}/wolbachia_followup/wolbachia_secretion_all.tsv",
        f"--de_table {RESCUE_DIR}/wolbachia_followup/wolbachia_de_sig_checked.tsv"),
    "fig7_titer_programs_bycelltype": (
        "plot_titer_program_correlations.py",
        f"{RESCUE_DIR}/titer_programs/titer_program_correlations.tsv", ""),
    "fig7_titer_program_tests": (
        "plot_titer_program_tests.py",
        f"{RESCUE_DIR}/titer_programs/titer_program_tests.tsv", ""),
    "fig7_titer_programs_per_sample": (
        "plot_titer_program_per_sample.py",
        f"{RESCUE_DIR}/titer_programs/titer_program_per_sample.tsv", ""),
    "supp_atlas_profile_correlations": (
        "plot_atlas_profile_correlations.py",
        f"{RESCUE_DIR}/atlas_fidelity/atlas_profile_correlations.tsv", ""),
    "supp_wolbachia_gsea": (
        "plot_wolbachia_gsea.py",
        f"{RESCUE_DIR}/wolbachia_followup/wolbachia_gsea.tsv", ""),
}

_S = {  # one figure: the table has nothing meaningful to filter on
    "fig1_wmel_breadth": (
        "plot_wmel_breadth.py",
        f"{RESCUE_DIR}/detection_bias/wmel_breadth.tsv", ""),
    # --total_loci lets the script tell complete recovery (the curve reached
    # every annotated locus) from an unsaturated one (it ran out of cells).
    "fig1_wmel_rarefaction": (
        "plot_wmel_rarefaction.py",
        f"{RESCUE_DIR}/detection_bias/wmel_rarefaction.tsv",
        "--total_loci {}".format(config.get("rescue_total_wmel_loci", 1056))),
    "fig2_wmel_feature_tests": (
        "plot_wmel_feature_tests.py",
        f"{RESCUE_DIR}/detection_bias/wmel_feature_tests.tsv", ""),
    "fig2_wmel_detection_model": (
        "plot_wmel_detection_model.py",
        f"{RESCUE_DIR}/detection_bias/wmel_detection_model.tsv", ""),
    "fig2_uninfected_floor": (
        "plot_uninfected_floor.py",
        f"{RESCUE_DIR}/uninfected_floor/uninfected_floor.tsv",
        "--floor_quantile {}".format(config.get("rescue_floor_quantile", 0.99))),
    "fig2_ambient_depth": (
        "plot_ambient_depth.py",
        f"{RESCUE_DIR}/uninfected_floor/ambient_depth.tsv", ""),
    "fig3_composition_cell_counts": (
        "plot_composition_cell_counts.py",
        f"{RESCUE_DIR}/composition/composition_cell_counts.tsv", ""),
    "supp_atlas_proportions": (
        "plot_atlas_proportions.py",
        f"{RESCUE_DIR}/atlas_fidelity/atlas_proportions.tsv", ""),
    "supp_atlas_correlation_matrix": (
        "plot_atlas_profile_correlation_matrix.py",
        f"{RESCUE_DIR}/atlas_fidelity/atlas_profile_correlation_matrix.tsv", ""),
    "supp_lineage_by_label": (
        "plot_lineage_by_label.py",
        f"{RESCUE_DIR}/label_qc/lineage_by_label.tsv", ""),
    "supp_cell_counts_debris": (
        "plot_cell_counts_by_sample_tissue.py",
        f"{RESCUE_DIR}/label_qc/cell_counts_by_sample_tissue.tsv", ""),
}

PLOTS = {**_V, **_S}


def _plot_field(stem, i):
    return PLOTS[stem][i]


rule rescue_plots_all:
    input:
        [f"{FIG_DIR}/{s}_{v}.pdf" for s in _V for v in ("robust", "all")],
        [f"{FIG_DIR}/{s}.pdf" for s in _S],


rule rescue_plot_versioned:
    """One table -> two figures: the filtered view and everything."""
    input:
        table = lambda w: _plot_field(w.stem, 1),
        script = lambda w: f"{PLOT_SCRIPTS}/{_plot_field(w.stem, 0)}",
        style = f"{PLOT_SCRIPTS}/figstyle.py",
    output:
        robust_pdf = f"{FIG_DIR}/{{stem}}_robust.pdf",
        robust_png = f"{FIG_DIR}/{{stem}}_robust.png",
        all_pdf = f"{FIG_DIR}/{{stem}}_all.pdf",
        all_png = f"{FIG_DIR}/{{stem}}_all.png",
    wildcard_constraints:
        stem = "|".join(_V),
    params:
        outdir = FIG_DIR,
        extra = lambda w: _plot_field(w.stem, 2),
    log: "logs/rescue/plots/{stem}.log"
    resources:
        slurm_partition = "short", mem_mb = 8000, runtime = hms_to_minutes("0:20:00")
    shell:
        """
        exec > {log} 2>&1
        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}
        {SCANPY_ENV}/bin/python {input.script} \
            --input {input.table} --outdir {params.outdir} \
            --stem {wildcards.stem} --version both {params.extra}
        """


rule rescue_plot_single:
    """One table -> one figure, where there is nothing to filter on."""
    input:
        table = lambda w: _plot_field(w.stem, 1),
        script = lambda w: f"{PLOT_SCRIPTS}/{_plot_field(w.stem, 0)}",
        style = f"{PLOT_SCRIPTS}/figstyle.py",
    output:
        pdf = f"{FIG_DIR}/{{stem}}.pdf",
        png = f"{FIG_DIR}/{{stem}}.png",
    wildcard_constraints:
        stem = "|".join(_S),
    params:
        outdir = FIG_DIR,
        extra = lambda w: _plot_field(w.stem, 2),
    log: "logs/rescue/plots/{stem}.log"
    resources:
        slurm_partition = "short", mem_mb = 8000, runtime = hms_to_minutes("0:20:00")
    shell:
        """
        exec > {log} 2>&1
        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {SCANPY_ENV}
        {SCANPY_ENV}/bin/python {input.script} \
            --input {input.table} --outdir {params.outdir} \
            --stem {wildcards.stem} {params.extra}
        """

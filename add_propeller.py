#!/usr/bin/env python3
"""Insert the propeller rule and its plot registry entry into the Snakefile.

Idempotent: run it twice and the second run reports that both are already
there. Edits in place after writing a .bak, because the Snakefile is large and
retyping it wholesale is a worse risk than a targeted insertion.
"""
import re
import shutil
import sys

PATH = sys.argv[1] if len(sys.argv) > 1 else "Snakefile"

RULE = '''
rule rescue_composition_propeller:
    # The moderated version of rescue_composition. 04_cell_composition.py fits
    # each cell type's logit proportion by OLS on 4 residual df; propeller fits
    # the same design and contrasts but borrows variance across cell types with
    # limma's empirical Bayes step, which is where the extra power comes from.
    # Both tables are kept so the difference between them stays visible.
    #
    # composition_cell_counts.tsv is written by rescue_composition but is not
    # one of its declared outputs, so the dependency is on the declared table
    # and the counts path is passed as a param. Declaring it instead would
    # change that rule's outputs and rerun everything downstream of it.
    input:
        tests = rules.rescue_composition.output,
    output:
        f"{RESCUE_DIR}/composition/composition_tests_propeller.tsv",
    params:
        counts = f"{RESCUE_DIR}/composition/composition_cell_counts.tsv",
        outdir = f"{RESCUE_DIR}/composition",
    log: "logs/rescue/composition_propeller.log"
    resources:
        slurm_partition = "medium", mem_mb = 32000, runtime = hms_to_minutes("0:20:00")
    shell:
        """
        exec > {log} 2>&1
        source /private/groups/russelllab/jodie/miniforge3/etc/profile.d/conda.sh
        conda activate {EDGER_ENV}
        {EDGER_ENV}/bin/Rscript {RESCUE_SCRIPTS}/04b_propeller.R \\
            --counts {params.counts} --outdir {params.outdir}
        """

'''

ENTRY = '''    "fig3_composition_tests_propeller": (
        "plot_composition_tests.py",
        f"{RESCUE_DIR}/composition/composition_tests_propeller.tsv", ""),
'''


def main():
    src = open(PATH).read()
    out = src
    done, skipped = [], []

    if "rule rescue_composition_propeller:" in out:
        skipped.append("rule")
    else:
        m = re.search(r"\nrule rescue_composition:\n", out)
        if not m:
            sys.exit("rule rescue_composition not found; nothing inserted")
        # End of that rule = the start of the next top-level rule after it.
        nxt = re.search(r"\nrule \w+:\n", out[m.end():])
        if not nxt:
            sys.exit("could not find the rule after rescue_composition")
        at = m.end() + nxt.start() + 1
        out = out[:at] + RULE.lstrip("\n") + out[at:]
        done.append("rule rescue_composition_propeller")

    if '"fig3_composition_tests_propeller":' in out:
        skipped.append("registry entry")
    else:
        # Insert straight after the existing fig3_composition_tests entry, which
        # ends at the first line closing with '"", ),' style after its key.
        k = out.find('    "fig3_composition_tests": (')
        if k == -1:
            sys.exit("fig3_composition_tests registry entry not found")
        end = out.find("),\n", k)
        if end == -1:
            sys.exit("could not find the end of the fig3_composition_tests entry")
        at = end + 3
        out = out[:at] + ENTRY + out[at:]
        done.append("fig3_composition_tests_propeller registry entry")

    if out == src:
        print("Already up to date: " + ", ".join(skipped))
        return
    shutil.copyfile(PATH, PATH + ".bak")
    open(PATH, "w").write(out)
    print("Added: " + "; ".join(done))
    if skipped:
        print("Already present: " + ", ".join(skipped))
    print(f"Backup at {PATH}.bak")


main()

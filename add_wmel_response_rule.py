#!/usr/bin/env python
"""Patch the Snakefile for the wMel response figure.

  - rule rescue_wolbachia_followup declares the two tables the figure reads:
    wolbachia_de_checked.tsv.gz (already written, not declared before) and
    wolbachia_logcpm_long.tsv.gz (new in 09_wolbachia_followup.R)
  - adds fig_wmel_response to the _S plot registry

Writes Snakefile.bak_<date>_wmel first; safe to run twice.

Usage (repo root):  python add_wmel_response_rule.py
"""
import datetime
import shutil
import sys

SF = "Snakefile"
s = open(SF).read()
orig = s

OUT_OLD = '''        f"{RESCUE_DIR}/wolbachia_followup/wolbachia_candidates.tsv",
    params:
        indir = f"{RESCUE_DIR}/pseudobulk", outdir = f"{RESCUE_DIR}/wolbachia_followup",'''
OUT_NEW = '''        f"{RESCUE_DIR}/wolbachia_followup/wolbachia_candidates.tsv",
        f"{RESCUE_DIR}/wolbachia_followup/wolbachia_de_checked.tsv.gz",
        f"{RESCUE_DIR}/wolbachia_followup/wolbachia_logcpm_long.tsv.gz",
    params:
        indir = f"{RESCUE_DIR}/pseudobulk", outdir = f"{RESCUE_DIR}/wolbachia_followup",'''
if "wolbachia_logcpm_long.tsv.gz\",\n    params:" in s:
    print("  already applied: followup outputs")
elif OUT_OLD in s:
    s = s.replace(OUT_OLD, OUT_NEW, 1)
else:
    sys.exit("rule rescue_wolbachia_followup outputs not in the expected form; "
             "nothing written")

REG = '''    # wMel expression by host genotype: four volcanoes, the module matrix and
    # the replicates behind the paired claims. Both tables come from
    # rescue_wolbachia_followup.
    "fig_wmel_response": (
        "plot_wmel_response.py",
        f"{RESCUE_DIR}/wolbachia_followup/wolbachia_de_checked.tsv.gz",
        f"--logcpm {RESCUE_DIR}/wolbachia_followup/wolbachia_logcpm_long.tsv.gz"),
}

PLOTS = {**_V, **_S}'''
if '"fig_wmel_response": (' in s:
    print("  already applied: registry entry")
elif "}\n\nPLOTS = {**_V, **_S}" in s:
    s = s.replace("}\n\nPLOTS = {**_V, **_S}", REG, 1)
else:
    sys.exit("plot registry not in the expected form; nothing written")

if s == orig:
    print("nothing to change")
    sys.exit(0)
bak = f"{SF}.bak_{datetime.date.today():%Y%m%d}_wmel"
shutil.copy(SF, bak)
open(SF, "w").write(s)
print(f"wrote {SF} (backup {bak})")

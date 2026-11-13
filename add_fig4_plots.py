#!/usr/bin/env python3
"""Add the two new Figure 4 plot entries to the Snakefile registry.

Idempotent. fig4_meip26_expression writes one figure, so it goes in the
single-version registry; fig4_interaction_volcano writes robust and all, so it
goes in the two-version one. Both need extra CLI args, which is what the third
tuple field is for.
"""
import re
import shutil
import sys

PATH = sys.argv[1] if len(sys.argv) > 1 else "Snakefile"

V_ENTRY = '''    "fig4_interaction_volcano": (
        "plot_interaction_volcano.py",
        f"{RESCUE_DIR}/de/host_de_all.tsv.gz",
        f"--null {RESCUE_DIR}/rescue_null/rescue_null_summary.tsv"),
'''

S_ENTRY = '''    "fig4_meip26_expression": (
        "plot_meip26_expression.py",
        f"{RESCUE_DIR}/de/host_de_all.tsv.gz",
        f"--logcpm_dir {RESCUE_DIR}/de"),
'''


def insert_after_entry(text, anchor_key, entry, label):
    if entry.split('"')[1] + '":' in text:
        return text, f"{label} already present"
    k = text.find(f'    "{anchor_key}": (')
    if k == -1:
        sys.exit(f"anchor entry {anchor_key} not found")
    end = text.find("),\n", k)
    if end == -1:
        sys.exit(f"could not find the end of the {anchor_key} entry")
    at = end + 3
    return text[:at] + entry + text[at:], f"added {label}"


def main():
    src = open(PATH).read()
    out = src
    msgs = []
    # The two-version dict is keyed by stems that get _robust and _all.
    out, m = insert_after_entry(out, "fig3_composition_tests_propeller",
                                V_ENTRY, "fig4_interaction_volcano")
    msgs.append(m)
    # The single-version dict. fig5_host_de_summary is NOT in it, it is in _V,
    # so anchor on fig1_wmel_breadth, which is the first entry of _S.
    out, m = insert_after_entry(out, "fig1_wmel_breadth",
                                S_ENTRY, "fig4_meip26_expression")
    msgs.append(m)

    if out == src:
        print("; ".join(msgs))
        return
    shutil.copyfile(PATH, PATH + ".bak2")
    open(PATH, "w").write(out)
    print("; ".join(msgs))
    print(f"Backup at {PATH}.bak2")


main()

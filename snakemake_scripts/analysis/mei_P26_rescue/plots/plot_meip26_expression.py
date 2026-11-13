#!/usr/bin/env python
"""Figure 4: wMel bypasses mei-P26 rather than restoring it.

Reads the per-cell-type pseudobulk tables de/host_logcpm_<tissue>.tsv.gz and
the contrast table de/host_de_all.tsv.gz.

This panel is a control, and it answers the first objection any reader will
raise: does wMel rescue simply by restoring mei-P26 expression? It does not,
and the figure has to show three things at once for that answer to hold.

  A  mei-P26 logCPM per cell type, one point per sample, four groups. The
     gene is expressed in every somatic cell type and not only in the
     germline, which matters because the mei-P26[1] hypomorph is whole-animal:
     if the gene were germline-restricted, the soma would be genetically
     unaffected and every somatic result would have to be non-cell-autonomous.
     The same panel shows the hypomorph depleting the gene in every cell type,
     and shows that infection does not bring it back.
  B  the genotype-by-infection interaction per cell type. wMel does nothing to
     mei-P26 in OreR and lowers it further in the mutant, so the interaction
     is negative where it is significant. Rescue happens while the gene's
     expression falls, which is what makes it a bypass.

Read panel B's direction with the caveat in --help for --note_mixture: within
a cell-type label, pseudobulk shifts if the substage mixture shifts, and the
composition results say it does. The claim that survives that caveat is the
one the figure is for, that wMel does not restore mei-P26 expression.

Non-ovarian labels are dropped by default, not only the low-confidence ones.
Tracheal cells and oviduct are real cells in a dissected ovary rather than
artefacts, but they are not the tissue the claim is about, and adult tracheal
cell otherwise carries the largest interaction count in the dataset.

Usage:
  python plot_meip26_expression.py \
      --logcpm_dir results/rescue/de \
      --de results/rescue/de/host_de_all.tsv.gz \
      --outdir results/rescue/figures
"""
import argparse
import glob
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (COL, GROUP_BASE, GROUP_ORDER, LOW_CONFIDENCE, WIDE,
                      apply_style, ct_short, fdr_stars, note, panel_letter,
                      sample_group, save)
import matplotlib.pyplot as plt

NON_OVARIAN = ("adult tracheal cell", "oviduct")
NICE_GROUP = {"OreR_uninf": "OreR", "OreR_wMel": "OreR + wMel",
              "MeiP26_uninf": "mei-P26", "MeiP26_wMel": "mei-P26 + wMel"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--logcpm_dir", required=True,
                   help="directory holding host_logcpm_<tissue>.tsv.gz")
    p.add_argument("--de", required=True, help="de/host_de_all.tsv.gz")
    p.add_argument("--gene", default="FBgn0026206",
                   help="FlyBase id of the gene to plot (default mei-P26)")
    p.add_argument("--gene_label", default="mei-P26")
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig4_meip26_expression")
    p.add_argument("--alpha", type=float, default=0.05)
    p.add_argument("--drop_labels",
                   default=",".join(tuple(LOW_CONFIDENCE) + NON_OVARIAN),
                   help="labels excluded; '' keeps every one of them")
    p.add_argument("--note_mixture", action="store_true", default=True,
                   help="print the substage-mixture caveat under panel B. "
                        "Pseudobulk within a cell-type label moves if the "
                        "mixture of substages inside it moves, and the "
                        "composition results say it does, so the magnitude of "
                        "a within-label change is softer evidence than its "
                        "presence. The claim the panel carries, that wMel does "
                        "not restore mei-P26, does not depend on the magnitude")
    return p.parse_args()


def load_logcpm(d, gene):
    """One row per sample per cell type for the gene, from the per-type tables.

    The cell type comes from the column header (<sample>|<cell type>) rather
    than from the filename, which has been through a path-safe mangling.
    """
    rows = []
    files = sorted(glob.glob(os.path.join(d, "host_logcpm_*.tsv.gz")))
    if not files:
        raise SystemExit(f"no host_logcpm_*.tsv.gz under {d}")
    for f in files:
        t = pd.read_csv(f, sep="\t", index_col=0)
        if gene not in t.index:
            continue
        v = t.loc[gene]
        for col, val in v.items():
            sample, _, tissue = str(col).partition("|")
            rows.append({"sample": sample, "tissue": tissue or "unknown",
                         "logcpm": float(val), "group": sample_group(sample)})
    if not rows:
        raise SystemExit(f"{gene} is in none of the {len(files)} logCPM tables")
    return pd.DataFrame(rows)


def main():
    args = parse_args()
    apply_style()
    drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}

    d = load_logcpm(args.logcpm_dir, args.gene)
    d = d[~d["tissue"].isin(drop)]
    de = pd.read_csv(args.de, sep="\t")
    de = de[(de["gene_id"] == args.gene) & (~de["tissue"].isin(drop))]
    for c in ("logFC", "FDR"):
        de[c] = pd.to_numeric(de[c], errors="coerce")
    if not len(d):
        raise SystemExit("every cell type was dropped; check --drop_labels")
    print(f"{args.gene_label}: {d['tissue'].nunique()} cell types, "
          f"{len(d)} sample by cell-type values; "
          f"dropped {sorted(drop) or 'none'}")

    # Order by wild-type expression, so the germline enrichment is the shape of
    # the panel rather than something the reader has to hunt for.
    order = (d[d["group"] == "OreR_uninf"].groupby("tissue")["logcpm"].mean()
             .sort_values().index.tolist())
    order += [t for t in d["tissue"].unique() if t not in order]
    y = {t: i for i, t in enumerate(order)}

    fig = plt.figure(figsize=(WIDE, max(2.8, 0.30 * len(order) + 1.5)))
    gs = fig.add_gridspec(1, 2, wspace=0.10, width_ratios=[1.5, 1.0])

    # ---- A: expression per cell type and group ----
    ax = fig.add_subplot(gs[0, 0])
    off = {g: (i - 1.5) * 0.19 for i, g in enumerate(GROUP_ORDER)}
    for g in GROUP_ORDER:
        s = d[d["group"] == g]
        if not len(s):
            continue
        ax.scatter(s["logcpm"], [y[t] + off[g] for t in s["tissue"]], s=15,
                   color=GROUP_BASE[g], edgecolors=COL["black"], linewidths=0.3,
                   zorder=3, label=NICE_GROUP.get(g, g))
    for t in order:                       # a rule per cell type, to read along
        ax.axhline(y[t] + 0.5, color=COL["lightgrey"], lw=0.4, zorder=1)
    ax.set_yticks(range(len(order)), ct_short(order), fontsize=6)
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.set_xlabel(f"{args.gene_label} expression (log$_2$ CPM)")
    ax.set_title("Expressed in the soma, not restored by wMel", pad=4,
                 fontsize=7.5)
    ax.legend(loc="lower right", fontsize=5.5, frameon=True, framealpha=0.92,
              borderpad=0.4).get_frame().set_edgecolor("none")
    note(ax, "one point per sample, two per group", y=-0.17)
    panel_letter(ax, "A", dx=-0.40)

    # ---- B: the interaction ----
    ax = fig.add_subplot(gs[0, 1])
    it = de[de["contrast"] == "interaction"].set_index("tissue").reindex(order)
    sig = it["FDR"] < args.alpha
    ax.scatter(it["logFC"], range(len(order)), s=26, zorder=3,
               color=[COL["vermillion"] if s else "white" for s in sig.fillna(False)],
               edgecolors=COL["vermillion"], linewidths=0.8)
    for i, (_, r) in enumerate(it.iterrows()):
        st = fdr_stars(r["FDR"])
        if st:
            ax.text(r["logFC"], i + 0.22, st, ha="center", fontsize=5.5,
                    color=COL["vermillion"])
    ax.axvline(0, color=COL["black"], lw=0.7)
    ax.set_yticks(range(len(order)), [])
    ax.set_ylim(-0.7, len(order) - 0.3)
    ax.set_xlabel("interaction log$_2$FC")
    ax.set_title("Rescue happens while it falls", pad=4, fontsize=7.5)
    cap = f"filled points FDR < {args.alpha}"
    if args.note_mixture:
        cap += ("\nwithin-label pseudobulk moves with the substage mixture,\n"
                "so read the sign, not the size")
    note(ax, cap, y=-0.17)
    panel_letter(ax, "B", dx=-0.06)

    save(fig, args.outdir, args.stem)

    piv = (d.groupby(["tissue", "group"])["logcpm"].mean().unstack()
           .reindex(order)[[g for g in GROUP_ORDER]])
    print("\nMean log2 CPM by cell type and group:")
    print(piv.round(2).to_string())
    for con in ("genotype_in_uninf", "infection_in_OreR", "infection_in_MeiP26",
                "interaction"):
        s = de[de["contrast"] == con]
        n = int((s["FDR"] < args.alpha).sum())
        print(f"{con:<22} {n} of {len(s)} cell types at FDR < {args.alpha}")
    o = de[(de["contrast"] == "infection_in_OreR") & (de["FDR"] < args.alpha)]
    if not len(o):
        print(f"\nNo cell type changes {args.gene_label} with wMel in OreR, so the "
              "mutant-specific decrease is not a general infection response. "
              f"wMel does not rescue by restoring {args.gene_label}; the rescue "
              "runs around it.")


if __name__ == "__main__":
    main()

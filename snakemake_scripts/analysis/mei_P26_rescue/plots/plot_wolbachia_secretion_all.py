#!/usr/bin/env python
"""Figure 6: which wMel loci are both secretion-competent and up in a mutant host.

Plots results/rescue/wolbachia_followup/wolbachia_secretion_all.tsv.

The question the panel answers is narrow: a protein that substitutes for a
missing host factor has to reach the host and has to be no lower in the mutant.
Panel A puts those two axes against each other so the reader can see how few
loci satisfy both. Panel B splits the candidates by where they can act, which
matters more than the score: a Sec signal peptide means the bacterial envelope,
so those proteins cannot reach cytoplasmic host targets whatever their
statistics, and only T4SS delivery can, which is not predictable from sequence.

Usage:
  python plot_wolbachia_secretion_all.py \
      --input results/rescue/wolbachia_followup/wolbachia_secretion_all.tsv \
      --outdir figures/fig6
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import COL, WIDE, apply_style, panel_letter, save
import matplotlib.pyplot as plt

ROUTE_COL = {"host_cytoplasm_possible": COL["vermillion"],
             "membrane_or_T4SS": COL["orange"],
             "bacterial_surface": COL["blue"],
             "unknown": COL["grey"]}
TIER_ORDER = ["published", "predicted", "domain", "surface", "uncharacterised", ""]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig6_secretion_candidates")
    p.add_argument("--label_top", type=int, default=8,
                   help="how many candidates to name in panel A")
    p.add_argument("--de_table", default=None,
                   help="wolbachia_de_sig_checked.tsv. Given it, contrast counts are "
                        "recomputed over real cell types only, because the counts in "
                        "the secretion table include the pooled contrast and the "
                        "low-confidence labels")
    p.add_argument("--drop_labels", default="artefact,unannotated,pooled_all_tissues")
    p.add_argument("--version", choices=("robust", "all", "both"), default="both",
                   help="robust = composite flag only; all = every significant row; "
                        "both = write one figure of each (default)")
    return p.parse_args()


def recount(de_path, drop, robust_only):
    """Count up-contrasts over real cell types, optionally robust rows only.

    The counts shipped in the secretion table include the pooled contrast and
    the low-confidence labels, and (before the fix) came from a weaker filter
    than the composite `robust` flag. Recomputing here keeps the figure honest
    whichever table version it is handed.
    """
    d = pd.read_csv(de_path, sep="\t")
    if robust_only and "robust" in d.columns:
        d = d[d["robust"].astype(str).str.upper().isin(("TRUE", "T", "1", "YES"))]
    d = d[(~d["tissue"].isin(drop)) & (d["logFC"] > 0)]
    g = d.groupby("gene_id")
    return g.size().rename("n_up_contrasts"), g["logFC"].max().rename("best_up_logFC")


def num(s):
    return pd.to_numeric(s, errors="coerce")


def build(d, args, suffix, note):
    d = d.copy()
    for c in ("n_up_contrasts", "n_down_contrasts", "best_up_logFC", "best_up_FDR",
              "germline_cpm", "evidence_score", "de_score", "priority_score"):
        if c in d.columns:
            d[c] = num(d[c])
    if "n_up_contrasts" not in d.columns:
        raise SystemExit("this table predates the rescoring; rerun rule rescue_secretion "
                         "so n_up_contrasts, route and evidence_score exist")
    d["evidence_tier"] = d["evidence_tier"].fillna("")
    d["route"] = d.get("route", pd.Series("unknown", index=d.index)).fillna("unknown")
    d["label"] = d["old"].fillna(d["locus_tag"]).astype(str).str.split(",").str[0]
    counted = "contrast counts as shipped in the table"
    if args.de_table and os.path.exists(args.de_table):
        drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}
        n, best = recount(args.de_table, drop, robust_only=(suffix == "robust"))
        d = d.set_index("locus_tag")
        d["n_up_contrasts"] = n.reindex(d.index).fillna(0)
        d["best_up_logFC"] = best.reindex(d.index)
        d["de_score"] = np.where(d["n_up_contrasts"] > 0,
                                 6 + 2 * np.minimum(d["n_up_contrasts"], 4),
                                 d["de_score"].fillna(0))
        d = d.reset_index()
        counted = (f"{'robust rows only, ' if suffix == 'robust' else ''}"
                   f"cell types only, excluding {sorted(drop)}")
    print(f"\n[{suffix}] {len(d)} loci; {(d['n_up_contrasts'] > 0).sum()} up in at "
          f"least one contrast ({counted}); {note}")

    fig = plt.figure(figsize=(WIDE, 2.9))
    gs = fig.add_gridspec(1, 2, wspace=0.34, width_ratios=[1.35, 1.0])

    # ---- A: evidence against differential expression ----
    ax = fig.add_subplot(gs[0, 0])
    jit = np.random.default_rng(0).normal(0, 0.09, len(d))
    bg = d[d["n_up_contrasts"].fillna(0) == 0]
    ax.scatter(bg["evidence_score"] + jit[bg.index.map(d.index.get_loc)],
               bg["de_score"], s=6, color=COL["lightgrey"], edgecolors="none",
               label="not up", zorder=1)
    up = d[d["n_up_contrasts"].fillna(0) > 0]
    ax.scatter(up["evidence_score"] + jit[up.index.map(d.index.get_loc)],
               up["de_score"], s=26, zorder=3, edgecolors=COL["black"], linewidths=0.4,
               color=[ROUTE_COL.get(r, COL["grey"]) for r in up["route"]])
    named = up[up["evidence_tier"] != ""].nlargest(args.label_top, "priority_score")
    for _, r in named.iterrows():
        ax.annotate(r["label"], (r["evidence_score"], r["de_score"]),
                    textcoords="offset points", xytext=(5, 3), fontsize=6.5)
    ax.set_xlabel("evidence score (none → published)")
    ax.set_ylabel("differential expression score")
    ax.set_title("Both axes must be satisfied", pad=4)
    ax.margins(y=0.14)
    h = [plt.Line2D([], [], marker="o", ls="", color=c, markersize=5,
                    markeredgecolor=COL["black"], markeredgewidth=0.4)
         for c in (ROUTE_COL["host_cytoplasm_possible"], ROUTE_COL["membrane_or_T4SS"],
                   ROUTE_COL["bacterial_surface"])]
    ax.legend(h, ["cytoplasm possible", "membrane / T4SS", "surface only"],
              title="route", loc="lower right", fontsize=6, title_fontsize=6.5)
    panel_letter(ax, "A")

    # ---- B: the candidates that satisfy both ----
    ax = fig.add_subplot(gs[0, 1])
    cand = (d[(d["n_up_contrasts"] > 0) & (d["evidence_tier"] != "")]
            .sort_values("best_up_logFC"))
    if len(cand):
        y = np.arange(len(cand))
        ax.barh(y, cand["best_up_logFC"], height=0.6,
                color=[ROUTE_COL.get(r, COL["grey"]) for r in cand["route"]])
        for i, (_, r) in enumerate(cand.iterrows()):
            ax.text(r["best_up_logFC"] + 0.06, i,
                    f"{int(r['n_up_contrasts'])}×", va="center", fontsize=6)
        ax.set_yticks(y, [f"{r['label']}" for _, r in cand.iterrows()])
        ax.set_xlabel("best log$_2$FC, mei-P26 / OreR")
        ax.set_title("Secretion-competent and up", pad=4)
        ax.margins(x=0.22)
        ax.text(0.02, -0.34, "n× = independent contrasts", transform=ax.transAxes,
                fontsize=6, color=COL["grey"])
    else:
        ax.text(0.5, 0.5, "no locus satisfies both", ha="center", va="center",
                transform=ax.transAxes, color=COL["grey"])
        ax.set_axis_off()
    panel_letter(ax, "B")

    save(fig, args.outdir, f"{args.stem}_{suffix}")
    if len(cand):
        print("Candidates:")
        cols = [c for c in ("label", "product", "route", "evidence_tier",
                            "n_up_contrasts", "best_up_logFC", "best_up_FDR",
                            "germline_cpm") if c in cand.columns]
        print(cand[cols].sort_values("best_up_FDR").to_string(index=False))
        surf = (cand["route"] == "bacterial_surface").sum()
        print(f"{surf} of {len(cand)} are surface-only and cannot reach cytoplasmic "
              "host targets, whatever their statistics.")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    versions = (["robust", "all"] if args.version == "both" else [args.version])
    notes = {"robust": "robust contrasts only (survives renormalisation without "
                       "the dominant rRNA loci)",
             "all": "every significant contrast, robust or not"}
    for v in versions:
        build(d, args, v, notes[v])


if __name__ == "__main__":
    main()

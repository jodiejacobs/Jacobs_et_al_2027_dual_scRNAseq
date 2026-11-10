#!/usr/bin/env python
"""Figure 6: the robust wMel response to a mei-P26 host, gene by gene.

Plots results/rescue/wolbachia_followup/wolbachia_de_sig_checked.tsv.

Four panels, in the order the argument runs:

  A  reproducibility. Loci ranked by how many independent cell-type contrasts
     they are significant in, which at two replicates per group is better
     evidence than any single p-value. GroES at six contexts and the P44/Msp2
     protein WD1063 at six are the only loci that replicate this widely.
  B  the GroES and GroEL dissociation. Both are plotted because they sit in one
     operon 101 bp apart, so the reader should see them together and judge the
     stoichiometry claim rather than read two independent findings.
  C  functional modules. Chaperone and stress loci against respiration and
     provisioning loci, which is the coherent part of the response.
  D  the characterised effector repertoire, where everything falls except TomO.

Only rows flagged robust (replicate-consistent and leave-one-out sign-stable)
are drawn, and --all_rows lifts that. Positive logFC means higher in the
mei-P26 host throughout.

Usage:
  python plot_wolbachia_de_sig_checked.py \
      --input results/rescue/wolbachia_followup/wolbachia_de_sig_checked.tsv \
      --outdir figures/fig6
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (empty_figure, COL, COLUMN, NON_CELLTYPE, RRNA_PATTERN, WIDE,
                      add_robust_args, apply_style, ct_short, de_versions,
                      panel_letter, save, shorten)
import matplotlib.pyplot as plt

# Modules named by product, so the panels do not depend on a hand-kept locus list.
MODULES = {
    "chaperone / stress": r"chaperon|heat shock|GroES|GroEL|DnaK|HtpG|cold-shock|"
                          r"thioredoxin|superoxide|peroxid",
    "respiration / provisioning": r"NADH|oxidoreductase|cytochrome|ATP synthase|"
                                  r"transporter|permease|symporter",
    "translation": r"ribosomal protein|elongation factor|tRNA|aminoacyl",
    "surface": r"outer membrane|P44|Msp2|surface|lipoprotein",
}
EFFECTORS = {"CifA": r"incompatibility factor CifA", "CifB": r"incompatibility factor CifB",
             "TomO": r"host RNA manipulator|TomO", "WalE1": r"actin bundler",
             "ankyrin": r"ankyrin repeat", "collagen-like": r"collagen-like"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig6_wolbachia_de")
    p.add_argument("--top_n", type=int, default=14, help="loci in panel A")
    add_robust_args(p)
    return p.parse_args()


def truthy(s):
    return s.astype(str).str.upper().isin(("TRUE", "T", "1", "YES"))


def gene_label(r):
    g = str(r.get("gene") or "").strip()
    if g and g.lower() not in ("na", "nan", ""):
        return g
    old = str(r.get("old_locus_tag") or "").split(",")[0].strip()
    return old or str(r.name)


def build(d, args, suffix, note):
    d = d.copy()
    d["label"] = d.apply(gene_label, axis=1)
    d["product"] = d.get("product", pd.Series("", index=d.index)).fillna("")
    # rRNA loci carry 85% of the read budget and are not genes of interest.
    is_rrna = d["product"].str.contains(RRNA_PATTERN, case=False, na=False)
    if is_rrna.any():
        print(f"Excluding {int(is_rrna.sum())} rRNA rows from the gene panels")
        d = d[~is_rrna]
    cell = d[~d["tissue"].isin(NON_CELLTYPE)]
    print(f"\n[{suffix}] {len(d)} rows, {d['gene_id'].nunique()} loci, "
          f"{cell['tissue'].nunique()} cell-type contrasts")
    print(f"  {note}")

    fig = plt.figure(figsize=(WIDE, 5.6))
    gs = fig.add_gridspec(2, 2, hspace=0.75, wspace=0.42,
                          height_ratios=[1.15, 1.0])

    # ---- A: how many cell types does each locus replicate in ----
    ax = fig.add_subplot(gs[0, 0])
    rep = (cell.groupby(["gene_id", "label"])
           .agg(n=("tissue", "nunique"), mean_lfc=("logFC", "mean")).reset_index())
    rep = rep.sort_values(["n", "mean_lfc"], ascending=[True, True]).tail(args.top_n)
    y = np.arange(len(rep))
    colours = [COL["vermillion"] if v > 0 else COL["blue"] for v in rep["mean_lfc"]]
    ax.barh(y, rep["n"], color=colours, height=0.72)
    ax.set_yticks(y, shorten(rep["label"], 16))
    ax.set_xlabel("independent cell-type contrasts")
    ax.set_xlim(0, max(rep["n"]) + 0.6)
    ax.set_xticks(range(0, int(max(rep["n"])) + 1))
    ax.set_title("Replication across cell types", pad=4)
    h = [plt.Line2D([], [], color=COL["vermillion"], lw=4),
         plt.Line2D([], [], color=COL["blue"], lw=4)]
    ax.legend(h, ["up in mei-P26", "down"], loc="lower right")
    panel_letter(ax, "A")

    # ---- B: GroES against GroEL, one operon ----
    ax = fig.add_subplot(gs[0, 1])
    pair = d[d["product"].str.contains("GroES|GroEL|chaperonin", case=False, na=False)]
    if len(pair):
        pair = pair.copy()
        pair["which"] = np.where(pair["product"].str.contains("GroES", case=False),
                                 "GroES", "GroEL")
        order = sorted(pair["tissue"].unique(),
                       key=lambda t: -pair.loc[pair["tissue"] == t, "logFC"].max())
        x = np.arange(len(order))
        for i, (nm, c) in enumerate((("GroES", COL["vermillion"]), ("GroEL", COL["blue"]))):
            sub = pair[pair["which"] == nm].set_index("tissue").reindex(order)
            ax.bar(x + (i - 0.5) * 0.36, sub["logFC"].fillna(0), width=0.34,
                   color=c, label=nm)
        ax.axhline(0, color=COL["black"], lw=0.6)
        ax.set_xticks(x, ct_short(order), rotation=30, ha="right")
        ax.set_ylabel("log$_2$FC, mei-P26 / OreR")
        ax.set_title("One operon, opposite directions", pad=4)
        ax.legend(loc="upper right")
    panel_letter(ax, "B")

    # ---- C: functional modules ----
    ax = fig.add_subplot(gs[1, 0])
    rows = []
    for name, pat in MODULES.items():
        m = d[d["product"].str.contains(pat, case=False, na=False)]
        for _, r in m.iterrows():
            rows.append({"module": name, "logFC": r["logFC"], "label": r["label"]})
    mod = pd.DataFrame(rows)
    if len(mod):
        names = [n for n in MODULES if n in set(mod["module"])]
        for i, n in enumerate(names):
            v = mod.loc[mod["module"] == n, "logFC"].to_numpy()
            jitter = (np.random.default_rng(0).random(len(v)) - 0.5) * 0.3
            ax.scatter(v, np.full(len(v), i) + jitter, s=12, alpha=0.85,
                       color=[COL["vermillion"] if x > 0 else COL["blue"] for x in v],
                       edgecolors="none")
            ax.plot([np.median(v)] * 2, [i - 0.26, i + 0.26],
                    color=COL["black"], lw=1.2)
        ax.axvline(0, color=COL["grey"], lw=0.6, ls="--")
        ax.set_yticks(range(len(names)), [n.replace(" / ", "\n") for n in names])
        ax.set_xlabel("log$_2$FC, mei-P26 / OreR")
        ax.set_title("Functional modules (bar = median)", pad=4)
    panel_letter(ax, "C")

    # ---- D: characterised effectors ----
    ax = fig.add_subplot(gs[1, 1])
    rows = []
    for name, pat in EFFECTORS.items():
        m = d[d["product"].str.contains(pat, case=False, na=False)]
        for _, r in m.iterrows():
            rows.append({"effector": name, "logFC": r["logFC"],
                         "tissue": r["tissue"], "label": r["label"]})
    eff = pd.DataFrame(rows)
    if len(eff):
        agg = (eff.groupby("effector")
               .agg(lo=("logFC", "min"), hi=("logFC", "max"),
                    mid=("logFC", "mean"), n=("logFC", "size"))
               .sort_values("mid"))
        y = np.arange(len(agg))
        for i, (_, r) in enumerate(agg.iterrows()):
            ax.plot([r["lo"], r["hi"]], [i, i], color=COL["grey"], lw=1.0, zorder=1)
        ax.scatter(agg["mid"], y, s=18, zorder=2,
                   color=[COL["vermillion"] if v > 0 else COL["blue"] for v in agg["mid"]])
        ax.axvline(0, color=COL["black"], lw=0.6)
        ax.set_yticks(y, [f"{n} ({int(r)})" for n, r in zip(agg.index, agg["n"])])
        ax.set_xlabel("log$_2$FC, mei-P26 / OreR")
        ax.set_title("Characterised effectors", pad=4)
        ax.margins(x=0.18)
    panel_letter(ax, "D")

    fig.suptitle(f"wMel response to a mei-P26 host — {note}", y=1.005,
                 fontsize=6.5, color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}
    d = d[~d["tissue"].isin(drop)].copy()
    d["product"] = d.get("product", pd.Series("", index=d.index)).fillna("")
    is_rrna = d["product"].str.contains(RRNA_PATTERN, case=False, na=False)
    if is_rrna.any():
        print(f"Excluding {int(is_rrna.sum())} rRNA rows from the gene panels")
        d = d[~is_rrna]
    if not len(d):
        raise SystemExit("no rows left after filtering; try --drop_labels ''")
    print(f"Labels excluded as low-confidence: {sorted(drop) or 'none'}")
    for suffix, sub, note in de_versions(d, args.version):
        if not len(sub):
            empty_figure(args.outdir, f"{args.stem}_{suffix}",
                         f"No rows survive the {suffix} filter.\n"
                         f"An empty set here is a result, not a failure:\n"
                         f"say so in the text rather than citing the\n"
                         f"unfiltered figure as though it passed.")
            continue
        build(sub, args, suffix, note)
    print("\nPanel A ranks by replication across cell types, not by p-value: at two "
          "replicates per group, concordance across independent cell types is the "
          "stronger evidence. The pooled contrast is not counted as a cell type.")
    print("Panel B: GroES and GroEL sit 101 bp apart in one operon, so read the "
          "ratio rather than two independent effects.")


if __name__ == "__main__":
    main()

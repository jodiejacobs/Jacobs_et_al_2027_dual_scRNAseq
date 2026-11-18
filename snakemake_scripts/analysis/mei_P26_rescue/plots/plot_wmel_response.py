#!/usr/bin/env python
"""wMel expression by host genotype: volcanoes, module matrix, replicates.

Plots results/rescue/wolbachia_followup/wolbachia_de_checked.tsv.gz (every
tested locus, per cell type and pooled) with wolbachia_logcpm_long.tsv.gz
(per-library logCPM behind each fit), both from 09_wolbachia_followup.R.

Every contrast is wMel in a mei-P26 host against wMel in an OreR host, among
infected libraries only (two per genotype). Positive logFC means higher in the
mutant host. Nothing here says what wMel does to the host; it says how wMel's
own expression depends on which host it is in.

  A-D  one volcano per cell type with enough infected cells to test, titled
       with the number of loci tested. A cell type with nothing significant is
       drawn anyway: the response is local, and the empty panel shows that.
  E    module matrix. Every locus significant (FDR < --fdr) in any panel
       context, grouped by function, against the cell types and the pooled
       test. Colour is logFC, size is -log10 FDR, a black ring marks a robust
       call (replicate-consistent, leave-one-out sign-stable, significant
       under both normalisations, no zero libraries), a grey ring a call that
       is significant but not robust. A blank cell was not tested there.
  F-H  replicates behind the three paired claims: the two P44/Msp2 paralogs,
       GroES against GroEL (one operon), CifA against TomO. Each library's
       logCPM minus the OreR mean for that gene and context, so OreR points sit
       around zero and mei-P26 points sit at the fold change.

The pooled test sums the confidently labelled ovarian cell types per library
(09 drops artefact, unannotated, adult tracheal cell and oviduct). Those
labels are also dropped here. rRNA loci are left out of every panel: they
carry most of the wMel read budget and are not genes of interest.

Usage:
  python plot_wmel_response.py \
      --input results/rescue/wolbachia_followup/wolbachia_de_checked.tsv.gz \
      --logcpm results/rescue/wolbachia_followup/wolbachia_logcpm_long.tsv.gz \
      --outdir results/rescue/figures
"""
import argparse
import os
import re
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (COL, GENOTYPE, RRNA_PATTERN, WIDE, apply_style,
                      ct_short, note, panel_letter, save)
import matplotlib.pyplot as plt
from matplotlib.colors import LinearSegmentedColormap

POOLED = "pooled_all_tissues"
DROP = ("artefact", "unannotated", "adult tracheal cell", "oviduct")
CELL_TYPES = ["16-cell germline cyst in germarium region 2a and 2b",
              "central main body follicle cell ca. St. 6-8",
              "choriogenic main body follicle cell and corpus luteum",
              "polar follicle cell"]
# Function by product text, first match wins, so a locus is in one module.
# Effectors go first: "incompatibility" contains "ompa".
MODULES = [
    ("characterised effectors", r"incompatibility factor|host RNA manipulator|"
                                r"actin bundler|ankyrin|collagen-like"),
    ("surface", r"outer membrane|P44|Msp2|\bOmpA\b|surface"),
    ("chaperone / stress", r"chaperon|GroES|GroEL|DnaK|HtpG|cold-shock|"
                           r"thioredoxin|superoxide|peroxid|heat shock"),
    ("respiration / provisioning", r"NADH|oxidoreductase|cytochrome|ATP synthase|"
                                   r"transporter|permease|symporter|biotin"),
    ("translation", r"ribosomal protein|elongation factor|tRNA|aminoacyl"),
    ("other annotated", r"."),
]
PAIRS = [("WD_RS04815", "WD_RS00060", "Two P44/Msp2 paralogs, opposite ways"),
         ("WD_RS01365", "WD_RS01360", "GroES against GroEL, one operon"),
         ("WD_RS02835", "WD_RS05785", "CifA down, TomO up")]
CMAP = LinearSegmentedColormap.from_list(
    "lfc", [COL["blue"], "#FFFFFF", COL["vermillion"]])


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True, help="wolbachia_de_checked.tsv.gz")
    p.add_argument("--logcpm", required=True, help="wolbachia_logcpm_long.tsv.gz")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--lfc_max", type=float, default=4.0,
                   help="colour scale limit for panel E")
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig_wmel_response")
    return p.parse_args()


def truthy(s):
    return s.astype(str).str.upper().isin(("TRUE", "T", "1", "YES"))


def label(r):
    g = str(r.get("gene", "") or "")
    if g and g.lower() not in ("nan", "na"):
        return g
    prod = str(r.get("product", "") or "")
    short = {"P44/Msp2": "P44/Msp2", "GroES": "GroES", "TomO": "TomO",
             "cold-shock": "cold-shock", "thioredoxin": "thioredoxin",
             "superoxide": "SOD", "ankyrin": "ankyrin", "collagen-like":
             "collagen-like", "HU family": "HU", "Rpn family": "Rpn",
             "OmpA": "OmpA", "hypothetical": "hyp."}
    for k, v in short.items():
        if k.lower() in prod.lower():
            old = str(r.get("old_locus_tag", "") or "").split(",")[0]
            return f"{v} {old}" if old and old != "nan" else f"{v} {r['gene_id']}"
    return prod[:22] or r["gene_id"]


def module_of(prod):
    for name, pat in MODULES:
        if re.search(pat, str(prod), flags=re.I):
            return name
    return "other annotated"


def volcano(ax, s, t, fdr):
    s = s.copy()
    s["y"] = -np.log10(s["PValue"].clip(lower=1e-30))
    sig = s["FDR"] < fdr
    rob = sig & s["robust"]
    ax.scatter(s.loc[~sig, "logFC"], s.loc[~sig, "y"], s=7, color=COL["lightgrey"],
               edgecolors="none", zorder=1)
    for m, filled in ((sig & ~rob, False), (rob, True)):
        c = np.where(s.loc[m, "logFC"] > 0, COL["vermillion"], COL["blue"])
        if filled:
            ax.scatter(s.loc[m, "logFC"], s.loc[m, "y"], s=14, color=c,
                       edgecolors=COL["black"], linewidths=0.3, zorder=3)
        else:
            ax.scatter(s.loc[m, "logFC"], s.loc[m, "y"], s=14, facecolors="none",
                       edgecolors=c, linewidths=0.8, zorder=2)
    # Label the three strongest robust calls, stepping each label down so
    # neighbours at similar p do not print over one another.
    for k, (_, r) in enumerate(s[rob].nlargest(3, "y").iterrows()):
        right = r["logFC"] < 0
        ax.annotate(r["label"], (r["logFC"], r["y"]), fontsize=4.6,
                    xytext=(10 if right else -12, 6 - 10 * k), textcoords="offset points",
                    ha="left" if right else "right", va="center", zorder=4,
                    arrowprops=dict(arrowstyle="-", lw=0.4, color=COL["black"]))
    ax.set_ylim(-0.3, float(s["y"].max()) * 1.30 + 0.3)  # room for labels
    lim = max(2.0, float(s["logFC"].abs().max()) * 1.15)
    ax.set_xlim(-lim, lim)
    ax.axvline(0, color=COL["grey"], lw=0.5, ls="--", zorder=0)
    ax.set_title(f"{ct_short([t])[0]}\n{len(s)} loci tested\n{int(sig.sum())} "
                 f"significant, {int(rob.sum())} robust", fontsize=5.8, pad=3)
    ax.set_xlabel("log$_2$FC, mei-P26 / OreR host", fontsize=5.8)
    ax.tick_params(labelsize=5.5)


def matrix(ax, d, contexts, fdr, lfc_max):
    sig = d[d["FDR"] < fdr]
    keep = sig["gene_id"].unique()
    m = d[d["gene_id"].isin(keep)].copy()
    info = (m.drop_duplicates("gene_id").set_index("gene_id")
            [["label", "product"]])
    info["module"] = info["product"].map(module_of)
    info.loc[info["product"].str.contains("hypothetical", case=False, na=False),
             "module"] = "hypothetical"
    mean_lfc = m.groupby("gene_id")["logFC"].mean()
    order_mod = [n for n, _ in MODULES] + ["hypothetical"]
    info["mi"] = info["module"].map(order_mod.index)
    info["lfc"] = mean_lfc
    info = info.sort_values(["mi", "lfc"], ascending=[True, False])

    y, ticks, seps, ylab = 0, [], [], []
    pos = {}
    for mod, g in info.groupby("mi", sort=True):
        if ticks:
            y += 0.6
            seps.append(y - 0.8)
        ax.text(len(contexts) - 0.35, y - 0.05, order_mod[mod], fontsize=5.5,
                ha="left", va="top", color=COL["grey"], style="italic")
        for gid in g.index:
            pos[gid] = y
            ticks.append(y)
            ylab.append(info.loc[gid, "label"])
            y += 1
    for _, r in m.iterrows():
        if r["tissue"] not in contexts:
            continue
        x, yy = contexts.index(r["tissue"]), pos[r["gene_id"]]
        v = float(np.clip(r["logFC"], -lfc_max, lfc_max))
        c = CMAP((v + lfc_max) / (2 * lfc_max))
        if r["FDR"] < fdr:
            size = 10 + 14 * min(-np.log10(max(r["FDR"], 1e-6)), 6) / 2
            ax.scatter(x, yy, s=size, color=c, zorder=3, linewidths=0.9,
                       edgecolors=COL["black"] if r["robust"] else COL["grey"])
        else:
            ax.scatter(x, yy, s=6, color=c, zorder=2, edgecolors="none")
    ax.set_yticks(ticks, ylab, fontsize=5.3)
    ax.set_ylim(y - 0.4, -0.7)
    ax.set_xticks(range(len(contexts)),
                  ["pooled" if c == POOLED else ct_short([c])[0] for c in contexts],
                  rotation=45, ha="right", fontsize=5.6)
    ax.set_xlim(-0.6, len(contexts) - 0.4)
    for s_ in seps:
        ax.axhline(s_, color=COL["lightgrey"], lw=0.5)
    ax.grid(axis="x", color=COL["lightgrey"], lw=0.4, zorder=0)
    sm = plt.cm.ScalarMappable(cmap=CMAP, norm=plt.Normalize(-lfc_max, lfc_max))
    cax = ax.inset_axes([1.06, -0.02, 0.42, 0.012])
    cb = plt.colorbar(sm, cax=cax, orientation="horizontal")
    cb.set_label("log$_2$FC, mei-P26 / OreR host", fontsize=5.3)
    cb.ax.tick_params(labelsize=5)
    cb.outline.set_visible(False)
    return info


def pair_panel(ax, lc, d, a, b, title, fdr):
    rows = d[d["gene_id"].isin([a, b])]
    ctx = [POOLED] + [t for t in CELL_TYPES
                      if ((rows["tissue"] == t) & (rows["FDR"] < fdr)).any()]
    names = {g: rows.loc[rows["gene_id"] == g, "label"].iloc[0]
             for g in (a, b) if (rows["gene_id"] == g).any()}
    rng = np.random.default_rng(0)
    for i, t in enumerate(ctx):
        for k, (g, mk) in enumerate(((a, "o"), (b, "s"))):
            s = lc[(lc["tissue"] == t) & (lc["gene_id"] == g)]
            if not len(s):
                continue
            ref = s.loc[s["geno"] == "OreR", "logCPM"].mean()
            x0 = i + (k - 0.5) * 0.34
            for geno in ("OreR", "MeiP26"):
                v = s.loc[s["geno"] == geno, "logCPM"].to_numpy() - ref
                ax.scatter(x0 + rng.uniform(-0.05, 0.05, len(v)), v, s=14,
                           marker=mk, color=GENOTYPE[geno], edgecolors=COL["black"],
                           linewidths=0.3, zorder=3)
    ax.axhline(0, color=COL["grey"], lw=0.5, ls="--", zorder=0)
    ax.set_xticks(range(len(ctx)),
                  ["pooled" if c == POOLED else
                   ct_short([c])[0].replace(" ", "\n", 1) for c in ctx],
                  fontsize=5.2)
    ax.set_xlim(-0.6, len(ctx) - 0.4)
    ax.set_ylabel("log$_2$CPM − OreR mean", fontsize=5.6)
    ax.tick_params(labelsize=5.4)
    ax.set_title(title, fontsize=6.3, pad=3)
    h = [plt.Line2D([], [], marker="o", ls="", color=COL["black"], ms=3.5,
                    mfc="none", label=names.get(a, a)),
         plt.Line2D([], [], marker="s", ls="", color=COL["black"], ms=3.5,
                    mfc="none", label=names.get(b, b))]
    ax.legend(handles=h, loc="best", fontsize=5, frameon=False,
              handletextpad=0.2, borderaxespad=0.2)


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    lc = pd.read_csv(args.logcpm, sep="\t")
    d = d[~d["tissue"].isin(DROP)].copy()
    d["product"] = d.get("product", pd.Series("", index=d.index)).fillna("")
    d = d[~d["product"].str.contains(RRNA_PATTERN, case=False, na=False)]
    d["robust"] = ((d["FDR"] < args.fdr) & truthy(d["rep_consistent"])
                   & truthy(d["loo_sign_stable"])
                   & (d["FDR_excl_dominant"] < args.fdr)
                   & ((d["n_zero_MeiP26"] + d["n_zero_OreR"]) == 0))
    d["label"] = d.apply(label, axis=1)
    cts = [t for t in CELL_TYPES if t in set(d["tissue"])]
    contexts = cts + ([POOLED] if POOLED in set(d["tissue"]) else [])

    fig = plt.figure(figsize=(WIDE, 8.4))
    outer = fig.add_gridspec(2, 1, height_ratios=[1.0, 3.1], hspace=0.30)
    top = outer[0].subgridspec(1, 4, wspace=0.38)
    for i, t in enumerate(cts):
        ax = fig.add_subplot(top[0, i])
        volcano(ax, d[d["tissue"] == t], t, args.fdr)
        if i == 0:
            ax.set_ylabel("−log$_{10}$ p", fontsize=5.8)
        panel_letter(ax, "ABCD"[i], dx=-0.30, dy=1.30)

    bot = outer[1].subgridspec(3, 2, width_ratios=[1.25, 1.0], wspace=0.78,
                               hspace=0.70)
    ax = fig.add_subplot(bot[:, 0])
    info = matrix(ax, d, contexts, args.fdr, args.lfc_max)
    ax.set_title("Every locus significant anywhere, by function", fontsize=6.8,
                 pad=4)
    note(ax, "filled: FDR < {:g}; black ring robust, grey ring not;\nsmall dot: "
             "tested, not significant; blank: not tested there".format(args.fdr),
         y=-0.16)
    panel_letter(ax, "E", dx=-0.62, dy=1.02)
    for k, (a, b, title) in enumerate(PAIRS):
        ax = fig.add_subplot(bot[k, 1])
        pair_panel(ax, lc, d, a, b, title, args.fdr)
        panel_letter(ax, "FGH"[k], dx=-0.30, dy=1.16)
    fig.text(0.99, 0.005, "blue: OreR host, vermillion: mei-P26 host; one point "
             "per infected library", fontsize=5.5, color=COL["grey"], ha="right")
    save(fig, args.outdir, args.stem)

    print(info[["module", "label", "lfc"]].round(2).to_string())
    print("\nPer context, loci tested / FDR < {:g} / robust:".format(args.fdr))
    for t in contexts:
        s = d[d["tissue"] == t]
        print(f"  {t[:50]:<52}{len(s):>5}{int((s['FDR'] < args.fdr).sum()):>5}"
              f"{int(s['robust'].sum()):>5}")


if __name__ == "__main__":
    main()

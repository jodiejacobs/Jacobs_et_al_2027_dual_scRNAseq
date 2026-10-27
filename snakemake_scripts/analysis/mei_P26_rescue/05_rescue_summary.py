#!/usr/bin/env python
"""Summarise host DE as a rescue analysis and plot the edgeR results.

Rescue classification (per tissue), for genes dysregulated in the uninfected
mutant (genotype_in_uninf FDR < --fdr):

    rescue_fraction = 1 - logFC(genotype_in_wMel) / logFC(genotype_in_uninf)

    1 = wMel brings Mei-P26 fully back to OreR; 0 = no change; < 0 = wMel makes
    the mutant phenotype worse.

    full          rescue_fraction >= 0.75 and genotype_in_wMel FDR >= --fdr
    partial       rescue_fraction >= 0.25 (and not full)
    none          -0.25 < rescue_fraction < 0.25
    exacerbated   rescue_fraction <= -0.25

The 0.25 / 0.75 cutoffs are descriptive, not statistical. interaction_fdr is
reported next to each gene as the formal test that the infection response
differs between genotypes; with 2 replicates many visibly rescued genes will
not reach significance on the interaction alone.

Outputs: rescue_genes.tsv, rescue_category_counts.tsv, rescue_scatter.pdf,
de_counts_heatmap.pdf, volcano_<contrast>.pdf, wolbachia_volcano.pdf
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns

CAT_COLORS = {"full": "#1a9641", "partial": "#a6d96a", "none": "#bdbdbd", "exacerbated": "#d7191c"}
CONTRAST_ORDER = ["genotype_in_uninf", "genotype_in_wMel", "infection_in_OreR",
                  "infection_in_MeiP26", "infection_main", "interaction"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--de_dir", required=True, help="output dir of 03_edger_de.R")
    p.add_argument("--outdir", required=True)
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--n_label", type=int, default=8, help="genes labeled per volcano panel")
    return p.parse_args()


def classify(rf, fdr_wmel, fdr):
    if rf >= 0.75 and fdr_wmel >= fdr:
        return "full"
    if rf >= 0.25:
        return "partial"
    if rf > -0.25:
        return "none"
    return "exacerbated"


def panel_grid(n, ncol=4, size=3.6):
    ncol = min(ncol, n)
    nrow = int(np.ceil(n / ncol))
    fig, axes = plt.subplots(nrow, ncol, figsize=(size * ncol, size * nrow), squeeze=False)
    for ax in axes.ravel()[n:]:
        ax.axis("off")
    return fig, axes.ravel()


def label_of(row):
    return row["symbol"] if isinstance(row.get("symbol"), str) and row["symbol"] else row["gene_id"]


def volcano(ax, d, fdr, n_label, title):
    sig = d["FDR"] < fdr
    ax.scatter(d["logFC"], -np.log10(d["PValue"]), s=4, c=np.where(sig, "#d7301f", "#bdbdbd"),
               rasterized=True)
    for _, r in d[sig].nsmallest(n_label, "PValue").iterrows():
        ax.text(r["logFC"], -np.log10(r["PValue"]), label_of(r), fontsize=6)
    ax.axvline(0, c="k", lw=0.5)
    ax.set_title(f"{title}\n{int(sig.sum())} genes FDR<{fdr}", fontsize=9)
    ax.set_xlabel("log2FC", fontsize=8)
    ax.set_ylabel("-log10 p", fontsize=8)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    host = pd.read_csv(os.path.join(args.de_dir, "host_de_all.tsv.gz"), sep="\t")
    if "symbol" not in host.columns:
        host["symbol"] = np.nan

    # ---- rescue classification ----
    wide = host.pivot_table(index=["tissue", "gene_id"], columns="contrast",
                            values=["logFC", "FDR"]).reset_index()
    wide.columns = ["_".join(c).strip("_") if c[1] else c[0] for c in wide.columns]
    wide = wide.merge(host.drop_duplicates("gene_id")[["gene_id", "symbol"]], on="gene_id", how="left")

    mut = wide[wide["FDR_genotype_in_uninf"] < args.fdr].copy()
    mut["rescue_fraction"] = 1 - mut["logFC_genotype_in_wMel"] / mut["logFC_genotype_in_uninf"]
    mut["category"] = [classify(rf, f, args.fdr) for rf, f in
                       zip(mut["rescue_fraction"], mut["FDR_genotype_in_wMel"])]
    mut["direction_in_mutant"] = np.where(mut["logFC_genotype_in_uninf"] > 0, "up", "down")
    keep_cols = ["tissue", "gene_id", "symbol", "category", "rescue_fraction", "direction_in_mutant",
                 "logFC_genotype_in_uninf", "FDR_genotype_in_uninf",
                 "logFC_genotype_in_wMel", "FDR_genotype_in_wMel",
                 "logFC_interaction", "FDR_interaction",
                 "logFC_infection_in_OreR", "logFC_infection_in_MeiP26"]
    mut = mut[[c for c in keep_cols if c in mut.columns]].sort_values(["tissue", "category", "rescue_fraction"],
                                                                      ascending=[True, True, False])
    mut.to_csv(os.path.join(args.outdir, "rescue_genes.tsv"), sep="\t", index=False)

    cat_counts = (mut.groupby(["tissue", "category"]).size().unstack(fill_value=0)
                  .reindex(columns=list(CAT_COLORS), fill_value=0))
    cat_counts["n_mutant_DE"] = cat_counts.sum(axis=1)
    cat_counts["median_rescue_fraction"] = mut.groupby("tissue")["rescue_fraction"].median()
    cat_counts.to_csv(os.path.join(args.outdir, "rescue_category_counts.tsv"), sep="\t")
    print("Rescue categories (genes DE in uninfected Mei-P26 vs OreR):")
    print(cat_counts.to_string())

    # ---- rescue scatter: mutant effect without vs with wMel ----
    tissues = sorted(mut["tissue"].unique())
    if tissues:
        fig, axes = panel_grid(len(tissues))
        for ax, t in zip(axes, tissues):
            d = mut[mut["tissue"] == t]
            ax.scatter(d["logFC_genotype_in_uninf"], d["logFC_genotype_in_wMel"], s=10,
                       c=d["category"].map(CAT_COLORS))
            lim = np.nanmax(np.abs(d[["logFC_genotype_in_uninf", "logFC_genotype_in_wMel"]].to_numpy())) * 1.1
            ax.plot([-lim, lim], [-lim, lim], "k--", lw=0.6)   # no rescue
            ax.axhline(0, c="k", lw=0.6)                         # full rescue
            ax.set(xlim=(-lim, lim), ylim=(-lim, lim))
            for _, r in d.nlargest(min(args.n_label, len(d)), "rescue_fraction").iterrows():
                ax.text(r["logFC_genotype_in_uninf"], r["logFC_genotype_in_wMel"], label_of(r), fontsize=6)
            ax.set_title(f"{t} (n={len(d)})", fontsize=9)
            ax.set_xlabel("Mei-P26 vs OreR, uninfected (log2FC)", fontsize=7)
            ax.set_ylabel("Mei-P26 vs OreR, wMel (log2FC)", fontsize=7)
        handles = [plt.Line2D([], [], marker="o", ls="", c=c, label=k) for k, c in CAT_COLORS.items()]
        fig.legend(handles=handles, loc="center left", bbox_to_anchor=(1.0, 0.5), fontsize=8)
        fig.suptitle("Dashed = no rescue; horizontal = full rescue", fontsize=10)
        fig.tight_layout()
        fig.savefig(os.path.join(args.outdir, "rescue_scatter.pdf"), bbox_inches="tight")
        fig.savefig(os.path.join(args.outdir, "rescue_scatter.png"), dpi=200, bbox_inches="tight")
        plt.close(fig)

    # ---- DE count heatmap (tissue x contrast) ----
    summ = pd.read_csv(os.path.join(args.de_dir, "host_de_summary.tsv"), sep="\t")
    hm = summ.pivot(index="tissue", columns="contrast", values="n_sig")
    hm = hm[[c for c in CONTRAST_ORDER if c in hm.columns]]
    fig, ax = plt.subplots(figsize=(1.1 * hm.shape[1] + 3, 0.45 * hm.shape[0] + 1.5))
    sns.heatmap(hm, annot=True, fmt="d", cmap="viridis", cbar_kws={"label": f"genes FDR<{args.fdr}"}, ax=ax)
    ax.set(xlabel="", ylabel="", title="Host DE genes per tissue and contrast")
    plt.setp(ax.get_xticklabels(), rotation=40, ha="right")
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "de_counts_heatmap.pdf"))
    fig.savefig(os.path.join(args.outdir, "de_counts_heatmap.png"), dpi=200)
    plt.close(fig)

    # ---- volcanoes: one file per contrast, one panel per tissue ----
    for con in [c for c in CONTRAST_ORDER if c in host["contrast"].unique()]:
        d_con = host[host["contrast"] == con]
        ts = sorted(d_con["tissue"].unique())
        fig, axes = panel_grid(len(ts))
        for ax, t in zip(axes, ts):
            volcano(ax, d_con[d_con["tissue"] == t], args.fdr, args.n_label, t)
        fig.suptitle(con)
        fig.tight_layout()
        fig.savefig(os.path.join(args.outdir, f"volcano_{con}.pdf"))
        plt.close(fig)

    # ---- Wolbachia per-tissue genotype volcano ----
    wpath = os.path.join(args.de_dir, "wolbachia_de_all.tsv.gz")
    if os.path.exists(wpath):
        wol = pd.read_csv(wpath, sep="\t")
        wol = wol[wol["tissue"] != "all"].assign(symbol=np.nan)
        ts = sorted(wol["tissue"].unique())
        if ts:
            fig, axes = panel_grid(len(ts))
            for ax, t in zip(axes, ts):
                volcano(ax, wol[wol["tissue"] == t], args.fdr, args.n_label, t)
            fig.suptitle(f"Wolbachia genes: {wol['contrast'].iloc[0]}")
            fig.tight_layout()
            fig.savefig(os.path.join(args.outdir, "wolbachia_volcano.pdf"))
            plt.close(fig)


if __name__ == "__main__":
    main()

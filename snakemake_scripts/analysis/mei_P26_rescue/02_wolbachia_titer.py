#!/usr/bin/env python
"""Wolbachia titer across tissues: Mei-P26_wMel vs OreR_wMel.

Titer per cell = wMel UMIs / (wMel + host UMIs), as in filter_h5ad.

Three readouts per sample x tissue:
  pb_titer          pooled titer (sum wMel / sum total) -- the modeled quantity
  mean_cell_titer   average of per-cell titers
  frac_above_floor  fraction of cells whose titer exceeds the uninfected floor

Uninfected floor: the --floor_quantile of per-cell titer across ALL uninfected
cells. wMel UMIs in uninfected samples come from ambient RNA / index hopping
(all libraries share a lane), so a cell must beat this floor before it is
called infected. The floor is pooled across tissues because some tissues have
too few uninfected cells for a stable per-tissue quantile.

Model: quasi-binomial GLM on infected pseudobulks,
    cbind(wmel_umis, host_umis) ~ tissue * genotype
with dispersion estimated from Pearson chi2 (one shared dispersion across
tissues, which gives far more residual df than per-tissue 2-vs-2 tests).
The genotype effect within each tissue is the odds ratio Mei-P26 / OreR,
which approximates the titer fold change when titer is small.

Caveat: each sample contributes one pseudobulk per tissue, so observations
from the same sample are not independent. With 2 replicates per genotype a
sample random effect is not estimable, so treat p-values as approximate and
check that both replicates agree in the per-sample plot.
"""
import argparse
import os

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from patsy import dmatrix
from statsmodels.stats.multitest import multipletests

REF_GENO = "OreR"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pseudobulk_dir", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--min_cells", type=int, default=20,
                   help="minimum cells for a sample x tissue pseudobulk to be modeled")
    p.add_argument("--floor_quantile", type=float, default=0.99)
    return p.parse_args()


def genotype_by_tissue(meta):
    """Quasi-binomial tissue * genotype model; per-tissue genotype odds ratios."""
    d = meta.copy()
    d["tissue"] = pd.Categorical(d["tissue"])
    others = [g for g in d["genotype"].unique() if g != REF_GENO]
    d["genotype"] = pd.Categorical(d["genotype"], categories=[REF_GENO] + others)
    alt = others[0]
    tissues = list(d["tissue"].cat.categories)

    X = dmatrix("C(tissue) * C(genotype)", d, return_type="dataframe")
    names = list(X.columns)
    total = (d["wmel_umis"] + d["host_umis"]).to_numpy(float)
    prop = d["wmel_umis"].to_numpy(float) / total
    fit = sm.GLM(prop, X, family=sm.families.Binomial(), var_weights=total).fit(scale="X2")

    rows = []
    geno_term = f"C(genotype)[T.{alt}]"
    for t in tissues:
        c = np.zeros(len(names))
        c[names.index(geno_term)] = 1
        inter = f"C(tissue)[T.{t}]:{geno_term}"
        if inter in names:
            c[names.index(inter)] = 1
        tt = fit.t_test(c)
        est, se = float(np.ravel(tt.effect)[0]), float(np.ravel(tt.sd)[0])
        rows.append({
            "tissue": t, "contrast": f"{alt}_vs_{REF_GENO}",
            "log2_titer_ratio": est / np.log(2),
            "titer_ratio": np.exp(est),
            "ci_low": np.exp(est - 1.96 * se), "ci_high": np.exp(est + 1.96 * se),
            "pvalue": float(np.ravel(tt.pvalue)[0]),
        })
    res = pd.DataFrame(rows)
    res["fdr"] = multipletests(res["pvalue"], method="fdr_bh")[1]

    inter_idx = [i for i, n in enumerate(names) if ":" in n]
    if inter_idx:
        R = np.zeros((len(inter_idx), len(names)))
        for r, i in enumerate(inter_idx):
            R[r, i] = 1
        w = fit.wald_test(R, scalar=True)
        interaction = {"test": "genotype effect differs across tissues (tissue:genotype)",
                       "F": float(w.statistic), "df_num": len(inter_idx),
                       "df_denom": float(fit.df_resid), "pvalue": float(w.pvalue)}
    else:
        interaction = {}
    return res, interaction, fit.scale


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    meta = pd.read_csv(os.path.join(args.pseudobulk_dir, "pseudobulk_meta.tsv"), sep="\t")
    cells = ad.read_h5ad(os.path.join(args.pseudobulk_dir, "counts_labeled.h5ad"),
                         backed="r").obs.copy()
    for c in ["sample", "genotype", "infection_status", "tissue"]:
        cells[c] = cells[c].astype(str)

    # ---- uninfected floor and per-cell infection call ----
    uninf = cells["infection_status"] == "uninf"
    floor = float(np.quantile(cells.loc[uninf, "titer"], args.floor_quantile))
    cells["above_floor"] = cells["titer"] > floor
    print(f"Uninfected floor ({args.floor_quantile:.0%} quantile of uninfected per-cell titer): "
          f"{floor:.2e}")

    frac = (cells.groupby(["sample", "tissue"], observed=True)["above_floor"].mean()
            .rename("frac_above_floor").reset_index())
    summary = meta.merge(frac, on=["sample", "tissue"], how="left")
    summary["passes_min_cells"] = summary["n_cells"] >= args.min_cells
    summary.to_csv(os.path.join(args.outdir, "titer_by_sample_tissue.tsv"), sep="\t", index=False)

    floor_tbl = (summary[summary["infection_status"] == "uninf"]
                 .groupby("tissue")[["wmel_umis", "host_umis"]].sum())
    floor_tbl["uninf_pb_titer"] = floor_tbl["wmel_umis"] / floor_tbl.sum(axis=1)

    # ---- model on infected pseudobulks ----
    inf = summary[(summary["infection_status"] == "wMel") & summary["passes_min_cells"]]
    keep = inf.groupby("tissue")["genotype"].nunique()
    keep = keep[keep == inf["genotype"].nunique()].index
    dropped = sorted(set(inf["tissue"]) - set(keep))
    if dropped:
        print(f"Tissues without >= {args.min_cells} cells in both genotypes (not modeled): {dropped}")
    inf = inf[inf["tissue"].isin(keep)]

    res, interaction, scale = genotype_by_tissue(inf)
    res = res.merge(floor_tbl[["uninf_pb_titer"]], left_on="tissue", right_index=True, how="left")
    means = (inf.pivot_table(index="tissue", columns="genotype", values="pb_titer", aggfunc="mean")
             .add_prefix("mean_pb_titer_"))
    res = res.merge(means, left_on="tissue", right_index=True)
    res.sort_values("pvalue").to_csv(os.path.join(args.outdir, "titer_genotype_by_tissue.tsv"),
                                     sep="\t", index=False)
    with open(os.path.join(args.outdir, "titer_model_summary.txt"), "w") as fh:
        fh.write(f"uninfected_floor\t{floor:.4e}\n")
        fh.write(f"quasibinomial_dispersion\t{scale:.3f}\n")
        for k, v in interaction.items():
            fh.write(f"interaction_{k}\t{v}\n")
    print(res[["tissue", "titer_ratio", "ci_low", "ci_high", "pvalue", "fdr"]].to_string(index=False))
    if interaction:
        print(f"\ntissue:genotype F={interaction['F']:.2f}, p={interaction['pvalue']:.3g}")

    # ---- plots ----
    order = (inf.groupby("tissue")["pb_titer"].mean().sort_values(ascending=False).index.tolist())
    pal = {g: c for g, c in zip(sorted(inf["genotype"].unique(), key=lambda g: g != REF_GENO),
                                ["#4C72B0", "#DD8452"])}

    fig, axes = plt.subplots(1, 3, figsize=(17, 5))
    ax = axes[0]
    pc = cells[(cells["infection_status"] == "wMel") & cells["tissue"].isin(order)].copy()
    pc["log10_titer"] = np.log10(pc["titer"] + 1e-5)
    sns.violinplot(data=pc, x="tissue", y="log10_titer", hue="genotype", order=order, hue_order=list(pal),
                   palette=pal, cut=0, inner=None, linewidth=0.5, ax=ax)
    ax.axhline(np.log10(floor + 1e-5), ls="--", c="grey", lw=1)
    ax.set(title="Per-cell titer (dashed = uninfected floor)", ylabel="log10(titer + 1e-5)", xlabel="")

    ax = axes[1]
    sns.stripplot(data=inf, x="tissue", y="pb_titer", hue="genotype", order=order, hue_order=list(pal),
                  palette=pal, dodge=True, size=7, ax=ax)
    ax.scatter(range(len(order)), floor_tbl.reindex(order)["uninf_pb_titer"],
               marker="_", s=300, c="grey", label="uninfected")
    ax.set_yscale("log")
    ax.set(title="Pseudobulk titer per sample", ylabel="wMel / total UMIs", xlabel="")

    ax = axes[2]
    sns.stripplot(data=inf, x="tissue", y="frac_above_floor", hue="genotype", order=order, hue_order=list(pal),
                  palette=pal, dodge=True, size=7, ax=ax)
    ax.set(title="Fraction of cells above floor", ylabel="fraction", xlabel="", ylim=(0, 1.02))

    for ax in axes:
        ax.tick_params(axis="x", rotation=60)
        for lab in ax.get_xticklabels():
            lab.set_ha("right")
    for ax in axes[1:]:
        h, l = ax.get_legend_handles_labels()
        seen = dict(zip(l, h))
        ax.legend(seen.values(), seen.keys(), fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "titer_by_tissue.pdf"))
    fig.savefig(os.path.join(args.outdir, "titer_by_tissue.png"), dpi=200)

    fig, ax = plt.subplots(figsize=(5, 0.45 * len(res) + 1.2))
    r = res.set_index("tissue").reindex(order[::-1])
    ax.errorbar(r["titer_ratio"], range(len(r)),
                xerr=[r["titer_ratio"] - r["ci_low"], r["ci_high"] - r["titer_ratio"]],
                fmt="o", c="k", capsize=3)
    ax.axvline(1, ls="--", c="grey")
    ax.set_xscale("log")
    ax.set_yticks(range(len(r)), r.index)
    for i, (t, row) in enumerate(r.iterrows()):
        if row["fdr"] < 0.05:
            ax.text(row["ci_high"] * 1.05, i, "*", va="center")
    ax.set(xlabel=f"Titer ratio ({res['contrast'].iloc[0].replace('_vs_', ' / ')}), 95% CI",
           title="Genotype effect on titer by tissue")
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "titer_ratio_forest.pdf"))
    fig.savefig(os.path.join(args.outdir, "titer_ratio_forest.png"), dpi=200)


if __name__ == "__main__":
    main()

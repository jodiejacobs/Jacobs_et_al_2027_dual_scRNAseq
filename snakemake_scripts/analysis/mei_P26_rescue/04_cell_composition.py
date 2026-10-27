#!/usr/bin/env python
"""Cell-type composition: ~ genotype * infection, per tissue (propeller-style).

mei-P26 loss changes germline differentiation, so shifts in cell-type
proportions are part of the rescue phenotype and also matter when reading the
DE results (a tissue that expands in one condition contributes more cells to
its pseudobulk).

Per sample, proportions are logit-transformed (0.5-cell pseudocount) and each
tissue is fit with OLS using the same contrasts as the edgeR host model:
  infection_in_OreR, infection_in_MeiP26, genotype_in_uninf,
  genotype_in_wMel, interaction
With 2 replicates per group there are 4 residual df per tissue.

Proportions use the cells that passed --min_confidence in 01_build_pseudobulk,
so low-confidence labels are excluded from both numerator and denominator.
"""
import argparse
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
import statsmodels.api as sm
from patsy import dmatrix
from statsmodels.stats.multitest import multipletests

REF_GENO, REF_INF = "OreR", "uninf"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--pseudobulk_dir", required=True)
    p.add_argument("--outdir", required=True)
    return p.parse_args()


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    meta = pd.read_csv(os.path.join(args.pseudobulk_dir, "pseudobulk_meta.tsv"), sep="\t")

    counts = meta.pivot_table(index="sample", columns="tissue", values="n_cells",
                              aggfunc="sum", fill_value=0)
    samples = meta.drop_duplicates("sample").set_index("sample").loc[
        counts.index, ["genotype", "infection_status", "replicate"]]
    totals = counts.sum(axis=1)
    props = counts.div(totals, axis=0)
    logit = np.log((counts + 0.5).div(totals + 1, axis=0) /
                   (1 - (counts + 0.5).div(totals + 1, axis=0)))

    alt = [g for g in samples["genotype"].unique() if g != REF_GENO][0]
    alt_inf = [i for i in samples["infection_status"].unique() if i != REF_INF][0]
    d = samples.copy()
    d["genotype"] = pd.Categorical(d["genotype"], [REF_GENO, alt])
    d["infection_status"] = pd.Categorical(d["infection_status"], [REF_INF, alt_inf])
    X = dmatrix("C(genotype) * C(infection_status)", d, return_type="dataframe")
    names = list(X.columns)
    gi = names.index(f"C(genotype)[T.{alt}]")
    ii = names.index(f"C(infection_status)[T.{alt_inf}]")
    xi = [i for i, n in enumerate(names) if ":" in n][0]
    contrasts = {
        "infection_in_OreR": {ii: 1},
        f"infection_in_{alt.replace('-', '')}": {ii: 1, xi: 1},
        "genotype_in_uninf": {gi: 1},
        "genotype_in_wMel": {gi: 1, xi: 1},
        "interaction": {xi: 1},
    }

    rows = []
    for t in counts.columns:
        fit = sm.OLS(logit[t].to_numpy(), X).fit()
        for name, cdict in contrasts.items():
            c = np.zeros(len(names))
            for k, v in cdict.items():
                c[k] = v
            tt = fit.t_test(c)
            rows.append({"tissue": t, "contrast": name,
                         "logit_diff": float(np.ravel(tt.effect)[0]),
                         "se": float(np.ravel(tt.sd)[0]),
                         "pvalue": float(np.ravel(tt.pvalue)[0])})
    res = pd.DataFrame(rows)
    res["fdr"] = res.groupby("contrast")["pvalue"].transform(
        lambda p: multipletests(p, method="fdr_bh")[1])
    mean_props = props.groupby(samples["genotype"] + "_" + samples["infection_status"]).mean().T
    res = res.merge(mean_props.add_prefix("mean_prop_"), left_on="tissue", right_index=True)
    res.to_csv(os.path.join(args.outdir, "composition_tests.tsv"), sep="\t", index=False)
    props.assign(**samples).to_csv(os.path.join(args.outdir, "composition_proportions.tsv"), sep="\t")
    counts.to_csv(os.path.join(args.outdir, "composition_cell_counts.tsv"), sep="\t")
    sig = res[res["fdr"] < 0.05]
    print(sig.to_string(index=False) if len(sig) else "No composition changes at FDR < 0.05")

    # ---- plots ----
    order = props.mean().sort_values(ascending=False).index
    cond_order = [f"{g}_{i}" for g in [REF_GENO, alt] for i in [REF_INF, alt_inf]]
    srt = samples.assign(cond=samples["genotype"] + "_" + samples["infection_status"])
    srt["cond"] = pd.Categorical(srt["cond"], cond_order)
    srt = srt.sort_values(["cond", "replicate"])

    fig, axes = plt.subplots(1, 2, figsize=(15, 5.5), gridspec_kw={"width_ratios": [1, 1.6]})
    colors = sns.color_palette("tab20", len(order))
    props.loc[srt.index, order].plot.bar(stacked=True, ax=axes[0], color=colors, width=0.85, legend=False)
    axes[0].set(ylabel="Fraction of cells", xlabel="", title="Composition per sample")
    axes[0].legend(order, bbox_to_anchor=(1.0, 1), loc="upper left", fontsize=7)

    long = props.assign(cond=srt["cond"]).melt(id_vars="cond", var_name="tissue", value_name="prop")
    sns.stripplot(data=long, x="tissue", y="prop", hue="cond", order=order, hue_order=cond_order,
                  dodge=True, size=6, palette=["#9ecae1", "#3182bd", "#fdae6b", "#e6550d"], ax=axes[1])
    axes[1].set(ylabel="Fraction of cells", xlabel="", title="Per-tissue proportions")
    axes[1].tick_params(axis="x", rotation=60)
    for lab in axes[1].get_xticklabels():
        lab.set_ha("right")
    axes[1].legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "composition.pdf"))
    fig.savefig(os.path.join(args.outdir, "composition.png"), dpi=200)


if __name__ == "__main__":
    main()

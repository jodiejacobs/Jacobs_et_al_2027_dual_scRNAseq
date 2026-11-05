#!/usr/bin/env python
"""The specificity control: what wMel signal looks like where there is none.

Uninfected OreR and Mei-P26 libraries were prepared and sequenced alongside the
infected ones, so every wMel read in them is ambient RNA or index hopping. That
makes them the internal negative control for the whole method, and the strongest
validation available without a comparison run. This script quantifies it.

Four things, in the order a reviewer will ask for them:

  floor         per-cell titer quantiles in uninfected cells, pooled and per
                cell type, which is the threshold a cell must beat to be called
                infected. Per-cell-type floors matter because ambient uptake
                is not uniform across cell sizes.
  separation    AUROC of per-cell titer, infected against uninfected, computed
                within cell type so composition cannot inflate it, plus the
                fraction of infected cells above the pooled floor. AUROC near 1
                means titer separates infection status cell by cell; AUROC near
                0.5 in a cell type means titer is uninformative there, and that
                cell type should carry no titer claim.
  depth         whether uninfected wMel UMIs scale with droplet depth. Ambient
                contamination does (it is a fixed fraction of the soup);
                genuine low-level infection does not. A strong correlation here
                is the signature that the floor is ambient.
  profile       the wMel gene composition of the ambient signal against the
                infected signal. If ambient is a scaled copy of real expression
                the two profiles correlate tightly, which is both the expected
                result and the reason an ambient-driven titer difference cannot
                be diagnosed from the gene list alone.

Outputs (--outdir):
  uninfected_floor.tsv        pooled and per-cell-type floor quantiles
  titer_separation.tsv        AUROC and above-floor fraction per cell type
  ambient_depth.tsv           wMel vs host UMIs in uninfected cells, per sample
  ambient_profile.tsv.gz      per-gene CPM, ambient against infected
  uninfected_floor.{pdf,png}  the four panels

Usage:
  python 16_uninfected_floor.py \
      --counts results/rescue/pseudobulk_all/counts_labeled.h5ad \
      --outdir results/rescue/uninfected_floor
"""
import argparse
import os

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import seaborn as sns
from scipy import stats

WMEL_PREFIX = ("WD", "GQX67")
QUANTILES = (0.5, 0.9, 0.95, 0.99, 0.999, 1.0)


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts", required=True, help="counts_labeled.h5ad (cells x genes)")
    p.add_argument("--outdir", required=True)
    p.add_argument("--floor_quantile", type=float, default=0.99,
                   help="quantile of uninfected per-cell titer used as the floor")
    p.add_argument("--min_cells", type=int, default=20,
                   help="minimum cells per group for a per-cell-type statistic")
    return p.parse_args()


def auroc(a, b):
    """P(random infected > random uninfected), from the Mann-Whitney U."""
    if len(a) < 2 or len(b) < 2:
        return np.nan, np.nan
    u = stats.mannwhitneyu(a, b, alternative="two-sided")
    return float(u.statistic) / (len(a) * len(b)), float(u.pvalue)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    a = ad.read_h5ad(args.counts)
    obs = a.obs.copy()
    for c in ("sample", "genotype", "infection_status", "tissue"):
        obs[c] = obs[c].astype(str)
    uninf = (obs["infection_status"] == "uninf").to_numpy()
    print(f"{a.n_obs:,} cells: {int(uninf.sum()):,} uninfected, "
          f"{int((~uninf).sum()):,} infected")
    print(obs.groupby(["infection_status", "sample"], observed=True)
          .agg(cells=("titer", "size"), median_titer=("titer", "median"))
          .round(6).to_string())

    # ---- floor ----
    pooled = {f"q{q}": float(np.quantile(obs.loc[uninf, "titer"], q)) for q in QUANTILES}
    floor = float(np.quantile(obs.loc[uninf, "titer"], args.floor_quantile))
    rows = [{"group": "pooled_uninfected", "n_cells": int(uninf.sum()), **pooled}]
    for ct, g in obs[uninf].groupby("tissue", observed=True):
        if len(g) < args.min_cells:
            continue
        rows.append({"group": ct, "n_cells": len(g),
                     **{f"q{q}": float(np.quantile(g["titer"], q)) for q in QUANTILES}})
    floor_tbl = pd.DataFrame(rows)
    floor_tbl.to_csv(os.path.join(args.outdir, "uninfected_floor.tsv"),
                     sep="\t", index=False)
    print(f"\nPooled uninfected floor at q{args.floor_quantile} = {floor:.3e}")
    print(floor_tbl.head(12).to_string(index=False))

    # ---- separation, within cell type ----
    obs["above_floor"] = obs["titer"] > floor
    rows = []
    for ct, g in obs.groupby("tissue", observed=True):
        inf_t = g.loc[g["infection_status"] != "uninf", "titer"].to_numpy()
        un_t = g.loc[g["infection_status"] == "uninf", "titer"].to_numpy()
        if len(inf_t) < args.min_cells or len(un_t) < args.min_cells:
            continue
        au, p = auroc(inf_t, un_t)
        rows.append({"cell_type": ct, "n_infected": len(inf_t), "n_uninf": len(un_t),
                     "median_titer_infected": float(np.median(inf_t)),
                     "median_titer_uninf": float(np.median(un_t)),
                     "auroc": au, "pvalue": p,
                     "frac_infected_above_floor":
                         float((inf_t > floor).mean()),
                     "frac_uninf_above_floor": float((un_t > floor).mean())})
    sep = pd.DataFrame(rows).sort_values("auroc", ascending=False)
    sep.to_csv(os.path.join(args.outdir, "titer_separation.tsv"), sep="\t", index=False)
    print("\nPer-cell separation of infection status by titer (within cell type):")
    print(sep.round(4).to_string(index=False))
    weak = sep[sep["auroc"] < 0.7]
    if len(weak):
        print(f"\n{len(weak)} cell type(s) where titer barely separates infection status "
              f"(AUROC < 0.7): {', '.join(weak['cell_type'].astype(str))}. "
              "Titer claims in these cell types are not supported at the per-cell level.")

    # ---- ambient scales with depth ----
    rows = []
    for samp, g in obs[uninf].groupby("sample", observed=True):
        if len(g) < args.min_cells:
            continue
        r = stats.spearmanr(g["host_umis"], g["wmel_umis"])
        rows.append({"sample": samp, "n_cells": len(g),
                     "median_host_umis": float(g["host_umis"].median()),
                     "median_wmel_umis": float(g["wmel_umis"].median()),
                     "pb_titer": float(g["wmel_umis"].sum()
                                       / (g["wmel_umis"].sum() + g["host_umis"].sum())),
                     "spearman_wmel_vs_host_umis": float(r.statistic),
                     "pvalue": float(r.pvalue)})
    depth = pd.DataFrame(rows)
    depth.to_csv(os.path.join(args.outdir, "ambient_depth.tsv"), sep="\t", index=False)
    print("\nDoes uninfected wMel signal scale with droplet depth (ambient signature)?")
    print(depth.round(4).to_string(index=False))

    # ---- ambient gene profile ----
    is_w = np.asarray(pd.Index(a.var_names).str.startswith(WMEL_PREFIX))
    W = sp.csr_matrix(a[:, is_w].X)
    names = pd.Index(a.var_names[is_w])
    u_sum = np.asarray(W[uninf].sum(axis=0)).ravel()
    i_sum = np.asarray(W[~uninf].sum(axis=0)).ravel()
    prof = pd.DataFrame({
        "uninf_umis": u_sum, "infected_umis": i_sum,
        "uninf_cpm": u_sum / max(u_sum.sum(), 1) * 1e6,
        "infected_cpm": i_sum / max(i_sum.sum(), 1) * 1e6}, index=names)
    both = prof[(prof["uninf_umis"] > 0) & (prof["infected_umis"] > 0)]
    r = stats.spearmanr(np.log1p(both["uninf_cpm"]), np.log1p(both["infected_cpm"]))
    prof["log2_ambient_enrichment"] = np.log2((prof["uninf_cpm"] + 1)
                                              / (prof["infected_cpm"] + 1))
    prof.sort_values("infected_cpm", ascending=False).to_csv(
        os.path.join(args.outdir, "ambient_profile.tsv.gz"), sep="\t")
    print(f"\nAmbient vs infected wMel gene profile over {len(both)} loci: "
          f"Spearman rho = {r.statistic:.3f} (p = {r.pvalue:.3g})")
    print(f"  uninfected wMel UMIs: {int(u_sum.sum()):,}; "
          f"infected: {int(i_sum.sum()):,} "
          f"({u_sum.sum() / max(i_sum.sum(), 1) * 100:.2f}% as many)")
    if r.statistic > 0.7:
        print("  the ambient profile is close to a scaled copy of the infected profile, "
              "which is what free RNA in shared lanes produces")

    # ---- figures ----
    fig, axes = plt.subplots(2, 2, figsize=(13, 9))

    ax = axes[0, 0]
    pc = obs.copy()
    pc["log10_titer"] = np.log10(pc["titer"] + 1e-6)
    pc["status"] = np.where(pc["infection_status"] == "uninf", "uninfected", "wMel")
    sns.violinplot(data=pc, x="status", y="log10_titer", hue="status",
                   palette={"uninfected": "grey", "wMel": "#4C72B0"}, cut=0,
                   inner="quartile", legend=False, ax=ax)
    ax.axhline(np.log10(floor + 1e-6), ls="--", c="crimson", lw=1)
    ax.set(title=f"Per-cell titer (dashed = q{args.floor_quantile} uninfected floor)",
           ylabel="log10(titer + 1e-6)", xlabel="")

    ax = axes[0, 1]
    if len(sep):
        s = sep.sort_values("auroc")
        ax.barh(range(len(s)), s["auroc"], color="#4C72B0")
        ax.axvline(0.5, ls="--", c="grey", lw=1)
        ax.set_yticks(range(len(s)), [str(c)[:28] for c in s["cell_type"]], fontsize=7)
        ax.set(xlim=(0.4, 1.0), xlabel="AUROC, infected vs uninfected per-cell titer",
               title="Titer separates infection status, by cell type")

    ax = axes[1, 0]
    u = obs[uninf]
    ax.scatter(u["host_umis"], u["wmel_umis"] + 0.5, s=3, alpha=0.25, c="grey",
               label="uninfected")
    i = obs[~uninf].sample(min(20000, int((~uninf).sum())), random_state=0)
    ax.scatter(i["host_umis"], i["wmel_umis"] + 0.5, s=3, alpha=0.25, c="#4C72B0",
               label="wMel")
    ax.set(xscale="log", yscale="log", xlabel="host UMIs per cell",
           ylabel="wMel UMIs per cell + 0.5",
           title="Ambient signal scales with droplet depth")
    ax.legend(fontsize=8, markerscale=3)

    ax = axes[1, 1]
    ax.scatter(both["infected_cpm"], both["uninf_cpm"], s=6, alpha=0.5, c="k")
    lim = [1, max(both[["infected_cpm", "uninf_cpm"]].to_numpy().max() * 1.5, 10)]
    ax.plot(lim, lim, ls="--", c="grey", lw=1)
    ax.set(xscale="log", yscale="log", xlim=lim, ylim=lim,
           xlabel="CPM in infected cells", ylabel="CPM in uninfected cells",
           title=f"Ambient is a scaled copy (rho = {r.statistic:.2f})")

    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "uninfected_floor.pdf"))
    fig.savefig(os.path.join(args.outdir, "uninfected_floor.png"), dpi=200)
    plt.close(fig)
    print(f"\nWrote {args.outdir}")


if __name__ == "__main__":
    main()

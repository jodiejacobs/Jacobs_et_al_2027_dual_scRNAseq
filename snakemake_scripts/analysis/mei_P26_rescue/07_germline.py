#!/usr/bin/env python
"""Germline-focused rescue analysis, independent of atlas labels.

Atlas labels are unreliable for germline here: in Mei-P26, 'young germ cell'
is mostly tj+/eya+ somatic cells with no vas/orb/osk. This script defines
germline cells from germline-specific markers instead and asks whether wMel
moves the Mei-P26 germline toward OreR.

1. Germline call: a cell is germline if >= --min_markers of --markers
   (default vas, aub, AGO3, orb, osk) have nonzero counts. mei-P26 itself is
   not used (it is the mutated gene).
2. Composition: per-sample germline fraction with three denominators (all
   cells; cells not labelled 'artefact'; germline + tj+ somatic cells). Each is
   tested on logit fractions with OLS ~ 0 + group over the 8 samples (residual
   df = 4): genotype_in_uninf, genotype_in_wMel, infection_in_OreR,
   infection_in_MeiP26, interaction. Samples, not cells, are the replicates.
   Both test tables carry a BH 'fdr' column, corrected within contrast across
   metrics, and a 'clipped_groups' column. A metric with a group entirely at 0
   or 1 cannot be put on a logit scale, so it is estimated but not tested, and
   its se, t and p are blank: see clipped_groups() for why a clipped group
   otherwise manufactures significance.
3. Differentiation: germline cells are re-embedded (HVG, PCA, Leiden, UMAP)
   and ordered by diffusion pseudotime rooted in the cell with the highest
   early-germline score (bam, bgcn, Dad). Per-sample summaries (median
   pseudotime, fraction early, subcluster proportions) are tested as above.
4. Germline pseudobulks (all germline cells, and each subcluster with enough
   cells) written in 01_build_pseudobulk.py format, so 03_edger_de.R and
   08_rescue_null.R run on them unchanged with --pseudobulk_dir <outdir>/pseudobulk.

Usage:
  python 07_germline.py --counts results/rescue/pseudobulk/counts_labeled.h5ad \
      --annotation reference/fbgn_annotation_ID_fb_2025_04.tsv.gz \
      --outdir results/rescue/germline
"""
import argparse
import gzip
import os

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scanpy as sc
import scipy.sparse as sp
from scipy import stats
from statsmodels.stats.multitest import multipletests

HOST_PREFIX = "FB"
WMEL_PREFIX = ("WD", "GQX67")
GROUPS = ["OreR_uninf", "OreR_wMel", "MeiP26_uninf", "MeiP26_wMel"]
CONTRASTS = {   # weights on GROUPS
    "genotype_in_uninf":   [-1, 0, 1, 0],
    "genotype_in_wMel":    [0, -1, 0, 1],
    "infection_in_OreR":   [-1, 1, 0, 0],
    "infection_in_MeiP26": [0, 0, -1, 1],
    "interaction":         [1, -1, -1, 1],
}
COLORS = {"OreR_uninf": "#9ecae1", "OreR_wMel": "#08519c",
          "MeiP26_uninf": "#fdae6b", "MeiP26_wMel": "#a63603"}
EARLY = ["bam", "bgcn", "Dad"]
LATE = ["orb", "osk", "grk", "BicD", "egl"]


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts", required=True)
    p.add_argument("--annotation", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--markers", nargs="+", default=["vas", "aub", "AGO3", "orb", "osk"])
    p.add_argument("--min_markers", type=int, default=2)
    p.add_argument("--resolution", type=float, default=0.5)
    p.add_argument("--min_subcluster_cells", type=int, default=20,
                   help="per sample; subclusters below this in any sample get no pseudobulk")
    p.add_argument("--seed", type=int, default=0)
    return p.parse_args()


def load_symbols(path):
    opener = gzip.open if path.endswith(".gz") else open
    sym = {}
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                sym[parts[0]] = parts[2]
    return sym          # symbol -> FBgn


def sample_test(values, groups):
    """OLS on 8 per-sample values ~ 0 + group; t tests for CONTRASTS (df = n - 4)."""
    X = pd.get_dummies(pd.Categorical(groups, categories=GROUPS)).to_numpy(float)
    y = np.asarray(values, float)
    beta, *_ = np.linalg.lstsq(X, y, rcond=None)
    df = len(y) - X.shape[1]
    s2 = np.sum((y - X @ beta) ** 2) / df
    xtx_inv = np.linalg.inv(X.T @ X)
    rows = []
    for name, w in CONTRASTS.items():
        w = np.asarray(w, float)
        est = w @ beta
        se = np.sqrt(s2 * w @ xtx_inv @ w)
        t = est / se if se > 0 else np.nan
        rows.append({"contrast": name, "estimate": est, "se": se, "t": t,
                     "pvalue": 2 * stats.t.sf(abs(t), df) if np.isfinite(t) else np.nan, "df": df})
    return pd.DataFrame(rows)


def logit(p, eps=1e-4):
    p = np.clip(np.asarray(p, float), eps, 1 - eps)
    return np.log(p / (1 - p))


def detected(a, genes):
    X = a[:, genes].X
    X = X.toarray() if sp.issparse(X) else X
    return pd.DataFrame(X > 0, columns=genes, index=a.obs_names)


def per_sample_table(obs, value_cols):
    s = obs.groupby("sample", observed=True).agg(
        group=("group", "first"), replicate=("replicate", "first"), **{c: (c, "mean") for c in value_cols})
    return s.reset_index()


def clipped_groups(values, groups, eps=1e-4):
    """Groups whose every sample sits at a bound the logit cannot represent.

    A proportion of exactly 0 or 1 has no logit, so logit() clips it to eps.
    Every sample in such a group then carries the same invented number: the
    group contributes no variance, its contrasts are set by eps rather than by
    data, and the pooled residual variance is deflated for every other
    contrast on that metric. prop_c4 is the worked example. All four OreR
    samples hold exactly zero cells in cluster 4, which made
    infection_in_OreR come out as exactly 0.000 with p = 1.0000, made the
    interaction numerically identical to infection_in_MeiP26, and gave
    genotype_in_uninf a t of 68 against a standard error no real variation
    produced.

    Such a metric is still estimated, because "no wild-type cell is in this
    cluster" is worth recording, but it is not tested. The descriptive version
    belongs in germline_subclusters.tsv.
    """
    v = np.asarray(values, float)
    g = np.asarray(groups)
    return [str(k) for k in pd.unique(g)
            if ((g == k).any()
                and (np.all(v[g == k] <= 0) or np.all(v[g == k] >= 1)))]


def add_fdr(df):
    """BH within contrast, across metrics, as 04_cell_composition.py does.

    Within contrast rather than over the whole table: the five contrasts
    answer different questions and are not a family. Suppressed rows carry no
    p-value and are left out of the correction rather than counted in it.
    """
    df = df.copy()
    df["fdr"] = np.nan
    for _, idx in df.groupby("contrast").groups.items():
        p = df.loc[idx, "pvalue"]
        ok = p.notna()
        if ok.any():
            df.loc[p.index[ok], "fdr"] = multipletests(p[ok], method="fdr_bh")[1]
    return df


def test_metrics(samp, metrics, transform):
    out = []
    for m in metrics:
        bad = clipped_groups(samp[m], samp["group"]) if transform is logit else []
        v = transform(samp[m]) if transform else samp[m]
        t = sample_test(v, samp["group"])
        if bad:
            t[["se", "t", "pvalue"]] = np.nan
            print(f"  {m}: {', '.join(bad)} entirely at the bound, so the "
                  "estimate is set by eps and not by data. Estimate kept, test "
                  "suppressed.")
        t.insert(0, "metric", m)
        t["clipped_groups"] = ",".join(bad)
        out.append(t)
    return pd.concat(out, ignore_index=True)


def write_pseudobulk(full, germ_obs, outdir):
    """01_build_pseudobulk.py-format pseudobulks for 'germline' and each subcluster."""
    os.makedirs(outdir, exist_ok=True)
    host = full.var_names.str.startswith(HOST_PREFIX)
    wmel = full.var_names.str.startswith(WMEL_PREFIX)
    X = full[germ_obs.index].X.tocsr()
    pbs = {}
    for samp, idx in germ_obs.groupby("sample", observed=True).groups.items():
        pbs[(samp, "germline")] = idx
    for (samp, cl), idx in germ_obs.groupby(["sample", "pb_tissue"], observed=True).groups.items():
        if cl:
            pbs[(samp, cl)] = idx
    pos = pd.Series(np.arange(len(germ_obs)), index=germ_obs.index)
    hcols, wcols, meta = {}, {}, []
    for (samp, tis), idx in pbs.items():
        rows = pos[idx].to_numpy()
        sub = X[rows]
        h = np.asarray(sub[:, host].sum(axis=0)).ravel()
        w = np.asarray(sub[:, wmel].sum(axis=0)).ravel()
        o = germ_obs.loc[idx]
        pb = f"{samp}|{tis}"
        hcols[pb], wcols[pb] = h, w
        meta.append({"pb_id": pb, "sample": samp, "genotype": o["genotype"].iloc[0],
                     "infection_status": o["infection_status"].iloc[0], "replicate": o["replicate"].iloc[0],
                     "tissue": tis, "n_cells": len(o), "host_umis": h.sum(), "wmel_umis": w.sum(),
                     "mean_cell_titer": o["titer"].mean(), "median_cell_titer": o["titer"].median(),
                     "pb_titer": w.sum() / max(w.sum() + h.sum(), 1)})
    pd.DataFrame(hcols, index=full.var_names[host]).rename_axis("gene_id").to_csv(
        os.path.join(outdir, "pseudobulk_host_counts.tsv.gz"), sep="\t")
    pd.DataFrame(wcols, index=full.var_names[wmel]).rename_axis("gene_id").to_csv(
        os.path.join(outdir, "pseudobulk_wolbachia_counts.tsv.gz"), sep="\t")
    pd.DataFrame(meta).to_csv(os.path.join(outdir, "pseudobulk_meta.tsv"), sep="\t", index=False)
    return pd.DataFrame(meta)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    sc.settings.verbosity = 1
    sym2fb = load_symbols(args.annotation)

    full = ad.read_h5ad(args.counts)
    full.obs["group"] = (full.obs["genotype"].astype(str).str.replace("-", "") + "_" +
                         full.obs["infection_status"].astype(str))
    fb = lambda genes: [sym2fb[g] for g in genes if g in sym2fb and sym2fb[g] in full.var_names]

    # ---- 1. germline call ----
    mk = fb(args.markers)
    if len(mk) < args.min_markers:
        raise ValueError(f"only {len(mk)} germline markers found in matrix")
    n_mk = detected(full, mk).sum(axis=1)
    full.obs["is_germline"] = (n_mk >= args.min_markers).values
    tj = sym2fb.get("tj")
    full.obs["is_tj_soma"] = (~full.obs["is_germline"]) & \
        (detected(full, [tj])[tj].values if tj in full.var_names else False)
    print(f"Germline cells (>= {args.min_markers} of {len(mk)} markers): {full.obs['is_germline'].sum():,}")

    # ---- 2. composition ----
    o = full.obs
    comp = o.groupby("sample", observed=True).agg(
        group=("group", "first"), replicate=("replicate", "first"),
        n_cells=("is_germline", "size"), n_germline=("is_germline", "sum"),
        n_not_artefact=("tissue", lambda t: int((t != "artefact").sum())),
        n_tj_soma=("is_tj_soma", "sum")).reset_index()
    comp["frac_of_all"] = comp["n_germline"] / comp["n_cells"]
    comp["frac_of_not_artefact"] = comp["n_germline"] / comp["n_not_artefact"]
    comp["frac_of_germ_plus_tj_soma"] = comp["n_germline"] / (comp["n_germline"] + comp["n_tj_soma"])
    comp = comp.sort_values(["group", "replicate"])
    comp.to_csv(os.path.join(args.outdir, "germline_composition_by_sample.tsv"), sep="\t", index=False)
    fmetrics = ["frac_of_all", "frac_of_not_artefact", "frac_of_germ_plus_tj_soma"]
    ctest = add_fdr(test_metrics(comp, fmetrics, logit))
    ctest.to_csv(os.path.join(args.outdir, "germline_composition_tests.tsv"), sep="\t", index=False)
    print("\nGermline fraction by sample:")
    print(comp[["sample", "group", "n_germline"] + fmetrics].round(4).to_string(index=False))
    print("\nTests on logit(fraction), OLS over samples (df=4), BH within contrast:")
    print(ctest.round(4).to_string(index=False))

    # ---- 3. differentiation within germline ----
    g = full[full.obs["is_germline"].values, full.var_names.str.startswith(HOST_PREFIX)].copy()
    g.X = g.X.astype(np.float32)
    sc.pp.normalize_total(g, target_sum=1e4)
    sc.pp.log1p(g)
    early, late = fb(EARLY), fb(LATE)
    sc.tl.score_genes(g, early, score_name="early_score", random_state=args.seed)
    sc.tl.score_genes(g, late, score_name="late_score", random_state=args.seed)
    sc.pp.highly_variable_genes(g, n_top_genes=2000, flavor="seurat", batch_key="sample")
    sc.pp.pca(g, n_comps=20, mask_var="highly_variable", random_state=args.seed)
    sc.pp.neighbors(g, n_neighbors=15, random_state=args.seed)
    sc.tl.leiden(g, resolution=args.resolution, key_added="germ_cluster", flavor="igraph",
                 n_iterations=2, random_state=args.seed)
    sc.tl.umap(g, random_state=args.seed)
    sc.tl.diffmap(g)
    g.uns["iroot"] = int(np.argmax(g.obs["early_score"].to_numpy()))
    sc.tl.dpt(g)
    pt = g.obs["dpt_pseudotime"].replace(np.inf, np.nan)
    cut = np.nanquantile(pt[g.obs["group"] == "OreR_uninf"], 0.25)
    g.obs["early_by_pt"] = pt < cut

    # subcluster markers (symbols)
    fb2sym = {v: k for k, v in sym2fb.items()}
    sc.tl.rank_genes_groups(g, "germ_cluster", method="wilcoxon", n_genes=12)
    names = g.uns["rank_genes_groups"]["names"]
    ctab = []
    for cl in sorted(g.obs["germ_cluster"].unique(), key=int):
        oc = g.obs[g.obs["germ_cluster"] == cl]
        ctab.append({"germ_cluster": cl, "n_cells": len(oc),
                     "median_pseudotime": oc["dpt_pseudotime"].median(),
                     "mean_early_score": oc["early_score"].mean(), "mean_late_score": oc["late_score"].mean(),
                     **{f"frac_{k}": (oc["group"] == k).mean() for k in GROUPS},
                     "top_markers": ", ".join(fb2sym.get(x, x) for x in names[cl])})
    ctab = pd.DataFrame(ctab)
    ctab.to_csv(os.path.join(args.outdir, "germline_subclusters.tsv"), sep="\t", index=False)
    print("\nGermline subclusters:")
    print(ctab.drop(columns="top_markers").round(3).to_string(index=False))

    props = pd.crosstab(g.obs["sample"], g.obs["germ_cluster"], normalize="index")
    props.columns = [f"prop_c{c}" for c in props.columns]
    gobs = g.obs.assign(pt=pt)
    samp = per_sample_table(gobs, ["pt", "early_by_pt", "early_score", "late_score", "titer"])
    samp["median_pt"] = gobs.groupby("sample", observed=True)["pt"].median().reindex(samp["sample"]).values
    samp = samp.merge(props, left_on="sample", right_index=True).sort_values(["group", "replicate"])
    samp.to_csv(os.path.join(args.outdir, "germline_state_by_sample.tsv"), sep="\t", index=False)
    stest = add_fdr(pd.concat(
        [test_metrics(samp, ["median_pt", "early_score", "late_score"], None),
         test_metrics(samp, ["early_by_pt"] + list(props.columns), logit)],
        ignore_index=True))
    stest.to_csv(os.path.join(args.outdir, "germline_state_tests.tsv"), sep="\t", index=False)
    print("\nGermline state by sample:")
    print(samp.drop(columns=["pt"]).round(3).to_string(index=False))
    print("\nState tests (interaction + infection_in_MeiP26 rows):")
    print(stest[stest["contrast"].isin(["infection_in_MeiP26", "infection_in_OreR", "interaction"])]
          .round(4).to_string(index=False))
    # The rescue claim stands or falls on the interaction, so say outright
    # whether it survives correction instead of leaving it to be read off.
    si = stest[(stest["contrast"] == "interaction") & (stest["fdr"] < 0.05)]
    if len(si):
        print(f"\nInteraction at FDR < 0.05: {', '.join(si['metric'])}")
    else:
        print("\nNo germline state metric has an interaction at FDR < 0.05. The "
              "differentiation block is real and large, but these data do not "
              "show wMel relieving it. Report that, and do not let a composition "
              "result stand in for it.")

    # ---- 4. germline pseudobulks ----
    keep = (g.obs.groupby(["germ_cluster", "sample"], observed=True).size().unstack(fill_value=0)
            .min(axis=1) >= args.min_subcluster_cells)
    gobs = g.obs[["sample", "genotype", "infection_status", "replicate", "titer"]].copy()
    gobs["pb_tissue"] = [f"germline_c{c}" if keep.get(c, False) else "" for c in g.obs["germ_cluster"]]
    meta = write_pseudobulk(full, gobs, os.path.join(args.outdir, "pseudobulk"))
    print(f"\nGermline pseudobulks: {len(meta)} ({meta['tissue'].nunique()} 'tissues')")
    g.obs[["sample", "group", "germ_cluster", "dpt_pseudotime", "early_score", "late_score", "titer"]] \
        .to_csv(os.path.join(args.outdir, "germline_cells.tsv.gz"), sep="\t")

    # ---- figures ----
    fig, axes = plt.subplots(1, 3, figsize=(13, 3.8))
    for ax, m in zip(axes, fmetrics):
        for i, k in enumerate(GROUPS):
            v = comp.loc[comp["group"] == k, m]
            ax.scatter(np.full(len(v), i) + np.array([-0.06, 0.06])[: len(v)], v * 100, c=COLORS[k], s=40)
        ax.set_xticks(range(4), GROUPS, rotation=30, fontsize=7)
        q = ctest.query("metric == @m and contrast == 'interaction'")["fdr"].iloc[0]
        lab = "not estimable" if not np.isfinite(q) else f"interaction FDR = {q:.3g}"
        ax.set_title(f"germline % ({m.replace('frac_of_', 'of ')})\n{lab}", fontsize=8)
        ax.set_ylim(bottom=0)
    fig.tight_layout()
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(args.outdir, f"germline_composition.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)

    g.obs["dpt_plot"] = pt
    fig, axes = plt.subplots(1, 5, figsize=(26, 4.8))
    sc.pl.umap(g, color="germ_cluster", ax=axes[0], show=False, legend_loc="on data", title="germline subcluster")
    sc.pl.umap(g, color="group", ax=axes[1], show=False, palette=COLORS, title="condition")
    sc.pl.umap(g, color="dpt_plot", ax=axes[2], show=False, title="diffusion pseudotime")
    sc.pl.umap(g, color="early_score", ax=axes[3], show=False, title="early score (bam, bgcn, Dad)")
    sc.pl.umap(g, color="late_score", ax=axes[4], show=False, title="late score (orb, osk, grk, BicD, egl)")
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(args.outdir, f"germline_umap.{ext}"), dpi=150, bbox_inches="tight")
    plt.close(fig)

    fig, ax = plt.subplots(figsize=(5, 4))
    for s_, d in gobs.assign(pt=pt, group=g.obs["group"]).groupby("sample", observed=True):
        x = np.sort(d["pt"].dropna())
        ax.plot(x, np.arange(1, len(x) + 1) / len(x), c=COLORS[d["group"].iloc[0]], lw=1.2)
    ax.set(xlabel="diffusion pseudotime (germline)", ylabel="cumulative fraction of cells")
    ax.legend(handles=[plt.Line2D([], [], c=c, label=k) for k, c in COLORS.items()], fontsize=7)
    for ext in ("pdf", "png"):
        fig.savefig(os.path.join(args.outdir, f"germline_pseudotime_ecdf.{ext}"), dpi=200, bbox_inches="tight")
    plt.close(fig)


if __name__ == "__main__":
    main()

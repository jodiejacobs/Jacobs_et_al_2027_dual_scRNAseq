#!/usr/bin/env python
"""Per-cell Wolbachia titer against host NMF program usage, within cell type.

This is the analysis only dual single-cell data can do: separate host and
symbiont experiments give titer per sample and host state per cell, never both
in the same cell, so they cannot ask whether a cell carrying more wMel is a
cell using a different host program.

The question is asked WITHIN cell type and sample, so a correlation cannot be
produced by composition (more germ cells, lower titer) or by batch.

Three confounds are controlled explicitly, because each on its own would
manufacture a correlation:

  depth      program usage and titer both track UMI depth. Both are rank
             transformed and residualised on rank(log10 n_counts) inside each
             group, so the reported value is a partial Spearman given depth.
  ambient    wMel reads appear in cells that cannot host wMel. The same
             statistic is computed on the UNINFECTED samples, where every wMel
             read is ambient; that is the null this is tested against, not zero.
  n          the unit of replication is the sample, not the cell. Group-level
             correlations are averaged across cell types within a sample, so
             each program is tested on 4 infected and 4 uninfected samples.
             With 2 replicates per genotype nothing better is available, and a
             cell-level p-value here would be a measure of cell count.

Primary statistic per program: Welch t-test of mean Fisher z, infected samples
vs uninfected samples. A program passes only if its titer association is
stronger in infected material than in ambient material.

Outputs (--outdir):
  titer_program_correlations.tsv   one row per sample x cell type x program
  titer_program_tests.tsv          one row per program: infected vs ambient
  titer_program_heatmap.{pdf,png}  partial Spearman per cell type x program
  titer_program_forest.{pdf,png}   the per-program test, infected vs ambient

Usage:
  python 13_titer_vs_programs.py \
      --programs results/nmf_programs/adata_with_programs.h5ad \
      --outdir results/rescue/titer_programs
"""
import argparse
import os
import re

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import seaborn as sns
from scipy import stats
from statsmodels.stats.multitest import multipletests

PROGRAM_RE = re.compile(r"^Program_\d+$")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--programs", required=True,
                   help="h5ad with Program_* columns in .obs (nmf_programs output)")
    p.add_argument("--outdir", required=True)
    p.add_argument("--label_col", default="cell_type_annotation")
    p.add_argument("--titer_col", default="wolbachia_titer")
    p.add_argument("--depth_col", default="n_counts")
    p.add_argument("--min_cells", type=int, default=50,
                   help="minimum cells for a sample x cell type group to be correlated")
    p.add_argument("--min_samples", type=int, default=3,
                   help="minimum samples per infection group for a program to be tested")
    return p.parse_args()


def partial_spearman(x, y, z):
    """Spearman of x and y with z partialled out, on ranks."""
    rx, ry, rz = (stats.rankdata(v).astype(float) for v in (x, y, z))
    if np.ptp(rz) == 0:
        return float(stats.spearmanr(rx, ry).statistic)
    Z = np.column_stack([np.ones_like(rz), rz])
    ex = rx - Z @ np.linalg.lstsq(Z, rx, rcond=None)[0]
    ey = ry - Z @ np.linalg.lstsq(Z, ry, rcond=None)[0]
    if np.std(ex) == 0 or np.std(ey) == 0:
        return np.nan
    return float(np.corrcoef(ex, ey)[0, 1])


def build_obs(path, args):
    obs = ad.read_h5ad(path, backed="r").obs.copy()
    progs = sorted([c for c in obs.columns if PROGRAM_RE.match(c)],
                   key=lambda c: int(c.split("_")[1]))
    if not progs:
        raise SystemExit("no Program_* columns in .obs -- point --programs at the "
                         "nmf_programs output (adata_with_programs.h5ad)")
    for c in (args.label_col, args.titer_col, args.depth_col):
        if c not in obs.columns:
            raise SystemExit(f"missing .obs column {c!r}; present: {sorted(obs.columns)}")

    if "sample" not in obs.columns:
        parts = [obs[c].astype(str) for c in ("genotype", "infection_status", "replicate")
                 if c in obs.columns]
        if len(parts) < 2:
            raise SystemExit("cannot build a sample column; need genotype / "
                             "infection_status / replicate or a 'sample' column")
        obs["sample"] = parts[0].str.cat(parts[1:], sep="_")
    for c in ("sample", "infection_status", args.label_col):
        obs[c] = obs[c].astype(str)
    obs["cell_type"] = obs[args.label_col]
    obs["log_depth"] = np.log10(obs[args.depth_col].astype(float) + 1)
    obs["titer"] = obs[args.titer_col].astype(float)
    return obs, progs


def correlate(obs, progs, min_cells):
    rows = []
    for (samp, ct), g in obs.groupby(["sample", "cell_type"], observed=True):
        if len(g) < min_cells:
            continue
        t, d = g["titer"].to_numpy(), g["log_depth"].to_numpy()
        if np.ptp(t) == 0:
            continue
        inf = g["infection_status"].iloc[0]
        for p in progs:
            u = g[p].to_numpy(float)
            if np.ptp(u) == 0:
                continue
            r = partial_spearman(t, u, d)
            if np.isnan(r):
                continue
            rows.append({"sample": samp, "infection_status": inf, "cell_type": ct,
                         "program": p, "n_cells": len(g), "partial_spearman": r,
                         "fisher_z": np.arctanh(np.clip(r, -0.999999, 0.999999))})
    return pd.DataFrame(rows)


def test_programs(cor, min_samples):
    """Average z across cell types within sample, then infected vs uninfected."""
    per_sample = (cor.groupby(["program", "infection_status", "sample"], observed=True)
                  ["fisher_z"].mean().reset_index())
    rows = []
    for p, g in per_sample.groupby("program", observed=True):
        a = g.loc[g["infection_status"] != "uninf", "fisher_z"].to_numpy()
        b = g.loc[g["infection_status"] == "uninf", "fisher_z"].to_numpy()
        if len(a) < min_samples:
            continue
        rec = {"program": p, "n_infected_samples": len(a), "n_uninf_samples": len(b),
               "mean_z_infected": a.mean(),
               "mean_r_infected": np.tanh(a.mean()),
               "mean_z_uninf": b.mean() if len(b) else np.nan,
               "mean_r_uninf": np.tanh(b.mean()) if len(b) else np.nan}
        if len(b) >= min_samples:
            t = stats.ttest_ind(a, b, equal_var=False)
            rec.update(test="infected_vs_ambient", statistic=float(t.statistic),
                       df=float(t.df), pvalue=float(t.pvalue))
        else:
            t = stats.ttest_1samp(a, 0.0)
            rec.update(test="infected_vs_zero", statistic=float(t.statistic),
                       df=float(len(a) - 1), pvalue=float(t.pvalue))
        rows.append(rec)
    res = pd.DataFrame(rows)
    if len(res):
        res["fdr"] = multipletests(res["pvalue"], method="fdr_bh")[1]
        res = res.sort_values("pvalue")
    return res, per_sample


def plots(cor, res, outdir):
    inf = cor[cor["infection_status"] != "uninf"]
    if len(inf):
        m = (inf.pivot_table(index="cell_type", columns="program",
                             values="partial_spearman", aggfunc="mean")
             .reindex(columns=sorted(inf["program"].unique(),
                                     key=lambda c: int(c.split("_")[1]))))
        lim = float(np.nanmax(np.abs(m.to_numpy()))) or 0.1
        fig, ax = plt.subplots(figsize=(1 + 0.45 * m.shape[1], 1 + 0.32 * m.shape[0]))
        sns.heatmap(m, cmap="RdBu_r", center=0, vmin=-lim, vmax=lim, ax=ax,
                    cbar_kws={"label": "partial Spearman (titer, usage | depth)"})
        ax.set(title="Titer against program usage, infected samples", xlabel="", ylabel="")
        fig.tight_layout()
        fig.savefig(os.path.join(outdir, "titer_program_heatmap.pdf"))
        fig.savefig(os.path.join(outdir, "titer_program_heatmap.png"), dpi=200)
        plt.close(fig)

    if not len(res):
        return
    r = res.sort_values("program", key=lambda s: s.map(lambda c: int(c.split("_")[1])))
    y = np.arange(len(r))
    fig, ax = plt.subplots(figsize=(6, 0.4 * len(r) + 1.5))
    ax.scatter(r["mean_r_infected"], y, c="#4C72B0", label="infected", zorder=3)
    if r["mean_r_uninf"].notna().any():
        ax.scatter(r["mean_r_uninf"], y, c="grey", marker="s", s=25,
                   label="uninfected (ambient)", zorder=3)
    ax.axvline(0, ls="--", c="grey", lw=1)
    for i, (_, row) in enumerate(r.iterrows()):
        if row.get("fdr", 1) < 0.05:
            ax.text(max(row["mean_r_infected"], row.get("mean_r_uninf", -np.inf)) + 0.01,
                    i, "*", va="center")
    ax.set_yticks(y, r["program"])
    ax.set(xlabel="mean partial Spearman across cell types",
           title="Per-cell titer vs host program usage")
    ax.legend(fontsize=8)
    fig.tight_layout()
    fig.savefig(os.path.join(outdir, "titer_program_forest.pdf"))
    fig.savefig(os.path.join(outdir, "titer_program_forest.png"), dpi=200)
    plt.close(fig)


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    obs, progs = build_obs(args.programs, args)
    print(f"{len(obs):,} cells, {len(progs)} programs, "
          f"{obs['cell_type'].nunique()} cell types, {obs['sample'].nunique()} samples")
    print(obs.groupby(["infection_status", "sample"], observed=True).size().to_string())

    cor = correlate(obs, progs, args.min_cells)
    if not len(cor):
        raise SystemExit(f"no sample x cell type group reached --min_cells {args.min_cells}")
    cor.to_csv(os.path.join(args.outdir, "titer_program_correlations.tsv"),
               sep="\t", index=False)
    print(f"\n{len(cor):,} sample x cell type x program correlations "
          f"({cor['cell_type'].nunique()} cell types passed --min_cells)")

    res, per_sample = test_programs(cor, args.min_samples)
    if not len(res):
        raise SystemExit("no program had enough infected samples to test")
    res.to_csv(os.path.join(args.outdir, "titer_program_tests.tsv"), sep="\t", index=False)

    cols = ["program", "mean_r_infected", "mean_r_uninf", "test", "statistic", "df",
            "pvalue", "fdr"]
    print("\n" + res[cols].round(4).to_string(index=False))
    hit = res[res["fdr"] < 0.05]
    if len(hit):
        print(f"\n{len(hit)} program(s) with FDR < 0.05: "
              f"{', '.join(hit['program'])}")
    else:
        print(f"\nNothing passes FDR 0.05; best is {res.iloc[0]['program']} "
              f"(p = {res.iloc[0]['pvalue']:.3g}). With "
              f"{int(res.iloc[0]['n_infected_samples'])} infected samples this is a "
              f"bounded negative, not evidence of no relationship -- report the effect "
              f"sizes and the sample count.")

    plots(cor, res, args.outdir)
    print(f"\nWrote {args.outdir}")


if __name__ == "__main__":
    main()

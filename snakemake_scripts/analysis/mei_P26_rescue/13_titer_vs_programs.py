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
             transformed and residualised on depth inside each group, so the
             reported value is a partial Spearman given depth. The covariate
             defaults to a HOST-ONLY count column (qc_host_counts) rather than
             total UMIs: titer is wMel / (wMel + host), so total depth shares
             the wMel term with the quantity being correlated and partialling
             it out removes part of the signal along with the confound.
  pairing    the primary null is a within-group permutation of titer across
             cells, which keeps every group's size, depth structure and program
             distribution and destroys only the cell-level pairing. That is the
             null that answers "does it matter which cell the wMel is in".
  n          the unit of replication is the sample, not the cell. Group-level
             correlations are averaged across cell types within a sample, and
             the permutation statistic is aggregated the same way, so the test
             is on the 4-sample mean rather than on 170,000 cells.

The uninfected samples are reported as a secondary comparison but are NOT the
primary null, because uninfected per-cell titer is almost entirely zeros: very
little correlation could arise there whatever the biology, so "stronger in
infected than in ambient" is close to guaranteed and says less than it looks.

Read n_samples_same_sign before any p-value. A program whose correlation holds
its sign in all 4 infected samples across two genotypes is a result; a small
permutation p-value on its own only means the cell-level pairing beats noise,
which thousands of cells per group make easy.

Outputs (--outdir):
  titer_program_correlations.tsv   one row per sample x cell type x program,
                                   with the per-group permutation null
  titer_program_per_sample.tsv     mean Fisher z per program x sample
  titer_program_tests.tsv          one row per program: permutation test first,
                                   ambient comparison as secondary columns
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
    p.add_argument("--depth_col", default="auto",
                   help="covariate partialled out of both titer and usage. 'auto' prefers "
                        "host-only depth (qc_host_counts) over total depth, because titer "
                        "has wMel UMIs in its numerator and total depth therefore shares "
                        "that term")
    p.add_argument("--min_cells", type=int, default=50,
                   help="minimum cells for a sample x cell type group to be correlated")
    p.add_argument("--min_samples", type=int, default=3,
                   help="minimum samples per infection group for a program to be tested")
    p.add_argument("--permutations", type=int, default=2000,
                   help="within-group permutations of titer for the primary null. The "
                        "empirical p-value cannot go below 1/(permutations+1), and BH "
                        "across the programs multiplies that floor, so keep this well "
                        "above n_programs / target FDR")
    p.add_argument("--seed", type=int, default=1)
    return p.parse_args()


DEPTH_PREFERENCE = ("qc_host_counts", "qc_log10_host_counts", "n_counts")


def pick_depth_col(obs, requested):
    if requested != "auto":
        if requested not in obs.columns:
            raise SystemExit(f"--depth_col {requested!r} not in .obs")
        return requested
    for c in DEPTH_PREFERENCE:
        if c in obs.columns:
            return c
    raise SystemExit("no depth column found; pass --depth_col explicitly "
                     f"(tried {DEPTH_PREFERENCE})")


def _residualise(v, Z, P):
    """Residual of v on the columns of Z, using the precomputed projector P = pinv(Z)."""
    return v - Z @ (P @ v)


def _unit(v):
    s = v.std()
    return None if s == 0 else (v - v.mean()) / (s * np.sqrt(len(v)))


def group_partials(titer, usage, depth, n_perm, rng):
    """Partial Spearman of titer against every program, given depth, plus the
    within-group permutation null.

    Residualising on ranks makes the partial correlation a dot product of unit
    vectors, so one permutation costs one matrix multiply across all programs
    at once instead of a loop. Returns (observed r per program, (n_perm, P)
    null r) with NaN where a program or titer has no variance.
    """
    n, P = usage.shape
    rz = stats.rankdata(depth).astype(float)
    Z = (np.column_stack([np.ones(n), rz]) if np.ptp(rz) > 0
         else np.ones((n, 1)))
    proj = np.linalg.pinv(Z)

    U = np.empty((n, P))
    ok = np.zeros(P, bool)
    for j in range(P):
        u = _unit(_residualise(stats.rankdata(usage[:, j]).astype(float), Z, proj))
        if u is not None:
            U[:, j], ok[j] = u, True
    U[:, ~ok] = 0.0

    rt = stats.rankdata(titer).astype(float)
    t = _unit(_residualise(rt, Z, proj))
    obs = np.full(P, np.nan)
    if t is None:
        return obs, np.full((n_perm, P), np.nan)
    obs[ok] = (U.T @ t)[ok]

    null = np.full((n_perm, P), np.nan)
    for b in range(n_perm):
        tb = _unit(_residualise(rng.permutation(rt), Z, proj))
        if tb is not None:
            null[b, ok] = (U.T @ tb)[ok]
    return obs, null


def build_obs(path, args):
    obs = ad.read_h5ad(path, backed="r").obs.copy()
    progs = sorted([c for c in obs.columns if PROGRAM_RE.match(c)],
                   key=lambda c: int(c.split("_")[1]))
    if not progs:
        raise SystemExit("no Program_* columns in .obs -- point --programs at the "
                         "nmf_programs output (adata_with_programs.h5ad)")
    for c in (args.label_col, args.titer_col):
        if c not in obs.columns:
            raise SystemExit(f"missing .obs column {c!r}; present: {sorted(obs.columns)}")
    depth_col = pick_depth_col(obs, args.depth_col)
    print(f"depth covariate: {depth_col}"
          + ("" if depth_col.startswith("qc_host")
             else "  (total depth shares the wMel term with titer; prefer a host-only "
                  "count column if one exists)"))

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
    obs["depth"] = obs[depth_col].astype(float)
    obs["titer"] = obs[args.titer_col].astype(float)
    return obs, progs, depth_col


def z(r):
    return np.arctanh(np.clip(r, -0.999999, 0.999999))


def correlate(obs, progs, min_cells, n_perm, seed):
    """Observed partials per sample x cell type x program, plus the permutation
    null kept as a (n_perm, P) array per group for the aggregate test."""
    rng = np.random.default_rng(seed)
    rows, nulls = [], {}
    for (samp, ct), g in obs.groupby(["sample", "cell_type"], observed=True):
        if len(g) < min_cells or np.ptp(g["titer"].to_numpy()) == 0:
            continue
        obs_r, null_r = group_partials(g["titer"].to_numpy(float),
                                       g[progs].to_numpy(float),
                                       g["depth"].to_numpy(float), n_perm, rng)
        inf = g["infection_status"].iloc[0]
        nulls[(samp, ct)] = (inf, null_r)
        for j, p in enumerate(progs):
            if np.isnan(obs_r[j]):
                continue
            nn = null_r[:, j]
            nn = nn[np.isfinite(nn)]
            rows.append({
                "sample": samp, "infection_status": inf, "cell_type": ct,
                "program": p, "n_cells": len(g), "partial_spearman": obs_r[j],
                "fisher_z": z(obs_r[j]),
                "perm_null_sd": float(nn.std()) if len(nn) else np.nan,
                "perm_p_group": (float((1 + (np.abs(nn) >= abs(obs_r[j])).sum())
                                       / (len(nn) + 1)) if len(nn) else np.nan)})
    return pd.DataFrame(rows), nulls


def permutation_stat(nulls, progs, infected_only=True):
    """Aggregate each permutation the same way as the observed statistic:
    mean Fisher z across cell types within sample, then mean across samples."""
    by_sample = {}
    for (samp, _ct), (inf, null_r) in nulls.items():
        if infected_only and inf == "uninf":
            continue
        by_sample.setdefault(samp, []).append(z(null_r))
    if not by_sample:
        return None
    per_sample = [np.nanmean(np.stack(v), axis=0) for v in by_sample.values()]
    return np.nanmean(np.stack(per_sample), axis=0)  # (n_perm, P)


def test_programs(cor, nulls, progs, min_samples):
    """Average z across cell types within sample, then test two ways.

    Primary: the within-sample permutation null, which breaks the cell-level
    pairing of titer and usage while keeping every group's size, depth
    structure and program distribution. Secondary: the uninfected samples.
    The uninfected comparison is a weak null here because uninfected titer is
    almost all zeros, so very little correlation could arise there whatever the
    biology; it is reported for context and is not what the FDR is computed on.
    """
    per_sample = (cor.groupby(["program", "infection_status", "sample"], observed=True)
                  ["fisher_z"].mean().reset_index())
    null_stat = permutation_stat(nulls, progs)       # (n_perm, P) or None
    pidx = {p: j for j, p in enumerate(progs)}
    rows = []
    for p, g in per_sample.groupby("program", observed=True):
        a = g.loc[g["infection_status"] != "uninf", "fisher_z"].to_numpy()
        b = g.loc[g["infection_status"] == "uninf", "fisher_z"].to_numpy()
        if len(a) < min_samples:
            continue
        obs = a.mean()
        rec = {"program": p, "n_infected_samples": len(a), "n_uninf_samples": len(b),
               "mean_z_infected": obs, "mean_r_infected": np.tanh(obs),
               "n_samples_same_sign": int((np.sign(a) == np.sign(obs)).sum()),
               "mean_z_uninf": b.mean() if len(b) else np.nan,
               "mean_r_uninf": np.tanh(b.mean()) if len(b) else np.nan}
        if null_stat is not None:
            nn = null_stat[:, pidx[p]]
            nn = nn[np.isfinite(nn)]
            rec["perm_null_sd"] = float(nn.std()) if len(nn) else np.nan
            rec["perm_z_score"] = (obs / nn.std() if len(nn) and nn.std() > 0 else np.nan)
            rec["pvalue"] = float((1 + (np.abs(nn - nn.mean()) >= abs(obs - nn.mean())).sum())
                                  / (len(nn) + 1)) if len(nn) else np.nan
            rec["perm_p_floor"] = 1.0 / (len(nn) + 1) if len(nn) else np.nan
            rec["test"] = "within_sample_permutation"
            rec["n_permutations"] = int(len(nn))
        else:
            t = stats.ttest_1samp(a, 0.0)
            rec.update(test="infected_vs_zero", pvalue=float(t.pvalue))
        if len(b) >= min_samples:
            t = stats.ttest_ind(a, b, equal_var=False)
            rec["ambient_t"] = float(t.statistic)
            rec["ambient_pvalue"] = float(t.pvalue)
        rows.append(rec)
    res = pd.DataFrame(rows)
    if len(res):
        res["fdr"] = multipletests(res["pvalue"], method="fdr_bh")[1]
        if "ambient_pvalue" in res.columns and res["ambient_pvalue"].notna().all():
            res["ambient_fdr"] = multipletests(res["ambient_pvalue"], method="fdr_bh")[1]
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

    obs, progs, depth_col = build_obs(args.programs, args)
    print(f"{len(obs):,} cells, {len(progs)} programs, "
          f"{obs['cell_type'].nunique()} cell types, {obs['sample'].nunique()} samples")
    print(obs.groupby(["infection_status", "sample"], observed=True).size().to_string())

    print(f"\nCorrelating with {args.permutations} within-group permutations")
    cor, nulls = correlate(obs, progs, args.min_cells, args.permutations, args.seed)
    if not len(cor):
        raise SystemExit(f"no sample x cell type group reached --min_cells {args.min_cells}")
    cor.to_csv(os.path.join(args.outdir, "titer_program_correlations.tsv"),
               sep="\t", index=False)
    print(f"{len(cor):,} sample x cell type x program correlations "
          f"({cor['cell_type'].nunique()} cell types passed --min_cells)")

    res, per_sample = test_programs(cor, nulls, progs, args.min_samples)
    if not len(res):
        raise SystemExit("no program had enough infected samples to test")
    res.to_csv(os.path.join(args.outdir, "titer_program_tests.tsv"), sep="\t", index=False)
    per_sample.to_csv(os.path.join(args.outdir, "titer_program_per_sample.tsv"),
                      sep="\t", index=False)

    cols = [c for c in ["program", "mean_r_infected", "n_samples_same_sign",
                        "perm_null_sd", "perm_z_score", "pvalue", "fdr",
                        "mean_r_uninf", "ambient_pvalue"] if c in res.columns]
    print("\nPrimary test: within-sample permutation of titer. mean_r_uninf and "
          "ambient_pvalue are the secondary (weak) ambient comparison.")
    print(res[cols].round(4).to_string(index=False))
    hit = res[res["fdr"] < 0.05]
    if len(hit):
        print(f"\n{len(hit)} program(s) with FDR < 0.05 against the permutation null: "
              f"{', '.join(hit['program'])}")
        n_samp = int(res["n_infected_samples"].iloc[0])
        solid = hit[hit["n_samples_same_sign"] == n_samp]
        print(f"of those, {len(solid)} agree in sign across all {n_samp} infected samples"
              + (f": {', '.join(solid['program'])}" if len(solid) else ""))
        print("Sign agreement across samples is the load-bearing evidence; the "
              "permutation p-value only says the cell-level pairing is not noise.")
    elif ("perm_p_floor" in res.columns
          and (res["pvalue"] <= res["perm_p_floor"] + 1e-12).any()):
        n_at = int((res["pvalue"] <= res["perm_p_floor"] + 1e-12).sum())
        print(f"\n{n_at} program(s) hit the permutation floor "
              f"p = {res['perm_p_floor'].iloc[0]:.2g}, so the FDR cannot resolve them. "
              f"Rerun with --permutations above {len(res)} / 0.05 "
              f"= {int(len(res) / 0.05)} to tell them apart. Their perm_z_score values "
              "are not floored and are usable meanwhile.")
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

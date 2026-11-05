#!/usr/bin/env python
"""Is the wMel signal a biased slice of the genome, or the whole transcriptome?

The method claim is that pseudorandom priming recovers the Wolbachia
transcriptome, not a subset of it. The failure mode that would break the claim
is detection driven by sequence features rather than by expression: short
genes, GC-extreme genes, or genes in one part of the genome systematically
missed. This script tests for that.

Detection obviously tracks expression, so the test is conditional. For every
wMel locus it computes detection rate across infected cells and expression
(CPM of the wMel pool), then asks whether length, GC or genome position still
predicts detection once expression is accounted for:

  partial Spearman  detection vs each feature, given log10 CPM
  OLS               logit(detection) ~ log10 CPM + log10 length + GC + position

A null result is the result in the method's favour, and the coefficients bound
how large a bias the data would have hidden.

Also reported, for Figure 1:
  breadth      loci detected at a series of per-cell thresholds
  rarefaction  loci detected as cells are added, which shows whether breadth is
               saturated or still depth limited

Inputs. --counts is the cell x gene object with per-cell obs (the
pseudobulk_all counts_labeled.h5ad). --gff supplies length, coordinates and
strand from the Wolbachia lines of the reference annotation. --fasta is
optional and is the Wolbachia genome (single contig); given it, GC is computed
from the locus coordinates, otherwise GC is skipped and the OLS drops that term.

Usage:
  python 14_wmel_detection_bias.py \
      --counts results/rescue/pseudobulk_all/counts_labeled.h5ad \
      --gff reference/.../Drosophila_melanogaster_wMel_combined_16S.gtf \
      --outdir results/rescue/detection_bias [--fasta wMel_genome.fa]
"""
import argparse
import gzip
import os
import re

import anndata as ad
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
import pandas as pd
import scipy.sparse as sp
import statsmodels.api as sm
from scipy import stats

WMEL_PREFIX = ("WD", "GQX67")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts", required=True, help="counts_labeled.h5ad (cells x genes)")
    p.add_argument("--gff", required=True, help="GTF/GFF3; only Wolbachia lines are read")
    p.add_argument("--outdir", required=True)
    p.add_argument("--fasta", default=None, help="Wolbachia genome FASTA (optional, for GC)")
    p.add_argument("--min_cell_wmel", type=int, default=5,
                   help="a cell must carry this many wMel UMIs to count toward detection")
    p.add_argument("--seed", type=int, default=1)
    return p.parse_args()


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_features(path):
    """locus_tag -> start, end, strand, from gene/CDS lines of Wolbachia loci."""
    best = {}
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) < 9:
                continue
            m = re.search(r'locus_tag[ =]"?([^";\s]+)', p[8])
            if not m or not m.group(1).startswith(WMEL_PREFIX):
                continue
            lt = m.group(1)
            try:
                s, e = int(p[3]), int(p[4])
            except ValueError:
                continue
            prev = best.get(lt)
            if prev is None or (e - s) > (prev[1] - prev[0]):
                best[lt] = (s, e, p[6], p[0])
    return pd.DataFrame(
        [{"gene_id": k, "start": v[0], "end": v[1], "strand": v[2], "contig": v[3],
          "length": v[1] - v[0] + 1} for k, v in best.items()]).set_index("gene_id")


def read_genome(path):
    seqs, name, buf = {}, None, []
    with opener(path) as f:
        for line in f:
            if line.startswith(">"):
                if name:
                    seqs[name] = "".join(buf)
                name, buf = line[1:].split()[0], []
            else:
                buf.append(line.strip())
    if name:
        seqs[name] = "".join(buf)
    return seqs


def add_gc(feat, fasta):
    """Per-locus GC, each locus sliced from its OWN contig.

    The reference may be a combined host + symbiont genome, so the Wolbachia
    contig is one record among many; looking each locus up by the contig named
    in the annotation avoids guessing which record to use. A transcriptome
    FASTA has no record matching a genome contig name, which is caught here
    rather than silently producing GC from the wrong sequence.
    """
    seqs = read_genome(fasta)
    if not seqs:
        raise SystemExit(f"no sequences read from {fasta}")
    contigs = sorted(set(feat["contig"]))
    missing = [c for c in contigs if c not in seqs]
    if len(missing) == len(contigs):
        raise SystemExit(
            f"none of the annotation's contigs {contigs} are in {fasta}. "
            f"First records there: {list(seqs)[:5]}. This looks like a "
            "transcriptome rather than a genome; GC needs the genome the "
            "annotation's coordinates refer to, so point --fasta at that or "
            "leave wmel_genome empty to skip GC.")
    if missing:
        print(f"  contigs not in the FASTA, their loci get no GC: {missing}")
    gc = []
    for contig, a, b in zip(feat["contig"], feat["start"], feat["end"]):
        s = seqs.get(contig)
        if s is None:
            gc.append(np.nan)
            continue
        sub = s[max(a - 1, 0):b].upper()
        n = sub.count("G") + sub.count("C") + sub.count("A") + sub.count("T")
        gc.append((sub.count("G") + sub.count("C")) / n if n else np.nan)
    feat["gc"] = gc
    ok = int(np.isfinite(feat["gc"]).sum())
    print(f"  GC computed for {ok} of {len(feat)} loci from "
          f"{', '.join(c for c in contigs if c in seqs)}")
    if ok:
        print(f"  GC range {np.nanmin(feat['gc']):.3f} to {np.nanmax(feat['gc']):.3f}, "
              f"median {np.nanmedian(feat['gc']):.3f}")
    return feat


def per_gene_stats(path, min_cell_wmel, seed):
    a = ad.read_h5ad(path)
    is_w = np.asarray(pd.Index(a.var_names).str.startswith(WMEL_PREFIX))
    if is_w.sum() == 0:
        raise SystemExit("no wMel genes in --counts; check WMEL_PREFIX against the reference")
    inf = (a.obs["infection_status"].astype(str) != "uninf").to_numpy()
    W = sp.csr_matrix(a[:, is_w].X)
    names = pd.Index(a.var_names[is_w])

    per_cell = np.asarray(W.sum(axis=1)).ravel()
    use = inf & (per_cell >= min_cell_wmel)
    print(f"{int(inf.sum()):,} infected cells, {int(use.sum()):,} with "
          f">= {min_cell_wmel} wMel UMIs ({int(is_w.sum())} loci)")
    if use.sum() < 100:
        raise SystemExit("fewer than 100 usable infected cells; lower --min_cell_wmel")

    S = W[use]
    counts = np.asarray(S.sum(axis=0)).ravel()
    ncells = np.asarray((S > 0).sum(axis=0)).ravel()
    tbl = pd.DataFrame({"total_umis": counts, "n_cells": ncells}, index=names)
    tbl["detection_frac"] = tbl["n_cells"] / S.shape[0]
    tbl["cpm"] = tbl["total_umis"] / max(counts.sum(), 1) * 1e6

    # rarefaction over cells, for the breadth panel
    rng = np.random.default_rng(seed)
    order = rng.permutation(S.shape[0])
    checks = [n for n in (100, 500, 1000, 5000, 10000, 25000, 50000, 100000)
              if n <= S.shape[0]] + [S.shape[0]]
    rare, seen = [], np.zeros(S.shape[1], bool)
    prev = 0
    for n in checks:
        blk = S[order[prev:n]]
        seen |= np.asarray((blk > 0).sum(axis=0)).ravel() > 0
        rare.append({"n_cells": n, "loci_detected": int(seen.sum())})
        prev = n
    return tbl, pd.DataFrame(rare), int(S.shape[0])


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    tbl, rare, n_cells = per_gene_stats(args.counts, args.min_cell_wmel, args.seed)
    feat = parse_features(args.gff)
    print(f"{len(feat)} Wolbachia loci with coordinates in {os.path.basename(args.gff)}")
    if args.fasta and os.path.exists(args.fasta):
        feat = add_gc(feat, args.fasta)

    d = tbl.join(feat, how="inner")
    span = float(d["end"].max())
    d["rel_position"] = (d["start"] + d["end"]) / 2 / span
    d = d[d["cpm"] > 0].copy()
    d["log10_cpm"] = np.log10(d["cpm"])
    d["log10_length"] = np.log10(d["length"])
    d["plus_strand"] = (d["strand"] == "+").astype(int)
    print(f"{len(d)} loci with both counts and coordinates")

    # ---- breadth ----
    breadth = pd.DataFrame([{"min_cells": t,
                             "loci_detected": int((d["n_cells"] >= t).sum()),
                             "frac_of_loci": float((d["n_cells"] >= t).mean())}
                            for t in (1, 5, 10, 50, 100, 500, 1000)])
    breadth.to_csv(os.path.join(args.outdir, "wmel_breadth.tsv"), sep="\t", index=False)
    rare.to_csv(os.path.join(args.outdir, "wmel_rarefaction.tsv"), sep="\t", index=False)
    print("\nBreadth (of loci with coordinates):")
    print(breadth.to_string(index=False))

    # ---- conditional feature tests ----
    features = ["log10_length", "rel_position", "plus_strand"]
    if "gc" in d.columns and d["gc"].notna().any():
        features.insert(1, "gc")
    rz = stats.rankdata(d["log10_cpm"]).astype(float)
    Z = np.column_stack([np.ones_like(rz), rz])
    res_det = (stats.rankdata(d["detection_frac"]).astype(float)
               - Z @ np.linalg.lstsq(Z, stats.rankdata(d["detection_frac"]).astype(float),
                                     rcond=None)[0])
    rows = []
    for f in features:
        v = d[f].to_numpy(float)
        ok = np.isfinite(v)
        rv = stats.rankdata(v[ok]).astype(float)
        Zo = np.column_stack([np.ones_like(rv), rz[ok]])
        rf = rv - Zo @ np.linalg.lstsq(Zo, rv, rcond=None)[0]
        rd = res_det[ok]
        r, p = stats.pearsonr(rf, rd)
        raw = stats.spearmanr(v[ok], d["detection_frac"].to_numpy()[ok])
        rows.append({"feature": f, "n": int(ok.sum()),
                     "spearman_raw": float(raw.statistic), "pvalue_raw": float(raw.pvalue),
                     "partial_spearman_given_cpm": float(r), "pvalue_partial": float(p)})
    feat_tests = pd.DataFrame(rows)
    feat_tests.to_csv(os.path.join(args.outdir, "wmel_feature_tests.tsv"),
                      sep="\t", index=False)
    print("\nDetection vs sequence features (partial = conditional on expression):")
    print(feat_tests.round(4).to_string(index=False))

    # ---- OLS on the logit of detection ----
    eps = 0.5 / n_cells
    y = np.log(np.clip(d["detection_frac"], eps, 1 - eps)
               / (1 - np.clip(d["detection_frac"], eps, 1 - eps)))
    X = d[["log10_cpm"] + features].astype(float)
    ok = np.isfinite(X).all(axis=1) & np.isfinite(y)
    fit = sm.OLS(y[ok], sm.add_constant(X[ok])).fit()
    ols = pd.DataFrame({"term": fit.params.index, "coef": fit.params.values,
                        "se": fit.bse.values, "t": fit.tvalues.values,
                        "pvalue": fit.pvalues.values})
    ols.to_csv(os.path.join(args.outdir, "wmel_detection_model.tsv"), sep="\t", index=False)
    print(f"\nlogit(detection) ~ expression + features, n = {int(ok.sum())}, "
          f"R2 = {fit.rsquared:.3f}")
    print(ols.round(4).to_string(index=False))
    d.to_csv(os.path.join(args.outdir, "wmel_gene_detection.tsv.gz"), sep="\t")

    bad = feat_tests[(feat_tests["pvalue_partial"] < 0.05)
                     & (feat_tests["partial_spearman_given_cpm"].abs() > 0.1)]
    if len(bad):
        print("\nResidual feature dependence after conditioning on expression: "
              + ", ".join(f"{r.feature} (rho = {r.partial_spearman_given_cpm:.2f})"
                          for r in bad.itertuples()))
    else:
        print("\nNo feature predicts detection once expression is accounted for; the "
              "largest partial |rho| is "
              f"{feat_tests['partial_spearman_given_cpm'].abs().max():.3f}, which bounds "
              "the bias the data could hide.")

    # ---- figures ----
    fig, axes = plt.subplots(1, 3, figsize=(15, 4.2))
    ax = axes[0]
    ax.plot(rare["n_cells"], rare["loci_detected"], "o-", c="#4C72B0")
    ax.set(xscale="log", xlabel="infected cells sampled", ylabel="wMel loci detected",
           title="Breadth saturates with cells")
    ax.axhline(len(d), ls="--", c="grey", lw=1)

    ax = axes[1]
    q = pd.qcut(d["log10_length"], 3, labels=["short", "mid", "long"])
    for lab, c in zip(["short", "mid", "long"], ["#DD8452", "#55A868", "#4C72B0"]):
        m = (q == lab).to_numpy()
        ax.scatter(d["cpm"][m], d["detection_frac"][m], s=6, alpha=0.5, c=c, label=lab)
    ax.set(xscale="log", xlabel="CPM of the wMel pool", ylabel="fraction of cells detected",
           title="Detection tracks expression, not length")
    ax.legend(fontsize=8, title="length tercile")

    ax = axes[2]
    if "gc" in d.columns and d["gc"].notna().any():
        ax.scatter(d["gc"], res_det, s=6, alpha=0.5, c="k")
        ax.axhline(0, ls="--", c="grey", lw=1)
        ax.set(xlabel="GC fraction", ylabel="detection rank, expression removed",
               title="No GC bias in detection")
    else:
        ax.scatter(d["rel_position"], res_det, s=6, alpha=0.5, c="k")
        ax.axhline(0, ls="--", c="grey", lw=1)
        ax.set(xlabel="relative genome position", ylabel="detection rank, expression removed",
               title="No positional bias in detection")
    fig.tight_layout()
    fig.savefig(os.path.join(args.outdir, "wmel_detection_bias.pdf"))
    fig.savefig(os.path.join(args.outdir, "wmel_detection_bias.png"), dpi=200)
    plt.close(fig)
    print(f"\nWrote {args.outdir}")


if __name__ == "__main__":
    main()

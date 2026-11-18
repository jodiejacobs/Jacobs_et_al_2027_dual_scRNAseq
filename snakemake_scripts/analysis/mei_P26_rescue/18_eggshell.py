#!/usr/bin/env python
"""Per-cell eggshell program in follicle cells, by stage class.

The pseudobulk DE says eggshell genes are down in the uninfected mutant and up
with wMel in the mutant only. A pseudobulk mean cannot say why: more follicle
cells reaching the secreting stages raises it, and so does each cell expressing
more. Those are different biological claims, and this script separates them.

Per cell, for each stage class in plots/gene_sets.py (vitelline membrane,
yolk, chorion) and for all eggshell genes together:

  score = log1p(class UMIs / host UMIs x 1e4)

A cell is "on" for a class when its score exceeds the ambient floor of its own
library. The floor is the --floor_quantile of the same score in that library's
marker-defined germline cells (>= --min_markers of vas, aub, AGO3, orb, osk
detected, as in 07_germline.py). Germ cells do not make eggshell, so their
eggshell reads are ambient, and ambient RNA is shared by every droplet in a
library. The atlas label "young germ cell" is not used for this because in the
mutant it is mostly somatic cells (07_germline.py).

Per sample and follicle-cell type (and all follicle cells pooled):
  frac_on         fraction of cells above the floor: more cells at the stage
  score_on        mean score among on cells: more expression per cell
  score_all       mean score over every cell, the pseudobulk-like quantity

Each is tested over the 8 samples with OLS ~ 0 + group and the five contrasts
used throughout (residual df = 4), frac_on on the logit scale. BH is applied
within contrast and metric, across cell types and classes. A metric with a
group entirely at 0 or 1 is estimated but not tested, as in 07_germline.py.

Two controls:
  ambient   per sample, eggshell score in germline cells against follicle
            cells. If germline tracks follicle across libraries, eggshell
            reads in germline droplets are carry-over, which is what the
            Figure 4 volcano should treat them as.
  titer     within infected libraries, Spearman of per-cell score against
            per-cell wMel titer, per sample and cell type. Exploratory: titer
            and stage are confounded, since both change along oogenesis.

Outputs, in --outdir:
  eggshell_genes.tsv       the gene set as resolved against the matrix
  eggshell_by_sample.tsv   per sample x cell type x class summaries
  eggshell_tests.tsv       contrasts on those summaries
  eggshell_floor.tsv       the ambient floor per sample and class
  eggshell_ambient.tsv     germline against follicle, per sample and class
  eggshell_titer.tsv       per-cell score against titer, infected samples

Usage:
  python 18_eggshell.py --counts results/rescue/pseudobulk/counts_labeled.h5ad \
      --annotation reference/fbgn_annotation_ID_fb_2025_04.tsv.gz \
      --outdir results/rescue/eggshell
"""
import argparse
import gzip
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats
from statsmodels.stats.multitest import multipletests

sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "plots"))
from gene_sets import (EGGSHELL, EGGSHELL_ORDER, FBGN, FOLLICLE_KEY,
                       GERMLINE_MARKERS)

GROUPS = ["OreR_uninf", "OreR_wMel", "MeiP26_uninf", "MeiP26_wMel"]
CONTRASTS = {   # weights on GROUPS, identical to 07_germline.py
    "genotype_in_uninf":   [-1, 0, 1, 0],
    "genotype_in_wMel":    [0, -1, 0, 1],
    "infection_in_OreR":   [-1, 1, 0, 0],
    "infection_in_MeiP26": [0, 0, -1, 1],
    "interaction":         [1, -1, -1, 1],
}
POOLED = "all follicle cells"


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--counts", required=True, help="counts_labeled.h5ad")
    p.add_argument("--annotation", default=None,
                   help="FlyBase fbgn_annotation_ID table, to resolve symbols. "
                        "Without it the ids stored in plots/gene_sets.py are "
                        "used")
    p.add_argument("--outdir", required=True)
    p.add_argument("--floor_quantile", type=float, default=0.99)
    p.add_argument("--min_markers", type=int, default=2)
    p.add_argument("--min_floor_cells", type=int, default=50,
                   help="germline cells needed in a library for its own "
                        "floor; below this the floor pools every library")
    p.add_argument("--min_cells", type=int, default=30,
                   help="per sample; a follicle type below this in any "
                        "sample is not summarised")
    return p.parse_args()


def load_symbols(path):
    """symbol -> FBgn, same reader as 07_germline.py."""
    opener = gzip.open if path.endswith(".gz") else open
    sym = {}
    with opener(path, "rt") as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3:
                sym[parts[0]] = parts[2]
    return sym


def group_of(genotype, infection):
    return f"{str(genotype).replace('-', '')}_{infection}"


def logit(p, eps=1e-4):
    p = np.clip(np.asarray(p, float), eps, 1 - eps)
    return np.log(p / (1 - p))


def sample_test(values, groups):
    """OLS ~ 0 + group on per-sample values; t tests on CONTRASTS."""
    X = pd.get_dummies(pd.Categorical(groups, categories=GROUPS)).to_numpy(float)
    y = np.asarray(values, float)
    if not np.all(np.isfinite(y)):
        return pd.DataFrame([{"contrast": c, "estimate": np.nan, "se": np.nan,
                              "t": np.nan, "pvalue": np.nan, "df": np.nan}
                             for c in CONTRASTS])
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
                     "pvalue": 2 * stats.t.sf(abs(t), df) if np.isfinite(t) else np.nan,
                     "df": df})
    return pd.DataFrame(rows)


def clipped_groups(values, groups):
    v = np.asarray(values, float)
    g = np.asarray(groups)
    return [str(k) for k in pd.unique(g)
            if (g == k).any() and (np.all(v[g == k] <= 0) or np.all(v[g == k] >= 1))]


def main():
    args = parse_args()
    import anndata as ad
    import scipy.sparse as sp
    os.makedirs(args.outdir, exist_ok=True)

    a = ad.read_h5ad(args.counts)
    o = a.obs
    print(f"{a.n_obs:,} cells x {a.n_vars:,} genes")
    sym2fb = load_symbols(args.annotation) if args.annotation else {}
    var = pd.Index(a.var_names.astype(str))

    def resolve(symbol, fallback):
        fb = sym2fb.get(symbol, fallback)
        return fb if fb in var else None

    rows = []
    for c in EGGSHELL_ORDER:
        for s in EGGSHELL[c]:
            fb = resolve(s, FBGN.get(s))
            rows.append({"symbol": s, "gene_id": fb or FBGN.get(s, ""),
                         "class": c, "in_matrix": fb is not None})
    genes = pd.DataFrame(rows)
    genes.to_csv(os.path.join(args.outdir, "eggshell_genes.tsv"), sep="\t",
                 index=False)
    miss = genes.loc[~genes["in_matrix"], "symbol"].tolist()
    print(f"eggshell genes: {int(genes['in_matrix'].sum())} of {len(genes)} "
          f"in the matrix" + (f"; missing {miss}" if miss else ""))

    X = a.X.tocsc() if sp.issparse(a.X) else sp.csc_matrix(a.X)

    def col_sum(ids):
        idx = [var.get_loc(i) for i in ids]
        return np.asarray(X[:, idx].sum(axis=1)).ravel() if idx else np.zeros(a.n_obs)

    depth = pd.to_numeric(o["host_umis"], errors="coerce").to_numpy(float)
    depth = np.where(np.isfinite(depth) & (depth > 0), depth, np.nan)
    classes = EGGSHELL_ORDER + ["all eggshell"]
    score = {}
    for c in EGGSHELL_ORDER:
        ids = genes.loc[(genes["class"] == c) & genes["in_matrix"], "gene_id"]
        score[c] = np.log1p(col_sum(ids) / depth * 1e4)
    ids = genes.loc[genes["in_matrix"], "gene_id"]
    score["all eggshell"] = np.log1p(col_sum(ids) / depth * 1e4)

    mk = [fb for s, fb in GERMLINE_MARKERS.items()
          if resolve(s, fb) is not None]
    nmk = np.zeros(a.n_obs, int)
    for fb in mk:
        col = X[:, var.get_loc(fb)]
        nmk += (np.asarray(col.todense()).ravel() > 0).astype(int)
    germ = nmk >= args.min_markers
    tissue = o["tissue"].astype(str).to_numpy()
    sample = o["sample"].astype(str).to_numpy()
    group = np.array([group_of(g, i) for g, i in
                      zip(o["genotype"].astype(str), o["infection_status"].astype(str))])
    titer = pd.to_numeric(o["titer"], errors="coerce").to_numpy(float) \
        if "titer" in o.columns else np.full(a.n_obs, np.nan)
    foll = np.char.find(np.char.lower(tissue.astype(str)), FOLLICLE_KEY) >= 0
    print(f"marker-defined germline cells: {int(germ.sum()):,} "
          f"(>= {args.min_markers} of {len(mk)} markers); follicle-labelled "
          f"cells: {int(foll.sum()):,}; overlap {int((germ & foll).sum()):,}")
    germ = germ & ~foll        # a follicle-labelled cell never sets the floor

    samples = sorted(pd.unique(sample))
    sgroup = {s: group[sample == s][0] for s in samples}
    srep = {s: (o.loc[sample == s, "replicate"].iloc[0]
                if "replicate" in o.columns else "") for s in samples}

    # ---- ambient floor per library ----
    floor_rows = []
    floor = {}
    for c in classes:
        pooled = float(np.nanquantile(score[c][germ], args.floor_quantile))
        for s in samples:
            m = germ & (sample == s)
            ok = int(m.sum()) >= args.min_floor_cells
            f = float(np.nanquantile(score[c][m], args.floor_quantile)) if ok else pooled
            floor[(s, c)] = f
            floor_rows.append({"sample": s, "group": sgroup[s], "class": c,
                               "n_germline": int(m.sum()), "floor": f,
                               "source": "own library" if ok else "pooled"})
    pd.DataFrame(floor_rows).to_csv(os.path.join(args.outdir, "eggshell_floor.tsv"),
                                    sep="\t", index=False)

    # ---- per sample x follicle type x class ----
    ftypes = sorted(pd.unique(tissue[foll]))
    keep_types = [t for t in ftypes
                  if all(int(((tissue == t) & (sample == s)).sum()) >= args.min_cells
                         for s in samples)]
    print(f"follicle types with >= {args.min_cells} cells in every sample: "
          f"{len(keep_types)} of {len(ftypes)}")
    units = [(t, (tissue == t)) for t in keep_types] + [(POOLED, foll)]
    by = []
    for t, tm in units:
        for c in classes:
            for s in samples:
                m = tm & (sample == s)
                v = score[c][m]
                v = v[np.isfinite(v)]
                on = v > floor[(s, c)]
                by.append({"sample": s, "group": sgroup[s], "replicate": srep[s],
                           "tissue": t, "class": c, "n_cells": len(v),
                           "floor": floor[(s, c)],
                           "frac_on": float(on.mean()) if len(v) else np.nan,
                           "score_on": float(v[on].mean()) if on.any() else np.nan,
                           "score_all": float(v.mean()) if len(v) else np.nan})
    by = pd.DataFrame(by)
    by.to_csv(os.path.join(args.outdir, "eggshell_by_sample.tsv"), sep="\t",
              index=False)

    # ---- tests ----
    tests = []
    for (t, c), s in by.groupby(["tissue", "class"], sort=False):
        s = s.set_index("sample").reindex(samples)
        for metric, tf in (("frac_on", logit), ("score_on", None),
                           ("score_all", None)):
            bad = clipped_groups(s[metric], s["group"]) if metric == "frac_on" else []
            v = tf(s[metric]) if tf else s[metric]
            r = sample_test(v, s["group"])
            if bad:
                r[["se", "t", "pvalue"]] = np.nan
            r.insert(0, "metric", metric)
            r.insert(0, "class", c)
            r.insert(0, "tissue", t)
            r["clipped_groups"] = ",".join(bad)
            tests.append(r)
    tests = pd.concat(tests, ignore_index=True)
    tests["fdr"] = np.nan
    for _, idx in tests.groupby(["contrast", "metric"]).groups.items():
        p = tests.loc[idx, "pvalue"]
        ok = p.notna()
        if ok.any():
            tests.loc[p.index[ok], "fdr"] = multipletests(p[ok], method="fdr_bh")[1]
    tests.to_csv(os.path.join(args.outdir, "eggshell_tests.tsv"), sep="\t",
                 index=False)

    # ---- ambient control ----
    amb = []
    for c in classes:
        for s in samples:
            gm = germ & (sample == s)
            fm = foll & (sample == s)
            amb.append({"sample": s, "group": sgroup[s], "class": c,
                        "germline_mean": float(np.nanmean(score[c][gm])),
                        "follicle_mean": float(np.nanmean(score[c][fm])),
                        "n_germline": int(gm.sum()), "n_follicle": int(fm.sum())})
    amb = pd.DataFrame(amb)
    amb.to_csv(os.path.join(args.outdir, "eggshell_ambient.tsv"), sep="\t",
               index=False)

    # ---- titer link, infected libraries only ----
    tit = []
    for t, tm in units:
        for c in classes:
            for s in samples:
                if not sgroup[s].endswith("wMel"):
                    continue
                m = tm & (sample == s) & np.isfinite(titer) & np.isfinite(score[c])
                if int(m.sum()) < args.min_cells:
                    continue
                r = stats.spearmanr(score[c][m], titer[m])
                tit.append({"sample": s, "group": sgroup[s], "tissue": t,
                            "class": c, "n_cells": int(m.sum()),
                            "rho": float(r.statistic), "pvalue": float(r.pvalue)})
    pd.DataFrame(tit).to_csv(os.path.join(args.outdir, "eggshell_titer.tsv"),
                             sep="\t", index=False)

    # ---- the numbers the text will quote ----
    pd.set_option("display.width", 200)
    print("\nPooled follicle cells, per group (mean of the two libraries):")
    pv = (by[by["tissue"] == POOLED].groupby(["class", "group"])
          [["frac_on", "score_on", "score_all"]].mean().round(3))
    print(pv.to_string())
    print("\nInteraction, pooled follicle cells:")
    it = tests[(tests["tissue"] == POOLED) & (tests["contrast"] == "interaction")]
    print(it[["class", "metric", "estimate", "pvalue", "fdr"]].round(4)
          .to_string(index=False))
    for c in classes:
        x = amb[amb["class"] == c]
        r = stats.spearmanr(x["follicle_mean"], x["germline_mean"])
        print(f"ambient check, {c:<20} germline vs follicle across 8 "
              f"libraries: rho = {r.statistic:+.2f}")
    print("\nRead frac_on against score_on. If frac_on carries the "
          "interaction and score_on does not, wMel puts more follicle cells "
          "into the secreting stages; if score_on carries it, each cell "
          "makes more.")


if __name__ == "__main__":
    main()

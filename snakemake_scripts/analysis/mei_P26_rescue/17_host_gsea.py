#!/usr/bin/env python
"""Host pathway enrichment on the rescue contrasts, per cell type.

Replaces the one-vs-rest analysis in snakemake_scripts/analysis/cluster_gsea.py
(rule cluster_gsea, now removed). That analysis asked what distinguishes each
cell type from the rest of the ovary, pooling all four groups, so it could say
nothing about genotype or infection. This asks the rescue question directly:
which pathways move with wMel differently in the mutant than in wild type.

For each ovarian cell type and each contrast in --contrasts, genes are ranked
by sign(logFC) x -log10(p) from the edgeR fits in de/host_de_all.tsv.gz and run
through gseapy prerank against each gene set library. The default contrasts:

  interaction           the rescue test; terms here respond to wMel
                        differently in the mutant
  infection_in_MeiP26   what wMel does in the mutant, for direction
  genotype_in_uninf     the defect, so a rescued term can be checked against
                        what the mutant did to it in the first place

FDR is gseapy's, computed within one library for one cell type and one
contrast. It is not corrected across cell types or contrasts, so a term at
FDR 0.04 in one cell type of seven is weaker evidence than its q-value looks.

Low-confidence and non-ovarian labels are dropped by default (artefact,
unannotated, adult tracheal cell, oviduct), as in the Figure 4 panels. In
germline cell types the eggshell genes of plots/gene_sets.py are removed from
the ranking before enrichment: they are follicle-cell transcripts carried into
germline droplets, and left in, the only interaction terms in 16-cell cysts are
eggshell and vitelline membrane formation. --keep_ambient restores them.

Gene sets are the Enrichr libraries the old analysis cached as .gmt files.
They are read from --gmt_dir and downloaded through gseapy only if missing, so
the run needs no network when the cache exists.

Outputs, in --outdir:
  host_gsea.tsv           every term tested, every cell type, every contrast
  host_gsea_summary.tsv   terms passing FDR < 0.05 and < 0.25, per cell type
                          and contrast

Usage:
  python 17_host_gsea.py --de results/rescue/de/host_de_all.tsv.gz \
      --gmt_dir results/cluster_gsea/gene_sets --outdir results/rescue/host_gsea
"""
import argparse
import os
import warnings

import sys

import numpy as np
import pandas as pd

warnings.filterwarnings("ignore")
sys.path.insert(0, os.path.join(os.path.dirname(os.path.abspath(__file__)),
                                "plots"))
from gene_sets import eggshell_symbols, is_germline

LOW_CONFIDENCE = ("artefact", "unannotated")
NON_OVARIAN = ("adult tracheal cell", "oviduct")


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--de", required=True, help="de/host_de_all.tsv.gz")
    p.add_argument("--gmt_dir", required=True,
                   help="directory of <library>.gmt files")
    p.add_argument("--outdir", required=True)
    p.add_argument("--libraries", nargs="+",
                   default=["KEGG_2019", "GO_Biological_Process_2018"])
    p.add_argument("--contrasts", nargs="+",
                   default=["interaction", "infection_in_MeiP26",
                            "genotype_in_uninf"])
    p.add_argument("--drop_labels",
                   default=",".join(LOW_CONFIDENCE + NON_OVARIAN))
    p.add_argument("--keep_ambient", action="store_true",
                   help="keep eggshell genes in germline rankings. By default "
                        "they are removed there: germ cells do not make "
                        "eggshell, the reads track follicle-cell levels across "
                        "libraries (18_eggshell.py), and left in they make "
                        "eggshell formation the top germline term")
    p.add_argument("--min_size", type=int, default=10)
    p.add_argument("--max_size", type=int, default=500)
    p.add_argument("--permutations", type=int, default=1000)
    p.add_argument("--seed", type=int, default=0)
    p.add_argument("--threads", type=int, default=4)
    p.add_argument("--organism", default="Fly",
                   help="passed to gseapy.get_library when a .gmt is missing")
    return p.parse_args()


def load_library(gmt_dir, name, organism):
    import gseapy as gp
    path = os.path.join(gmt_dir, f"{name}.gmt")
    if not os.path.exists(path):
        print(f"  {path} missing, downloading {name} through gseapy")
        os.makedirs(gmt_dir, exist_ok=True)
        gs = gp.get_library(name=name, organism=organism)
        with open(path, "w") as f:
            for term, genes in gs.items():
                f.write("\t".join([term, "NA"] + list(genes)) + "\n")
    sets = {}
    with open(path) as f:
        for line in f:
            p = line.rstrip("\n").split("\t")
            if len(p) > 2:
                sets[p[0]] = [g for g in p[2:] if g]
    print(f"  {name}: {len(sets):,} sets")
    return sets


def ranking(s):
    s = s.dropna(subset=["symbol", "PValue", "logFC"])
    s = s[s["symbol"].astype(str).str.len() > 0]
    s = s[s["symbol"].astype(str) != "nan"]
    r = (np.sign(s["logFC"]) * -np.log10(s["PValue"].clip(lower=1e-300))
         # A tiny logFC term orders genes that share a p-value, so the
         # ranking does not depend on input order between runs.
         + s["logFC"] * 1e-8)
    return (pd.Series(r.to_numpy(), index=s["symbol"].astype(str))
            .groupby(level=0).mean().sort_values(ascending=False))


def main():
    args = parse_args()
    import gseapy as gp
    os.makedirs(args.outdir, exist_ok=True)
    d = pd.read_csv(args.de, sep="\t")
    drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}
    d = d[~d["tissue"].isin(drop) & d["contrast"].isin(args.contrasts)]
    print(f"{d['tissue'].nunique()} cell types x {d['contrast'].nunique()} "
          f"contrasts; dropped {sorted(drop) or 'none'}")
    libs = {n: load_library(args.gmt_dir, n, args.organism) for n in args.libraries}

    out = []
    egg = set(eggshell_symbols())
    for (con, t), s in d.groupby(["contrast", "tissue"]):
        if is_germline(t) and not args.keep_ambient:
            n0 = len(s)
            s = s[~s["symbol"].isin(egg)]
            if n0 - len(s) and con == args.contrasts[0]:
                print(f"  {t[:46]}: {n0 - len(s)} eggshell genes removed as "
                      "ambient")
        r = ranking(s)
        for ln, sets in libs.items():
            try:
                res = gp.prerank(rnk=r, gene_sets=sets, min_size=args.min_size,
                                 max_size=args.max_size,
                                 permutation_num=args.permutations,
                                 seed=args.seed, threads=args.threads,
                                 no_plot=True, outdir=None, verbose=False).res2d
            except Exception as e:      # a cell type with too few genes
                print(f"  skipped {con} / {t} / {ln}: {e}")
                continue
            res = res.rename(columns={"Term": "term", "ES": "es", "NES": "nes",
                                      "NOM p-val": "pval", "FDR q-val": "fdr",
                                      "FWER p-val": "fwer", "Tag %": "tag",
                                      "Gene %": "gene_pct",
                                      "Lead_genes": "lead_genes"})
            res = res.drop(columns=[c for c in ("Name",) if c in res.columns])
            res.insert(0, "library", ln)
            res.insert(0, "tissue", t)
            res.insert(0, "contrast", con)
            res["n_ranked"] = len(r)
            out.append(res)
        print(f"  {con:<20} {t[:46]:<48} {len(r):>6} genes ranked")

    o = pd.concat(out, ignore_index=True)
    for c in ("es", "nes", "pval", "fdr", "fwer"):
        o[c] = pd.to_numeric(o[c], errors="coerce")
    o = o.sort_values(["contrast", "tissue", "fdr", "pval"])
    f = os.path.join(args.outdir, "host_gsea.tsv")
    o.to_csv(f, sep="\t", index=False)

    summ = (o.groupby(["contrast", "tissue"])
            .agg(n_terms=("term", "size"),
                 fdr_05=("fdr", lambda v: int((v < 0.05).sum())),
                 fdr_25=("fdr", lambda v: int((v < 0.25).sum())),
                 min_fdr=("fdr", "min"))
            .reset_index())
    summ.to_csv(os.path.join(args.outdir, "host_gsea_summary.tsv"),
                sep="\t", index=False)
    print(f"\nwrote {f} ({len(o):,} rows)")
    print(summ.to_string(index=False))

    top = o[(o["contrast"] == "interaction") & (o["fdr"] < 0.05)]
    if len(top):
        print("\nInteraction terms at FDR < 0.05, top 5 per cell type by FDR:")
        for t, s in top.groupby("tissue"):
            print(f"  {t}")
            for _, r in s.head(5).iterrows():
                print(f"    {r['nes']:+.2f}  {r['fdr']:.3g}  {r['term'][:70]}")
    print("\nFDR is per library, cell type and contrast, not across them. "
          "Ribosome and translation terms deserve the most suspicion: they "
          "move with depth and with stage mixture inside a label, and depth "
          "is confounded with genotype in this design.")


if __name__ == "__main__":
    main()

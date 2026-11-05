#!/usr/bin/env python
"""Gene-set enrichment on the Wolbachia DE ranking.

Asks whether the wMel response to a Mei-P26 host is a coherent functional
shift rather than a list of individual genes. Three set families are tested
against the same ranking (pooled Mei-P26 vs OreR logFC):

  GO          parsed from the annotation's go_process / go_function /
              go_component attributes (no network needed)
  KEGG        pathways for the wMel organism code, fetched from the KEGG REST
              API or read from --kegg_dir if already downloaded
  category    product-name categories (ribosomal, chaperone, effector-like,
              hypothetical, type IV secretion, mobile element, ...), which
              cover the genes GO and KEGG leave unannotated

GO and KEGG are scored by preranked GSEA (gseapy); the product categories by a
two-sided Mann-Whitney U of in-set vs out-of-set logFC, which behaves better
for a handful of small sets. Both are competitive tests: they ask whether a
set sits unusually high or low in the ranking, not whether its genes are
significant on their own.

WHAT TO EXPECT. On the Nov 2026 data every family is null (best nominal
p = 0.08, all FDR > 0.14). Three things limit it, and all are structural:

  - only ~210 of 1,154 wMel genes survive filterByExpr, so the ranking is
    short and biased toward highly expressed genes
  - GO and KEGG annotate housekeeping, not host interaction: 271 of 273
    hypothetical proteins have no GO term, and neither do CifA, TomO or WalE1
  - only ~14 sets per family reach 5 tested genes

So a null here is weak evidence of no pathway-level shift, and the gene-level
results stay the stronger read. Rerun if the DE gains power or the annotation
improves.

Usage:
  python 12_wolbachia_gsea.py \
      --gff reference/.../Drosophila_melanogaster_wMel_combined_16S.gtf \
      --secretion results/rescue/wolbachia_followup/wolbachia_secretion_all.tsv \
      --de results/rescue/wolbachia_followup/wolbachia_de_checked.tsv.gz \
      --outdir results/rescue/wolbachia_followup [--kegg_org wol] [--kegg_dir DIR]
"""
import argparse
import collections
import gzip
import os
import re
import urllib.request

import numpy as np
import pandas as pd
from scipy import stats
from scipy.stats import false_discovery_control

CATEGORIES = {
    "ribosomal protein": "ribosomal protein",
    "chaperone / heat shock": "chaperon|heat shock|groes|groel|dnak|clp|hsp",
    "effector-like (ankyrin/TPR/LRR/collagen)": "ankyrin|tetratricopeptide|leucine-rich|collagen",
    "hypothetical / DUF": "hypothetical|duf|uncharacteri",
    "outer membrane / surface": "outer membrane|p44|msp2|surface|lipoprotein",
    "type IV secretion": "type iv secretion|virb|trb",
    "mobile element": "transposase|phage|integrase|recombinase",
    "energy (NADH/ATP synthase/cytochrome)": "nadh|atp synthase|cytochrome|ubiquinol",
    "tRNA / aminoacyl": "trna",
    "transporter": "transporter|permease|symporter|abc ",
}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gff", required=True, help="GTF/GFF3 with go_* attributes")
    p.add_argument("--secretion", required=True, help="wolbachia_secretion_all.tsv (for products, old tags)")
    p.add_argument("--de", required=True, help="wolbachia_de_checked.tsv.gz")
    p.add_argument("--outdir", required=True)
    p.add_argument("--tissue", default="pooled_all_tissues", help="which DE contrast to rank by")
    p.add_argument("--kegg_org", default="wol", help="KEGG organism code ('' to skip KEGG)")
    p.add_argument("--kegg_dir", default=None,
                   help="directory holding kegg_links.tsv and kegg_names.tsv; "
                        "downloaded there when absent")
    p.add_argument("--min_set", type=int, default=5, help="min tested genes per set")
    p.add_argument("--permutations", type=int, default=2000)
    return p.parse_args()


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_go(path):
    """locus_tag -> GO terms, from go_process / go_function / go_component."""
    sets, names = collections.defaultdict(set), {}
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            lt = re.search(r'locus_tag[ =]"?([^";\s]+)', line)
            if not lt:
                continue
            for key in ("go_process", "go_function", "go_component"):
                m = re.search(key + r'[ =]"?([^";]*)', line)
                if not m:
                    continue
                for item in m.group(1).split(","):
                    parts = item.split("|")
                    if len(parts) >= 2 and parts[1].strip():
                        gid = "GO:" + parts[1].strip()
                        sets[gid].add(lt.group(1))
                        names[gid] = parts[0].strip()
    return sets, names


def fetch_kegg(org, cache_dir):
    """pathway -> KEGG gene ids, and pathway -> name. Cached as TSV."""
    out = {}
    for kind, url in (("links", f"https://rest.kegg.jp/link/{org}/pathway"),
                      ("names", f"https://rest.kegg.jp/list/pathway/{org}")):
        path = os.path.join(cache_dir, f"kegg_{kind}.tsv") if cache_dir else None
        if path and os.path.exists(path):
            out[kind] = open(path).read()
            continue
        try:
            out[kind] = urllib.request.urlopen(url, timeout=60).read().decode()
        except Exception as e:
            print(f"  KEGG {kind} unavailable ({e}); skipping KEGG")
            return None, None
        if path:
            os.makedirs(cache_dir, exist_ok=True)
            open(path, "w").write(out[kind])
    links = collections.defaultdict(set)
    for line in out["links"].splitlines():
        if "\t" not in line:
            continue
        p, g = line.split("\t")
        links[p.replace("path:", "")].add(g.split(":", 1)[1])
    names = {}
    for line in out["names"].splitlines():
        if "\t" in line:
            k, v = line.split("\t", 1)
            names[k] = v.split(" - ")[0]
    return links, names


def prerank(rnk, sets, label, n_perm):
    try:
        import gseapy as gp
    except ImportError:
        print(f"  gseapy not installed; skipping {label}")
        return pd.DataFrame()
    if not sets:
        print(f"  no {label} sets reach the size cutoff")
        return pd.DataFrame()
    r = gp.prerank(rnk=rnk, gene_sets=sets, permutation_num=n_perm, min_size=2,
                   max_size=1000, threads=2, seed=1, no_plot=True, outdir=None).res2d
    r = r.rename(columns={"NOM p-val": "p", "FDR q-val": "FDR_gseapy"})
    r["family"] = label
    r["n_genes"] = [len(sets[t]) for t in r["Term"]]
    for c in ("NES", "p", "FDR_gseapy"):
        r[c] = r[c].astype(float)
    return r[["family", "Term", "n_genes", "NES", "p", "FDR_gseapy"]]


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)

    d = pd.read_csv(args.secretion, sep="\t", index_col=0)
    de = pd.read_csv(args.de, sep="\t")
    de = de[de["tissue"] == args.tissue].dropna(subset=["logFC"]).set_index("gene_id")
    rnk = de["logFC"].sort_values(ascending=False)
    tested = set(rnk.index)
    print(f"Ranking: {len(rnk)} genes from contrast '{args.tissue}'")

    results = []

    go_sets, go_names = parse_go(args.gff)
    go_use = {f"{k} {go_names.get(k, '')}".strip(): sorted(v & tested)
              for k, v in go_sets.items()}
    go_use = {k: v for k, v in go_use.items() if len(v) >= args.min_set}
    print(f"GO: {len(go_sets)} terms, {len(go_use)} with >= {args.min_set} tested genes")
    results.append(prerank(rnk, go_use, "GO", args.permutations))

    if args.kegg_org:
        links, knames = fetch_kegg(args.kegg_org, args.kegg_dir)
        if links:
            old2rs = {str(r["old"]).replace("_", ""): rs
                      for rs, r in d.iterrows() if pd.notna(r.get("old"))}
            kegg_use = {}
            for p, genes in links.items():
                mapped = {old2rs[g.replace("_", "")] for g in genes
                          if g.replace("_", "") in old2rs} & tested
                if len(mapped) >= args.min_set:
                    kegg_use[f"{p} {knames.get(p, '')}".strip()] = sorted(mapped)
            print(f"KEGG: {len(links)} pathways, {len(kegg_use)} with >= {args.min_set} tested genes")
            results.append(prerank(rnk, kegg_use, "KEGG", args.permutations))

    prod = d["product"].fillna("").str.lower()
    lfc = d.index.map(lambda g: rnk.get(g, np.nan))
    rows = []
    for name, pat in CATEGORIES.items():
        m = prod.str.contains(pat).to_numpy() & ~pd.isna(lfc)
        inside = pd.Series(lfc)[m].astype(float)
        outside = pd.Series(lfc)[~m & ~pd.isna(lfc)].astype(float)
        if len(inside) < args.min_set - 1:
            continue
        p = stats.mannwhitneyu(inside, outside, alternative="two-sided")[1]
        rows.append({"family": "category", "Term": name, "n_genes": len(inside),
                     "NES": float(inside.median() - outside.median()), "p": p,
                     "FDR_gseapy": np.nan})
    print(f"category: {len(rows)} product categories with >= {args.min_set - 1} tested genes")
    results.append(pd.DataFrame(rows))

    r = pd.concat([x for x in results if len(x)], ignore_index=True)
    r["FDR"] = np.nan
    for fam, sub in r.groupby("family"):
        r.loc[sub.index, "FDR"] = false_discovery_control(sub["p"].values)
    r = r.sort_values(["family", "p"])
    r.to_csv(os.path.join(args.outdir, "wolbachia_gsea.tsv"), sep="\t", index=False)

    print("\nFor the category family NES is the median logFC difference, in-set minus rest.")
    print(r.round(4).to_string(index=False))
    best = r.loc[r["p"].idxmin()]
    if (r["FDR"] >= 0.05).all():
        print(f"\nNothing significant: best is {best['Term']} (p = {best['p']:.3g}, "
              f"FDR = {best['FDR']:.2g}). See the docstring for why this is expected.")


if __name__ == "__main__":
    main()

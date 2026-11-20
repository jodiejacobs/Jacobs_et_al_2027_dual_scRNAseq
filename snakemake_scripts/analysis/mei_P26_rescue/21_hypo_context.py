#!/usr/bin/env python
"""Genomic neighbourhood of each hit locus.

Lists the --window genes on either side (by position along the wMel genome), the
number of neighbours with mobile-element / prophage-like products, and the
number of neighbours that are also significantly up (up_loci.txt). A hit whose
neighbours are mostly up is a region-level effect (prophage module or operon),
not a lone gene. The product keywords are a hint, so read the neighbour lists.
"""
import argparse
import gzip
import re

import pandas as pd

PREFIXES = ("WD", "GQX67")
MOBILE = re.compile(r"phage|baseplate|terminase|capsid|portal|tail|integrase|transposase|"
                    r"recombinase|resolvase|Rpn|\bWO\b|insertion|IS\d", re.I)


def parse_genes(path):
    op = gzip.open if path.endswith(".gz") else open
    genes = {}
    with op(path, "rt") as f:
        for line in f:
            if line.startswith("#"):
                continue
            p = line.rstrip("\n").split("\t")
            if len(p) != 9 or p[2] not in ("gene", "CDS"):
                continue
            m = re.search(r'locus_tag[ =]"?([^";]+)', p[8])
            if not m or not m.group(1).startswith(PREFIXES):
                continue
            r = genes.setdefault(m.group(1), {"locus_tag": m.group(1), "start": int(p[3]),
                                              "end": int(p[4]), "strand": p[6], "product": ""})
            pr = re.search(r'product[ =]"?([^";]+)', p[8])
            if pr and not r["product"]:
                r["product"] = pr.group(1)
    return pd.DataFrame(genes.values()).sort_values("start").reset_index(drop=True)


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--gtf", required=True)
    p.add_argument("--hits", required=True, help="hits.tsv")
    p.add_argument("--up_loci", required=True, help="up_loci.txt")
    p.add_argument("--window", type=int, default=5)
    p.add_argument("--out", required=True)
    a = p.parse_args()

    g = parse_genes(a.gtf)
    up = set(open(a.up_loci).read().split())
    hits = pd.read_csv(a.hits, sep="\t")["locus_tag"]
    rows = []
    for lt in hits:
        i = g.index[g["locus_tag"] == lt]
        if len(i) == 0:
            rows.append({"locus_tag": lt})
            continue
        i = i[0]
        me = g.loc[i]
        nb = pd.concat([g.loc[max(0, i - a.window): i - 1], g.loc[i + 1: i + a.window]])
        rows.append({
            "locus_tag": lt, "start": me.start, "end": me.end, "strand": me.strand,
            "n_mobile_neighbours": int(nb["product"].str.contains(MOBILE).sum()),
            "n_up_neighbours": int(nb["locus_tag"].isin(up).sum()),
            "n_neighbours": len(nb),
            "neighbours": " | ".join(f"{r.locus_tag}({r.strand}){'*' if r.locus_tag in up else ''}:{r.product}"
                                     for r in nb.itertuples()),
        })
    pd.DataFrame(rows).to_csv(a.out, sep="\t", index=False)
    print(f"context for {len(rows)} hits ({len(g)} wMel genes in the GTF); * = also significant up")


if __name__ == "__main__":
    main()

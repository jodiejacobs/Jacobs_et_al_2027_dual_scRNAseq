#!/usr/bin/env python
"""Select hypothetical wMel loci that are significantly up in the mutant.

A locus is a hit when, in at least one tissue of --contrast, FDR < --fdr and
logFC > 0, and its product matches --regex (default: hypothetical / DUF /
uncharacterised). Products and the replicate / leave-one-out robustness flags
come from 09_wolbachia_followup.R (wolbachia_de_checked.tsv.gz).

Writes
  hits.tsv            one row per hit locus
  hits_by_tissue.tsv  every significant-up (locus, tissue) row, all products
  up_loci.txt         all significant-up loci, any product (for neighbour checks)
  hits.faa            hit proteins, taken from --proteome_faa (header = locus tag)
"""
import argparse
import os

import pandas as pd


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as f:
        for line in f:
            line = line.rstrip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name:
                seqs[name].append(line)
    return {k: "".join(v) for k, v in seqs.items()}


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--de", required=True, help="wolbachia_de_all.tsv.gz")
    p.add_argument("--checked", required=True, help="wolbachia_de_checked.tsv.gz")
    p.add_argument("--proteome_faa", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--contrast", default="MeiP26_vs_OreR")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--regex", default="hypothetical|DUF|uncharacteri[sz]ed")
    a = p.parse_args()
    os.makedirs(a.outdir, exist_ok=True)

    de = pd.read_csv(a.de, sep="\t")
    chk = pd.read_csv(a.checked, sep="\t")
    flags = [c for c in ("rep_consistent", "loo_sign_stable", "FDR_excl_dominant")
             if c in chk.columns]
    prod = chk.drop_duplicates("gene_id")[["gene_id", "product", "old_locus_tag"]]

    up = de[(de["contrast"] == a.contrast) & (de["FDR"] < a.fdr) & (de["logFC"] > 0)]
    up = (up.merge(prod, on="gene_id", how="left")
            .merge(chk[["gene_id", "tissue"] + flags], on=["gene_id", "tissue"], how="left"))
    up["hypothetical"] = up["product"].fillna("").str.contains(a.regex, case=False, regex=True)
    up.sort_values(["gene_id", "FDR"]).to_csv(f"{a.outdir}/hits_by_tissue.tsv", sep="\t", index=False)
    with open(f"{a.outdir}/up_loci.txt", "w") as f:
        f.write("\n".join(sorted(up["gene_id"].unique())) + "\n")

    rows = []
    for gid, g in up[up["hypothetical"]].groupby("gene_id"):
        g = g.sort_values("FDR")
        robust = (g["rep_consistent"].eq(True) & g["loo_sign_stable"].eq(True)
                  if {"rep_consistent", "loo_sign_stable"} <= set(g.columns) else pd.Series(False, index=g.index))
        rows.append({"locus_tag": gid, "old_locus_tag": g["old_locus_tag"].iloc[0],
                     "product": g["product"].iloc[0], "n_tissues_up": len(g),
                     "tissues": "; ".join(g["tissue"]), "best_logFC": g["logFC"].iloc[0],
                     "best_FDR": g["FDR"].iloc[0], "n_robust_tissues": int(robust.sum())})
    hits = pd.DataFrame(rows, columns=["locus_tag", "old_locus_tag", "product", "n_tissues_up",
                                       "tissues", "best_logFC", "best_FDR", "n_robust_tissues"])
    hits = hits.sort_values("best_FDR")
    hits.to_csv(f"{a.outdir}/hits.tsv", sep="\t", index=False)

    seqs = read_fasta(a.proteome_faa)
    missing = [g for g in hits["locus_tag"] if g not in seqs]
    if missing:
        print(f"WARNING: no protein sequence for {missing}")
    with open(f"{a.outdir}/hits.faa", "w") as f:
        for g in hits["locus_tag"]:
            if g in seqs:
                f.write(f">{g}\n{seqs[g]}\n")
    print(f"{len(hits)} hypothetical hit loci of {up['gene_id'].nunique()} significant-up loci "
          f"({a.contrast}, FDR < {a.fdr})")


if __name__ == "__main__":
    main()

#!/usr/bin/env python
"""One evidence table for the hypothetical hit loci.

Joins DE (hits.tsv), secretion class (11_wolbachia_secretion.py), InterProScan
domains, eggNOG-mapper annotation, ColabFold mean pLDDT, the top Foldseek hit
per database, and genomic context. Missing inputs leave their columns empty.
"""
import argparse
import glob
import json
import os

import pandas as pd


def interpro(path):
    cols = ["locus_tag", "md5", "len", "analysis", "acc", "desc", "start", "end", "score",
            "status", "date", "ipr", "ipr_desc", "go", "pathways"]
    df = pd.read_csv(path, sep="\t", header=None, names=cols, dtype=str, on_bad_lines="skip")
    sc = pd.to_numeric(df["score"], errors="coerce")
    df = df[sc.isna() | (sc <= 1e-3)]          # keep unscored (-) and E <= 1e-3
    df = df[df["acc"].notna() & (df["acc"] != "-")]
    out = {}
    for lt, g in df.groupby("locus_tag"):
        out[lt] = {"interpro_domains": "; ".join(sorted({f"{r.analysis}:{r.acc}({r.desc})" for r in g.itertuples()})),
                   "interpro_entries": "; ".join(sorted({f"{r.ipr}({r.ipr_desc})" for r in g.itertuples() if r.ipr != "-"})),
                   "interpro_go": "; ".join(sorted({x for v in g["go"].dropna() for x in v.split("|") if x != "-"}))}
    return pd.DataFrame.from_dict(out, orient="index").rename_axis("locus_tag").reset_index()


def eggnog(path):
    lines = [l.rstrip("\n") for l in open(path)]
    head = next(i for i, l in enumerate(lines) if l.startswith("#query"))
    cols = lines[head].lstrip("#").split("\t")
    rows = [l.split("\t") for l in lines[head + 1:] if l and not l.startswith("#")]
    df = pd.DataFrame(rows, columns=cols).rename(columns={"query": "locus_tag"})
    keep = [c for c in ("locus_tag", "Description", "Preferred_name", "COG_category", "max_annot_lvl") if c in df]
    return df[keep].rename(columns=lambda c: c if c == "locus_tag" else f"eggnog_{c}")


def plddt(dirpath):
    rows = []
    for f in glob.glob(os.path.join(dirpath, "*.json")):
        vals = json.load(open(f)).get("plddt")
        if vals:
            rows.append({"locus_tag": os.path.basename(f)[:-len(".json")], "mean_plddt": sum(vals) / len(vals)})
    return pd.DataFrame(rows, columns=["locus_tag", "mean_plddt"])


def foldseek(name, path):
    cols = ["query", "target", "fident", "alnlen", "evalue", "bits", "alntmscore", "prob", "theader"]
    if not os.path.exists(path) or os.path.getsize(path) == 0:
        return pd.DataFrame(columns=["locus_tag"])
    df = pd.read_csv(path, sep="\t", header=None, names=cols)
    df["locus_tag"] = df["query"].str.replace(r"\.pdb.*$", "", regex=True)
    top = df.sort_values("evalue").drop_duplicates("locus_tag")
    top = top[["locus_tag", "target", "theader", "evalue", "alntmscore", "prob"]]
    return top.rename(columns=lambda c: c if c == "locus_tag" else f"fs_{name}_{c}")


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument("--hits", required=True)
    p.add_argument("--secretion", required=True)
    p.add_argument("--interpro", required=True)
    p.add_argument("--eggnog", required=True)
    p.add_argument("--scores_dir", required=True)
    p.add_argument("--context", required=True)
    p.add_argument("--foldseek", nargs="*", default=[], help="name=path pairs")
    p.add_argument("--out", required=True)
    a = p.parse_args()

    out = pd.read_csv(a.hits, sep="\t")
    sec = pd.read_csv(a.secretion, sep="\t")
    sec = sec[["locus_tag"] + [c for c in ("aa_length", "signal_peptide", "tm_helices_kd",
                                           "secretion_class", "euk_domain") if c in sec]]
    parts = [sec, interpro(a.interpro), eggnog(a.eggnog), plddt(a.scores_dir),
             pd.read_csv(a.context, sep="\t")]
    parts += [foldseek(*kv.split("=", 1)) for kv in a.foldseek]
    for d in parts:
        if len(d) or "locus_tag" in d:
            out = out.merge(d, on="locus_tag", how="left")
    out.to_csv(a.out, sep="\t", index=False)
    print(f"{len(out)} hit loci, {out.shape[1]} columns -> {a.out}")


if __name__ == "__main__":
    main()

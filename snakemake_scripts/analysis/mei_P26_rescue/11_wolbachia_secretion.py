#!/usr/bin/env python
"""Rank wMel genes as candidate host-acting (secreted) rescue factors.

Wolbachia cannot be transformed, so a candidate is tested by expressing it in
the fly (UAS/Gal4), the way TomO and WalE1 were. This script builds the
shortlist: every wMel gene with its secretion class, its host-acting evidence,
and its behaviour in this dataset.

SECRETION CLASS (what sequence can and cannot tell you)

  A Sec signal peptide sends a protein to the periplasm, outer membrane or
  bacterial surface -- NOT into the host cytoplasm. Proteins that reach the
  host cytoplasm go through the type IV secretion system (T4SS), and T4SS
  substrates have no Sec signal and no transmembrane helices; their secretion
  signal is C-terminal and poorly defined, so sequence-based T4SS prediction is
  weak. The classes here are therefore a filter, not a prediction of host
  translocation:

    surface_periplasmic   Sec signal peptide -> bacterial envelope/surface
    membrane              no signal peptide, >= 1 predicted TM helix
    cytoplasmic_or_T4SS   no signal peptide, no TM helix -> can be a T4SS
                          substrate (also most housekeeping proteins)

  Host-acting evidence, not the class, is what ranks a candidate:

    published    shown to act in host cells (TomO, WalE1, CifA, CifB)
    predicted    predicted T4SS effector with experimental support
                 (Rice, Sheehan & Newton 2017 yeast screen)
    domain       eukaryote-like domain (ankyrin, TPR, leucine-rich,
                 collagen-like) -- the commonest effector domains in Wolbachia
    uncharacterised  large (>= --min_aa) soluble protein of unknown function:
                 no annotation to rule it in or out, so expression decides

  CAVEAT on tm_helices_kd: the Kyte-Doolittle screen calls TomO and WalE1
  "membrane" even though both act in host cells, so a membrane call must not be
  used to drop a candidate. Use it to prioritise, never to exclude, and replace
  it with DeepTMHMM output when you have it.

INPUTS: the pipeline's own GTF (locus tags, products) and Wolbachia DE tables.
--proteins and --signalp are optional; without them the secretion class is
"unknown" and ranking falls back to evidence plus expression.

Signal peptides can be predicted with any tool that writes one line per
protein; pass its output with --signalp (DeepSig GFF, or a two-column TSV of
protein_id and Y/N). Without a local tool, upload --write_faa output to the
DeepTMHMM or SignalP web servers and pass the result back.

Usage:
  python 11_wolbachia_secretion.py \
      --gtf reference/KB_.../Drosophila_melanogaster_wMel_combined_16S.gtf \
      --de_dir results/rescue --outdir results/rescue/wolbachia_followup \
      [--proteins wmel_proteins.faa] [--signalp wmel_deepsig.gff]
"""
import argparse
import gzip
import os
import re

import numpy as np
import pandas as pd

# Published host-acting wMel proteins, keyed by old locus tag.
PUBLISHED = {
    "WD1278": "TomO: germline expression enhances nanos translation and prevents "
              "Sxl-mutant GSC loss (Ote et al. 2016)",
    "WD0830": "WalE1: predicted T4SS effector, bundles actin, transgenic protein "
              "localises to the oocyte and raises oocyte titer (Sheehan et al. 2016)",
    "WD0631": "CifA: acts in host cells to rescue cytoplasmic incompatibility",
    "WD0632": "CifB: acts in host cells to induce cytoplasmic incompatibility",
}
# Predicted T4SS effectors with a yeast growth defect (Rice, Sheehan & Newton 2017, GBE)
PREDICTED = ["WD1171", "WD0338", "WD0033", "WD0830", "WD1223", "WD0290", "WD0811",
             "WD0462", "WD0438", "WD1321", "WD0385", "WD0292", "WD0353", "WD0465"]
# Eukaryote-like domains, matched against the product name
EUK_DOMAINS = {
    "ankyrin repeat": "ankyrin", "tetratricopeptide": "TPR", "leucine-rich": "LRR",
    "collagen-like": "collagen-like", "patatin": "patatin", "f-box": "F-box",
    "ubiquitin": "ubiquitin", "deubiquitinase": "deubiquitinase",
}
# Kyte-Doolittle hydropathy
KD = {"A": 1.8, "R": -4.5, "N": -3.5, "D": -3.5, "C": 2.5, "Q": -3.5, "E": -3.5,
      "G": -0.4, "H": -3.2, "I": 4.5, "L": 3.8, "K": -3.9, "M": 1.9, "F": 2.8,
      "P": -1.6, "S": -0.8, "T": -0.7, "W": -0.9, "Y": -1.3, "V": 4.2}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--gtf", required=True, help="combined fly + wMel GTF (or a wMel GFF3)")
    p.add_argument("--de_dir", required=True, help="results/rescue")
    p.add_argument("--outdir", required=True)
    p.add_argument("--proteins", default=None, help="wMel protein FASTA (for TM helices)")
    p.add_argument("--signalp", default=None,
                   help="signal-peptide predictions: DeepSig GFF, or TSV of protein_id + Y/N")
    p.add_argument("--write_faa", default=None,
                   help="write the wMel proteins named in the GTF to this FASTA and stop "
                        "(needs --proteins); use it to feed a web predictor")
    p.add_argument("--min_aa", type=int, default=400,
                   help="minimum length for the 'uncharacterised' tier")
    p.add_argument("--min_cpm", type=float, default=100,
                   help="candidates must reach this mean CPM of the wMel transcriptome "
                        "in infected germline pseudobulks")
    return p.parse_args()


def opener(path):
    return gzip.open(path, "rt") if path.endswith(".gz") else open(path)


def parse_annotation(path):
    """locus_tag -> old_locus_tag, gene, product, protein_id. Accepts GTF or GFF3."""
    pat_gtf = re.compile(r'(?:^|; *)%s "([^"]*)"')
    pat_gff = re.compile(r"(?:^|;)%s=([^;]*)")

    def attr(a, key):
        m = pat_gtf.pattern % key
        hit = re.search(m, a)
        if hit:
            return hit.group(1)
        hit = re.search(pat_gff.pattern % key, a)
        return hit.group(1) if hit else None

    rows = {}
    with opener(path) as f:
        for line in f:
            if line.startswith("#"):
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) != 9:
                continue
            lt = attr(parts[8], "locus_tag")
            if not lt or not lt.startswith(("WD", "GQX67")):
                continue
            r = rows.setdefault(lt, {"locus_tag": lt, "old_locus_tag": None, "gene": None,
                                     "product": None, "protein_id": None, "length_nt": 0})
            for key, col in (("old_locus_tag", "old_locus_tag"), ("gene", "gene"),
                             ("product", "product"), ("protein_id", "protein_id")):
                v = attr(parts[8], key)
                if v and not r[col]:
                    r[col] = v
            if parts[2] == "gene":
                r["length_nt"] = int(parts[4]) - int(parts[3]) + 1
    df = pd.DataFrame(rows.values())
    df["old"] = df["old_locus_tag"].fillna("").str.split(",").str[0]
    return df


def read_fasta(path):
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


def count_tm(seq, window=19, cutoff=1.6):
    """Kyte-Doolittle heuristic: windows of mean hydropathy above cutoff.

    A rough stand-in for TMHMM/DeepTMHMM, used only to separate membrane
    proteins from soluble ones. Reported as tm_helices_kd so it is never
    mistaken for a trained predictor's output.
    """
    if len(seq) < window:
        return 0
    h = np.array([KD.get(c, 0.0) for c in seq])
    means = np.convolve(h, np.ones(window) / window, mode="valid")
    above = means >= cutoff
    n, prev = 0, False
    for a in above:
        if a and not prev:
            n += 1
        prev = a
    return n


def read_signalp(path):
    """DeepSig GFF ('Signal peptide' feature) or a two-column TSV."""
    sp = {}
    with opener(path) as f:
        for line in f:
            if line.startswith("#") or not line.strip():
                continue
            parts = line.rstrip("\n").split("\t")
            if len(parts) >= 3 and parts[1] in ("DeepSig", "SignalP"):
                sp.setdefault(parts[0], False)
                if "signal peptide" in parts[2].lower():
                    sp[parts[0]] = True
            elif len(parts) >= 2:
                sp[parts[0]] = parts[1].strip().upper() in ("Y", "YES", "SP", "TRUE", "1")
    return sp


def wmel_cpm(pb_dir, tissue=None):
    """Mean CPM of the wMel transcriptome across infected pseudobulks."""
    counts = pd.read_csv(os.path.join(pb_dir, "pseudobulk_wolbachia_counts.tsv.gz"),
                         sep="\t", index_col=0)
    meta = pd.read_csv(os.path.join(pb_dir, "pseudobulk_meta.tsv"), sep="\t")
    keep = meta[meta["infection_status"] == "wMel"]
    if tissue:
        keep = keep[keep["tissue"] == tissue]
    cols = [c for c in keep["pb_id"] if c in counts.columns]
    if not cols:
        return pd.Series(dtype=float)
    sub = counts[cols]
    return (sub.div(sub.sum(axis=0), axis=1) * 1e6).mean(axis=1)


def de_table(path, tissue):
    if not os.path.exists(path):
        return pd.DataFrame()
    d = pd.read_csv(path, sep="\t")
    d = d[d["tissue"] == tissue]
    return d.set_index("gene_id")[["logFC", "FDR"]]


def de_across_tissues(path, fdr_max=0.05):
    """Per-locus summary of EVERY cell-type contrast, not just germline and pooled.

    The germline and pooled contrasts miss candidates whose response is confined
    to one cell type, which is how WD1063 (up in six contexts, FDR 2.9e-6) came
    to be scored as tier-less and dropped from the shortlist. Robust means
    replicate-consistent and leave-one-out sign-stable where those columns exist.
    """
    if not os.path.exists(path):
        return pd.DataFrame()
    d = pd.read_csv(path, sep="\t")
    d = d[d["FDR"] < fdr_max].copy()
    for c in ("rep_consistent", "loo_sign_stable"):
        if c in d.columns:
            d = d[d[c].astype(str).str.upper().isin(("TRUE", "T", "1", "YES"))]
    if not len(d):
        return pd.DataFrame()
    up, dn = d[d["logFC"] > 0], d[d["logFC"] < 0]
    out = pd.DataFrame(index=sorted(set(d["gene_id"])))
    out["n_up_contrasts"] = up.groupby("gene_id").size()
    out["n_down_contrasts"] = dn.groupby("gene_id").size()
    out[["n_up_contrasts", "n_down_contrasts"]] = \
        out[["n_up_contrasts", "n_down_contrasts"]].fillna(0).astype(int)
    out["best_up_logFC"] = up.groupby("gene_id")["logFC"].max()
    out["best_up_FDR"] = up.groupby("gene_id")["FDR"].min()
    out["up_contrast_tissues"] = (up.sort_values("FDR").groupby("gene_id")["tissue"]
                                  .apply(lambda s: "; ".join(s.head(3))))
    return out


def main():
    args = parse_args()
    os.makedirs(args.outdir, exist_ok=True)
    ann = parse_annotation(args.gtf)
    print(f"{len(ann):,} Wolbachia loci in {os.path.basename(args.gtf)}")

    seqs = read_fasta(args.proteins) if args.proteins else {}
    if args.write_faa:
        n = 0
        with open(args.write_faa, "w") as out:
            for _, r in ann.iterrows():
                s = seqs.get(r["protein_id"])
                if s:
                    out.write(f">{r['locus_tag']}|{r['protein_id']}\n{s}\n")
                    n += 1
        print(f"Wrote {n:,} proteins to {args.write_faa}")
        return

    sp = read_signalp(args.signalp) if args.signalp else {}
    ann["signal_peptide"] = ann["protein_id"].map(sp) if sp else np.nan
    ann["tm_helices_kd"] = [count_tm(seqs[p]) if p in seqs else np.nan
                            for p in ann["protein_id"]]
    ann["aa_length"] = [len(seqs[p]) if p in seqs else np.nan for p in ann["protein_id"]]

    def klass(r):
        if pd.isna(r["signal_peptide"]) and pd.isna(r["tm_helices_kd"]):
            return "unknown"
        if r["signal_peptide"] is True:
            return "surface_periplasmic"
        if r["tm_helices_kd"] and r["tm_helices_kd"] > 0:
            return "membrane"
        return "cytoplasmic_or_T4SS"
    ann["secretion_class"] = ann.apply(klass, axis=1)

    prod = ann["product"].fillna("").str.lower()
    ann["euk_domain"] = [next((v for k, v in EUK_DOMAINS.items() if k in p), "") for p in prod]
    unchar = (prod.str.contains("hypothetical|DUF|uncharacteri") &
              (ann["secretion_class"] == "cytoplasmic_or_T4SS") &
              (ann["aa_length"].fillna(0) >= args.min_aa))
    # A Sec signal peptide is host-acting evidence of its own kind: the protein
    # reaches the bacterial surface, which is the host interface. Without its own
    # tier these loci scored 0 and the shortlist filter dropped them, however
    # strong their differential expression.
    ann["evidence_tier"] = np.where(ann["old"].isin(list(PUBLISHED)), "published",
                            np.where(ann["old"].isin(PREDICTED), "predicted",
                             np.where(ann["euk_domain"] != "", "domain",
                              np.where(ann["signal_peptide"] == True, "surface",
                               np.where(unchar, "uncharacterised", "")))))
    ann["evidence_note"] = ann["old"].map(PUBLISHED).fillna("")
    ann.loc[(ann["evidence_tier"] == "predicted") & (ann["evidence_note"] == ""),
            "evidence_note"] = "Predicted T4SS effector, yeast growth defect (Rice et al. 2017)"
    ann.loc[ann["evidence_tier"] == "domain", "evidence_note"] = \
        "Eukaryote-like domain: " + ann.loc[ann["evidence_tier"] == "domain", "euk_domain"]
    ann.loc[ann["evidence_tier"] == "uncharacterised", "evidence_note"] = \
        "Large soluble protein of unknown function (no signal peptide, no TM helix)"
    ann.loc[ann["evidence_tier"] == "surface", "evidence_note"] = \
        "Sec signal peptide: exported to the envelope / bacterial surface"

    ann = ann.set_index("locus_tag")
    ann["germline_cpm"] = wmel_cpm(os.path.join(args.de_dir, "germline", "pseudobulk"), "germline")
    ann["all_cells_cpm"] = wmel_cpm(os.path.join(args.de_dir, "pseudobulk"))
    g = de_table(os.path.join(args.de_dir, "germline", "de", "wolbachia_de_all.tsv.gz"), "germline")
    p = de_table(os.path.join(args.de_dir, "wolbachia_followup", "wolbachia_de_checked.tsv.gz"),
                 "pooled_all_tissues")
    ann = ann.join(g.add_prefix("germline_")).join(p.add_prefix("pooled_"))
    across = de_across_tissues(os.path.join(args.de_dir, "wolbachia_followup",
                                            "wolbachia_de_sig_checked.tsv"))
    if len(across):
        ann = ann.join(across)
        print(f"Per-cell-type DE: {int((ann['n_up_contrasts'].fillna(0) > 0).sum())} loci up "
              f"in at least one contrast, "
              f"{int((ann['n_down_contrasts'].fillna(0) > 0).sum())} down")
    else:
        for c in ("n_up_contrasts", "n_down_contrasts", "best_up_logFC",
                  "best_up_FDR", "up_contrast_tissues"):
            ann[c] = np.nan
        print("No per-cell-type DE table found; scoring on germline and pooled only")

    # Direction: a factor that substitutes for Mei-P26 is expected to be no
    # lower in the mutant host, so "up" is the a priori interesting direction.
    lfc = ann["germline_logFC"].fillna(ann["pooled_logFC"])
    fdr = ann[["germline_FDR", "pooled_FDR"]].min(axis=1)
    ann["de_direction"] = np.select(
        [fdr.isna(), fdr >= 0.05, lfc > 0], ["not_tested", "n.s.", "up"], default="down")
    ann["best_FDR"] = fdr

    # Where the protein can act, stated rather than buried in the score. A Sec
    # signal peptide means the envelope, so those loci cannot reach cytoplasmic
    # host targets; only T4SS delivery can, and that signal is not predictable
    # from sequence. Filter on this column rather than inferring it from a rank.
    ann["route"] = np.where(ann["secretion_class"] == "surface_periplasmic",
                            "bacterial_surface",
                            np.where(ann["secretion_class"] == "cytoplasmic_or_T4SS",
                                     "host_cytoplasm_possible",
                                     np.where(ann["secretion_class"] == "membrane",
                                              "membrane_or_T4SS", "unknown")))

    # Evidence and differential expression are scored on comparable scales. The
    # previous version weighted tier at 10 per rank against at most 3 for DE, so
    # a domain-only locus with no DE outranked one that is up in six contrasts.
    tier_rank = {"published": 12, "predicted": 9, "domain": 6, "surface": 5,
                 "uncharacterised": 3, "": 0}
    n_up = ann["n_up_contrasts"].fillna(0)
    n_dn = ann["n_down_contrasts"].fillna(0)
    de_score = np.where(n_up > 0, 6 + 2 * np.minimum(n_up, 4),
                 np.where(n_dn > 0, 3 + np.minimum(n_dn, 3),
                  np.where(ann["de_direction"] == "n.s.", 1, 0)))
    ann["evidence_score"] = ann["evidence_tier"].map(tier_rank)
    ann["de_score"] = de_score
    ann["priority_score"] = (ann["evidence_score"] + de_score
                             + np.where(ann["germline_cpm"].fillna(0) >= args.min_cpm, 2, 0)
                             + np.where(ann["route"] != "bacterial_surface", 2, 0))

    for c in ("n_up_contrasts", "n_down_contrasts"):
        ann[c] = ann[c].fillna(0).astype(int)
    cols = ["old", "gene", "product", "aa_length", "secretion_class", "route",
            "signal_peptide", "tm_helices_kd", "euk_domain", "evidence_tier",
            "evidence_note", "germline_cpm", "all_cells_cpm",
            "germline_logFC", "germline_FDR", "pooled_logFC", "pooled_FDR",
            "n_up_contrasts", "n_down_contrasts", "best_up_logFC", "best_up_FDR",
            "up_contrast_tissues", "de_direction", "best_FDR",
            "evidence_score", "de_score", "priority_score"]
    cols = [c for c in cols if c in ann.columns]
    out = ann[cols].sort_values(["priority_score", "germline_cpm"], ascending=False)
    out.to_csv(os.path.join(args.outdir, "wolbachia_secretion_all.tsv"), sep="\t")

    # Admit a locus on EITHER host-acting evidence OR robust differential
    # expression. The old filter required a tier, so a strongly up-regulated
    # locus with no eukaryote-like domain was excluded no matter its statistics.
    short = out[((out["evidence_tier"] != "") | (out["n_up_contrasts"] > 0)
                 | (out["n_down_contrasts"] >= 2))
                & (out["germline_cpm"].fillna(0) >= args.min_cpm)]
    short.to_csv(os.path.join(args.outdir, "wolbachia_secretion_shortlist.tsv"), sep="\t")

    ups = out[(out["n_up_contrasts"] > 0) & (out["evidence_tier"] != "")]
    if len(ups):
        print(f"\nSecretion-competent and UP in a Mei-P26 host ({len(ups)} loci), "
              "the set a substitute-for-Mei-P26 hypothesis needs:")
        print(ups[["old", "product", "route", "evidence_tier", "n_up_contrasts",
                   "best_up_logFC", "best_up_FDR"]].head(15).round(3).to_string())
        print("  route = bacterial_surface cannot reach cytoplasmic host targets; "
              "only T4SS delivery can, and that signal is not sequence-predictable.")

    if "signal_peptide" in out and out["signal_peptide"].notna().any():
        print("\nSecretion class, all loci:")
        print(out["secretion_class"].value_counts().to_string())
    print(f"\nShortlist ({len(short)} loci: host-acting evidence and "
          f">= {args.min_cpm:g} CPM in infected germline)")
    show = ["old", "product", "secretion_class", "evidence_tier", "germline_cpm",
            "germline_logFC", "pooled_logFC", "de_direction", "priority_score"]
    print(short[show].round(3).to_string())


if __name__ == "__main__":
    main()

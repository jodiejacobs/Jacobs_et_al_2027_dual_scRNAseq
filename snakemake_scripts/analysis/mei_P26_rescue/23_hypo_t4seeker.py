#!/usr/bin/env python
"""Per-protein T4SE predictions from T4Seeker (Li lab, github.com/lijingtju/T4Seeker).

The T4Seeker scripts do not run as shipped: extract_features.py reads two
arguments it never defines, the checkpoint is named model_gpu.pth but loaded as
model.pth, Pse-in-One-2.0 needs Python 2 and sits in a different folder than the
script expects, and T4SeekerPredict.py prints only accuracy metrics (it needs
labels) and writes no per-protein output. This wrapper reuses T4Seeker's own
dataset class, model and weights and replaces the feature step:

  DR features   20 amino-acid counts + 400 residue-pair counts at gaps 1, 2, 3
                (1220 values), a Python 3 version of Pse-in-One-2.0 `nac.py ... DR`
  ESM features  ESM-2 8M, layer 6, first 33 residues x 320, flattened, then the
                1024 columns in models/selected_feature_indices.npy

Features are standardised with a scaler fitted on the proteins in the input
(T4Seeker does the same), so a score depends on the other proteins in the file.
Run it on a whole proteome, never on a handful of candidates.

Proteins shorter than --min_len (33) have no complete ESM block and get no score.
Sequences longer than 1022 residues are truncated for ESM only.

Output TSV: locus_tag, n_aa, p_t4se (softmax, not calibrated), p_t4se_pctile
(rank among scored proteins), call (Y/N/NA), note.
"""
import argparse
import importlib
import os
import sys
from collections import Counter

import numpy as np
import pandas as pd
import torch

# fair-esm checkpoints store an argparse.Namespace; torch >= 2.6 refuses it unless allowed
if hasattr(torch.serialization, "add_safe_globals"):
    torch.serialization.add_safe_globals([argparse.Namespace])

AA = "ACDEFGHIKLMNPQRSTVWY"
PAIRS = [a + b for a in AA for b in AA]      # first letter varies slowest, as in nac.py
ESM_RES, ESM_DIM, ESM_MAXLEN = 33, 320, 1022
VALID = set("ACDEFGHIKLMNPQRSTVWYXBUZO")


def read_fasta(path):
    seqs, name = {}, None
    with open(path) as f:
        for line in f:
            line = line.strip()
            if line.startswith(">"):
                name = line[1:].split()[0]
                seqs[name] = []
            elif name and line:
                seqs[name].append(line)
    out = {}
    for k, v in seqs.items():
        s = "".join(v).upper().rstrip("*")
        out[k] = "".join(c if c in VALID else "X" for c in s)
    return out


def dr_features(seq, max_dis=3):
    """Pse-in-One-2.0 DR for one protein: residue counts, then pair counts per gap."""
    vec = [seq.count(a) for a in AA]
    for d in range(1, max_dis + 1):
        c = Counter(seq[i] + seq[i + d] for i in range(len(seq) - d))
        vec.extend(c.get(p, 0) for p in PAIRS)
    return vec


def esm_features(seqs, repo, device="cpu"):
    """First 33 residue embeddings (ESM-2 8M, layer 6), flattened, one protein at a time."""
    import esm
    path = os.path.join(repo, "ESM_models", "esm2_t6_8M_UR50D.pt")
    model, alphabet = esm.pretrained.load_model_and_alphabet(path)
    model.eval().to(device)
    conv = alphabet.get_batch_converter()
    rows = []
    for i, s in enumerate(seqs):
        _, _, tok = conv([(f"p{i}", s[:ESM_MAXLEN])])
        with torch.no_grad():
            rep = model(tok.to(device), repr_layers=[6])["representations"][6]
        rows.append(rep[0, 1:len(s[:ESM_MAXLEN]) + 1][:ESM_RES].reshape(-1).cpu().numpy())
    return np.vstack(rows)


def load_t4seeker(repo):
    sys.path.insert(0, repo)
    os.chdir(repo)                         # T4Seeker loads models/... relative to cwd
    return importlib.import_module("T4SeekerPredict")


def predict(seqs, dr, esm_sel, repo, weights=None, device="cpu", batch=32):
    """Return softmax P(T4SE) for each row, using T4Seeker's dataset, model and weights."""
    t4 = load_t4seeker(repo)

    class DS(t4.ProteinDataset):           # T4Seeker's tokenising and standardising, arrays in
        def __init__(self, seqs, dr, esm_sel):
            self.kmers, self.max_len = 1, 1000
            self.seqs = np.array(seqs)
            self.labels = np.zeros(len(seqs), dtype=int)
            self.vocab = self.create_vocab_by_amino_acids()
            x = np.concatenate([dr, esm_sel], axis=1).astype(np.float32)
            self.features = self.standardize(x)

    ds = DS(seqs, dr, esm_sel)
    dl = torch.utils.data.DataLoader(ds, batch_size=batch, shuffle=False)
    model = t4.T4Seeker().to(device)
    w = weights or next(p for p in (os.path.join(repo, "models", "model.pth"),
                                    os.path.join(repo, "models", "model_gpu.pth"))
                        if os.path.exists(p))
    ckpt = torch.load(w, map_location=device)
    model.load_state_dict(ckpt["model_state_dict"])
    model.eval()
    probs = []
    with torch.no_grad():
        for x, feats, _ in dl:
            out, _ = model(x.to(device), feats.to(device))
            probs.append(torch.softmax(out, dim=1)[:, 1].cpu().numpy())
    return np.concatenate(probs)


def main():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--repo", required=True, help="clone of lijingtju/T4Seeker")
    p.add_argument("--fasta", required=True, help="protein FASTA, ideally a whole proteome")
    p.add_argument("--out", required=True)
    p.add_argument("--weights", default=None)
    p.add_argument("--min_len", type=int, default=33)
    p.add_argument("--min_proteins", type=int, default=200,
                   help="refuse to score fewer proteins than this (scaler is fitted on the input)")
    a = p.parse_args()

    fasta, out, repo = os.path.abspath(a.fasta), os.path.abspath(a.out), os.path.abspath(a.repo)
    seqs = read_fasta(fasta)
    ids = [k for k, s in seqs.items() if len(s) >= a.min_len]
    if len(ids) < a.min_proteins:
        sys.exit(f"only {len(ids)} proteins long enough to score; T4Seeker standardises features "
                 f"on its input, so give it a proteome (--min_proteins {a.min_proteins})")
    device = "cuda" if torch.cuda.is_available() else "cpu"
    print(f"{len(seqs)} proteins, {len(ids)} scored on {device}")

    s_list = [seqs[k] for k in ids]
    dr = np.array([dr_features(s) for s in s_list], dtype=np.float32)
    esm_all = esm_features(s_list, repo, device)
    idx = np.load(os.path.join(repo, "models", "selected_feature_indices.npy"))
    prob = predict(s_list, dr, esm_all[:, idx], repo, a.weights, device)

    res = pd.DataFrame({"locus_tag": ids, "n_aa": [len(s) for s in s_list], "p_t4se": prob})
    res["p_t4se_pctile"] = res["p_t4se"].rank(pct=True)
    res["call"] = np.where(res["p_t4se"] >= 0.5, "Y", "N")
    res["note"] = np.where(res["n_aa"] > ESM_MAXLEN, "truncated_for_esm", "")
    skipped = pd.DataFrame({"locus_tag": [k for k in seqs if k not in set(ids)]})
    skipped["n_aa"] = skipped["locus_tag"].map(lambda k: len(seqs[k]))
    skipped["call"], skipped["note"] = "NA", f"shorter_than_{a.min_len}_aa"
    res = pd.concat([res, skipped], ignore_index=True)
    os.makedirs(os.path.dirname(out), exist_ok=True)
    res.to_csv(out, sep="\t", index=False)
    print(f"{(res['call'] == 'Y').sum()} of {len(ids)} scored proteins called T4SE (p >= 0.5)")


if __name__ == "__main__":
    main()

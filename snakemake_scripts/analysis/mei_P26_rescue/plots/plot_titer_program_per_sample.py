#!/usr/bin/env python
"""Figure 7: the titer-program effect in each sample, infected against uninfected.

Plots results/rescue/titer_programs/titer_program_per_sample.tsv.

With two replicates per genotype no single p-value carries much weight, so the
evidence is that the same sign appears in every infected sample and not in the
uninfected ones. This figure shows that directly: every sample is a point, so
a reader can count the replication rather than take a summary statistic on
trust. Fisher z is plotted rather than the correlation because z is what the
test averages, and the two are monotonic, so nothing is hidden by the choice.

  A  per program, every sample as a point, infected and uninfected offset
  B  the paired view for the programs with the largest infected effect, each
     sample labelled, so an outlying replicate is visible rather than averaged
     away

Two versions are written. "robust" keeps programs where every infected sample
agrees in sign; "all" keeps every program.

Usage:
  python plot_titer_program_per_sample.py \
      --input results/rescue/titer_programs/titer_program_per_sample.tsv \
      --outdir figures/fig7
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import empty_figure, COL, INFECTION, WIDE, apply_style, panel_letter, save
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig7_titer_programs_per_sample")
    p.add_argument("--top_n", type=int, default=5,
                   help="programs shown sample by sample in panel B")
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def short(p):
    return str(p).replace("Program_", "P")


def is_inf(s):
    return s.astype(str) != "uninf"


def build(d, args, suffix, note):
    d = d.copy()
    d["fisher_z"] = pd.to_numeric(d["fisher_z"], errors="coerce")
    d["infected"] = is_inf(d["infection_status"])
    order = (d[d["infected"]].groupby("program")["fisher_z"].mean()
             .sort_values().index.tolist())
    print(f"\n[{suffix}] {len(order)} programs, "
          f"{d.loc[d['infected'], 'sample'].nunique()} infected and "
          f"{d.loc[~d['infected'], 'sample'].nunique()} uninfected samples; {note}")

    fig = plt.figure(figsize=(WIDE, 2.7))
    gs = fig.add_gridspec(1, 2, wspace=0.26, width_ratios=[1.6, 1.0])

    # ---- A: every program, every sample ----
    ax = fig.add_subplot(gs[0, 0])
    x = np.arange(len(order))
    for grp, off, c, lab in ((True, -0.17, INFECTION["wMel"], "infected"),
                             (False, 0.17, INFECTION["uninf"], "uninfected")):
        sub = d[d["infected"] == grp]
        for i, p in enumerate(order):
            v = sub.loc[sub["program"] == p, "fisher_z"].to_numpy()
            if not len(v):
                continue
            j = (np.random.default_rng(i).random(len(v)) - 0.5) * 0.12
            ax.scatter(np.full(len(v), i + off) + j, v, s=11, color=c,
                       edgecolors="none", zorder=3,
                       label=lab if i == 0 else None)
            ax.plot([i + off - 0.11, i + off + 0.11], [v.mean()] * 2,
                    color=COL["black"], lw=1.0, zorder=4)
    ax.axhline(0, color=COL["grey"], lw=0.6, ls="--")
    ax.set_xticks(x, [short(p) for p in order])
    ax.set_ylabel("Fisher z (titer vs usage | depth)")
    ax.set_title("Each point is one sample, bar is the mean", pad=4)
    ax.legend(loc="upper left", fontsize=6, ncol=2)
    ax.margins(x=0.02)
    panel_letter(ax, "A", dx=-0.09)

    # ---- B: the largest effects, sample by sample ----
    ax = fig.add_subplot(gs[0, 1])
    top = (d[d["infected"]].groupby("program")["fisher_z"].mean().abs()
           .sort_values().tail(args.top_n).index.tolist())
    samples = sorted(d["sample"].unique())
    marks = ["o", "s", "^", "D", "v", "P", "X", "*"]
    for k, sm in enumerate(samples):
        sub = d[(d["sample"] == sm) & (d["program"].isin(top))]
        if not len(sub):
            continue
        inf = bool(sub["infected"].iloc[0])
        xs = [top.index(p) for p in sub["program"]]
        ax.scatter(xs, sub["fisher_z"], s=22, marker=marks[k % len(marks)],
                   color=INFECTION["wMel"] if inf else INFECTION["uninf"],
                   edgecolors=COL["black"], linewidths=0.3, zorder=3,
                   label=sm.replace("_", " "))
    ax.axhline(0, color=COL["grey"], lw=0.6, ls="--")
    ax.set_xticks(range(len(top)), [short(p) for p in top])
    ax.set_ylabel("Fisher z")
    ax.set_title("Largest effects, by sample", pad=4)
    ax.legend(loc="center left", bbox_to_anchor=(1.01, 0.5), fontsize=5.5,
              handletextpad=0.3, labelspacing=0.3)
    ax.margins(x=0.12)
    panel_letter(ax, "B")

    fig.suptitle(note, y=1.04, fontsize=6.5, color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")

    summ = (d.groupby(["program", "infected"])["fisher_z"]
            .agg(["mean", "min", "max", "size"]).reset_index())
    summ = summ[summ["program"].isin(top)].sort_values(["program", "infected"])
    summ["infected"] = np.where(summ["infected"], "infected", "uninfected")
    print(summ.to_string(index=False))


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    d["fisher_z"] = pd.to_numeric(d["fisher_z"], errors="coerce")
    inf = d[is_inf(d["infection_status"])]
    sign_ok = inf.groupby("program")["fisher_z"].apply(
        lambda s: len(set(np.sign(s))) == 1 and len(s) > 1)
    keep = set(sign_ok[sign_ok].index)
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        if v == "robust":
            sub = d[d["program"].isin(keep)]
            note = "programs where every infected sample agrees in sign"
        else:
            sub = d
            note = "every program tested, sign-consistent or not"
        if not len(sub):
            empty_figure(args.outdir, f"{args.stem}_{v}",
                         f"No rows survive the {v} filter.\n"
                         f"An empty set here is a result, not a failure:\n"
                         f"say so in the text rather than citing the\n"
                         f"unfiltered figure as though it passed.")
            continue
        build(sub, args, v, note)


if __name__ == "__main__":
    main()

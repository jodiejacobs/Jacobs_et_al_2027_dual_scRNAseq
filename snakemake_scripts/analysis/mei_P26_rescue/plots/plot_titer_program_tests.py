#!/usr/bin/env python
"""Figure 7: significance of the titer-program relationship, and its controls.

Plots results/rescue/titer_programs/titer_program_tests.tsv.

Two things have to be shown together or the result is overstated. The first is
the permutation p-value floor: with N permutations the smallest attainable
p-value is 1/(N+1), and several programs sit exactly on it, so their FDR is a
bound rather than an estimate. The panel marks those rather than printing a
number that looks sharper than the test can be. The second is the uninfected
control, which is what separates a titer effect from a depth or ambient
artifact, so it is plotted beside the effect and not relegated to a table.

  A  effect against significance. Mean within-sample partial Spearman in
     infected samples against -log10(FDR), with floored programs ringed
  B  infected against uninfected mean effect, the control that matters. A
     program on the diagonal is explained by something other than titer
  C  replication: how many infected samples agree in sign, per program

Two versions are written. "robust" keeps programs where every infected sample
agrees in sign, which at two replicates per genotype is the standard the paper
leans on; "all" keeps every program tested.

Usage:
  python plot_titer_program_tests.py \
      --input results/rescue/titer_programs/titer_program_tests.tsv \
      --outdir figures/fig7
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import empty_figure, COL, WIDE, apply_style, panel_letter, save
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig7_titer_program_tests")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--label_top", type=int, default=6,
                   help="how many programs to name, by absolute effect")
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def short(p):
    return str(p).replace("Program_", "P")


def build(d, args, suffix, note):
    d = d.copy()
    for c in ("mean_r_infected", "mean_r_uninf", "fdr", "pvalue", "perm_p_floor",
              "ambient_fdr", "n_samples_same_sign", "n_infected_samples",
              "perm_z_score"):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    d["label"] = d["program"].map(short)
    # A p-value equal to the floor is a bound, not an estimate.
    floored = np.isclose(d["pvalue"], d["perm_p_floor"], rtol=1e-6)
    d["floored"] = floored
    sig = d["fdr"] < args.fdr
    # ambient_* is a Welch test of infected against uninfected sample effects.
    # Significant there means the effect is NOT reproduced in uninfected samples,
    # which is the control passing. A program that fails it is explained by
    # something other than titer and should not be read as a rescue signal.
    d["ambient_ok"] = (d["ambient_fdr"] < args.fdr if "ambient_fdr" in d.columns
                       else True)
    print(f"\n[{suffix}] {len(d)} programs, {int(sig.sum())} at FDR < {args.fdr}, "
          f"{int((sig & ~d['ambient_ok']).sum())} of those fail the uninfected "
          f"control, {int(floored.sum())} on the permutation floor; {note}")

    fig = plt.figure(figsize=(WIDE, 2.5))
    gs = fig.add_gridspec(1, 3, wspace=0.42, width_ratios=[1.15, 1.0, 1.0])

    # ---- A: effect against strength relative to the null ----
    # Not -log10(FDR): most programs sit on the permutation floor, so that axis
    # would collapse to one horizontal line and imply equal evidence. The
    # permutation z-score is continuous and is what the floored p came from.
    ax = fig.add_subplot(gs[0, 0])
    if "perm_z_score" in d.columns and d["perm_z_score"].notna().any():
        y = d["perm_z_score"].abs()
        ylab, ytitle = "|permutation z|", "Effect and strength vs null"
    else:
        y = -np.log10(d["fdr"].clip(lower=1e-12))
        ylab, ytitle = "$-$log$_{10}$ FDR", "Effect and significance"
    ax.scatter(d.loc[~sig, "mean_r_infected"], y[~sig], s=16,
               color=COL["lightgrey"], edgecolors="none", zorder=1)
    s = d[sig]
    cols = [COL["vermillion"] if v > 0 else COL["blue"] for v in s["mean_r_infected"]]
    ok = s["ambient_ok"].to_numpy()
    ax.scatter(s.loc[ok, "mean_r_infected"], y[sig][ok], s=26, zorder=3,
               color=[c for c, k in zip(cols, ok) if k],
               edgecolors=COL["black"], linewidths=0.4, label="passes control")
    if (~ok).any():
        ax.scatter(s.loc[~ok, "mean_r_infected"], y[sig][~ok], s=26, zorder=3,
                   facecolors="white", linewidths=0.8,
                   edgecolors=[c for c, k in zip(cols, ok) if not k],
                   label="same effect uninfected")
    ax.legend(loc="lower right", fontsize=5.5)
    for i in d["mean_r_infected"].abs().sort_values().tail(args.label_top).index:
        ax.annotate(d.at[i, "label"], (d.at[i, "mean_r_infected"], y[i]),
                    textcoords="offset points", xytext=(4, 3), fontsize=6.5)
    ax.axvline(0, color=COL["grey"], lw=0.6, ls="--")
    ax.set_xlabel("mean partial Spearman, infected")
    ax.set_ylabel(ylab)
    ax.set_title(ytitle, pad=4)
    ax.margins(x=0.16, y=0.16)
    panel_letter(ax, "A")

    # ---- B: the uninfected control ----
    ax = fig.add_subplot(gs[0, 1])
    lim = float(np.nanmax(np.abs(
        np.r_[d["mean_r_infected"], d["mean_r_uninf"]]))) * 1.18
    ax.plot([-lim, lim], [-lim, lim], color=COL["grey"], lw=0.6, ls="--", zorder=1)
    ax.axhline(0, color=COL["lightgrey"], lw=0.6, zorder=0)
    ax.axvline(0, color=COL["lightgrey"], lw=0.6, zorder=0)
    ax.scatter(d.loc[~sig, "mean_r_uninf"], d.loc[~sig, "mean_r_infected"],
               s=16, color=COL["lightgrey"], edgecolors="none", zorder=2)
    ax.scatter(s.loc[ok, "mean_r_uninf"], s.loc[ok, "mean_r_infected"], s=26,
               zorder=3, color=[c for c, k in zip(cols, ok) if k],
               edgecolors=COL["black"], linewidths=0.4)
    if (~ok).any():
        ax.scatter(s.loc[~ok, "mean_r_uninf"], s.loc[~ok, "mean_r_infected"], s=26,
                   zorder=3, facecolors="white", linewidths=0.8,
                   edgecolors=[c for c, k in zip(cols, ok) if not k])
    for _, r in s.reindex(s["mean_r_infected"].abs().sort_values()
                          .tail(args.label_top).index).iterrows():
        ax.annotate(r["label"], (r["mean_r_uninf"], r["mean_r_infected"]),
                    textcoords="offset points", xytext=(4, 3), fontsize=6.5)
    ax.set_xlim(-lim, lim)
    ax.set_ylim(-lim, lim)
    ax.set_xlabel("mean effect, uninfected")
    ax.set_ylabel("mean effect, infected")
    ax.set_title("Off the diagonal is the result", pad=4)
    panel_letter(ax, "B")

    # ---- C: replication across infected samples ----
    ax = fig.add_subplot(gs[0, 2])
    r = d.sort_values(["n_samples_same_sign", "mean_r_infected"])
    yy = np.arange(len(r))
    ax.barh(yy, r["n_samples_same_sign"], height=0.7,
            color=[COL["vermillion"] if v > 0 else COL["blue"]
                   for v in r["mean_r_infected"]])
    n_inf = int(r["n_infected_samples"].max())
    ax.axvline(n_inf, color=COL["black"], lw=0.7, ls="--")
    ax.text(n_inf, len(r) - 0.2, " all", fontsize=6, color=COL["grey"], va="top")
    ax.set_yticks(yy, r["label"], fontsize=6)
    ax.set_xlabel("infected samples agreeing in sign")
    ax.set_xlim(0, n_inf + 0.6)
    ax.set_xticks(range(n_inf + 1))
    ax.set_title("Replication", pad=4)
    panel_letter(ax, "C")

    fig.suptitle(note, y=1.04, fontsize=6.5, color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")

    show = [c for c in ("program", "mean_r_infected", "mean_r_uninf",
                        "n_samples_same_sign", "pvalue", "fdr", "ambient_fdr",
                        "ambient_ok") if c in d.columns]
    print(d[sig].sort_values("mean_r_infected")[show].to_string(index=False))
    bad = d.loc[sig & ~d["ambient_ok"], "label"].tolist()
    if bad:
        print(f"{', '.join(bad)} show the same effect in uninfected samples, so "
              "titer does not explain them. Do not cite them as rescue signal.")
    if floored.any():
        print(f"{int((sig & floored).sum())} significant programs sit on the "
              f"permutation floor (p = {d['perm_p_floor'].min():.4g}), so their "
              "FDR is an upper bound rather than an estimate. Raise "
              "--permutations in rule rescue_titer_programs to sharpen them.")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        if v == "robust":
            sub = d[pd.to_numeric(d["n_samples_same_sign"], errors="coerce")
                    == pd.to_numeric(d["n_infected_samples"], errors="coerce")]
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

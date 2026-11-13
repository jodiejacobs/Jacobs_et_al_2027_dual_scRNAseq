#!/usr/bin/env python
"""Figure 3: tests for shifts in cell-type composition.

Plots results/rescue/composition/composition_tests.tsv.

Composition is fitted as a two-by-two on logit proportions, which gives four
contrasts plus their interaction, and reading any one of them alone invites a
wrong conclusion. The genotype effect in uninfected animals is the defect; the
genotype effect with wMel is what remains of it; the interaction is the formal
rescue test. The panel lays the contrasts out together per cell type so the
reader can see that a cell type whose defect shrinks also has an interaction
term, rather than taking the rescue claim on the strength of two separately
significant tests.

A significant interaction is the claim. A shrinking genotype effect without
one is a description, and the panel marks which is which.

Two versions are written. "robust" keeps cell types with a significant
interaction, where the rescue test itself passed; "all" keeps every cell type.

Usage:
  python plot_composition_tests.py \
      --input results/rescue/composition/composition_tests.tsv \
      --outdir figures/fig3
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import (empty_figure, COL, LOW_CONFIDENCE, WIDE, apply_style, ct_short,
                      fdr_stars, note, panel_letter, save)
import matplotlib.pyplot as plt

ORDER = ["genotype_in_uninf", "genotype_in_wMel", "interaction",
         "infection_in_OreR", "infection_in_MeiP26"]
NICE = {"genotype_in_uninf": "the defect\n(mei-P26 vs OreR, uninfected)",
        "genotype_in_wMel": "what remains\n(mei-P26 vs OreR, wMel)",
        "interaction": "rescue test\n(interaction)",
        "infection_in_OreR": "wMel effect\nin OreR",
        "infection_in_MeiP26": "wMel effect\nin mei-P26"}


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig3_composition_tests")
    p.add_argument("--fdr", type=float, default=0.05)
    p.add_argument("--drop_labels", default=",".join(LOW_CONFIDENCE))
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    p.add_argument("--fit_label", default=None,
                   help="named in the figure's subtitle, so the moderated and "
                        "unmoderated fits cannot be mistaken for each other. "
                        "Inferred from the input filename when not given")
    return p.parse_args()


def build(d, args, suffix, note_txt):
    cons = [c for c in ORDER if c in set(d["contrast"])]
    eff = d.pivot_table(index="tissue", columns="contrast", values="logit_diff")
    fdr = d.pivot_table(index="tissue", columns="contrast", values="fdr")
    se = d.pivot_table(index="tissue", columns="contrast", values="se")
    key = "genotype_in_uninf" if "genotype_in_uninf" in eff.columns else cons[0]
    order = eff[key].sort_values().index
    eff, fdr, se = eff.loc[order, cons], fdr.loc[order, cons], se.loc[order, cons]
    print(f"\n[{suffix}] {len(eff)} cell types x {len(cons)} contrasts; {note_txt}")

    fig = plt.figure(figsize=(WIDE, 3.2))
    gs = fig.add_gridspec(1, 2, wspace=0.08, width_ratios=[1.0, 1.35])
    y = np.arange(len(eff))

    # ---- A: effect sizes with their uncertainty ----
    ax = fig.add_subplot(gs[0, 0])
    for i, c in enumerate(("genotype_in_uninf", "genotype_in_wMel")):
        if c not in eff.columns:
            continue
        off = (i - 0.5) * 0.30
        col = COL["vermillion"] if i == 0 else COL["blue"]
        ax.errorbar(eff[c], y + off, xerr=1.96 * se[c], fmt="none",
                    ecolor=COL["grey"], elinewidth=0.8, capsize=1.3, zorder=2)
        ax.scatter(eff[c], y + off, s=20, color=col, zorder=3,
                   edgecolors=COL["black"], linewidths=0.3,
                   label="uninfected" if i == 0 else "with wMel")
    ax.axvline(0, color=COL["black"], lw=0.7)
    ax.set_yticks(y, ct_short(eff.index), fontsize=6)
    ax.set_xlabel("mei-P26 vs OreR, logit proportion")
    ax.set_title("Does wMel shrink the defect?", pad=4)
    ax.legend(loc="lower left", fontsize=6)
    note(ax, "bars are 95% CI", y=-0.14)
    panel_letter(ax, "A", dx=-0.60)

    # ---- B: every contrast, with the rescue test marked ----
    ax = fig.add_subplot(gs[0, 1])
    lim = float(np.nanmax(np.abs(eff.to_numpy()))) or 1.0
    im = ax.imshow(eff.to_numpy(), cmap="RdBu_r", vmin=-lim, vmax=lim, aspect="auto")
    for i in range(eff.shape[0]):
        for j in range(eff.shape[1]):
            s = fdr_stars(fdr.iat[i, j])
            if s:
                ax.text(j, i, s, ha="center", va="center", fontsize=6.5,
                        color=COL["black"])
    if "interaction" in cons:
        k = cons.index("interaction")
        ax.axvline(k - 0.5, color=COL["black"], lw=0.9)
        ax.axvline(k + 0.5, color=COL["black"], lw=0.9)
    ax.set_xticks(range(len(cons)), [NICE.get(c, c) for c in cons],
                  fontsize=5.5, rotation=0)
    ax.set_yticks(y, [])
    cb = fig.colorbar(im, ax=ax, fraction=0.03, pad=0.015)
    cb.set_label("logit difference", fontsize=6)
    cb.ax.tick_params(labelsize=5.5)
    ax.set_title("All contrasts (boxed = the rescue test)", pad=4)
    panel_letter(ax, "B", dx=-0.05)

    fig.suptitle(f"{note_txt} · {args.fit_label}", y=1.03, fontsize=6.5,
                 color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")

    if "interaction" in fdr.columns:
        hit = fdr.index[fdr["interaction"] < args.fdr]
        print(f"Significant interaction (the rescue test) in {len(hit)} cell types: "
              f"{', '.join(ct_short(hit)) or 'none'}")
        shrink = eff.index[(eff.get("genotype_in_wMel", eff[key]).abs()
                            < eff[key].abs() * 0.5) & (~eff.index.isin(hit))]
        if len(shrink):
            print(f"{', '.join(ct_short(shrink))} show the defect halving with wMel "
                  "but no significant interaction. Describe these, do not call "
                  "them rescued.")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t")
    for c in ("logit_diff", "se", "pvalue", "fdr"):
        d[c] = pd.to_numeric(d[c], errors="coerce")
    drop = {s.strip() for s in args.drop_labels.split(",") if s.strip()}
    d = d[~d["tissue"].isin(drop)]
    print(f"Labels excluded as low-confidence: {sorted(drop) or 'none'}")
    # The two tables have identical columns and very different standard errors,
    # so the figure has to say which fit produced it.
    if args.fit_label is None:
        args.fit_label = ("propeller, variance moderated across cell types"
                          if "propeller" in os.path.basename(args.input)
                          else "OLS per cell type, 4 residual df, unmoderated")
    print(f"Fit: {args.fit_label}")
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        if v == "robust":
            keep = set(d.loc[(d["contrast"] == "interaction")
                             & (d["fdr"] < args.fdr), "tissue"])
            sub = d[d["tissue"].isin(keep)]
            txt = "cell types with a significant interaction, the formal rescue test"
        else:
            sub = d
            txt = "every cell type tested"
        if not len(sub):
            empty_figure(args.outdir, f"{args.stem}_{v}",
                         f"No rows survive the {v} filter.\n"
                         f"An empty set here is a result, not a failure:\n"
                         f"say so in the text rather than citing the\n"
                         f"unfiltered figure as though it passed.")
            continue
        build(sub, args, v, txt)


if __name__ == "__main__":
    main()

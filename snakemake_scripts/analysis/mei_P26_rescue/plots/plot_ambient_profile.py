#!/usr/bin/env python
"""Figure 2: does the wMel signal in uninfected cells look like real infection?

Plots results/rescue/uninfected_floor/ambient_profile.tsv.gz.

The depth panel asks where the uninfected signal comes from. This one asks what
it is made of. If the uninfected signal is ambient RNA spilled from infected
droplets, its gene composition should resemble the infected profile, because it
is the same transcripts. If instead it is dominated by a handful of
high-abundance loci (rRNA above all), it is not a faint copy of infection but a
different thing, and treating it as a scaled-down infection profile would be
wrong.

The panel also states how many loci the comparison actually rests on. After
debris removal the shared set is small, and a correlation computed on a dozen
loci is not evidence however strong it looks.

Two versions are written. "robust" keeps loci seen in both uninfected and
infected cells, which is the only set the comparison is defined on; "all" keeps
every locus with any signal, so the reader sees how much is dropped.

Usage:
  python plot_ambient_profile.py \
      --input results/rescue/uninfected_floor/ambient_profile.tsv.gz \
      --outdir figures/fig2
"""
import argparse
import os
import sys

import numpy as np
import pandas as pd
from scipy import stats

sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
from figstyle import empty_figure, COL, RRNA_PATTERN, WIDE, apply_style, note, panel_letter, save
import matplotlib.pyplot as plt


def parse_args():
    p = argparse.ArgumentParser(description=__doc__,
                                formatter_class=argparse.RawDescriptionHelpFormatter)
    p.add_argument("--input", required=True)
    p.add_argument("--outdir", required=True)
    p.add_argument("--stem", default="fig2_ambient_profile")
    p.add_argument("--min_shared", type=int, default=30,
                   help="below this many shared loci the correlation is not reported "
                        "as evidence")
    p.add_argument("--top_n", type=int, default=12, help="loci named in panel B")
    p.add_argument("--version", choices=("robust", "all", "both"), default="both")
    return p.parse_args()


def build(d, args, suffix, note_txt):
    both = d[(d["uninf_umis"] > 0) & (d["infected_umis"] > 0)]
    enough = len(both) >= args.min_shared
    if len(both) > 2:
        r = stats.spearmanr(np.log1p(both["uninf_cpm"]), np.log1p(both["infected_cpm"]))
        rho, pv = r.statistic, r.pvalue
    else:
        rho, pv = np.nan, np.nan
    print(f"\n[{suffix}] {len(d)} loci, {len(both)} seen in both; "
          f"Spearman ρ = {rho:.3f} (p = {pv:.3g}); {note_txt}")
    if not enough:
        print(f"  Only {len(both)} shared loci, below --min_shared "
              f"{args.min_shared}. Do not cite this correlation as evidence; "
              "say the comparison is underpowered after debris removal.")

    fig = plt.figure(figsize=(WIDE, 2.6))
    gs = fig.add_gridspec(1, 2, wspace=0.30, width_ratios=[1.0, 1.2])

    # ---- A: ambient against infected composition ----
    ax = fig.add_subplot(gs[0, 0])
    is_rrna = (d["product"].astype(str).str.contains(RRNA_PATTERN, case=False, na=False)
               if "product" in d.columns else pd.Series(False, index=d.index))
    x, y = np.log10(d["infected_cpm"] + 1), np.log10(d["uninf_cpm"] + 1)
    ax.scatter(x[~is_rrna], y[~is_rrna], s=6, alpha=0.5, color=COL["blue"],
               edgecolors="none", zorder=2)
    if is_rrna.any():
        ax.scatter(x[is_rrna], y[is_rrna], s=26, color=COL["vermillion"],
                   edgecolors=COL["black"], linewidths=0.4, zorder=3,
                   label=f"rRNA ({int(is_rrna.sum())})")
        ax.legend(loc="upper left", fontsize=6)
    lim = float(max(x.max(), y.max())) * 1.05
    ax.plot([0, lim], [0, lim], color=COL["grey"], lw=0.6, ls="--", zorder=1)
    ax.set_xlim(0, lim)
    ax.set_ylim(0, lim)
    ax.set_xlabel("infected, log$_{10}$(CPM + 1)")
    ax.set_ylabel("uninfected, log$_{10}$(CPM + 1)")
    title = (f"Composition agreement (ρ = {rho:.2f})" if enough
             else f"Underpowered: only {len(both)} shared loci")
    ax.set_title(title, pad=4, fontsize=7.5)
    panel_letter(ax, "A", dx=-0.26)

    # ---- B: what the uninfected signal is actually made of ----
    ax = fig.add_subplot(gs[0, 1])
    top = d.nlargest(args.top_n, "uninf_cpm").iloc[::-1]
    yy = np.arange(len(top))
    share = top["uninf_cpm"] / 1e6
    ax.barh(yy, share, height=0.68, zorder=3,
            color=[COL["vermillion"] if b else COL["blue"]
                   for b in (is_rrna.reindex(top.index).fillna(False))])
    ax.set_yticks(yy, [str(i) for i in top.index], fontsize=6)
    ax.set_xlabel("share of the uninfected wMel signal")
    ax.set_title("What the noise is made of", pad=4)
    ax.margins(x=0.16)
    carried = float(share.sum())
    note(ax, f"top {len(top)} loci carry {carried:.0%} of it", y=-0.24)
    panel_letter(ax, "B", dx=-0.38)

    fig.suptitle(note_txt, y=1.05, fontsize=6.5, color=COL["grey"])
    save(fig, args.outdir, f"{args.stem}_{suffix}")
    print(f"  top {len(top)} loci carry {carried:.1%} of the uninfected signal")


def main():
    args = parse_args()
    apply_style()
    d = pd.read_csv(args.input, sep="\t", index_col=0)
    for c in ("uninf_umis", "infected_umis", "uninf_cpm", "infected_cpm",
              "log2_ambient_enrichment"):
        if c in d.columns:
            d[c] = pd.to_numeric(d[c], errors="coerce")
    versions = ["robust", "all"] if args.version == "both" else [args.version]
    for v in versions:
        if v == "robust":
            sub = d[(d["uninf_umis"] > 0) & (d["infected_umis"] > 0)]
            txt = "loci seen in both uninfected and infected cells"
        else:
            sub = d[(d["uninf_umis"] > 0) | (d["infected_umis"] > 0)]
            txt = "every locus with signal in either group"
        if len(sub) < 3:
            empty_figure(args.outdir, f"{args.stem}_{v}",
                         f"No rows survive the {v} filter.\n"
                         f"An empty set here is a result, not a failure:\n"
                         f"say so in the text rather than citing the\n"
                         f"unfiltered figure as though it passed.")
            continue
        build(sub, args, v, txt)


if __name__ == "__main__":
    main()

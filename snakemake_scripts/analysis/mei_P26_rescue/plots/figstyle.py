#!/usr/bin/env python
"""Shared publication style for the mei-P26 rescue figure panels.

Every plot_*.py script in this directory imports this module, so font sizes,
colours and output conventions stay identical across panels that will sit in
one figure. Import it before creating any figure:

    from figstyle import apply_style, COL, PAL, save, panel_letter, COLUMN, WIDE
    apply_style()

Design choices and why:

  colours      Okabe-Ito, which is distinguishable under all common forms of
               colour vision deficiency and prints legibly in greyscale. Never
               encode a variable by colour alone; pair it with marker or
               position.
  sizes        7 pt body text, 8 pt axis labels. Journals typeset single-column
               figures at 85 mm, so anything smaller than 7 pt is unreadable
               after reduction. COLUMN and WIDE give the two standard widths in
               inches; set height per panel.
  no grid      Gridlines compete with data at figure scale. Where a reference
               value matters, draw one line for it instead.
  spines       Top and right removed. Left and bottom kept, because a panel
               with no frame is hard to read when several sit side by side.
  output       PDF for submission (vector, editable in Illustrator) and PNG at
               400 dpi for slides and drafts. save() writes both and prints the
               paths so a run log records what it produced.
"""
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# Okabe-Ito, colour-vision-deficiency safe.
COL = {
    "black": "#000000",
    "orange": "#E69F00",
    "skyblue": "#56B4E9",
    "green": "#009E73",
    "yellow": "#F0E442",
    "blue": "#0072B2",
    "vermillion": "#D55E00",
    "purple": "#CC79A7",
    "grey": "#999999",
    "lightgrey": "#D9D9D9",
}
PAL = [COL["blue"], COL["vermillion"], COL["green"], COL["orange"],
       COL["purple"], COL["skyblue"], COL["yellow"], COL["black"]]

# Fixed meanings, so the same thing is the same colour in every panel.
GENOTYPE = {"OreR": COL["blue"], "Mei-P26": COL["vermillion"],
            "MeiP26": COL["vermillion"], "mei-P26": COL["vermillion"]}
INFECTION = {"wMel": COL["green"], "uninf": COL["grey"],
             "uninfected": COL["grey"], "infected": COL["green"]}
DIRECTION = {"up": COL["vermillion"], "down": COL["blue"], "ns": COL["lightgrey"]}

# Per-sample colour. Eight libraries need eight colours, and eight unrelated
# hues read as noise, so each genotype-by-infection group keeps its own hue and
# its replicates are separated by shade. The four groups then stay legible at a
# glance while individual libraries remain distinguishable.
GROUP_BASE = {"OreR_uninf": COL["skyblue"], "OreR_wMel": COL["blue"],
              "MeiP26_uninf": COL["orange"], "MeiP26_wMel": COL["vermillion"]}
GROUP_ORDER = ["OreR_uninf", "OreR_wMel", "MeiP26_uninf", "MeiP26_wMel"]


def _lighten(hexcol, frac):
    c = hexcol.lstrip("#")
    rgb = [int(c[i:i + 2], 16) for i in (0, 2, 4)]
    return "#%02X%02X%02X" % tuple(int(v + (255 - v) * frac) for v in rgb)


def sample_group(sample):
    """Genotype-by-infection group of a sample id, from the id alone."""
    g = "MeiP26" if "Mei" in str(sample) else "OreR"
    i = "uninf" if "uninf" in str(sample) else "wMel"
    return f"{g}_{i}"


def sample_colours(samples):
    """Return (colour_by_sample, samples_ordered_by_group).

    Use the returned order for axes and legends so the same library is the
    same colour and sits in the same place in every figure.
    """
    uniq = list(dict.fromkeys(str(s) for s in samples))
    order = [s for g in GROUP_ORDER for s in sorted(
        {x for x in uniq if sample_group(x) == g})]
    out = {}
    for g in GROUP_ORDER:
        mem = [s for s in order if sample_group(s) == g]
        for k, s in enumerate(mem):
            out[s] = _lighten(GROUP_BASE[g], 0.46 * k / max(1, len(mem) - 1))
    return out, order


def sample_label(sample):
    return str(sample).replace("_10x", "").replace("_", " ")

COLUMN = 3.35   # inches, 85 mm single column
WIDE = 6.85     # inches, 174 mm double column


def apply_style():
    plt.rcParams.update({
        "figure.dpi": 150,
        "savefig.dpi": 400,
        "savefig.bbox": "tight",
        "savefig.pad_inches": 0.02,
        "pdf.fonttype": 42,      # TrueType, so text stays editable in Illustrator
        "ps.fonttype": 42,
        "font.family": "sans-serif",
        "font.sans-serif": ["Helvetica", "Arial", "DejaVu Sans"],
        "font.size": 7,
        "axes.titlesize": 8,
        "axes.labelsize": 8,
        "axes.linewidth": 0.6,
        "axes.spines.top": False,
        "axes.spines.right": False,
        "axes.grid": False,
        "axes.prop_cycle": plt.cycler(color=PAL),
        "xtick.labelsize": 7,
        "ytick.labelsize": 7,
        "xtick.major.width": 0.6,
        "ytick.major.width": 0.6,
        "xtick.major.size": 2.5,
        "ytick.major.size": 2.5,
        "legend.fontsize": 7,
        "legend.frameon": False,
        "legend.handlelength": 1.2,
        "legend.borderaxespad": 0.2,
        "lines.linewidth": 1.0,
        "lines.markersize": 3.5,
        "patch.linewidth": 0.5,
        "errorbar.capsize": 1.5,
    })


def panel_letter(ax, letter, dx=-0.28, dy=1.04):
    """Bold panel letter outside the axes, positioned as journals expect."""
    ax.text(dx, dy, letter, transform=ax.transAxes, fontsize=9,
            fontweight="bold", va="top", ha="left")


def save(fig, outdir, stem, formats=("pdf", "png")):
    os.makedirs(outdir, exist_ok=True)
    for ext in formats:
        path = os.path.join(outdir, f"{stem}.{ext}")
        fig.savefig(path)
        print(f"  wrote {path}")
    plt.close(fig)


def empty_figure(outdir, stem, message):
    """Write a placeholder when a version's filter leaves nothing.

    An empty robust set is a result, not a failure: it means no row in the
    table survives the standard the paper's claims rest on. The pipeline
    declares both the robust and the all figure as outputs, so writing this
    keeps the rule honest instead of letting a missing file look like a crash,
    and the message says on the page itself what was filtered away.
    """
    fig, ax = plt.subplots(figsize=(COLUMN, 1.5))
    ax.text(0.5, 0.5, message, ha="center", va="center", wrap=True,
            transform=ax.transAxes, fontsize=7.5, color=COL["black"])
    ax.set_axis_off()
    save(fig, outdir, stem)
    print(f"  (empty: {message})")


def wrap_labels(labels, width=28):
    """Cell-type names are long; wrap rather than rotate where there is room."""
    import textwrap
    return ["\n".join(textwrap.wrap(str(s), width)) for s in labels]


def shorten(labels, width=30):
    return [str(s) if len(str(s)) <= width else str(s)[:width - 1] + "…"
            for s in labels]


# Atlas cell-type names run to 50 characters, which no axis can carry. These are
# the forms to use on tick labels; spell the full name in the legend once.
CELLTYPE_SHORT = {
    "16-cell germline cyst in germarium region 2a and 2b": "16-cell cyst",
    "choriogenic main body follicle cell and corpus luteum": "choriogenic FC",
    "choriogenic main body follicle cell St. 12": "choriogenic St.12",
    "choriogenic main body follicle cell St. 14": "choriogenic St.14",
    "central main body follicle cell ca. St. 6-8": "main body FC St.6-8",
    "main body follicle cell ca. until St. 5": "main body FC <St.5",
    "posterior terminal follicle cell ca. St. 5-8": "post. terminal FC",
    "dorsal appendage forming follicle cell": "dorsal appendage FC",
    "follicle stem cell and prefollicle cell": "FSC / prefollicle",
    "germ cell stage 4 and later": "germ cell St.4+",
    "post-mitotic germ cell early 16-cell cyst": "post-mitotic cyst",
    "polar follicle cell": "polar FC",
    "stalk follicle cell": "stalk FC",
    "stretch follicle cell": "stretch FC",
    "ovarian sheath muscle": "sheath muscle",
    "adult tracheal cell": "tracheal",
    "young germ cell": "young germ cell",
    "escort cell": "escort",
    "oviduct": "oviduct",
    "unannotated": "unannotated",
    "artefact": "artefact",
    "pooled_all_tissues": "pooled",
    "germline": "germline",
}


def ct_short(labels):
    return [CELLTYPE_SHORT.get(str(s), str(s)) for s in labels]


# Loci that dominate the wMel read budget but are not genes of interest; keep
# them out of gene-level panels and report them separately.
RRNA_PATTERN = r"ribosomal RNA|rRNA|23S|16S|5S"

# Labels that are not cell types. Counting them inflates every "replicated in
# N cell types" statement, which is how "GroES in six contexts" became three.
NON_CELLTYPE = ("pooled_all_tissues", "germline", "all")
LOW_CONFIDENCE = ("artefact", "unannotated")


def truthy(s):
    import pandas as pd
    return pd.Series(s).astype(str).str.upper().isin(("TRUE", "T", "1", "YES"))


def add_robust_args(p):
    """Every DE-reading plot script offers the same two versions.

    `robust` is the composite flag from 09_wolbachia_followup.R: replicate
    consistent, leave-one-out sign stable, significant after renormalising
    without the dominant rRNA loci, and no zero-count samples. The last two
    matter because rRNA is 85% of the read budget. Plotting only the robust set
    hides what was filtered, and plotting everything overstates the evidence,
    so each script writes both and the stem says which is which.
    """
    p.add_argument("--version", choices=("robust", "all", "both"), default="both",
                   help="robust = composite flag only; all = every significant row; "
                        "both = write one figure of each (default)")
    p.add_argument("--drop_labels", default=",".join(LOW_CONFIDENCE),
                   help="labels excluded as low-confidence; '' keeps them")
    return p


def de_versions(d, version, robust_col="robust"):
    """Yield (suffix, frame, note) for each requested version."""
    has = robust_col in d.columns
    out = []
    if version in ("robust", "both"):
        sub = d[truthy(d[robust_col]).to_numpy()] if has else d
        note = ("robust only: replicate-consistent, leave-one-out stable, "
                "significant without the dominant rRNA loci, no zero counts"
                if has else "no robust column in this table")
        out.append(("robust", sub, note))
    if version in ("all", "both"):
        out.append(("all", d, "every significant row, robust or not"))
    return out


def fdr_stars(p):
    if p is None:
        return ""
    try:
        p = float(p)
    except (TypeError, ValueError):
        return ""
    return "***" if p < 0.001 else "**" if p < 0.01 else "*" if p < 0.05 else ""


def note(ax, text, y=-0.22):
    """One line of reading under a panel, for what the eye should take from it."""
    ax.text(0, y, text, transform=ax.transAxes, fontsize=6,
            color=COL["grey"], va="top", ha="left")

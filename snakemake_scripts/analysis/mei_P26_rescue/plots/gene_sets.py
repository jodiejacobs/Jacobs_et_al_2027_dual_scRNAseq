"""Gene sets shared by the rescue analysis and its plots.

One definition, imported everywhere, so the eggshell genes in a heatmap, a
volcano and a per-cell score cannot drift apart. Plain Python with no
plotting imports, so analysis scripts can use it too.

EGGSHELL groups the follicle-cell secretory genes by when they come on, since
the stage at which a program is restored says how far oogenesis gets:

  vitelline membrane   stages 8-10, the first eggshell layer
  yolk                 stages 8-10, vitellogenesis; not eggshell, but made by
                       the same follicle cells in the same window, so an
                       independent marker of reaching that stage
  chorion              stages 11-14, choriogenesis

dec (defective chorion 1) is expressed from stage 9 and its products go into
both the vitelline membrane and the chorion; it sits with chorion here because
that is where its loss-of-function phenotype is. Femcoat and Fcp3C are
eggshell-associated but their onset is not well enough established to place,
so they are left out rather than guessed.

Symbols are as they appear in the FlyBase 2025_04 annotation and in
de/host_de_all.tsv.gz. FBGN holds the ids resolved from that table on
2026-11-16, so code that only has FlyBase ids can use the same set.
"""

EGGSHELL = {
    "vitelline membrane": ["Vm26Aa", "Vm26Ab", "Vm26Ac", "Vm32E", "Vm34Ca",
                           "Vml", "psd"],
    "yolk": ["Yp1", "Yp2", "Yp3"],
    "chorion": ["Cp15", "Cp16", "Cp18", "Cp19", "Cp36", "Cp38", "Cp7Fa",
                "Cp7Fb", "Cp7Fc", "dec", "yellow-g", "yellow-g2"],
}
EGGSHELL_ORDER = ["vitelline membrane", "yolk", "chorion"]

FBGN = {
    "Vm26Aa": "FBgn0003979", "Vm26Ab": "FBgn0003980", "Vm26Ac": "FBgn0086266",
    "Vm32E": "FBgn0014076", "Vm34Ca": "FBgn0003983", "Vml": "FBgn0085362",
    "psd": "FBgn0086265",
    "Yp1": "FBgn0004045", "Yp2": "FBgn0005391", "Yp3": "FBgn0004047",
    "Cp15": "FBgn0000355", "Cp16": "FBgn0000356", "Cp18": "FBgn0000357",
    "Cp19": "FBgn0000358", "Cp36": "FBgn0000359", "Cp38": "FBgn0000360",
    "Cp7Fa": "FBgn0014464", "Cp7Fb": "FBgn0014465", "Cp7Fc": "FBgn0014466",
    "dec": "FBgn0000427", "yellow-g": "FBgn0041709", "yellow-g2": "FBgn0035328",
}

# Germline by markers, as 07_germline.py defines it: a cell with at least two
# of these detected. Used as the ambient reference for eggshell transcripts,
# since germ cells do not make eggshell. The atlas label "young germ cell" is
# not usable for that in the mutant, where 07 found it is mostly tj+/eya+
# somatic cells.
GERMLINE_MARKERS = {"vas": "FBgn0283442", "aub": "FBgn0000146",
                    "AGO3": "FBgn0250816", "orb": "FBgn0004882",
                    "osk": "FBgn0003015"}

# Genes drawn on the UMAP and tabulated per cell type (19_marker_expression.py,
# rule rescue_plot_umap_markers). Label -> id. WD_ ids are wMel loci and are
# normalised to the cell's wMel UMIs rather than its host UMIs.
MARKERS = {
    "mei-P26": "FBgn0026206",   # the mutated gene
    "Vm26Aa": "FBgn0003979",    # vitelline membrane, St.8-10 follicle cells
    "bam": "FBgn0000158",       # germline differentiation, cystoblast onward
    "orb": "FBgn0004882",       # oocyte determination
    "Sxl": "FBgn0264270",       # germline sex determination; Wolbachia rescues
                                # Sxl hypomorphs (Starr & Cline 2002)
    "TomO": "WD_RS05785",       # wMel host RNA manipulator
    "CifA": "WD_RS02835",       # wMel cytoplasmic incompatibility factor A
}

# Follicle-cell labels with enough cells per sample to summarise. Matching is
# by substring so the broad and fine annotations both work.
FOLLICLE_KEY = "follicle cell"
GERMLINE_LABELS = ("young germ cell",
                   "16-cell germline cyst in germarium region 2a and 2b")


def eggshell_symbols():
    return [s for c in EGGSHELL_ORDER for s in EGGSHELL[c]]


def eggshell_class():
    return {s: c for c, ss in EGGSHELL.items() for s in ss}


def is_germline(label):
    l = str(label).lower()
    return "germ" in l or "cyst" in l

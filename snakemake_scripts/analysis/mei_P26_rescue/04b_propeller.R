#!/usr/bin/env Rscript
# Cell-type composition with propeller (speckle), the moderated version of the
# test in 04_cell_composition.py.
#
# 04_cell_composition.py fits each cell type's logit proportion independently by
# OLS over the 8 samples, which leaves 4 residual df per cell type and a
# variance estimated from those 4 numbers alone. propeller fits the same design
# but moderates each cell type's variance toward the trend across all cell
# types with limma's empirical Bayes step, which is where its extra power comes
# from. With 19 cell types and 4 df each, that borrowing is worth a good deal.
#
# This does not depend on speckle. propeller's whole method is its variance
# stabilising transform followed by lmFit + contrasts.fit + eBayes, and limma
# is already in the edger env, so pulling in speckle would mean installing
# Seurat and SingleCellExperiment to gain five lines of arithmetic. Both the
# transform and the fit below are copied from speckle's own source
# (R/getTransformedProps.R and R/propeller.ttest.R, Oshlack/speckle), including
# its defaults transform = "logit", robust = TRUE, trend = FALSE. speckle's
# propeller() takes only a single grouping factor anyway, so a 2x2 with custom
# contrasts has to go through lmFit directly whichever way round it is done.
#
# Note the pseudocount differs from 04_cell_composition.py. propeller uses
# (count + 0.5) / rowSums(count + 0.5), so the denominator grows by half a cell
# per cell type; 04_cell_composition.py uses (count + 0.5) / (total + 1). The
# difference is small and matters only for rare cell types, but it is a second
# reason the two tables will not agree exactly.
#
# Design and contrasts match 04_cell_composition.py and 03_edger_de.R exactly:
#   ~ 0 + group   (group = genotype x infection, cell-means coding)
#     infection_in_OreR     OreR_wMel   - OreR_uninf
#     infection_in_MeiP26   MeiP26_wMel - MeiP26_uninf
#     genotype_in_uninf     MeiP26_uninf - OreR_uninf    (the defect)
#     genotype_in_wMel      MeiP26_wMel  - OreR_wMel     (what is left of it)
#     interaction           infection_in_MeiP26 - infection_in_OreR
#   The interaction is the rescue test.
#
# Output columns match composition_tests.tsv, so plot_composition_tests.py runs
# on either table unchanged. "logit_diff" is the contrast estimate on the
# transformed scale, as there.
#
# Usage:
#   Rscript 04b_propeller.R --counts results/rescue/composition/composition_cell_counts.tsv \
#     --outdir results/rescue/composition [--transform logit] [--no_robust]

suppressPackageStartupMessages({
  library(limma)
})

REF_GENO <- "OreR"
REF_INF  <- "uninf"

parse_args <- function(argv) {
  o <- list(counts = NA, outdir = NA, transform = "logit", robust = TRUE,
            stem = "composition_tests_propeller")
  i <- 1
  while (i <= length(argv)) {
    a <- argv[i]
    if (a == "--counts")         { o$counts <- argv[i + 1]; i <- i + 2 }
    else if (a == "--outdir")    { o$outdir <- argv[i + 1]; i <- i + 2 }
    else if (a == "--transform") { o$transform <- argv[i + 1]; i <- i + 2 }
    else if (a == "--stem")      { o$stem <- argv[i + 1]; i <- i + 2 }
    else if (a == "--no_robust") { o$robust <- FALSE; i <- i + 1 }
    else stop("unknown argument: ", a, call. = FALSE)
  }
  if (is.na(o$counts) || is.na(o$outdir))
    stop("--counts and --outdir are required", call. = FALSE)
  if (!o$transform %in% c("logit", "asin"))
    stop("--transform must be logit or asin", call. = FALSE)
  o
}

# Sample ids carry the design, as everywhere else in this pipeline.
geno_of <- function(s) ifelse(grepl("Mei", s), "MeiP26", "OreR")
inf_of  <- function(s) ifelse(grepl("uninf", s), "uninf", "wMel")

main <- function() {
  opt <- parse_args(commandArgs(trailingOnly = TRUE))
  dir.create(opt$outdir, showWarnings = FALSE, recursive = TRUE)

  # composition_cell_counts.tsv is samples x cell types, written by
  # 04_cell_composition.py from the same confident-label cells.
  m <- read.delim(opt$counts, row.names = 1, check.names = FALSE)
  m <- as.matrix(m[, vapply(m, is.numeric, logical(1)), drop = FALSE])
  storage.mode(m) <- "double"
  if (nrow(m) < 4 || ncol(m) < 2)
    stop("expected a samples x cell types count matrix, got ",
         nrow(m), " x ", ncol(m), call. = FALSE)
  cat(sprintf("%d samples x %d cell types, %s cells\n",
              nrow(m), ncol(m), format(sum(m), big.mark = ",")))

  # speckle::getTransformedProps, verbatim, for a samples x clusters table.
  # Transposed at the end because limma wants features in rows.
  proportions <- m / rowSums(m)
  if (opt$transform == "asin") {
    prop_trans <- asin(sqrt(proportions))
  } else {
    pseudo <- (m + 0.5) / rowSums(m + 0.5)
    prop_trans <- log(pseudo / (1 - pseudo))
  }
  props <- list(TransformedProps = t(prop_trans), Proportions = t(proportions))
  tp <- props$TransformedProps                 # cell types x samples

  sid <- colnames(tp)
  group <- factor(paste0(geno_of(sid), "_", inf_of(sid)))
  need <- c("OreR_uninf", "OreR_wMel", "MeiP26_uninf", "MeiP26_wMel")
  if (!all(need %in% levels(group)))
    stop("expected all four groups, found: ",
         paste(levels(group), collapse = ", "), call. = FALSE)
  group <- factor(group, levels = need)
  cat("Samples per group: ",
      paste(sprintf("%s=%d", levels(group), table(group)), collapse = ", "),
      "\n", sep = "")

  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  cm <- makeContrasts(
    infection_in_OreR   = OreR_wMel   - OreR_uninf,
    infection_in_MeiP26 = MeiP26_wMel - MeiP26_uninf,
    genotype_in_uninf   = MeiP26_uninf - OreR_uninf,
    genotype_in_wMel    = MeiP26_wMel  - OreR_wMel,
    interaction         = (MeiP26_wMel - MeiP26_uninf) - (OreR_wMel - OreR_uninf),
    levels = design)

  fit <- lmFit(tp, design)
  fit <- contrasts.fit(fit, cm)
  # robust = TRUE down-weights cell types whose variance is far off the trend,
  # so one erratic cell type cannot drag the moderated variance for the rest.
  fit <- eBayes(fit, robust = opt$robust)

  mean_prop <- t(apply(props$Proportions, 1, function(p)
    tapply(p, group, mean)))
  colnames(mean_prop) <- paste0("mean_prop_", colnames(mean_prop))

  out <- do.call(rbind, lapply(colnames(cm), function(cn) {
    tt <- topTable(fit, coef = cn, number = Inf, sort.by = "none")
    data.frame(tissue = rownames(tt), contrast = cn,
               logit_diff = tt$logFC,
               se = fit$coefficients[rownames(tt), cn] / fit$t[rownames(tt), cn],
               t = tt$t, pvalue = tt$P.Value, fdr = tt$adj.P.Val,
               row.names = NULL, stringsAsFactors = FALSE)
  }))
  out <- cbind(out, mean_prop[out$tissue, , drop = FALSE])
  # With robust = TRUE the prior df is per cell type, not a single number: a
  # cell type whose variance is an outlier gets less of the shared estimate.
  dp <- fit$df.prior
  dres <- as.integer(unique(fit$df.residual)[1])
  cat(sprintf(
    "Prior df from the empirical Bayes step: %s, added to the %d residual df each cell type has on its own\n",
    if (length(dp) > 1)
      sprintf("%.1f to %.1f across cell types (median %.1f)",
              min(dp), max(dp), stats::median(dp))
    else sprintf("%.1f", dp),
    dres))

  f <- file.path(opt$outdir, paste0(opt$stem, ".tsv"))
  write.table(out, f, sep = "\t", quote = FALSE, row.names = FALSE)
  cat("wrote ", f, "\n", sep = "")

  sig <- out[out$fdr < 0.05, ]
  if (nrow(sig)) {
    cat("\nFDR < 0.05:\n")
    print(sig[order(sig$contrast, sig$fdr),
              c("tissue", "contrast", "logit_diff", "pvalue", "fdr")],
          row.names = FALSE)
  } else {
    cat("\nNothing at FDR < 0.05. Moderation raises power but cannot ",
        "manufacture a result from 2 replicates per group.\n", sep = "")
  }
  ni <- sum(out$contrast == "interaction" & out$fdr < 0.05)
  cat(sprintf("\n%d cell types have a significant interaction, which is the ",
              ni),
      "rescue test. Compare this against composition_tests.tsv, the ",
      "unmoderated OLS fit: a difference between the two is the empirical ",
      "Bayes step, not new data.\n", sep = "")
}

main()

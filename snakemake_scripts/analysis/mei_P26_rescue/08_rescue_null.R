#!/usr/bin/env Rscript
# Is the host "rescue" larger than selection alone would produce?
#
# 05_rescue_summary.py selects genes on genotype_in_uninf (FDR < fdr) and
# reports rescue_fraction = 1 - logFC(genotype_in_wMel) / logFC(genotype_in_uninf).
# Selecting on significance inflates |logFC(genotype_in_uninf)| (winner's
# curse) while genotype_in_wMel comes from independent samples, so the ratio
# shrinks and rescue_fraction comes out > 0 even with no real rescue.
#
# This script does two things per tissue (same min_cells / min_reps filters
# as 03_edger_de.R):
#
# 1. NULL: relabel infection status at the SAMPLE level within each genotype.
#    Every way to call 2 of a genotype's 4 samples "uninfected" is used,
#    choose(4, 2) = 6 per genotype, so 6 x 6 = 36 global relabellings; the
#    observed one is dropped, leaving 35 null draws per tissue. Refit the same
#    edgeR QL model and recompute the median rescue_fraction and the n of genes
#    called full/partial, then compare the observed value with those 35. The
#    smallest attainable empirical p is 1/36 = 0.028.
#
#    This replaces an earlier scheme that used only the 4 "mixed" assignments
#    per genotype (16 null draws, floor 1/17 = 0.059). That floor sits above
#    0.05, so nothing could reach significance however large the effect, and
#    every p-value it reported was the floor rather than a measurement. Any
#    p-value from a run before 2026-11-10 is not interpretable; regenerate it.
#
# 2. SUPPORTED RESCUE: genes DE in the uninfected mutant (genotype_in_uninf
#    FDR < fdr) with a significant interaction (FDR < fdr) whose sign opposes
#    the mutant effect (wMel moves the mutant toward OreR). Each gene also gets
#    a replicate check: both Mei-P26 wMel samples lie beyond both Mei-P26
#    uninfected samples (logCPM) in the direction of the shift, i.e. the
#    infection effect in the mutant holds in each replicate.
#
# Usage:
#   Rscript 08_rescue_null.R --pseudobulk_dir results/rescue/pseudobulk \
#     --de_dir results/rescue/de --outdir results/rescue/rescue_null \
#     [--min_cells 20] [--min_reps 2] [--fdr 0.05]

suppressPackageStartupMessages({ library(edgeR); library(limma) })

parse_args <- function(argv) {
  opts <- list(min_cells = 20, min_reps = 2, fdr = 0.05)
  for (i in seq(1, length(argv), by = 2)) opts[[sub("^--", "", argv[i])]] <- argv[i + 1]
  for (k in c("min_cells", "min_reps", "fdr")) opts[[k]] <- as.numeric(opts[[k]])
  for (k in c("pseudobulk_dir", "de_dir", "outdir")) if (is.null(opts[[k]])) stop("need --", k)
  opts
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

meta <- read.delim(file.path(opt$pseudobulk_dir, "pseudobulk_meta.tsv"), check.names = FALSE,
                   stringsAsFactors = FALSE)
rownames(meta) <- meta$pb_id
meta$geno <- gsub("[^A-Za-z0-9]", "", meta$genotype)
meta$group <- paste(meta$geno, meta$infection_status, sep = "_")
host <- as.matrix(read.delim(gzfile(file.path(opt$pseudobulk_dir, "pseudobulk_host_counts.tsv.gz")),
                             row.names = 1, check.names = FALSE))[, meta$pb_id]
groups <- c("OreR_uninf", "OreR_wMel", "MeiP26_uninf", "MeiP26_wMel")

# genotype_in_uninf and genotype_in_wMel from one tissue's counts and a group vector
fit_genotype <- function(counts, grp) {
  grp <- factor(grp, levels = groups)
  design <- model.matrix(~ 0 + grp); colnames(design) <- groups
  y <- DGEList(counts)
  y <- y[filterByExpr(y, group = grp), , keep.lib.sizes = FALSE]
  y <- estimateDisp(normLibSizes(y), design)
  fit <- glmQLFit(y, design, robust = TRUE)
  con <- makeContrasts(u = MeiP26_uninf - OreR_uninf, w = MeiP26_wMel - OreR_wMel,
                       i = (MeiP26_wMel - MeiP26_uninf) - (OreR_wMel - OreR_uninf), levels = design)
  tt <- function(k) topTags(glmQLFTest(fit, contrast = con[, k]), n = Inf, sort.by = "none")$table
  u <- tt("u"); w <- tt("w"); i <- tt("i")
  list(y = y, tab = data.frame(gene_id = rownames(u), logFC_u = u$logFC, FDR_u = u$FDR,
                               logFC_w = w$logFC, FDR_w = w$FDR, logFC_i = i$logFC, FDR_i = i$FDR))
}

rescue_stats <- function(tab, fdr) {
  sel <- tab[tab$FDR_u < fdr, ]
  rf <- 1 - sel$logFC_w / sel$logFC_u
  c(n_mutant_DE = nrow(sel), median_rf = if (nrow(sel)) median(rf) else NA,
    frac_partial_or_full = if (nrow(sel)) mean(rf >= 0.25) else NA,
    n_interaction = sum(tab$FDR_i < fdr))
}

# balanced relabellings of infection within each genotype: pick one true-uninf
# and one true-wMel sample to call "uninf"
# Every way to call 2 of a genotype's 4 samples "uninfected": choose(4, 2) = 6
# per genotype, so 6 x 6 = 36 global relabellings. One of those is the observed
# labelling, leaving 35 null draws and a smallest attainable p of 1/36 = 0.028.
#
# The earlier version used only the 4 "mixed" assignments per genotype, each
# pseudo-group holding one true-uninfected and one true-wMel sample. That gave
# 16 null draws and a floor of 1/17 = 0.059, which sits ABOVE 0.05, so no
# tissue could reach significance however large its effect and every reported
# p-value was the floor rather than a measurement. The two assignments it left
# out per genotype, calling both true-uninfected or both true-wMel samples
# "uninfected", are exchangeable under the null of no infection effect and are
# exactly the draws that make the test able to reject.
label_sets <- function(samp_meta, geno) {
  s <- samp_meta[samp_meta$geno == geno, ]
  u <- s$sample[s$infection_status == "uninf"]; w <- s$sample[s$infection_status == "wMel"]
  if (length(u) != 2 || length(w) != 2) stop("expected 2 uninf + 2 wMel samples for ", geno)
  sets <- combn(c(u, w), 2)
  observed <- which(apply(sets, 2, function(p) setequal(p, u)))
  stopifnot(length(observed) == 1)
  list(sets = sets, observed = observed)
}
samp_meta <- unique(meta[, c("sample", "geno", "infection_status")])
ls_ore <- label_sets(samp_meta, "OreR")
ls_mei <- label_sets(samp_meta, "MeiP26")
perm_ore <- ls_ore$sets
perm_mei <- ls_mei$sets
n_null <- ncol(perm_ore) * ncol(perm_mei) - 1
cat(sprintf("Permutation set: %d x %d relabellings minus the observed = %d null draws; smallest attainable p = %.4f\n",
            ncol(perm_ore), ncol(perm_mei), n_null, 1 / (n_null + 1)))

null_rows <- list(); obs_rows <- list(); supported <- list()
for (tis in sort(unique(meta$tissue))) {
  m <- meta[meta$tissue == tis & meta$n_cells >= opt$min_cells, ]
  reps <- table(factor(m$group, levels = groups))
  if (any(reps < opt$min_reps) || nrow(m) != 8) next
  cat(sprintf("%-55s", tis))
  obs <- fit_genotype(host[, rownames(m)], m$group)
  obs_rows[[tis]] <- data.frame(tissue = tis, t(rescue_stats(obs$tab, opt$fdr)))

  k <- 0
  for (a in seq_len(ncol(perm_ore))) for (b in seq_len(ncol(perm_mei))) {
    if (a == ls_ore$observed && b == ls_mei$observed) next   # that one is the data
    as_uninf <- c(perm_ore[, a], perm_mei[, b])
    g <- paste(m$geno, ifelse(m$sample %in% as_uninf, "uninf", "wMel"), sep = "_")
    k <- k + 1
    null_rows[[paste(tis, k)]] <- data.frame(tissue = tis, perm = k,
                                             t(rescue_stats(fit_genotype(host[, rownames(m)], g)$tab, opt$fdr)))
  }
  cat(" done\n")

  # supported rescue genes + replicate consistency on observed fit
  tab <- obs$tab
  s <- tab[tab$FDR_u < opt$fdr & tab$FDR_i < opt$fdr, ]
  if (nrow(s)) {
    lc <- cpm(obs$y, log = TRUE, prior.count = 2)[s$gene_id, , drop = FALSE]
    grp <- m[colnames(lc), "group"]
    mu <- lc[, grp == "MeiP26_uninf", drop = FALSE]; mw <- lc[, grp == "MeiP26_wMel", drop = FALSE]
    # both Mei-P26 wMel samples beyond both Mei-P26 uninfected, in the direction of the shift
    up <- rowMeans(mw) > rowMeans(mu)
    s$rep_consistent <- ifelse(up, apply(mw, 1, min) > apply(mu, 1, max),
                                   apply(mw, 1, max) < apply(mu, 1, min))
    s$call <- ifelse(sign(s$logFC_i) == -sign(s$logFC_u), "rescue", "exacerbated")
    s$rescue_fraction <- 1 - s$logFC_w / s$logFC_u
    supported[[tis]] <- data.frame(tissue = tis, s)
  }
}

null <- do.call(rbind, null_rows); obs <- do.call(rbind, obs_rows)
summ <- do.call(rbind, lapply(split(null, null$tissue), function(d) {
  o <- obs[obs$tissue == d$tissue[1], ]
  data.frame(tissue = d$tissue[1],
             obs_median_rf = o$median_rf, null_median_rf_mean = mean(d$median_rf, na.rm = TRUE),
             null_median_rf_max = max(d$median_rf, na.rm = TRUE),
             p_median_rf = (1 + sum(d$median_rf >= o$median_rf, na.rm = TRUE)) / (1 + nrow(d)),
             obs_n_interaction = o$n_interaction,
             null_n_interaction_mean = mean(d$n_interaction),
             null_n_interaction_max = max(d$n_interaction),
             p_n_interaction = (1 + sum(d$n_interaction >= o$n_interaction)) / (1 + nrow(d)))
}))
write.table(null, file.path(opt$outdir, "rescue_null_permutations.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
write.table(summ, file.path(opt$outdir, "rescue_null_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ, row.names = FALSE, digits = 3)

sup <- do.call(rbind, supported)
ann_file <- file.path(opt$de_dir, "host_de_all.tsv.gz")
sym <- unique(read.delim(gzfile(ann_file))[, c("gene_id", "symbol")])
sup <- merge(sup, sym, by = "gene_id", all.x = TRUE)
sup <- sup[order(sup$tissue, sup$call, sup$FDR_i), c("tissue", "gene_id", "symbol", "call", "rep_consistent",
            "rescue_fraction", "logFC_u", "FDR_u", "logFC_w", "FDR_w", "logFC_i", "FDR_i")]
write.table(sup, file.path(opt$outdir, "rescue_supported_genes.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nInteraction-supported genes (mutant DE + interaction FDR <", opt$fdr, "):\n")
print(as.data.frame.matrix(table(sup$tissue, paste(sup$call, ifelse(sup$rep_consistent, "rep_ok", "rep_mixed")))))

# figure: observed vs null median rescue fraction and n interaction genes
pdf(file.path(opt$outdir, "rescue_null.pdf"), width = 10, height = 0.35 * nrow(summ) + 2.5)
par(mfrow = c(1, 2), mar = c(4, 22, 2, 1))
ord <- order(summ$obs_median_rf); ti <- summ$tissue[ord]
for (what in c("median_rf", "n_interaction")) {
  ov <- obs[match(ti, obs$tissue), what]
  rng <- range(c(ov, null[[what]]), na.rm = TRUE)
  plot(NA, xlim = rng, ylim = c(0.5, length(ti) + 0.5), yaxt = "n", ylab = "",
       xlab = if (what == "median_rf") "median rescue fraction" else "n interaction genes (FDR<0.05)",
       main = sprintf("grey = %d relabellings; red = observed", n_null), cex.main = 0.8)
  axis(2, at = seq_along(ti), labels = ti, las = 1, cex.axis = 0.7)
  # Take the y length from the data rather than a constant: the number of null
  # draws per tissue is set by the permutation set, so hardcoding it breaks the
  # plot whenever that set changes.
  for (j in seq_along(ti)) {
    v <- null[null$tissue == ti[j], what]
    points(v, rep(j, length(v)), pch = 16, col = "#00000040")
  }
  points(ov, seq_along(ti), pch = 18, col = "#d7301f", cex = 1.6)
}
invisible(dev.off())

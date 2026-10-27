#!/usr/bin/env Rscript
# Pseudobulk differential expression with edgeR quasi-likelihood.
#
# HOST genes (FB*), one model per tissue, all 8 samples:
#   design ~ 0 + group   (group = genotype x infection; cell-means coding)
#   [+ replicate if --block_replicate]
#   contrasts:
#     infection_in_OreR     OreR_wMel   - OreR_uninf
#     infection_in_MeiP26   MeiP26_wMel - MeiP26_uninf
#     infection_main        average of the two above
#     genotype_in_uninf     MeiP26_uninf - OreR_uninf     (mutant phenotype)
#     genotype_in_wMel      MeiP26_wMel  - OreR_wMel      (residual phenotype with wMel)
#     interaction           infection_in_MeiP26 - infection_in_OreR
#   The interaction is "~infection + genotype + infection:genotype": genes whose
#   infection response depends on genotype.
#
# WOLBACHIA genes (GQ*), infected samples only. Library sizes and TMM factors
# come from wMel genes alone, so logFCs describe the composition of the
# Wolbachia transcriptome (per-bacterium expression), not titer -- titer is
# handled by 02_wolbachia_titer.py.
#   per tissue:     ~ 0 + genotype, contrast MeiP26 - OreR
#   across tissues: ~ sample + tissue + tissue:genotype
#                   tissue_effect      genes whose expression varies by tissue
#                   tissue_x_genotype  genes whose tissue pattern differs by genotype
#
# Usage:
#   Rscript 03_edger_de.R --pseudobulk_dir DIR --outdir DIR [--min_cells 20]
#     [--min_reps 2] [--min_wmel_umis 2000] [--block_replicate]
#     [--annotation reference/fbgn_annotation_ID_fb_2025_04.tsv.gz]

suppressPackageStartupMessages({
  library(edgeR)
  library(limma)
})

# ---- args (base R, no optparse dependency) ----
parse_args <- function(argv) {
  opts <- list(min_cells = 20, min_reps = 2, min_wmel_umis = 2000,
               block_replicate = FALSE, annotation = NA, fdr = 0.05)
  i <- 1
  while (i <= length(argv)) {
    key <- sub("^--", "", argv[i])
    if (key == "block_replicate") { opts$block_replicate <- TRUE; i <- i + 1; next }
    opts[[key]] <- argv[i + 1]
    i <- i + 2
  }
  for (k in c("min_cells", "min_reps", "min_wmel_umis", "fdr")) opts[[k]] <- as.numeric(opts[[k]])
  if (is.null(opts$pseudobulk_dir) || is.null(opts$outdir)) stop("need --pseudobulk_dir and --outdir")
  opts
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

slug <- function(x) gsub("^_|_$", "", gsub("[^A-Za-z0-9]+", "_", x))

read_counts <- function(f) {
  m <- read.delim(gzfile(f), row.names = 1, check.names = FALSE)
  as.matrix(m)
}
meta <- read.delim(file.path(opt$pseudobulk_dir, "pseudobulk_meta.tsv"),
                   check.names = FALSE, stringsAsFactors = FALSE)
rownames(meta) <- meta$pb_id
meta$geno <- gsub("[^A-Za-z0-9]", "", meta$genotype)       # "Mei-P26" -> "MeiP26"
meta$group <- paste(meta$geno, meta$infection_status, sep = "_")
meta$replicate <- factor(meta$replicate)
host <- read_counts(file.path(opt$pseudobulk_dir, "pseudobulk_host_counts.tsv.gz"))
wmel <- read_counts(file.path(opt$pseudobulk_dir, "pseudobulk_wolbachia_counts.tsv.gz"))
host <- host[, meta$pb_id]; wmel <- wmel[, meta$pb_id]

ref_geno <- "OreR"
alt_geno <- setdiff(unique(meta$geno), ref_geno)
if (length(alt_geno) != 1) stop("expected exactly two genotypes, got: ", paste(unique(meta$geno), collapse = ", "))
groups <- c(paste0(ref_geno, "_uninf"), paste0(ref_geno, "_wMel"),
            paste0(alt_geno, "_uninf"), paste0(alt_geno, "_wMel"))

symbols <- NULL
if (!is.na(opt$annotation) && file.exists(opt$annotation)) {
  ann <- read.delim(gzfile(opt$annotation), header = FALSE, comment.char = "#",
                    stringsAsFactors = FALSE, quote = "")
  symbols <- setNames(ann$V1, ann$V3)                     # FBgn -> symbol
}
add_symbol <- function(df) {
  df$symbol <- if (is.null(symbols)) NA else unname(symbols[df$gene_id])
  df
}

fit_ql <- function(counts, group_or_design, design) {
  y <- DGEList(counts)
  keep <- filterByExpr(y, group = group_or_design)
  y <- y[keep, , keep.lib.sizes = FALSE]
  y <- normLibSizes(y)
  y <- estimateDisp(y, design)
  list(y = y, fit = glmQLFit(y, design, robust = TRUE))
}

tidy_test <- function(fit, contrast = NULL, coef = NULL, tissue, name) {
  res <- if (is.null(coef)) glmQLFTest(fit, contrast = contrast) else glmQLFTest(fit, coef = coef)
  tt <- topTags(res, n = Inf, sort.by = "none")$table
  out <- data.frame(tissue = tissue, contrast = name, gene_id = rownames(tt),
                    tt, row.names = NULL, check.names = FALSE)
  # coef-based F tests (across-tissue Wolbachia model): logFC sign would be
  # relative to an arbitrary reference tissue, so report the F test only
  if (!is.null(coef)) out$logFC <- NA
  out
}

summarise_de <- function(tab, fdr) {
  do.call(rbind, lapply(split(tab, list(tab$tissue, tab$contrast), drop = TRUE), function(d) {
    sig <- d$FDR < fdr
    data.frame(tissue = d$tissue[1], contrast = d$contrast[1], n_tested = nrow(d),
               n_sig = sum(sig), n_up = sum(sig & d$logFC > 0, na.rm = TRUE),
               n_down = sum(sig & d$logFC < 0, na.rm = TRUE))
  }))
}

# =====================================================================
# HOST
# =====================================================================
cat("\n=== Host DE per tissue ===\n")
host_res <- list(); host_skipped <- c()
for (tis in sort(unique(meta$tissue))) {
  m <- meta[meta$tissue == tis & meta$n_cells >= opt$min_cells, ]
  reps <- table(factor(m$group, levels = groups))
  if (any(reps < opt$min_reps)) {
    host_skipped <- c(host_skipped, sprintf("%s (reps: %s)", tis, paste(names(reps), reps, sep = "=", collapse = ", ")))
    next
  }
  grp <- factor(m$group, levels = groups)
  design <- if (opt$block_replicate) model.matrix(~ 0 + grp + m$replicate) else model.matrix(~ 0 + grp)
  colnames(design)[seq_along(groups)] <- groups
  colnames(design) <- make.names(colnames(design))
  g <- make.names(groups)   # OreR_uninf, OreR_wMel, MeiP26_uninf, MeiP26_wMel
  con_def <- c(
    infection_in_OreR   = paste(g[2], "-", g[1]),
    infection_in_MeiP26 = paste(g[4], "-", g[3]),
    infection_main      = sprintf("(%s - %s + %s - %s)/2", g[2], g[1], g[4], g[3]),
    genotype_in_uninf   = paste(g[3], "-", g[1]),
    genotype_in_wMel    = paste(g[4], "-", g[2]),
    interaction         = sprintf("(%s - %s) - (%s - %s)", g[4], g[3], g[2], g[1]))
  con <- makeContrasts(contrasts = con_def, levels = design)
  colnames(con) <- names(con_def)

  q <- fit_ql(host[, rownames(m)], grp, design)
  cat(sprintf("  %-30s %d samples, %d genes kept, residual df %d\n",
              tis, nrow(m), nrow(q$y), q$fit$df.residual[1]))
  for (cn in colnames(con)) {
    host_res[[paste(tis, cn)]] <- tidy_test(q$fit, contrast = con[, cn], tissue = tis, name = cn)
  }
  lc <- cpm(q$y, log = TRUE, prior.count = 2)
  write.table(data.frame(gene_id = rownames(lc), lc, check.names = FALSE),
              gzfile(file.path(opt$outdir, sprintf("host_logcpm_%s.tsv.gz", slug(tis)))),
              sep = "\t", quote = FALSE, row.names = FALSE)
}
if (length(host_skipped)) cat("  Skipped (< min_reps pseudobulks with >= min_cells):\n   ", paste(host_skipped, collapse = "\n    "), "\n")
host_tab <- add_symbol(do.call(rbind, host_res))
write.table(host_tab, gzfile(file.path(opt$outdir, "host_de_all.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
host_sum <- summarise_de(host_tab, opt$fdr)
write.table(host_sum, file.path(opt$outdir, "host_de_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(host_sum, row.names = FALSE)

# =====================================================================
# WOLBACHIA
# =====================================================================
cat("\n=== Wolbachia DE (infected samples) ===\n")
# ambient check: infected pseudobulk titer vs uninfected titer in the same tissue
amb <- aggregate(cbind(wmel_umis, host_umis) ~ tissue + infection_status, data = meta, FUN = sum)
amb$titer <- amb$wmel_umis / (amb$wmel_umis + amb$host_umis)
amb_w <- reshape(amb[, c("tissue", "infection_status", "titer")], idvar = "tissue",
                 timevar = "infection_status", direction = "wide")
amb_w$infected_over_uninf <- amb_w$titer.wMel / amb_w$titer.uninf
write.table(amb_w, file.path(opt$outdir, "wolbachia_ambient_check.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

mw <- meta[meta$infection_status == "wMel" & meta$n_cells >= opt$min_cells &
             meta$wmel_umis >= opt$min_wmel_umis, ]
wol_res <- list(); wol_skipped <- c()
for (tis in sort(unique(mw$tissue))) {
  m <- mw[mw$tissue == tis, ]
  reps <- table(factor(m$geno, levels = c(ref_geno, alt_geno)))
  if (any(reps < opt$min_reps)) { wol_skipped <- c(wol_skipped, tis); next }
  geno <- factor(m$geno, levels = c(ref_geno, alt_geno))
  design <- model.matrix(~ 0 + geno); colnames(design) <- levels(geno)
  q <- fit_ql(wmel[, rownames(m)], geno, design)
  con <- makeContrasts(contrasts = paste(alt_geno, "-", ref_geno), levels = design)
  wol_res[[tis]] <- tidy_test(q$fit, contrast = con[, 1], tissue = tis,
                              name = paste0(alt_geno, "_vs_", ref_geno))
  cat(sprintf("  %-30s %d samples, %d wMel genes kept\n", tis, nrow(m), nrow(q$y)))
}
if (length(wol_skipped)) cat("  Skipped:", paste(wol_skipped, collapse = ", "), "\n")

# across-tissue model: needs >= 2 tissues passing in every infected sample's genotype
tis_ok <- names(which(tapply(mw$geno, mw$tissue, function(x) all(table(factor(x, levels = c(ref_geno, alt_geno))) >= opt$min_reps))))
mt <- mw[mw$tissue %in% tis_ok, ]
if (length(tis_ok) >= 2) {
  tissue <- factor(mt$tissue); sample <- factor(mt$sample)
  geno <- factor(mt$geno, levels = c(ref_geno, alt_geno))
  design <- model.matrix(~ sample + tissue + tissue:geno)
  ne <- nonEstimable(design)
  if (!is.null(ne)) design <- design[, !colnames(design) %in% ne, drop = FALSE]
  q <- fit_ql(wmel[, rownames(mt)], tissue, design)
  tcoef <- grep("^tissue[^:]*$", colnames(design))
  icoef <- grep(":geno", colnames(design))
  wol_res[["across_tissue_effect"]] <- tidy_test(q$fit, coef = tcoef, tissue = "all", name = "tissue_effect")
  if (length(icoef)) {
    wol_res[["across_tissue_x_genotype"]] <- tidy_test(q$fit, coef = icoef, tissue = "all", name = "tissue_x_genotype")
  }
  lc <- cpm(q$y, log = TRUE, prior.count = 2)
  write.table(data.frame(gene_id = rownames(lc), lc, check.names = FALSE),
              gzfile(file.path(opt$outdir, "wolbachia_logcpm.tsv.gz")),
              sep = "\t", quote = FALSE, row.names = FALSE)
  cat(sprintf("  across-tissue model: %d tissues, %d pseudobulks, residual df %d\n",
              length(tis_ok), nrow(mt), q$fit$df.residual[1]))
}
# multi-coef F tests return logFC.<coef> columns; keep a common column set
common <- c("tissue", "contrast", "gene_id", "logFC", "logCPM", "F", "PValue", "FDR")
wol_tab <- do.call(rbind, lapply(wol_res, function(d) d[, common]))
write.table(wol_tab, gzfile(file.path(opt$outdir, "wolbachia_de_all.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)
wol_sum <- summarise_de(wol_tab, opt$fdr)
write.table(wol_sum, file.path(opt$outdir, "wolbachia_de_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(wol_sum, row.names = FALSE)

writeLines(capture.output(sessionInfo()), file.path(opt$outdir, "sessionInfo.txt"))

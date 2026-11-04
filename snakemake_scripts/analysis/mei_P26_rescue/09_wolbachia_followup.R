#!/usr/bin/env Rscript
# Robustness checks for the per-tissue Wolbachia DE (Mei-P26 vs OreR) from
# 03_edger_de.R, given only 2 infected samples per genotype.
#
# A. Replicate checks per tested gene and tissue:
#      rep_consistent  both Mei-P26 pseudobulks lie beyond both OreR
#                      pseudobulks in the direction of logFC (logCPM)
#      loo_sign_stable logFC keeps its sign when any one of the 4 infected
#                      samples is dropped (mean logCPM difference)
#      n_zero_*        pseudobulks with zero counts per genotype; logFC ~ +-10
#                      with zeros on one side is an on/off call, not a fold change
# B. Normalisation sensitivity: genes taking > --dominant_frac of wMel counts
#    (e.g. highly expressed ncRNAs) can dominate library size and TMM.
#    Refit with library sizes and TMM factors computed WITHOUT those genes,
#    still testing every gene, and compare logFC/FDR with the original fit.
# C. Pooled sample-level test: sum wMel counts over all tissues per infected
#    sample (4 libraries, 2 vs 2). Each gene gets one test on independent
#    samples, unlike the per-tissue tests, which reuse the same 4 samples.
# D. Optional --gff (RefSeq GFF3 for the wMel assembly): adds product and
#    old_locus_tag for each WD_* locus tag.
#
# Usage:
#   Rscript 09_wolbachia_followup.R --pseudobulk_dir results/rescue/pseudobulk \
#     --outdir results/rescue/wolbachia_followup [--min_cells 20] \
#     [--min_wmel_umis 2000] [--min_reps 2] [--fdr 0.05] [--dominant_frac 0.05] \
#     [--gff reference/wMel_genomic.gff.gz]

suppressPackageStartupMessages({ library(edgeR); library(limma) })

parse_args <- function(argv) {
  opts <- list(min_cells = 20, min_wmel_umis = 2000, min_reps = 2, fdr = 0.05,
               dominant_frac = 0.05, gff = NA)
  for (i in seq(1, length(argv), by = 2)) opts[[sub("^--", "", argv[i])]] <- argv[i + 1]
  for (k in c("min_cells", "min_wmel_umis", "min_reps", "fdr", "dominant_frac")) opts[[k]] <- as.numeric(opts[[k]])
  for (k in c("pseudobulk_dir", "outdir")) if (is.null(opts[[k]])) stop("need --", k)
  opts
}
opt <- parse_args(commandArgs(trailingOnly = TRUE))
dir.create(opt$outdir, recursive = TRUE, showWarnings = FALSE)

meta <- read.delim(file.path(opt$pseudobulk_dir, "pseudobulk_meta.tsv"), check.names = FALSE,
                   stringsAsFactors = FALSE)
rownames(meta) <- meta$pb_id
meta$geno <- factor(gsub("[^A-Za-z0-9]", "", meta$genotype), levels = c("OreR", "MeiP26"))
wmel <- as.matrix(read.delim(gzfile(file.path(opt$pseudobulk_dir, "pseudobulk_wolbachia_counts.tsv.gz")),
                             row.names = 1, check.names = FALSE))[, meta$pb_id]
mw <- meta[meta$infection_status == "wMel" & meta$n_cells >= opt$min_cells &
             meta$wmel_umis >= opt$min_wmel_umis, ]

# dominant genes: mean share of wMel counts across passing infected pseudobulks
share <- rowMeans(sweep(wmel[, rownames(mw)], 2, colSums(wmel[, rownames(mw)]), "/"))
dominant <- names(share)[share > opt$dominant_frac]
cat("Dominant wMel genes (mean share >", opt$dominant_frac, "):\n")
print(round(sort(share[dominant], decreasing = TRUE), 3))

fit_geno <- function(counts, geno, norm_genes = NULL) {
  design <- model.matrix(~ 0 + geno); colnames(design) <- levels(geno)
  y <- DGEList(counts)
  y <- y[filterByExpr(y, group = geno), , keep.lib.sizes = FALSE]
  if (is.null(norm_genes)) {
    y <- normLibSizes(y)
  } else {
    ref <- DGEList(counts[setdiff(rownames(counts), norm_genes), , drop = FALSE])
    ref <- normLibSizes(ref)
    y$samples$lib.size <- ref$samples$lib.size
    y$samples$norm.factors <- ref$samples$norm.factors
  }
  y <- estimateDisp(y, design)
  fit <- glmQLFit(y, design, robust = TRUE)
  tt <- topTags(glmQLFTest(fit, contrast = c(-1, 1)), n = Inf, sort.by = "none")$table
  list(y = y, tt = data.frame(gene_id = rownames(tt), logFC = tt$logFC, logCPM = tt$logCPM,
                              PValue = tt$PValue, FDR = tt$FDR))
}

rep_checks <- function(y, tt, geno, samp) {
  lc <- cpm(y, log = TRUE, prior.count = 2)[tt$gene_id, , drop = FALSE]
  raw <- y$counts[tt$gene_id, , drop = FALSE]
  mei <- geno == "MeiP26"; ore <- geno == "OreR"
  up <- tt$logFC > 0
  tt$rep_consistent <- ifelse(up, apply(lc[, mei, drop = FALSE], 1, min) > apply(lc[, ore, drop = FALSE], 1, max),
                                  apply(lc[, mei, drop = FALSE], 1, max) < apply(lc[, ore, drop = FALSE], 1, min))
  loo <- sapply(seq_along(samp), function(j) {
    k <- seq_along(samp) != j
    rowMeans(lc[, k & mei, drop = FALSE]) - rowMeans(lc[, k & ore, drop = FALSE])
  })
  tt$loo_sign_stable <- apply(sign(loo) == sign(tt$logFC), 1, all)
  tt$loo_logFC_min <- apply(loo, 1, min); tt$loo_logFC_max <- apply(loo, 1, max)
  tt$n_zero_MeiP26 <- rowSums(raw[, mei, drop = FALSE] == 0)
  tt$n_zero_OreR <- rowSums(raw[, ore, drop = FALSE] == 0)
  tt
}

res <- list()
for (tis in sort(unique(mw$tissue))) {
  m <- mw[mw$tissue == tis, ]
  if (any(table(m$geno) < opt$min_reps)) next
  o <- fit_geno(wmel[, rownames(m)], m$geno)
  s <- fit_geno(wmel[, rownames(m)], m$geno, norm_genes = dominant)
  tt <- rep_checks(o$y, o$tt, m$geno, m$sample)
  tt$logFC_excl_dominant <- s$tt$logFC[match(tt$gene_id, s$tt$gene_id)]
  tt$FDR_excl_dominant <- s$tt$FDR[match(tt$gene_id, s$tt$gene_id)]
  tt$dominant <- tt$gene_id %in% dominant
  res[[tis]] <- data.frame(tissue = tis, tt)
}
tab <- do.call(rbind, res)

# C. pooled per-sample test
inf <- meta[meta$infection_status == "wMel", ]
pooled <- sapply(split(inf$pb_id, inf$sample), function(ids) rowSums(wmel[, ids, drop = FALSE]))
pmeta <- unique(inf[, c("sample", "geno")]); pmeta <- pmeta[match(colnames(pooled), pmeta$sample), ]
po <- fit_geno(pooled, pmeta$geno)
pt <- rep_checks(po$y, po$tt, pmeta$geno, pmeta$sample)
ps <- fit_geno(pooled, pmeta$geno, norm_genes = dominant)
pt$logFC_excl_dominant <- ps$tt$logFC[match(pt$gene_id, ps$tt$gene_id)]
pt$FDR_excl_dominant <- ps$tt$FDR[match(pt$gene_id, ps$tt$gene_id)]
pt$dominant <- pt$gene_id %in% dominant
tab <- rbind(tab, data.frame(tissue = "pooled_all_tissues", pt))

# D. annotation
if (!is.na(opt$gff) && file.exists(opt$gff)) {
  g <- readLines(if (grepl("\\.gz$", opt$gff)) gzfile(opt$gff) else opt$gff)
  g <- strsplit(g[!startsWith(g, "#")], "\t")
  g <- do.call(rbind, lapply(g[lengths(g) == 9], function(x) c(x[3], x[9])))
  attr <- function(a, key) ifelse(grepl(paste0(key, "="), a),
                                  sub(paste0(".*", key, "=([^;]*).*"), "\\1", a), NA)
  ann <- data.frame(type = g[, 1], locus_tag = attr(g[, 2], "locus_tag"),
                    old_locus_tag = attr(g[, 2], "old_locus_tag"),
                    product = utils::URLdecode(as.character(attr(g[, 2], "product"))))
  gene <- ann[ann$type == "gene" & !is.na(ann$locus_tag), c("locus_tag", "old_locus_tag")]
  prod <- ann[!is.na(ann$product) & !is.na(ann$locus_tag), c("locus_tag", "product")]
  prod <- prod[!duplicated(prod$locus_tag), ]
  ann <- merge(gene[!duplicated(gene$locus_tag), ], prod, by = "locus_tag", all.x = TRUE)
  tab <- merge(tab, ann, by.x = "gene_id", by.y = "locus_tag", all.x = TRUE)
}

tab <- tab[order(tab$tissue, tab$FDR), ]
write.table(tab, gzfile(file.path(opt$outdir, "wolbachia_de_checked.tsv.gz")), sep = "\t", quote = FALSE, row.names = FALSE)

sig <- tab[tab$FDR < opt$fdr, ]
sig$robust <- sig$rep_consistent & sig$loo_sign_stable & !is.na(sig$FDR_excl_dominant) &
  sig$FDR_excl_dominant < opt$fdr & (sig$n_zero_MeiP26 + sig$n_zero_OreR) == 0
write.table(sig, file.path(opt$outdir, "wolbachia_de_sig_checked.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)

summ <- do.call(rbind, lapply(split(tab, tab$tissue), function(d) {
  s <- d[d$FDR < opt$fdr, ]
  data.frame(tissue = d$tissue[1], n_tested = nrow(d), n_sig = nrow(s),
             n_rep_consistent = sum(s$rep_consistent), n_loo_stable = sum(s$loo_sign_stable),
             n_sig_excl_dominant = sum(d$FDR_excl_dominant < opt$fdr, na.rm = TRUE),
             n_sig_both_norms = sum(s$FDR_excl_dominant < opt$fdr, na.rm = TRUE),
             n_onoff = sum(s$n_zero_MeiP26 > 0 | s$n_zero_OreR > 0),
             cor_logFC_norms = round(cor(d$logFC, d$logFC_excl_dominant, use = "complete.obs"), 3))
}))
write.table(summ, file.path(opt$outdir, "wolbachia_followup_summary.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
print(summ, row.names = FALSE)

# genes significant in >= 2 tissues or pooled: candidate list for follow-up
rec <- aggregate(cbind(n_tissues_sig = FDR < opt$fdr, n_up = FDR < opt$fdr & logFC > 0) ~ gene_id,
                 data = tab[tab$tissue != "pooled_all_tissues", ], FUN = sum)
pool <- tab[tab$tissue == "pooled_all_tissues", c("gene_id", "logFC", "FDR", "rep_consistent",
                                                  "FDR_excl_dominant")]
names(pool)[-1] <- paste0("pooled_", names(pool)[-1])
rec <- merge(rec, pool, by = "gene_id", all = TRUE)
extra <- intersect(c("old_locus_tag", "product"), names(tab))
if (length(extra)) rec <- merge(rec, unique(tab[, c("gene_id", extra)]), by = "gene_id", all.x = TRUE)
rec <- rec[order(rec$pooled_FDR), ]
write.table(rec, file.path(opt$outdir, "wolbachia_candidates.tsv"), sep = "\t", quote = FALSE, row.names = FALSE)
cat("\nPooled sample-level test (2 vs 2), genes FDR <", opt$fdr, ":\n")
print(head(rec[!is.na(rec$pooled_FDR) & rec$pooled_FDR < opt$fdr, ], 40), row.names = FALSE)

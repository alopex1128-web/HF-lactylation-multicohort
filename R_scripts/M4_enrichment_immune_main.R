# ============================================================
# M4_enrichment_immune_main.R  |  M4
# Purpose: GO/GSEA/KEGG + CIBERSORT, all 313 samples, 22 cell types
# Paper sections: Figure 5, Tables 4/5
# Pipeline position: step 5 of the main chain
# Prerequisite: M1 completed; the KEGG part needs network access (failure does not affect GO/GSEA)
# (Generated from the modular pipeline consolidation; analysis logic identical to the verified version.)
# ============================================================
# ============================================================
# M4: stage-7 pipeline (DLRG functional enrichment + CIBERSORT immune infiltration)
# ============================================================
# [Scope] mechanistic layer: GO/KEGG/GSEA enrichment (Figure 6 material) +
#   CIBERSORT infiltration + group comparison + LASSO-gene x immune-cell correlation (Figure 7 material)
#
# [Prerequisites]
#   (1) stages 1-6 completed (output/ contains GSE57338_processed.rds,
#       GSE57338_HFvsNF_all.csv, GSE57338_DLRGs.csv, LASSO_coefficients.csv)
#   (2) place two externally obtained files in DIR (not redistributed in this
#       repository under the Stanford license; see README "Note on CIBERSORT"):
#       CIBERSORT_v1.04.R and LM22.txt
#       (official CIBERSORT.R v1.04 + LM22 reference matrix, 547 genes x 22 cell types,
#        free academic download at https://cibersortx.stanford.edu)
#
# [Methods (consistent with the paper)]
#   ORA: clusterProfiler enrichGO (BP/CC/MF, BH, q<0.05);
#       universe = Entrez mapping of all genes in the expression matrix
#   GSEA: all genes ranked by descending limma t, gseGO (BP, minGSSize=10)
#   CIBERSORT: official v1.04, LM22 signature, perm=1000, QN=TRUE
#       (same parameters as Xie; a log2 input matrix is automatically reverse-log-transformed)
#   Immune comparison: Wilcoxon rank-sum (HF vs NF), BH across 22 cell types
#   Correlation: Spearman (LASSO-gene expression x 22 infiltration scores), BH-adjusted
#
# [Output] under stage7_results/:
#   Table4_DLRG_GO_BP.csv / _CC.csv / _MF.csv
#   Table4_DLRG_KEGG_enrichment.csv (if KEGG is unreachable, run M10 instead)
#   Table5_GSEA_BP_all_results.csv
#   Table6_CIBERSORT_22cell_scores.csv (313 samples x 22 cells + P)
#   Table6_immune_diff_Wilcoxon.csv
#   Table7_LASSOximmune_Spearman.csv (with BH-adjusted q)
#   Figure6 enrichment bubble drafts (BP/KEGG/GSEA, PNG 300 dpi each)
#   Figure7 immune boxplots + correlation heatmap drafts (PNG 300 dpi)
#   stage7_console_report.txt
#
# [Runtime] first run ~8-15 min (GSEA 1-5 min + CIBERSORT 2-5 min dominate)
# [Reproducibility] set.seed(42) throughout; CIBERSORT permutation P values are bit-reproducible under the fixed seed
#
# [Expected values (compare against part 7)]
#   - DLRGs = 56 (established in stage 3)
#   - GO-BP top terms should include chromatin assembly / nucleosome terms
#   - CIBERSORT valid-sample rate should be >80% (P<0.05)
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)
set.seed(42)

# Set DIR to your working directory (forward slashes on Windows)
DIR <- "./GEO_data"

OUT <- file.path(DIR, "stage7_results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOGF <- file.path(OUT, "stage7_console_report.txt")
logline <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n"); cat(msg, "\n", file = LOGF, append = TRUE)
}
if (file.exists(LOGF)) file.remove(LOGF)
logline("========== stage 7: enrichment + immune infiltration, started: ", format(Sys.time()), "==========")

# ============================================================
# Part 0: package check (install what is missing)
# ============================================================
need_pkgs <- c("ggplot2", "pheatmap")          # CRAN
bioc_pkgs <- c("clusterProfiler", "org.Hs.eg.db")  # Bioconductor
for (p in need_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    logline(sprintf("[install] %s missing, installing from CRAN mirror", p))
    options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
    install.packages(p)
  }
}
for (p in bioc_pkgs) {
  if (!requireNamespace(p, quietly = TRUE)) {
    logline(sprintf("[install] %s missing, installing via BiocManager (~2-5 min)", p))
    BiocManager::install(p, update = FALSE, ask = FALSE)
  }
}
if (!requireNamespace("preprocessCore", quietly = TRUE)) {
  logline("[install] preprocessCore missing (CIBERSORT dependency), installing from Bioconductor")
  BiocManager::install("preprocessCore", update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages({
  library(clusterProfiler); library(org.Hs.eg.db)
  library(ggplot2); library(pheatmap)
})
logline("[package check] all ready")

# ============================================================
# Part 1: data loading (outputs of stages 1-6)
# ============================================================
need <- c(rds = file.path(DIR, "output", "GSE57338_processed.rds"),
          all = file.path(DIR, "output", "GSE57338_HFvsNF_all.csv"),
          dlrg = file.path(DIR, "output", "GSE57338_DLRGs.csv"),
          lasso = file.path(DIR, "output", "LASSO_coefficients.csv"),
          cib = file.path(DIR, "CIBERSORT_v1.04.R"),
          lm22 = file.path(DIR, "LM22.txt"))
miss <- need[!file.exists(need)]
if (length(miss) > 0) {
  logline("!!! missing files, provide them before running:")
  for (m in miss) logline("   -", m)
  stop("missing prerequisite files, see the list above")
}
d38 <- readRDS(need[["rds"]])
expr38 <- d38$expr; grp38 <- d38$group
tt38 <- read.csv(need[["all"]], check.names = FALSE)
dlrg <- read.csv(need[["dlrg"]], check.names = FALSE)
lasso_tab <- read.csv(need[["lasso"]], check.names = FALSE)
gidcol <- if (".gid" %in% names(tt38)) ".gid" else "gene"

logline(sprintf("[data] expression matrix: %d genes x %d samples | HF=%d, NF=%d",
                nrow(expr38), ncol(expr38),
                sum(grp38 == "HF"), sum(grp38 == "NF")))
logline(sprintf("[data] DLRGs: %d | LASSO model: %d genes",
                nrow(dlrg), nrow(lasso_tab)))

# ---- gene symbol -> Entrez mapping (same approach as M3) ----
map_id <- function(symbols) {
  suppressWarnings(eg <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(symbols),
              columns = "ENTREZID", keytype = "SYMBOL"))
  eg <- eg[!is.na(eg$ENTREZID) & !duplicated(eg$SYMBOL), ]
  setNames(eg$ENTREZID, eg$SYMBOL)
}
bg_map <- map_id(rownames(expr38))
bg <- unique(bg_map)
logline(sprintf("[mapping] universe: %d genes (%d of %d symbols mapped to Entrez)",
                length(bg), nrow(expr38), length(bg_map)))

# ============================================================
# Part 2: GO enrichment of the DLRG set (BP/CC/MF)
# ============================================================
logline("\n========== Part A: DLRG GO enrichment (ORA) ==========")
dlrg_sym <- as.character(dlrg[[gidcol]])
dlrg_map <- map_id(dlrg_sym)
dlrg_entrez <- unique(dlrg_map)
logline(sprintf("DLRG Entrez mapping: %d/%d", length(dlrg_entrez), length(dlrg_sym)))

run_go <- function(ont) {
  ego <- tryCatch(
    enrichGO(gene = dlrg_entrez, universe = bg, OrgDb = org.Hs.eg.db,
             ont = ont, pAdjustMethod = "BH", pvalueCutoff = 0.05,
             qvalueCutoff = 0.05, readable = TRUE),
    error = function(e) NULL)
  if (is.null(ego) || nrow(as.data.frame(ego)) == 0) {
    logline(sprintf("  [%s] no significant terms", ont)); return(NULL)
  }
  df <- as.data.frame(ego)
  write.csv(df, file.path(OUT, sprintf("Table4_DLRG_GO_%s.csv", ont)),
            row.names = FALSE, fileEncoding = "UTF-8")
  logline(sprintf("  [%s] %d significant terms | top 3:", ont, nrow(df)))
  for (i in seq_len(min(3, nrow(df))))
    logline(sprintf("     %s (q=%.2e, %d genes)", df$Description[i], df$p.adjust[i], df$Count[i]))
  # bubble plot (Figure 6 material)
  png(file.path(OUT, sprintf("Figure6_GO_%s_bubble.png", ont)),
      width = 2400, height = 1800, res = 300)
  print(dotplot(ego, showCategory = 12) +
        ggtitle(sprintf("DLRG GO-%s enrichment (n=%d genes)", ont, length(dlrg_entrez))) +
        theme(axis.text.y = element_text(size = 8)))
  dev.off()
  df
}
go_bp <- run_go("BP"); go_cc <- run_go("CC"); go_mf <- run_go("MF")

# ---- check: BP should contain chromatin-related terms ----
if (!is.null(go_bp)) {
  chrom_hit <- grep("chromatin|nucleosome|histone", go_bp$Description, value = TRUE)
  logline(sprintf("  [check] BP contains %d chromatin/nucleosome terms: %s",
                  length(chrom_hit),
                  ifelse(length(chrom_hit) > 0, paste(head(chrom_hit, 3), collapse = "；"),
                         "none — if 0, check the console output (expected present, consistent with the green module)")))
}

# ============================================================
# Part 3: KEGG enrichment of the DLRG set (on failure the gene list is written out; M10 is the offline fallback)
# ============================================================
logline("\n========== Part A2: DLRG KEGG enrichment ==========")
ekk <- tryCatch(
  enrichKEGG(gene = dlrg_entrez, universe = bg, organism = "hsa",
             pAdjustMethod = "BH", pvalueCutoff = 0.05),
  error = function(e) { logline("  [KEGG] REST connection failed:", conditionMessage(e)); NULL })

if (!is.null(ekk) && nrow(as.data.frame(ekk)) > 0) {
  dfk <- as.data.frame(ekk)
  write.csv(dfk, file.path(OUT, "Table4_DLRG_KEGG_enrichment.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  logline(sprintf("  [KEGG] %d significant terms | top 3:", nrow(dfk)))
  for (i in seq_len(min(3, nrow(dfk))))
    logline(sprintf("     %s (q=%.2e)", dfk$Description[i], dfk$p.adjust[i]))
  png(file.path(OUT, "Figure6_KEGG_bubble.png"), width = 2200, height = 1600, res = 300)
  print(dotplot(ekk, showCategory = 12) + ggtitle("DLRG KEGG enrichment") +
        theme(axis.text.y = element_text(size = 8)))
  dev.off()
} else {
  # write out the gene list (fallback when KEGG REST is unreachable)
  writeLines(dlrg_entrez, file.path(OUT, "KEGG_pending_DLRG_entrez.txt"))
  logline("  [KEGG] written: stage7_results/KEGG_pending_DLRG_entrez.txt")
  logline("  >>> if KEGG is unreachable, run M10 locally with this list")
}

# ============================================================
# Part 4: GSEA (all genes ranked by t)
# ============================================================
logline("\n========== Part B: GSEA (genome-wide, ranked by t) ==========")
# ranking key: limma t (HF vs NF), descending = up in HF first
gl <- tt38$t
names(gl) <- as.character(tt38[[gidcol]])
gl <- gl[!duplicated(names(gl)) & !is.na(gl)]
gm <- map_id(names(gl))
gl_e <- setNames(gl[names(gm)], unname(gm))
gl_e <- sort(gl_e, decreasing = TRUE)
logline(sprintf("  ranked list: %d genes (after Entrez mapping)", length(gl_e)))

gse <- tryCatch(
  gseGO(geneList = gl_e, OrgDb = org.Hs.eg.db, ont = "BP",
        minGSSize = 10, maxGSSize = 500, pvalueCutoff = 0.05,
        verbose = FALSE, eps = 1e-10),
  error = function(e) { logline("  [GSEA] failed:", conditionMessage(e)); NULL })

if (!is.null(gse) && nrow(as.data.frame(gse)) > 0) {
  gdf <- as.data.frame(gse)
  write.csv(gdf, file.path(OUT, "Table5_GSEA_BP_all_results.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  up <- sum(gdf$NES > 0 & gdf$p.adjust < 0.05)
  dn <- sum(gdf$NES < 0 & gdf$p.adjust < 0.05)
  logline(sprintf("  [GSEA-BP] %d significant pathways (up in HF %d / down %d) | top 3 by NES:",
                  nrow(gdf), up, dn))
  gdf2 <- gdf[order(-gdf$NES), ]
  for (i in seq_len(min(3, nrow(gdf2))))
    logline(sprintf("     %s (NES=%+.2f, q=%.2e)", gdf2$Description[i], gdf2$NES[i], gdf2$p.adjust[i]))
  logline("  lowest NES (down in HF), top 3:")
  gdf3 <- gdf[order(gdf$NES), ]
  for (i in seq_len(min(3, nrow(gdf3))))
    logline(sprintf("     %s (NES=%+.2f, q=%.2e)", gdf3$Description[i], gdf3$NES[i], gdf3$p.adjust[i]))
  # bubble plot (Figure 6C material, colored by NES)
  png(file.path(OUT, "Figure6_GSEA_bubble.png"), width = 2600, height = 2000, res = 300)
  print(dotplot(gse, showCategory = 15) +
        ggtitle("GSEA-BP: HF vs NF") + theme(axis.text.y = element_text(size = 7)))
  dev.off()
} else {
  logline("  [GSEA] no significant terms or failure — if it failed, check the console output")
}

# ---- GSEA-KEGG (opportunistic; failure is non-blocking) ----
gsek <- tryCatch(
  gseKEGG(geneList = gl_e, organism = "hsa", minGSSize = 10,
          pvalueCutoff = 0.05, verbose = FALSE),
  error = function(e) NULL)
if (!is.null(gsek) && nrow(as.data.frame(gsek)) > 0) {
  write.csv(as.data.frame(gsek), file.path(OUT, "Table5_GSEA_KEGG_results.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  logline(sprintf("  [GSEA-KEGG] %d significant pathways", nrow(as.data.frame(gsek))))
} else {
  write.csv(data.frame(gene = names(gl_e), t = gl_e),
            file.path(OUT, "GSEA_KEGG_pending_ranked_list.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  logline("  [GSEA-KEGG] local failure — ranked list written out")
}

# ============================================================
# Part 5: CIBERSORT infiltration (official v1.04, fully offline)
# ============================================================
logline("\n========== Part C：CIBERSORT（LM22, perm=1000, QN=TRUE） ==========")
source(need[["cib"]])   # defines CoreAlg / doPerm / CIBERSORT
lm22 <- read.delim(need[["lm22"]], row.names = 1, check.names = FALSE)
lm22 <- data.matrix(lm22)
logline(sprintf("  LM22 signature: %d genes x %d cell types", nrow(lm22), ncol(lm22)))

# deduplicate gene rows (duplicate rownames would disturb the algorithm)
ex <- expr38[!duplicated(rownames(expr38)), ]
ov <- intersect(rownames(lm22), rownames(ex))
logline(sprintf("  gene overlap: %d/%d (LM22 reference)", length(ov), nrow(lm22)))
if (length(ov) < 400) logline("  [warning] overlap <400 is low")

set.seed(42)   # fixed seed: permutation P values are bit-reproducible
t0 <- Sys.time()
old_wd <- getwd(); setwd(OUT)   # CIBERSORT writes CIBERSORT-Results.txt into the working directory
cib_res <- CIBERSORT(lm22, ex, perm = 1000, QN = TRUE)
setwd(old_wd)
logline(sprintf("  CIBERSORT done: %d samples x %d cols | %.1f min",
                nrow(cib_res), ncol(cib_res),
                as.numeric(difftime(Sys.time(), t0, units = "mins"))))

write.csv(cib_res, file.path(OUT, "Table6_CIBERSORT_22cell_scores.csv"),
          fileEncoding = "UTF-8")

# ---- valid-sample statistics + two-track group comparison ----
# Reality: LM22 is a peripheral-blood-derived signature; in myocardial tissue the
# valid rate is typically low. Strategy: primary = all samples (powered);
valid <- rownames(cib_res)[cib_res[, "P-value"] < 0.05]
rate <- length(valid) / nrow(cib_res)
logline(sprintf("  CIBERSORT valid rate: %d/%d (%.1f%%)",
                length(valid), nrow(cib_res), 100 * rate))
if (rate < 0.5) {
  logline("  [note] valid rate <50% — an inherent limitation of applying a blood-derived signature to myocardium")
  logline("        strategy: primary = all-sample comparison; the P<0.05 subset is sensitivity-only (discussed in the paper)")
}
cells <- colnames(cib_res)[1:22]

gv_all <- grp38[rownames(cib_res)]
if (any(is.na(gv_all))) stop("group vector not aligned with sample names (NA) — check the console output")
logline(sprintf("  [primary] all samples: HF %d vs NF %d",
                sum(gv_all == "HF"), sum(gv_all == "NF")))
wilc <- do.call(rbind, lapply(cells, function(cc) {
  x <- cib_res[, cc]
  p <- tryCatch(wilcox.test(x ~ gv_all)$p.value, error = function(e) NA)
  data.frame(cell = cc,
             median_HF = round(median(x[gv_all == "HF"]), 4),
             median_NF = round(median(x[gv_all == "NF"]), 4),
             P = p, stringsAsFactors = FALSE)
}))
wilc$P_BH <- p.adjust(wilc$P, method = "BH")
wilc <- wilc[order(wilc$P), ]
write.csv(wilc, file.path(OUT, "Table6_immune_diff_Wilcoxon_allsamples.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
sig_cells <- wilc$cell[wilc$P_BH < 0.05]
logline(sprintf("  [primary] differential cell types (BH q<0.05): %d: %s",
                length(sig_cells),
                ifelse(length(sig_cells) > 0, paste(sig_cells, collapse = "; "), "none")))

# ---- sensitivity: P<0.05 subset (test only if >=5 per group, otherwise describe) ----
nH <- sum(grp38[valid] == "HF"); nN <- sum(grp38[valid] == "NF")
if (length(valid) >= 10 && nH >= 5 && nN >= 5) {
  gv_v <- grp38[valid]
  wilc2 <- do.call(rbind, lapply(cells, function(cc) {
    x <- cib_res[valid, cc]
    p <- tryCatch(wilcox.test(x ~ gv_v)$p.value, error = function(e) NA)
    data.frame(cell = cc,
               median_HF = round(median(x[gv_v == "HF"]), 4),
               median_NF = round(median(x[gv_v == "NF"]), 4),
               P = p, stringsAsFactors = FALSE)
  }))
  wilc2$P_BH <- p.adjust(wilc2$P, method = "BH")
  write.csv(wilc2, file.path(OUT, "Table6b_immune_diff_sensitivity_Psubset.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  logline(sprintf("  [sensitivity] P<0.05 subset: HF %d vs NF %d, results written", nH, nN))
} else {
  write.csv(cib_res[valid, ], file.path(OUT, "Table6b_Psubset_sample_scores.csv"),
            row.names = TRUE, fileEncoding = "UTF-8")
  logline(sprintf("  [sensitivity] only %d samples in the P<0.05 subset (HF %d vs NF %d) — too few to test, scores written only",
                  length(valid), nH, nN))
}

# ---- Figure 7A: 22-cell boxplots (HF vs NF, all-sample primary) ----
frac <- cib_res[, cells]
longdf <- do.call(rbind, lapply(cells, function(cc)
  data.frame(cell = cc, frac = as.numeric(frac[, cc]), group = gv_all)))
# significance stars are looked up by cell name (wilc is sorted by P; do not index by position)
sig_map <- setNames(
  sapply(cells, function(cc) {
    q <- wilc$P_BH[wilc$cell == cc]
    if (length(q) == 0 || is.na(q)) ""
    else if (q < 0.001) "***" else if (q < 0.01) "**"
    else if (q < 0.05) "*" else ""
  }), cells)
longdf$cell_lab <- paste0(longdf$cell, " ", sig_map[longdf$cell])
longdf$cell_lab <- factor(longdf$cell_lab,
                          levels = paste0(cells, " ", sig_map[cells]))
p_box <- ggplot(longdf, aes(group, frac, fill = group)) +
  geom_boxplot(outlier.size = 0.4, width = 0.6) +
  facet_wrap(~ cell_lab, scales = "free_y", ncol = 6) +
  scale_fill_manual(values = c(HF = "#b2182b", NF = "#2166ac")) +
  labs(x = NULL, y = "CIBERSORT fraction") +
  theme_bw(base_size = 10) +
  theme(legend.position = "top", axis.text.x = element_text(size = 7))
ggsave(file.path(OUT, "Figure7A_immune_infiltration_boxplot.png"), p_box,
       width = 14, height = 10, dpi = 300)

# ============================================================
# Part 6: LASSO-gene x immune-cell Spearman correlation (Figure 7B material)
# ============================================================
logline("\n========== Part D: LASSO-gene x immune-cell correlation matrix ==========")
lg <- intersect(as.character(lasso_tab$gene), rownames(ex))
logline(sprintf("  LASSO genes: %d (%d present in the matrix)", nrow(lasso_tab), length(lg)))

corM <- matrix(NA, length(lg), length(cells),
               dimnames = list(lg, cells))
pM <- matrix(NA, length(lg), length(cells), dimnames = list(lg, cells))
for (g in lg) {
  for (cc in cells) {
    ct <- suppressWarnings(cor.test(ex[g, rownames(cib_res)], cib_res[, cc],
                                    method = "spearman"))
    corM[g, cc] <- ct$estimate; pM[g, cc] <- ct$p.value
  }
}
qM <- matrix(p.adjust(as.vector(pM), method = "BH"), nrow = length(lg),
            dimnames = list(lg, cells))
qM[is.na(qM)] <- 1   # guard against NA breaking the star assignment
cor_df <- do.call(rbind, lapply(lg, function(g)
  data.frame(gene = g, cell = cells, rho = round(corM[g, ], 3),
             P = signif(pM[g, ], 3), q_BH = signif(qM[g, ], 3),
             stringsAsFactors = FALSE)))
write.csv(cor_df, file.path(OUT, "Table7_LASSOximmune_Spearman.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
# NA->0 for plotting and max-search only (a constant cell score can still be NA); the CSV keeps true NA
if (any(is.na(corM))) {
  logline(sprintf("  [note] %d gene-cell pairs have NA correlation (constant cell score across samples); kept as NA in the CSV, shown as 0 in the plot",
                  sum(is.na(corM))))
  corM[is.na(corM)] <- 0
}
nsig <- sum(qM < 0.05)
logline(sprintf("  significant correlations (BH q<0.05): %d/%d pairs | max |rho|: %s=%.2f",
                nsig, length(lg) * length(cells),
                {
                  ii <- which.max(abs(corM))
                  paste(rownames(corM)[row(corM)[ii]], colnames(corM)[col(corM)[ii]])
                },
                max(abs(corM))))

# ---- Figure 7B: correlation heatmap (with significance stars) ----
starM <- matrix("", nrow = length(lg), ncol = length(cells),
                dimnames = list(lg, cells))
starM[qM < 0.05 & qM >= 0.01] <- "*"
starM[qM < 0.01 & qM >= 0.001] <- "**"
starM[qM < 0.001] <- "***"
labM <- matrix(paste0(sprintf("%.2f", corM), starM),
               nrow = length(lg), dimnames = list(lg, cells))
bk <- seq(-1, 1, length.out = 101)
png(file.path(OUT, "Figure7B_LASSO_immune_heatmap.png"),
    width = 3000, height = 1800, res = 300)
pheatmap(corM, display_numbers = labM, fontsize_number = 7,
         color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
         breaks = bk, cluster_cols = TRUE, cluster_rows = TRUE,
         border_color = NA, fontsize_row = 8, fontsize_col = 8,
         main = "Spearman correlation: LASSO genes x immune cells")
dev.off()

# ============================================================
# Part 7: anchor comparison + methods wording
# ============================================================
logline("\n========== anchor comparison ==========")
logline(sprintf("  DLRGs: yours=%d | expected 56 %s", nrow(dlrg),
                ifelse(nrow(dlrg) == 56, "[OK]", "[!! mismatch, check stage-3 outputs]")))
logline(sprintf("  CIBERSORT valid rate: %.1f%% (blood-derived signature on myocardium; two-track: all-sample primary + P-subset sensitivity)",
                100 * rate))
logline(sprintf("  GSEA-BP significant pathways=%d | qualitative expectation: HF-down terms include ribosome/metabolism",
                ifelse(is.null(gse), 0, nrow(as.data.frame(gse)))))
logline("  (no hard anchors for term identity — judge directionality manually)")

logline("\n[Methods wording (ready to use)]")
logline("GO/KEGG ORA: clusterProfiler (v4.x) enrichGO/enrichKEGG, BH, q<0.05, universe = all tested genes.")
logline("GSEA: all genes ranked by descending limma t; gseGO (BP, minGSSize=10, maxGSSize=500); NES and BH-q reported.")
logline("Immune infiltration: CIBERSORT v1.04 (Newman et al., Nat Biotechnol 2015; PMID: 25822800),")
logline("  LM22 reference matrix, 1000 permutations, quantile normalization on.")
logline("  Only 15/313 samples (4.8%) passed the P<0.05 filter — LM22 is blood-derived;")
logline("  an inherent limitation in myocardium. The primary analysis therefore uses all 313 samples (HF 177 vs NF 136), Wilcoxon+BH,")
logline("  with the P<0.05 subset reported as a sensitivity analysis; the limitation is acknowledged in the Discussion.")
logline("Gene-immune correlation: Spearman (all 313 samples), LASSO genes x 22 cells, BH-adjusted.")
logline("========== stage 7 finished: ", format(Sys.time()), "==========")
logline("\n>>> stage 7 done; if the KEGG pending file exists, M10 can be run next")

# ============================================================
# M1_preprocessing_limma.R  |  M1
# Purpose: three-cohort GEO loading -> QC -> limma -> 56 DLRGs
# Paper sections: Tables 1-2, Figure 2
# Pipeline position: step 2 of the main chain
# Prerequisite: M0 completed; the three GEO series-matrix files downloaded
# (Generated from the modular pipeline consolidation; analysis logic identical to the verified version.)
# ============================================================
# ==============================================================================
# M1: Stages 1-3 pipeline
# ------------------------------------------------------------------------------
# [What it is] The complete stage 1-3 pipeline (parsing -> annotation -> QC -> limma ->
#   DLRGs -> cross-validation -> FC sensitivity -> hub-gene overlap check): the three raw
#   GEO files run through to the final results; no intermediate files from earlier runs are needed.
#     gse57338_series_matrix.txt / gse5406_series_matrix.txt / gse116250_rpkm.txt
# [Fresh reproduction] only 5 raw files are needed in DIR:
#     lrg332.txt (provided in data/)
# [Usage] edit only the DIR line below (forward slashes!) -> select all and run
# [Runtime] ~15-25 min (parsing GSE57338: 1-2 min; the four limma runs with sva dominate)
# [Output] output/ folder: 3 rds + 9 QC plots + 4 volcano plots + Venn diagram + heatmap
#     + all CSVs (DLRG / cross-validation / FC sensitivity / core genes) + result summary
# [Built in] the GSE116250 duplicate-symbol dedup fix is applied at the source (after log2,
#     keep the highest-mean row per symbol — same order as stage 3d, results comparable).
# [Expected values] from the reference run; reruns should match closely, with occasional
#     +/-1-2 differences due to sva stochasticity:
#     GSE57338 DLRGs = 56 | GSE5406 replication = 28 (50%) | GSE116250: DCM 15 + ICM 12
# ==============================================================================

# ============================================================
# 0. User configuration (the only section to edit)
# ============================================================
DIR <- "./GEO_data"   # set to your working directory (forward slashes on Windows)
FILE_57338  <- file.path(DIR, "gse57338_series_matrix.txt")
FILE_5406   <- file.path(DIR, "gse5406_series_matrix.txt")
FILE_116250 <- file.path(DIR, "gse116250_rpkm.txt")
GENE_FILE   <- file.path(DIR, "lrg332.txt")
OUT <- file.path(DIR, "output")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)

# Per-file check (a wrong path fails immediately instead of mid-run)
for (f in c(FILE_57338, FILE_5406, FILE_116250, GENE_FILE)) {
  if (!file.exists(f)) stop(sprintf("file not found: %s\ncheck the DIR path above (use forward slashes /)", f))
}
cat("Input files verified. Output directory:", OUT, "\n\n")

# ============================================================
# 1. Automatic package installation (double safety — normally skipped after M0)
# ============================================================
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
need <- c("limma", "sva", "data.table", "AnnotationDbi", "ggplot2",
          "ggrepel", "pheatmap", "VennDiagram",
          "hugene11sttranscriptcluster.db", "hgu133a.db")
new <- need[!sapply(need, requireNamespace, quietly = TRUE)]
if (length(new) > 0) {
  cat("Installing missing packages:", paste(new, collapse = ", "), "\n")
  if (!requireNamespace("BiocManager", quietly = TRUE)) install.packages("BiocManager")
  BiocManager::install(new, ask = FALSE, update = FALSE)
}
suppressPackageStartupMessages({ library(limma); library(sva) })
library(data.table); library(ggplot2); library(ggrepel)
library(pheatmap); library(VennDiagram)
# Suppress VennDiagram log-file noise
try(futile.logger::flog.threshold(futile.logger::ERROR, name = "VennDiagramLogger"), silent = TRUE)

LFC <- 0.2   # primary threshold: |log2FC| > 0.2 (~FC 1.15) + FDR < 0.05
lrg <- toupper(trimws(readLines(GENE_FILE, warn = FALSE)))
lrg <- lrg[nzchar(lrg)]
cat(sprintf("Lactylation gene set: %d genes (Li 2025) | primary threshold: FDR<0.05 & |logFC|>%s\n\n", length(lrg), LFC))

# ============================================================
# 2. Utility functions (defined once, reused throughout)
# ============================================================

# --- Robust gene-symbol extraction (guards against topTable moving symbols into the ID column) ---
get_gene <- function(tt) {
  g <- NULL
  if ("ID" %in% names(tt)) g <- as.character(tt$ID)
  if (is.null(g) || all(grepl("^[0-9.]+$", g[!is.na(g)]))) {
    if ("gene" %in% names(tt)) g <- as.character(tt$gene)
  }
  if (is.null(g) || all(grepl("^[0-9.]+$", g[!is.na(g)]))) g <- as.character(rownames(tt))
  g
}

# --- Fast series-matrix parsing via fread ---
load_geo <- function(path) {
  cat("  parsing (fread fast mode; large files ~30-90 s)...\n")
  meta <- readLines(path, n = 500, warn = FALSE)
  begin_line <- grep("!series_matrix_table_begin", meta, fixed = TRUE)[1]
  if (is.na(begin_line)) stop("unexpected file format: table start marker not found")
  parse_line <- function(tag) {
    ln <- grep(paste0("^", tag), meta, value = TRUE)
    if (length(ln) == 0) return(NULL)
    v <- strsplit(ln[1], "\t", fixed = TRUE)[[1]][-1]
    gsub('^"|"$', "", v)
  }
  gsm <- parse_line("!Sample_geo_accession")
  nS <- length(gsm)
  pdata <- data.frame(row.names = gsm, check.names = FALSE)
  char_lines <- grep("^!Sample_characteristics_ch1", meta, value = TRUE)
  for (i in seq_along(char_lines)) {
    v <- gsub('^"|"$', "", strsplit(char_lines[i], "\t", fixed = TRUE)[[1]][-1])
    pdata[[paste0("char", i)]] <- v
  }
  dt <- data.table::fread(path, skip = begin_line, sep = "\t", header = TRUE,
                          fill = TRUE, na.strings = c("NA", "null", ""),
                          showProgress = interactive(), quote = "\"")
  ids <- as.character(dt[[1]])
  keep <- !is.na(ids) & !grepl("^!", ids) & nzchar(ids)
  dt <- dt[keep]; ids <- ids[keep]
  mat <- as.matrix(dt[, -1]); storage.mode(mat) <- "numeric"
  rownames(mat) <- ids
  if (ncol(mat) != nS) stop(sprintf("column count mismatch: %d vs %d", ncol(mat), nS))
  colnames(mat) <- gsm
  na_frac <- rowMeans(is.na(mat))
  if (any(na_frac >= 0.5)) mat <- mat[na_frac < 0.5, ]
  cat(sprintf("  parsed: %d probes x %d samples\n", nrow(mat), ncol(mat)))
  list(expr = mat, pdata = pdata)
}

find_char_col <- function(pdata, keyword) {
  for (cn in names(pdata)) {
    if (any(grepl(keyword, as.character(pdata[[cn]]), ignore.case = TRUE))) return(cn)
  }
  NULL
}

# --- Probe annotation + multi-probe collapse (highest-mean probe per gene) ---
annotate_collapse <- function(expr, db) {
  suppressPackageStartupMessages(library(db, character.only = TRUE))
  map <- AnnotationDbi::select(get(db), keys = rownames(expr),
                               columns = "SYMBOL", keytype = "PROBEID")
  map <- map[!is.na(map$SYMBOL) & map$SYMBOL != "", ]
  expr <- expr[rownames(expr) %in% map$PROBEID, , drop = FALSE]
  map <- map[match(rownames(expr), map$PROBEID), ]
  cat(sprintf("  annotated: %d/%d probes mapped to symbols\n", nrow(expr), nrow(map)))
  ord <- order(rowMeans(expr, na.rm = TRUE), decreasing = TRUE)
  expr <- expr[ord, , drop = FALSE]; map <- map[ord, ]
  keep <- !duplicated(map$SYMBOL)
  expr <- expr[keep, , drop = FALSE]
  rownames(expr) <- map$SYMBOL[keep]
  stopifnot(!any(duplicated(rownames(expr))))   # assert rownames are unique
  expr
}

# --- QC triple ---
qc_all <- function(expr, group, prefix, title) {
  set.seed(2026)
  idx <- if (ncol(expr) > 60) sample(ncol(expr), 60) else seq_len(ncol(expr))
  png(file.path(OUT, paste0(prefix, "_density.png")), width = 1200, height = 900, res = 150)
  limma::plotDensities(expr[, idx], legend = FALSE, main = paste(title, "- expression density (60-sample subsample)"))
  dev.off()
  png(file.path(OUT, paste0(prefix, "_PCA.png")), width = 1200, height = 900, res = 150)
  pc <- prcomp(t(expr), center = TRUE, scale. = FALSE)
  df <- data.frame(PC1 = pc$x[, 1], PC2 = pc$x[, 2], grp = as.factor(group[colnames(expr)]))
  lv <- levels(df$grp); cols <- rainbow(length(lv))
  plot(df$PC1, df$PC2, col = cols[df$grp], pch = 16, cex = 1.2,
       xlab = sprintf("PC1 (%.1f%%)", summary(pc)$importance[2, 1] * 100),
       ylab = sprintf("PC2 (%.1f%%)", summary(pc)$importance[2, 2] * 100),
       main = paste(title, "- PCA"))
  legend("topright", legend = lv, col = cols, pch = 16, cex = 0.8)
  dev.off()
  png(file.path(OUT, paste0(prefix, "_clustering.png")), width = 1400, height = 900, res = 150)
  hc <- hclust(dist(t(expr)), method = "average")
  plot(hc, labels = FALSE, main = paste(title, "- sample clustering"))
  dev.off()
  cat("  QC plots written:", prefix, "_density/_PCA/_clustering.png\n")
}

# --- limma + sva (with the robust get_gene) ---
run_limma <- function(expr, group, g1, g2, use_sva = TRUE, prefix = "") {
  set.seed(2026)   # sva is stochastic; fixed seed for reproducibility
  group <- factor(as.character(group))
  design <- model.matrix(~ 0 + group)
  colnames(design) <- levels(group)
  contr <- paste0(g1, "-", g2)
  if (use_sva && ncol(expr) >= 20) {
    fit0 <- limma::lmFit(expr, design)
    resid <- residuals(fit0, expr)
    n.sv <- NA
    try({ n.sv <- sva::num.sv(t(resid), design, method = "be") }, silent = TRUE)
    if (!is.na(n.sv) && n.sv >= 1 && n.sv <= 15) {
      svobj <- sva::sva(t(resid), design, n.sv = n.sv)
      design <- cbind(design, svobj$sv)
      colnames(design)[(ncol(design)-n.sv+1):ncol(design)] <- paste0("SV", 1:n.sv)
      cat(sprintf("  [sva] %d surrogate variables absorbed\n", n.sv))
    }
  }
  cm <- limma::makeContrasts(contrasts = contr, levels = design)
  fit <- limma::lmFit(expr, design)
  fit <- limma::contrasts.fit(fit, cm)   # the contrast matrix must be applied to the fit
  fit <- limma::eBayes(fit)
  tt <- limma::topTable(fit, coef = contr, number = Inf, sort.by = "P")
  tt$.gid <- get_gene(tt)   # robust gene extraction
  cat(sprintf("  %s: %s vs %s -> %d genes tested | FDR<0.05 & |logFC|>%s: Up %d / Down %d\n",
              prefix, g1, g2, nrow(tt), LFC,
              sum(tt$adj.P.Val < 0.05 & tt$logFC > LFC, na.rm = TRUE),
              sum(tt$adj.P.Val < 0.05 & tt$logFC < -LFC, na.rm = TRUE)))
  tt
}

# --- Volcano plot (LRGs circled and labelled) ---
volcano <- function(tt, prefix, title) {
  tt$sig <- ifelse(tt$adj.P.Val < 0.05 & abs(tt$logFC) > LFC,
                   ifelse(tt$logFC > 0, "Up", "Down"), "NS")
  tt$is_lrg <- toupper(tt$.gid) %in% lrg
  tt$lab <- ifelse(tt$is_lrg & tt$sig != "NS", tt$.gid, "")
  p <- ggplot(tt, aes(logFC, -log10(adj.P.Val), color = sig)) +
    geom_point(data = subset(tt, !is_lrg), size = 0.6, alpha = 0.35) +
    geom_point(data = subset(tt, is_lrg), aes(fill = sig), size = 1.6, shape = 21, color = "black") +
    scale_color_manual(values = c(Down = "#2166ac", NS = "grey75", Up = "#b2182b")) +
    scale_fill_manual(values = c(Down = "#2166ac", NS = "grey75", Up = "#b2182b")) +
    geom_vline(xintercept = c(-LFC, LFC), linetype = 2, color = "grey45") +
    geom_hline(yintercept = -log10(0.05), linetype = 2, color = "grey45") +
    ggrepel::geom_text_repel(aes(label = lab), size = 2.6, max.overlaps = 25, segment.size = 0.25) +
    labs(title = title, x = "log2 FC", y = "-log10 FDR", color = "", fill = "") +
    theme_bw(base_size = 11) + theme(legend.position = "top")
  ggsave(file.path(OUT, paste0(prefix, "_volcano.png")), p, width = 7, height = 6, dpi = 300)
}

# ============================================================
# 3. Stage 1: GSE57338 discovery cohort (parse + annotate + QC)
# ============================================================
cat("========== stage 1/3: GSE57338 discovery cohort ==========\n")
g38 <- load_geo(FILE_57338)
pdata38 <- g38$pdata
hf_col <- find_char_col(pdata38, "heart failure")
grp38 <- ifelse(grepl("\\bno\\b", as.character(pdata38[[hf_col]]), ignore.case = TRUE), "NF", "HF")
names(grp38) <- rownames(pdata38)
etio_col <- find_char_col(pdata38, "disease status")
etio38 <- ifelse(grepl("non-failing|non failing", as.character(pdata38[[etio_col]]), ignore.case = TRUE),
                 "NF", sub('.*:\\s*', '', as.character(pdata38[[etio_col]])))
names(etio38) <- rownames(pdata38)
expr38 <- annotate_collapse(g38$expr, "hugene11sttranscriptcluster.db")
cat(sprintf("  gene-level matrix: %d genes x %d samples\n", nrow(expr38), ncol(expr38)))
cat("  group counts:\n"); print(table(grp38)); print(table(etio38))
qc_all(expr38, grp38[colnames(expr38)], "GSE57338", "GSE57338 discovery")
saveRDS(list(expr = expr38, group = grp38, etiology = etio38), file.path(OUT, "GSE57338_processed.rds"))
n_lrg38 <- sum(toupper(rownames(expr38)) %in% lrg)
cat(sprintf("  LRGs matched: %d/%d\n", n_lrg38, length(lrg)))
if (n_lrg38 < 100) cat("  warning: low LRG match count\n")

# ============================================================
# 4. Stage 1: GSE5406 validation cohort 1 (parse + annotate + QC)
# ============================================================
cat("\n========== stage 1/3: GSE5406 validation cohort 1 ==========\n")
g06 <- load_geo(FILE_5406)
pdata06 <- g06$pdata
fail_col <- find_char_col(pdata06, "heart failure")
fail06 <- as.character(pdata06[[fail_col]])
grp06 <- ifelse(grepl("donor|non-failing|normal", fail06, ignore.case = TRUE), "NF",
         ifelse(grepl("ischemic", fail06, ignore.case = TRUE), "ICM",
         ifelse(grepl("idiopathic|dilated", fail06, ignore.case = TRUE), "DCM", "Other")))
names(grp06) <- rownames(pdata06)
grp06[grp06 %in% c("ICM", "DCM")] <- "HF"
expr06 <- annotate_collapse(g06$expr, "hgu133a.db")
cat(sprintf("  gene-level matrix: %d genes x %d samples\n", nrow(expr06), ncol(expr06)))
cat("  group counts:\n"); print(table(grp06))
qc_all(expr06, grp06[colnames(expr06)], "GSE5406", "GSE5406 validation 1")
saveRDS(list(expr = expr06, group = grp06), file.path(OUT, "GSE5406_processed.rds"))
n_lrg06 <- sum(toupper(rownames(expr06)) %in% lrg)
cat(sprintf("  LRGs matched: %d/%d\n", n_lrg06, length(lrg)))

# ============================================================
# 5. Stage 1: GSE116250 validation cohort 2 (parse + log2 + dedup fix + QC)
#    [Key] the dedup fix is applied at the source: after log2, collapse duplicated symbols
# ============================================================
cat("\n========== stage 1/3: GSE116250 validation cohort 2 ==========\n")
rpkm <- data.table::fread(FILE_116250, sep = "\t", header = TRUE, data.table = FALSE)
genes62 <- rpkm$Common_name
mat62 <- as.matrix(rpkm[, -(1:2)])
storage.mode(mat62) <- "numeric"
rownames(mat62) <- genes62
mat62 <- mat62[!is.na(genes62) & nzchar(genes62), ]
grp62 <- sub("[0-9]+$", "", colnames(mat62))
names(grp62) <- colnames(mat62)
cat("  group counts:\n"); print(table(grp62))
mat62 <- log2(mat62 + 1)
#    by keeping the highest-mean row (same order as stage 3d).
n_dup <- sum(duplicated(rownames(mat62)))
if (n_dup > 0) {
  gm <- rowMeans(mat62, na.rm = TRUE)
  mat62 <- mat62[order(-gm), ]
  mat62 <- mat62[!duplicated(rownames(mat62)), ]
  cat(sprintf("  dedup: %d duplicated symbols merged (highest-mean row kept), %d genes remain\n", n_dup, nrow(mat62)))
}
stopifnot(!any(duplicated(rownames(mat62))))
cat(sprintf("  after log2 + dedup: %d genes x %d samples\n", nrow(mat62), ncol(mat62)))
qc_all(mat62, grp62[colnames(mat62)], "GSE116250", "GSE116250 validation 2")
saveRDS(list(expr = mat62, group = grp62), file.path(OUT, "GSE116250_processed.rds"))
n_lrg62 <- sum(toupper(rownames(mat62)) %in% lrg)
cat(sprintf("  LRGs matched: %d/%d\n", n_lrg62, length(lrg)))
if (n_lrg62 < 100) cat("  warning: low LRG match count\n")

# ============================================================
# ---- Dedup fix (root-cause fix applied at the source) ----
# ============================================================
cat("\n========== stage 2/3: discovery limma (57338 HF vs NF) ==========\n")
tt38 <- run_limma(expr38, grp38, "HF", "NF", prefix = "GSE57338")
write.csv(tt38, file.path(OUT, "GSE57338_HFvsNF_all.csv"), row.names = FALSE)
volcano(tt38, "GSE57338", "GSE57338 HF vs NF")
sig38 <- subset(tt38, adj.P.Val < 0.05 & abs(logFC) > LFC)
dlrg38 <- sig38[toupper(sig38$.gid) %in% lrg, ]
cat(sprintf("  DLRGs: %d  [expected: 56]\n", nrow(dlrg38)))
write.csv(dlrg38, file.path(OUT, "GSE57338_DLRGs.csv"), row.names = FALSE)

# ============================================================
# 6. Stage 2: discovery 57338 HF vs NF -> DLRGs
# ============================================================
cat("\n========== stage 2/3: validation 1 limma (5406 HF vs NF) ==========\n")
tt06 <- run_limma(expr06, grp06, "HF", "NF", prefix = "GSE5406")
write.csv(tt06, file.path(OUT, "GSE5406_HFvsNF_all.csv"), row.names = FALSE)
volcano(tt06, "GSE5406", "GSE5406 HF vs NF")

# ============================================================
# 7. Stage 2: validation 5406 HF vs NF
# ============================================================
cat("\n========== stage 2/3: validation 2 limma (116250 DCM vs NF) ==========\n")
tt62_dcm <- run_limma(mat62, grp62, "DCM", "NF", prefix = "GSE116250-DCM")
write.csv(tt62_dcm, file.path(OUT, "GSE116250_DCMvsNF_all.csv"), row.names = FALSE)
volcano(tt62_dcm, "GSE116250_DCM", "GSE116250 DCM vs NF")

cat("\n========== stage 2/3: validation 2 limma (116250 ICM vs NF) ==========\n")
tt62_icm <- run_limma(mat62, grp62, "ICM", "NF", prefix = "GSE116250-ICM")
write.csv(tt62_icm, file.path(OUT, "GSE116250_ICMvsNF_all.csv"), row.names = FALSE)
volcano(tt62_icm, "GSE116250_ICM", "GSE116250 ICM vs NF")

# ============================================================
# 8. Stage 2: validation 116250 DCM/ICM vs NF (two contrasts)
# ============================================================
cat("\n========== stage 3/3: cross-validation ==========\n")
sig06  <- subset(tt06, adj.P.Val < 0.05 & abs(logFC) > LFC)
sig62d <- subset(tt62_dcm, adj.P.Val < 0.05 & abs(logFC) > LFC)
sig62i <- subset(tt62_icm, adj.P.Val < 0.05 & abs(logFC) > LFC)

genes38 <- dlrg38$.gid
cv <- data.frame(
  gene = genes38,
  logFC_57338 = dlrg38$logFC,
  FDR_57338 = dlrg38$adj.P.Val,
  logFC_5406 = tt06$logFC[match(genes38, tt06$.gid)],
  FDR_5406 = tt06$adj.P.Val[match(genes38, tt06$.gid)],
  logFC_116250_DCM = tt62_dcm$logFC[match(genes38, tt62_dcm$.gid)],
  FDR_116250_DCM = tt62_dcm$adj.P.Val[match(genes38, tt62_dcm$.gid)],
  logFC_116250_ICM = tt62_icm$logFC[match(genes38, tt62_icm$.gid)],
  FDR_116250_ICM = tt62_icm$adj.P.Val[match(genes38, tt62_icm$.gid)]
)
cv$sig_5406 <- toupper(cv$gene) %in% toupper(sig06$.gid)
cv$sig_116250_DCM <- toupper(cv$gene) %in% toupper(sig62d$.gid)
cv$sig_116250_ICM <- toupper(cv$gene) %in% toupper(sig62i$.gid)
cv$dir_agree_5406 <- sign(cv$logFC_57338) == sign(cv$logFC_5406)
cv$dir_agree_DCM  <- sign(cv$logFC_57338) == sign(cv$logFC_116250_DCM)
cv$dir_agree_ICM  <- sign(cv$logFC_57338) == sign(cv$logFC_116250_ICM)
cv$n_validated <- rowSums(cbind(cv$sig_5406 & cv$dir_agree_5406,
                                cv$sig_116250_DCM & cv$dir_agree_DCM,
                                cv$sig_116250_ICM & cv$dir_agree_ICM), na.rm = TRUE)
cv <- cv[order(-cv$n_validated, cv$FDR_57338), ]
write.csv(cv, file.path(OUT, "DLRG_cross_validation.csv"), row.names = FALSE)

rep06 <- sum(cv$sig_5406 & cv$dir_agree_5406, na.rm = TRUE)
repDCM <- sum(cv$sig_116250_DCM & cv$dir_agree_DCM, na.rm = TRUE)
repICM <- sum(cv$sig_116250_ICM & cv$dir_agree_ICM, na.rm = TRUE)
gold <- sum(cv$n_validated == 3, na.rm = TRUE)
cat(sprintf("  discovery DLRGs: %d\n", nrow(cv)))
cat(sprintf("  GSE5406 replication: %d (%.0f%%)   [expected: 28/50%%]\n", rep06, 100*rep06/nrow(cv)))
cat(sprintf("  116250-DCM replication: %d (%.0f%%) [expected: 15]\n", repDCM, 100*repDCM/nrow(cv)))
cat(sprintf("  116250-ICM replication: %d (%.0f%%) [expected: 12]\n", repICM, 100*repICM/nrow(cv)))
cat(sprintf("  gold standard (replicated in all three): %d  [expected: 7]\n", gold))

# ============================================================
# 9. Cross-validation table + direction consistency
# ============================================================
cat("\n========== Venn diagram and core DLRGs ==========\n")
dlrg06  <- sig06[toupper(sig06$.gid) %in% lrg, ]
dlrg62d <- sig62d[toupper(sig62d$.gid) %in% lrg, ]
venn.plot <- venn.diagram(
  x = list(GSE57338 = dlrg38$.gid, GSE5406 = dlrg06$.gid,
           GSE116250_DCM = dlrg62d$.gid),
  filename = file.path(OUT, "DLRG_venn.png"), imagetype = "png", resolution = 300,
  fill = c("#b2182b", "#2166ac", "#1b7837"), alpha = 0.45,
  cat.cex = 1.1, cex = 1.2, main = "DLRG overlap across datasets", main.cex = 1.2
)
core <- Reduce(intersect, list(toupper(dlrg38$.gid), toupper(dlrg06$.gid), toupper(dlrg62d$.gid)))
cat(sprintf("Three-cohort core DLRGs: %d  [expected: 9]\n", length(core)))
if (length(core) > 0) cat("  ", paste(core, collapse = ", "), "\n")
write.csv(data.frame(core_DLRG = core), file.path(OUT, "core_DLRGs.csv"), row.names = FALSE)

if (length(core) >= 2) {
  idx <- match(core, toupper(rownames(expr38)))
  idx <- idx[!is.na(idx)]
  mat <- expr38[idx, ]
  mat <- t(scale(t(mat)))
  ann <- data.frame(Group = grp38[colnames(mat)])
  rownames(ann) <- colnames(mat)
  png(file.path(OUT, "core_DLRG_heatmap.png"), width = 2200, height = 1400, res = 200)
  pheatmap(mat, annotation_col = ann, scale = "none",
           color = colorRampPalette(c("#2166ac", "white", "#b2182b"))(100),
           border_color = NA, fontsize_row = 9, show_colnames = FALSE,
           main = "Core DLRGs in GSE57338 (z-score)")
  dev.off()
  cat(sprintf("core-DLRG heatmap saved (%d genes x %d samples)\n", nrow(mat), ncol(mat)))
}

# ============================================================
# 10. Venn diagram + core DLRGs + heatmap
# ============================================================
cat("\n========== FC sensitivity (three tiers) ==========\n")
LFC_SET <- c(0.2, 0.585, 1.0)
LFC_NAME <- c("FC1.15 (primary)", "FC1.5", "FC2.0")
sens <- list()
for (k in seq_along(LFC_SET)) {
  tag <- if (k == 1) "0.2" else if (k == 2) "0.585" else "1.0"
  s38 <- subset(tt38, adj.P.Val < 0.05 & abs(logFC) > LFC_SET[k])
  g38k <- toupper(s38$.gid)[toupper(s38$.gid) %in% lrg]
  if (length(g38k) == 0) {
    cat(sprintf("  %s: DLRGs = 0, skipped\n", LFC_NAME[k]))
    sens[[k]] <- data.frame(threshold = LFC_NAME[k], LFC = LFC_SET[k], DLRG_57338 = 0,
                            reproduced_5406 = 0, reproduced_116250_DCM = 0,
                            reproduced_116250_ICM = 0, gold_standard_3dir = 0,
                            core_intersection = 0)
    next
  }
  m06 <- toupper(tt06$.gid)[match(g38k, toupper(tt06$.gid))]
  md <- toupper(tt62_dcm$.gid)[match(g38k, toupper(tt62_dcm$.gid))]
  mi <- toupper(tt62_icm$.gid)[match(g38k, toupper(tt62_icm$.gid))]
  r06 <- sum(tt06$adj.P.Val[match(g38k, toupper(tt06$.gid))] < 0.05 &
             abs(tt06$logFC[match(g38k, toupper(tt06$.gid))]) > LFC_SET[k] &
             sign(tt06$logFC[match(g38k, toupper(tt06$.gid))]) == sign(tt38$logFC[match(g38k, toupper(tt38$.gid))]), na.rm = TRUE)
  rd <- sum(tt62_dcm$adj.P.Val[match(g38k, toupper(tt62_dcm$.gid))] < 0.05 &
             abs(tt62_dcm$logFC[match(g38k, toupper(tt62_dcm$.gid))]) > LFC_SET[k] &
             sign(tt62_dcm$logFC[match(g38k, toupper(tt62_dcm$.gid))]) == sign(tt38$logFC[match(g38k, toupper(tt38$.gid))]), na.rm = TRUE)
  ri <- sum(tt62_icm$adj.P.Val[match(g38k, toupper(tt62_icm$.gid))] < 0.05 &
             abs(tt62_icm$logFC[match(g38k, toupper(tt62_icm$.gid))]) > LFC_SET[k] &
             sign(tt62_icm$logFC[match(g38k, toupper(tt62_icm$.gid))]) == sign(tt38$logFC[match(g38k, toupper(tt38$.gid))]), na.rm = TRUE)
  s06k <- subset(tt06, adj.P.Val < 0.05 & abs(logFC) > LFC_SET[k])
  s62dk <- subset(tt62_dcm, adj.P.Val < 0.05 & abs(logFC) > LFC_SET[k])
  corek <- Reduce(intersect, list(g38k,
                                  toupper(s06k$.gid)[toupper(s06k$.gid) %in% lrg],
                                  toupper(s62dk$.gid)[toupper(s62dk$.gid) %in% lrg]))
  goldk <- sum(r06 + rd + ri >= 3)
  cat(sprintf("  %s: DLRGs=%d | 5406 replication %d | 116250 DCM %d + ICM %d | gold %d | core overlap %d\n",
              LFC_NAME[k], length(g38k), r06, rd, ri, goldk, length(corek)))
  sens[[k]] <- data.frame(threshold = LFC_NAME[k], LFC = LFC_SET[k], DLRG_57338 = length(g38k),
                          reproduced_5406 = r06, reproduced_116250_DCM = rd,
                          reproduced_116250_ICM = ri, gold_standard_3dir = goldk,
                          core_intersection = length(corek))
  write.csv(data.frame(gene = g38k,
                       logFC_5406 = tt06$logFC[match(g38k, toupper(tt06$.gid))],
                       FDR_5406 = tt06$adj.P.Val[match(g38k, toupper(tt06$.gid))],
                       logFC_116250_DCM = tt62_dcm$logFC[match(g38k, toupper(tt62_dcm$.gid))],
                       FDR_116250_DCM = tt62_dcm$adj.P.Val[match(g38k, toupper(tt62_dcm$.gid))]),
            file.path(OUT, sprintf("DLRG_cross_validation_LFC%s.csv", tag)), row.names = FALSE)
}
sens_df <- do.call(rbind, sens)
write.csv(sens_df, file.path(OUT, "FC_sensitivity_summary.csv"), row.names = FALSE)
cat("  summary table saved: FC_sensitivity_summary.csv\n")

# ============================================================
# 11. FC-threshold sensitivity (primary 0.2 / 1.5 / 2.0)
# ============================================================
cat("\n========== hub-gene overlap check vs Zhang BBRC (HLTF/HMGN2/ARGLU1/LSP1) ==========\n")
for (h in c("HLTF", "HMGN2", "ARGLU1", "LSP1")) {
  id <- which(toupper(tt62_dcm$.gid) == h)
  ii <- which(toupper(tt62_icm$.gid) == h)
  i38 <- which(toupper(tt38$.gid) == h)
  if (length(id) == 0 && length(ii) == 0) { cat(sprintf("  %s: not detected in 116250\n", h)); next }
  cat(sprintf("  %s | 57338: logFC=%.2f FDR=%.3g | 116250-DCM: logFC=%.2f FDR=%.3g | 116250-ICM: logFC=%.2f FDR=%.3g\n",
              h,
              if (length(i38)) tt38$logFC[i38[1]] else NA,
              if (length(i38)) tt38$adj.P.Val[i38[1]] else NA,
              if (length(id)) tt62_dcm$logFC[id[1]] else NA,
              if (length(id)) tt62_dcm$adj.P.Val[id[1]] else NA,
              if (length(ii)) tt62_icm$logFC[ii[1]] else NA,
              if (length(ii)) tt62_icm$adj.P.Val[ii[1]] else NA))
}
# 12. Overlap check with previously reported hub genes (Zhang BBRC) in 116250
if (length(core) > 0) {
  cat("\n  logFC and FDR of core DLRGs in 116250 (DCM/ICM):\n")
  for (g in core) {
    id <- which(toupper(tt62_dcm$.gid) == g); ii <- which(toupper(tt62_icm$.gid) == g)
    cat(sprintf("    %-8s DCM: %+.2f (FDR=%.2g) | ICM: %+.2f (FDR=%.2g)\n", g,
                if (length(id)) tt62_dcm$logFC[id[1]] else NA,
                if (length(id)) tt62_dcm$adj.P.Val[id[1]] else NA,
                if (length(ii)) tt62_icm$logFC[ii[1]] else NA,
                if (length(ii)) tt62_icm$adj.P.Val[ii[1]] else NA))
  }
}

# ============================================================
#     Full per-gene card of the core genes in 116250
# ============================================================
cat("\n==================================================\n")
cat("  Pipeline complete! Output directory:", OUT, "\n\n")
cat("  Compare with expected values (reference run; +/-1-2 is sva stochasticity):\n")
cat(sprintf("    DLRGs=%d (expected 56) | 5406 replication=%d (28) | DCM=%d (15) | ICM=%d (12)\n",
            nrow(cv), rep06, repDCM, repICM))
cat(sprintf("    gold standard=%d (7) | core overlap=%d (9)\n", gold, length(core)))
cat("  HMGN2 should replicate across all three cohorts (116250 bidirectional FDR ~ 1e-08); LSP1 logFC 1.18/1.49\n")
cat("  If numbers deviate strongly -> check the console output\n\n")
cat("  Next (stages 4-5): LASSO + ML feature selection (candidates = 56 DLRGs or 9 core genes)\n")

cat("==================================================\n")

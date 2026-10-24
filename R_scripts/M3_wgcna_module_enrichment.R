# ============================================================
# M3_wgcna_module_enrichment.R  |  M3
# Purpose: co-expression modules + green/yellow module enrichment
# Paper section: Figure 4
# Pipeline position: step 4 of the main chain
# Prerequisite: M1 completed
# ============================================================
# ==============================================================================
# M3: stage-5 pipeline (WGCNA co-expression network + module-HF association + DLRG enrichment + GO/KEGG functional enrichment)
# ------------------------------------------------------------------------------
# [Scope] complete stage 5 (WGCNA + module enrichment in one script):
#   WGCNA network -> module-trait association -> DLRG Fisher enrichment -> LASSO-gene kME
#   -> green/yellow module GO-BP enrichment (universe = all network genes) -> core-gene / LASSO-13 annotation -> KEGG enrichment
# [Input] M1 outputs (GSE57338_processed.rds + GSE57338_HFvsNF_all.csv) +
#         lrg332.txt + M2 output (LASSO_coefficients.csv)
# [Usage] set DIR, then select all and run
# [Runtime] first run ~8-12 min (WGCNA install 3-5 min + network build 5-8 min); afterwards 5-8 min
# [Output] output/ (nine WGCNA artifacts + seven enrichment artifacts):
#   four WGCNA plots + module_assignment_all_genes.csv + module_trait_correlation.csv
#   + DLRG_module_enrichment.csv + LASSO_genes_kME.csv + key_module_genes.csv
#   + enrich_green/yellow_module_GO.csv + enrich_core9/lasso13_GO.csv
#   + GO_green/yellow_module_top10.png + enrich_*_KEGG.csv (if the KEGG REST call fails it is skipped —
#   reference results are in results/ (paper Table 4))
# [Runtime] first run ~15-20 min; with WGCNA preinstalled ~10 min
# [Expected values] (reference run, numerically audited): soft-threshold power = 13 (R2 = 0.85) | 18 modules (grey excluded; 19 with grey) | green cor = -0.809, p = 1.35e-73
#   | DLRG enrichment in green OR = 7.47, p = 6.5e-07 (12 genes) | blue contains 22 DLRGs, p = 1.4e-05
#   | LASSO-13: 16 significant GO-BP terms, top = chromatin organization, FE = 15.2, padj = 3.2e-05
#   | yellow GO = ECM organization, padj = 0.023 | green module GO/KEGG not significant (genuine result)
# [Interpretation] 1-2 modules should correlate significantly with HF (|cor|>0.5, p<0.05); DLRGs should
#   enrich in the HF-related module (p<0.05); LASSO genes should cluster in the key module with kME>0.6
# ==============================================================================
rm(list = ls())
cat("== M3 stage-5 WGCNA started ==\n")

# ---- WGCNA auto-install (incl. the Bioconductor dependency impute) ----
if (!requireNamespace("BiocManager", quietly = TRUE))
  install.packages("BiocManager", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
if (!requireNamespace("impute", quietly = TRUE))
  BiocManager::install("impute", ask = FALSE, update = FALSE)
if (!requireNamespace("WGCNA", quietly = TRUE))
  install.packages("WGCNA", repos = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/")
suppressPackageStartupMessages(library(WGCNA))
cat("WGCNA ready, version:", as.character(packageVersion("WGCNA")), "\n")

DIR <- "./GEO_data"   # set to your working directory
OUT <- file.path(DIR, "output")

need <- c(file.path(OUT, "GSE57338_processed.rds"),
          file.path(OUT, "GSE57338_HFvsNF_all.csv"),
          file.path(OUT, "LASSO_coefficients.csv"),
          file.path(DIR, "lrg332.txt"))
missf <- need[!file.exists(need)]
if (length(missf) > 0) stop("file(s) not found: ", paste(basename(missf), collapse = ", "),
                            "\nrun M1 and M2 first")

# ------------------------------------------------------------------
# 1. Data preparation (group/etiology aligned by colnames as a safety guard)
# ------------------------------------------------------------------
d38  <- readRDS(need[1])
tt38 <- read.csv(need[2], check.names = FALSE)
co   <- read.csv(need[3])
sel  <- as.character(co$gene)
lrg  <- unique(toupper(trimws(readLines(need[4], warn = FALSE))))
lrg  <- lrg[nzchar(lrg)]

align_vec <- function(v, cn, tag) {
  if (is.null(v)) return(NULL)
  if (length(v) != length(cn)) { cat(sprintf("[%s] length %d != column count %d, trait ignored\n", tag, length(v), length(cn))); return(NULL) }
  if (!is.null(names(v))) {
    v2 <- v[cn]
    if (any(is.na(v2))) { cat(sprintf("[%s] name mismatch, ignored\n", tag)); return(NULL) }
    return(v2)
  }
  v
}
grp <- align_vec(d38$group, colnames(d38$expr), "group")
et  <- align_vec(d38$etiology, colnames(d38$expr), "etiology")
cat("group distribution:\n"); print(table(grp))
if (!is.null(et)) { cat("etiology distribution:\n"); print(table(et)) }

LFC <- 0.2
gidcol <- if (".gid" %in% names(tt38)) ".gid" else names(tt38)[ncol(tt38)]
sig38 <- subset(tt38, adj.P.Val < 0.05 & abs(logFC) > LFC)
genes56 <- as.character(sig38[[gidcol]])
genes56 <- genes56[toupper(genes56) %in% lrg]
cat(sprintf("DLRGs recomputed=%d | LASSO genes=%d\n", length(genes56), length(sel)))

# ------------------------------------------------------------------
# 2. Gene filtering: MAD top 5000 + forced inclusion of DLRGs and LASSO genes
# ------------------------------------------------------------------
expr <- d38$expr
expr[expr < 1e-3 & expr > -1e-3] <- 0
madv <- apply(expr, 1, mad)
top5k <- names(sort(madv, decreasing = TRUE))[1:min(5000, nrow(expr))]
useGenes <- unique(c(top5k, genes56, sel))
useGenes <- useGenes[useGenes %in% rownames(expr)]
datExpr <- t(expr[useGenes, ])          # WGCNA expects samples x genes
cat(sprintf("genes entering the network: %d (MAD top5000 + %d DLRGs + %d LASSO)\n",
            length(useGenes), sum(genes56 %in% useGenes), sum(sel %in% useGenes)))
gsg <- goodSamplesGenes(datExpr, verbose = 0)
if (!gsg$allOK) {
  datExpr <- datExpr[gsg$goodSamples, gsg$goodGenes]
  cat("problematic genes/samples removed\n")
}
cat(sprintf("network size: %d samples x %d genes\n", nrow(datExpr), ncol(datExpr)))

# ---- Sample outlier check (report only, no deletion) ----
sampleTree <- hclust(dist(datExpr), method = "average")
png(file.path(OUT, "WGCNA_sample_dendrogram.png"), 1200, 700, res = 120)
par(mar = c(1, 5, 3, 1))
plot(sampleTree, main = "Sample clustering (outlier check)", sub = "", xlab = "",
     cex = 0.35, cex.lab = 1.2)
abline(h = quantile(sampleTree$height, 0.995), col = "red", lty = 2)
dev.off()

# ------------------------------------------------------------------
# 3. Soft-threshold selection (signed network, R2 >= 0.85 automated)
# ------------------------------------------------------------------
cat("running pickSoftThreshold (~1-2 min)...\n")
sft <- pickSoftThreshold(datExpr, powerVector = 1:20, networkType = "signed",
                         verbose = 0)
pt <- sft$fitIndices
fitOk <- pt$SFT.R.sq >= 0.85
softPower <- if (any(fitOk)) pt$Power[which(fitOk)[1]] else
             { fitOk2 <- pt$SFT.R.sq >= 0.8
               if (any(fitOk2)) pt$Power[which(fitOk2)[1]] else 6 }
cat(sprintf("selected soft threshold power=%d (R2=%.2f)\n", softPower,
            pt$SFT.R.sq[pt$Power == softPower]))
png(file.path(OUT, "WGCNA_softthreshold.png"), 1300, 600, res = 120)
par(mfrow = c(1, 2), mar = c(4.5, 4.5, 3, 1))
plot(pt$Power, -sign(pt$slope) * pt$SFT.R.sq, type = "n",
     xlab = "Soft threshold power", ylab = "Scale-free topology R^2",
     main = "Scale-free topology fit")
text(pt$Power, -sign(pt$slope) * pt$SFT.R.sq, labels = pt$Power, cex = 0.8)
abline(v = softPower, col = "red", lty = 2)
plot(pt$Power, pt$mean.k., type = "n", xlab = "Soft threshold power",
     ylab = "Mean connectivity", main = "Mean connectivity")
text(pt$Power, pt$mean.k., labels = pt$Power, cex = 0.8)
abline(v = softPower, col = "red", lty = 2)
dev.off()

# ------------------------------------------------------------------
# 4. Network construction (single block, deterministic)
# ------------------------------------------------------------------
cat("building the network with blockwiseModules (3-6 min)...\n")
net <- blockwiseModules(datExpr, power = softPower, networkType = "signed",
                        TOMType = "signed", corType = "bicor",
                        maxBlockSize = 6000, minModuleSize = 30,
                        mergeCutHeight = 0.25, deepSplit = 2,
                        numericLabels = TRUE, pamRespectsDendro = FALSE,
                        saveTOMs = FALSE, verbose = 0)
moduleColors <- labels2colors(net$colors)
names(moduleColors) <- colnames(datExpr)
cat(sprintf("number of modules=%d (incl. unassigned grey)\n", length(unique(moduleColors))))
print(table(moduleColors))

png(file.path(OUT, "WGCNA_module_dendrogram.png"), 1400, 700, res = 120)
plotDendroAndColors(net$dendrograms[[1]],
                    moduleColors[net$blockGenes[[1]]],
                    "Module", dendroLabels = FALSE, hang = 0.03,
                    addGuide = TRUE, guideHang = 0.05, main = "Gene dendrogram and modules")
dev.off()

# Module assignment for all genes
write.csv(data.frame(gene = colnames(datExpr), module = moduleColors),
          file.path(OUT, "module_assignment_all_genes.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# 5. Module-trait association (HF + etiology)
# ------------------------------------------------------------------
traits <- data.frame(HF = as.integer(grp == "HF"),
                     NF = as.integer(grp == "NF"))
rownames(traits) <- colnames(d38$expr)
if (!is.null(et)) {
  for (lv in setdiff(unique(et), "NF")) {
    traits[[paste0("etio_", lv)]] <- as.integer(et == lv)
  }
}
traits <- traits[rownames(datExpr), , drop = FALSE]

MEs <- orderMEs(moduleEigengenes(datExpr, moduleColors)$eigengenes)
mtCor <- cor(MEs, traits, use = "p")
mtP <- corPvalueStudent(mtCor, nrow(datExpr))

png(file.path(OUT, "WGCNA_module_trait_heatmap.png"), 1500, 800, res = 120)
textMat <- paste0(signif(mtCor, 2), "\n(", signif(mtP, 1), ")")
dim(textMat) <- dim(mtCor)
labeledHeatmap(Matrix = mtCor, xLabels = colnames(traits),
               yLabels = colnames(MEs), ySymbols = colnames(MEs),
               colorLabels = FALSE, colors = blueWhiteRed(50),
               textMatrix = textMat, setStdMargins = FALSE, cex.text = 0.55,
               zlim = c(-1, 1), main = "Module-trait relationships")
dev.off()

mtOut <- data.frame(module = rownames(mtCor))
for (tl in colnames(traits)) {
  mtOut[[paste0(tl, "_cor")]] <- round(mtCor[, tl], 3)
  mtOut[[paste0(tl, "_p")]]   <- signif(mtP[, tl], 3)
}
write.csv(mtOut, file.path(OUT, "module_trait_correlation.csv"), row.names = FALSE)

hfCor <- mtCor[, "HF"]
names(hfCor) <- rownames(mtCor)
keyME <- names(sort(hfCor, decreasing = TRUE))[1]     # module most positively correlated with HF
negME <- names(sort(hfCor))[1]                        # module most negatively correlated with HF
keyMod <- sub("^ME", "", keyME)
cat(sprintf("\nmost positively HF-correlated module: %s (cor=%.2f, p=%.2g)\n", keyMod, hfCor[keyME], mtP[keyME, "HF"]))
cat(sprintf("most negatively HF-correlated module: %s (cor=%.2f, p=%.2g)\n", sub("^ME", "", negME), hfCor[negME], mtP[negME, "HF"]))

# ------------------------------------------------------------------
# 6. DLRG module enrichment (Fisher exact test)
# ------------------------------------------------------------------
inNet <- colnames(datExpr)
d <- genes56[genes56 %in% inNet]
cat(sprintf("\n%d of the 56 DLRGs entered the network\n", length(d)))
mods <- unique(moduleColors)
enr <- do.call(rbind, lapply(setdiff(mods, "grey"), function(m) {
  a <- sum(moduleColors[inNet] == m & inNet %in% d)          # DLRG in module
  b <- sum(inNet %in% d) - a                                  # DLRG outside
  ccc <- sum(moduleColors[inNet] == m) - a                    # non-DLRG in module
  dd <- length(inNet) - a - b - ccc
  ft <- fisher.test(matrix(c(a, b, ccc, dd), 2, 2), alternative = "greater")
  data.frame(module = m, module_size = a + ccc, DLRG_in_module = a,
             DLRG_total = a + b, odds_ratio = round(unname(ft$estimate), 2),
             p_value = signif(ft$p.value, 3))
}))
enr <- enr[order(enr$p_value), ]
cat("DLRG module enrichment (Fisher, top 5):\n"); print(head(enr, 5), row.names = FALSE)
write.csv(enr, file.path(OUT, "DLRG_module_enrichment.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# 7. LASSO-gene module assignment + kME + GS
# ------------------------------------------------------------------
kME <- as.data.frame(cor(datExpr, MEs, use = "p"))
GS_HF <- as.numeric(cor(datExpr, traits$HF, use = "p"))
names(GS_HF) <- colnames(datExpr)
lassoTab <- do.call(rbind, lapply(sel, function(g) {
  if (!g %in% colnames(datExpr)) return(data.frame(gene = g, module = NA, kME = NA, GS_HF = NA))
  m <- moduleColors[g]
  me <- paste0("ME", m)
  data.frame(gene = g, module = m,
             kME = round(if (me %in% names(kME)) kME[g, me] else NA, 3),
             GS_HF = round(GS_HF[g], 3))
}))
cat("\nmodule + kME of the 13 LASSO genes:\n"); print(lassoTab, row.names = FALSE)
write.csv(lassoTab, file.path(OUT, "LASSO_genes_kME.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# 8. Key-module gene list + hub genes (kME top 10)
# ------------------------------------------------------------------
keyGenes <- names(moduleColors)[moduleColors == keyMod]
me <- paste0("ME", keyMod)
kMEv <- setNames(kME[[me]], rownames(kME))   # use a named vector to preserve names (a data.frame subset would drop them)
hub <- names(sort(kMEv[keyGenes], decreasing = TRUE))[1:min(10, length(keyGenes))]
cat(sprintf("\nkey module %s (%d genes) hub top10: %s\n", keyMod, length(keyGenes), paste(hub, collapse = ", ")))
write.csv(data.frame(gene = keyGenes, kME = round(kME[keyGenes, me], 4),
                     GS_HF = round(GS_HF[keyGenes], 4)),
          file.path(OUT, "key_module_genes.csv"), row.names = FALSE)


# ================================================================
# 10. Module functional enrichment (GO-BP + KEGG)
# ================================================================
suppressPackageStartupMessages({
  library(clusterProfiler)
  library(org.Hs.eg.db)
})
options(timeout = 600)   # KEGG REST can be slow; relaxed timeout
greenG  <- names(moduleColors)[moduleColors == "green"]
yellowG <- names(moduleColors)[moduleColors == "yellow"]
cat(sprintf("enrichment targets — green module: %d genes | yellow module: %d genes\n", length(greenG), length(yellowG)))

# ---- correct universe: Entrez mapping of all network genes ----
map_id <- function(symbols) {
  suppressWarnings(eg <- AnnotationDbi::select(org.Hs.eg.db, keys = symbols,
              columns = "ENTREZID", keytype = "SYMBOL"))
  eg <- eg[!is.na(eg$ENTREZID) & !duplicated(eg$SYMBOL), ]
  setNames(eg$ENTREZID, eg$SYMBOL)
}
egAll <- map_id(colnames(datExpr))
bg <- unique(egAll)
cat(sprintf("universe: %d network genes -> %d mapped to Entrez\n",
            ncol(datExpr), length(bg)))
egG <- na.omit(egAll[greenG]); egY <- na.omit(egAll[yellowG])
cat(sprintf("green module mapped %d/%d | yellow module mapped %d/%d\n", length(egG), length(greenG), length(egY), length(yellowG)))

# ---- enrichment function (GO restricted to BP; KEGG fault-tolerant) ----
run_go <- function(entrez, tag) {
  kk <- tryCatch(enrichGO(gene = entrez, OrgDb = org.Hs.eg.db, ont = "BP",
                          universe = bg, pAdjustMethod = "BH",
                          pvalueCutoff = 0.05, qvalueCutoff = 0.1, readable = TRUE),
                 error = function(e) { cat(sprintf("[%s-GO] error: %s\n", tag, conditionMessage(e))); NULL })
  if (is.null(kk) || nrow(as.data.frame(kk)) == 0) {
    cat(sprintf("[%s-GO] no significant terms\n", tag)); return(NULL) }
  df <- as.data.frame(kk)
  write.csv(df, file.path(OUT, sprintf("enrich_%s_GO.csv", tag)), row.names = FALSE)
  cat(sprintf("[%s-GO] %d significant terms; top 3:\n", tag, nrow(df)))
  print(head(df[, c("Description", "p.adjust", "Count")], 3), row.names = FALSE)
  df
}
run_kegg <- function(entrez, tag) {
  kk <- tryCatch(enrichKEGG(gene = entrez, organism = "hsa", universe = bg,
                            pAdjustMethod = "BH", pvalueCutoff = 0.05, qvalueCutoff = 0.1),
                 error = function(e) { cat(sprintf("[%s-KEGG] failed (if unreachable, run the M10 local rerun instead): %s\n", tag, conditionMessage(e))); NULL })
  if (is.null(kk) || nrow(as.data.frame(kk)) == 0) {
    cat(sprintf("[%s-KEGG] no significant terms or network failure\n", tag)); return(NULL) }
  df <- as.data.frame(kk)
  write.csv(df, file.path(OUT, sprintf("enrich_%s_KEGG.csv", tag)), row.names = FALSE)
  cat(sprintf("[%s-KEGG] %d significant terms; top 3:\n", tag, nrow(df)))
  print(head(df[, c("Description", "p.adjust", "Count")], 3), row.names = FALSE)
  df
}

goG <- run_go(egG, "green_module")
kgG <- run_kegg(egG, "green_module")
goY <- run_go(egY, "yellow_module")
kgY <- run_kegg(egY, "yellow_module")

# ---- core 9 genes (three-cohort definition) + LASSO-13 annotation ----
lrg <- unique(toupper(trimws(readLines(file.path(DIR, "lrg332.txt"), warn = FALSE))))
lrg <- lrg[nzchar(lrg)]
gidcol <- function(df) if (".gid" %in% names(df)) ".gid" else names(df)[ncol(df)]
sig_of <- function(path) {
  df <- read.csv(path, check.names = FALSE)
  s <- subset(df, adj.P.Val < 0.05 & abs(logFC) > 0.2)
  g <- as.character(s[[gidcol(df)]]); g[toupper(g) %in% lrg]
}
g38 <- sig_of(file.path(OUT, "GSE57338_HFvsNF_all.csv"))
g06 <- sig_of(file.path(OUT, "GSE5406_HFvsNF_all.csv"))
g62d <- sig_of(file.path(OUT, "GSE116250_DCMvsNF_all.csv"))
g62i <- sig_of(file.path(OUT, "GSE116250_ICMvsNF_all.csv"))
gold7 <- Reduce(intersect, list(g38, g06, g62d, g62i))          # four-way gold standard
core9 <- intersect(intersect(g38, g06), union(g62d, g62i))      # three-cohort core
cat(sprintf("gold standard (four-way)=%d [expected 7] | core overlap (three cohorts)=%d [expected 9]\n",
            length(gold7), length(core9)))
cat(sprintf("core genes: %s\n", paste(core9, collapse = ", ")))

if (length(core9) >= 5) {
  egC <- na.omit(map_id(core9))
  goC <- run_go(egC, "core9")
  kgC <- run_kegg(egC, "core9")
}
sel13 <- as.character(read.csv(file.path(OUT, "LASSO_coefficients.csv"))$gene)
egL <- na.omit(map_id(sel13))
goL <- run_go(egL, "lasso13")

# ---- bar plot ----
bar_top <- function(df, tag, n = 10) {
  if (is.null(df)) return(invisible(NULL))
  d <- head(df[order(df$p.adjust), ], n)
  d$Description <- factor(d$Description, levels = rev(d$Description))
  p <- ggplot(d, aes(x = Count, y = Description, fill = p.adjust)) +
    geom_col() + scale_fill_gradient(low = "#b45309", high = "#fde68a") +
    labs(title = sprintf("GO-BP enrichment: %s", tag), x = "Gene count", y = NULL) +
    theme_bw(base_size = 11) + theme(plot.title = element_text(face = "bold"))
  ggsave(file.path(OUT, sprintf("GO_%s_top10.png", tag)), p, width = 9, height = 6, dpi = 300)
  cat(sprintf("written: GO_%s_top10.png\n", tag))
}
bar_top(goG, "green_module")
bar_top(goY, "yellow_module")

cat("\nstage 5 (WGCNA + enrichment) complete\n")
# ------------------------------------------------------------------
# 9. Interpretation
# ------------------------------------------------------------------
cat("\n========== stage-5 interpretation ==========\n")
cat(sprintf("  * key module %s correlates with HF, cor=%.2f (should be >0.5 with p<0.05)\n", keyMod, hfCor[keyME]))
cat(sprintf("  * most DLRG-enriched module: %s (p=%.2g, should be <0.05 and directionally consistent with the key module)\n", enr$module[1], enr$p_value[1]))
cat("  * most LASSO genes should cluster in the key module with kME>0.6 (check the table above)\n")
cat("  * self-check items: module counts / enrichment top 5 / LASSO-kME table\n")
cat("all done!\n")

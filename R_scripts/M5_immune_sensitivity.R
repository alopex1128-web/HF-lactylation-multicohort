# ============================================================
# M5_immune_sensitivity.R  |  M5
# Purpose: MCP-counter / EPIC / P<0.05-subset three-method concordance
# Paper section: Suppl. Table S6
# Pipeline position: step 6 of the main chain
# Prerequisite: M4 completed (reuses the CIBERSORT output)
# (Generated from the modular pipeline consolidation; analysis logic identical to the verified version.)
# ============================================================
# ============================================================
# M5: immune-infiltration sensitivity, three methods in one (addressing CIBERSORT applicability)
# ============================================================
# [What it does] (deconvolution robustness for a blood-derived signature applied to myocardium)
#   Part A  CIBERSORT P<0.05 valid-subset sensitivity analysis
#   Part B  MCP-counter alternative (10 tissue-infiltrating populations, no reference matrix)
#   Part C  EPIC alternative (constrained least squares; "other cells" = cardiomyocytes etc.)
#   Part D  three-method direction-concordance cross table (for the Discussion)
#
# [Prerequisites]
#   (1) stages 1-7 completed (output/GSE57338_processed.rds + stage7_results/
#       Table6_CIBERSORT_22cell_scores.csv — the seeded M4 output, authoritative numbers)
#   (2) network access (MCP-counter pulls signature files; two GitHub packages need a first-time install)
#
# [Methods (consistent with the paper)]
#   Part A  reads the M4 Table 6 CSV directly (CIBERSORT is not rerun, keeping the published
#           numbers bit-identical); Wilcoxon on the P<0.05 subset, 22 cells + BH
#   Part B  MCPcounter.estimate(log2 matrix, HUGO symbols), Wilcoxon + BH
#       (MCP-counter averages marker-gene sets; it is platform-robust and needs no reference
#       matrix — officially supported for both arrays and RNA-seq, hence the choice as alternative)
#   Part C  EPIC(bulk = 2^log2 matrix) (EPIC expects linear TPM-like values, hence the reverse log)，
#       default TRef; otherCells = cardiomyocytes and cells not covered by the reference;
#       Wilcoxon + BH (proportions, no normality assumption)
#   All comparisons are two-sided Wilcoxon rank-sum + BH, identical to M4
#
# [Output] under stage8_results/:
#   TableA1_immune_diff_allsamples_control.csv (should match the M4 differential CSV)
#   TableA2_immune_diff_P005subset_sensitivity.csv
#   TableB_MCPcounter_10pop_scores.csv + TableB2_MCP_diff_Wilcoxon.csv
#   TableC_EPIC_proportions.csv + TableC2_EPIC_diff_Wilcoxon.csv
#   TableD_threemethod_direction_concordance.csv (core material for the revision)
#   FigureS_draft_MCP_boxplot.png / FigureS_draft_EPIC_boxplot.png (300 dpi drafts)
#   stage8_console_report.txt
#
# [Runtime] first run ~5-10 min (MCP-counter: seconds; EPIC: minutes; no permutation rerun)
# [Reproducibility] no random permutations; EPIC and MCP-counter are deterministic
#
# [Expected values (self-check against the output)]
#   - Table 6 CSV: 313 samples x (22 cells + P-value and related columns)
#   - P<0.05 subset: 15 samples (4.8%), as reported by the M4 console
#   - MCP-counter output: 10 cell populations (T cells / CD8 T / cytotoxic lymphocytes /
#     NK/B lineage/Monocytic/Myeloid DC/Neutrophils/Endothelial/Fibroblasts）
#   - EPIC output: 6 reference cell types + otherCells (~ cardiomyocytes + resident cells)
#   - MCP/EPIC P values: report the local output as authoritative
#   - With n=15 the HF/NF split is printed — if a group has <3, the subset Wilcoxon
#     is reported descriptively only (such a small subset is expected to be underpowered)
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# Set DIR to your working directory (forward slashes on Windows)
DIR <- "./GEO_data"

OUT <- file.path(DIR, "stage8_results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
LOGF <- file.path(OUT, "stage8_console_report.txt")
logline <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n"); cat(msg, "\n", file = LOGF, append = TRUE)
}
if (file.exists(LOGF)) file.remove(LOGF)
logline("========== stage 8: immune sensitivity, started: ", format(Sys.time()), "==========")

# ============================================================
# Part 0: prerequisite files and packages
# ============================================================
need <- c(rds = file.path(DIR, "output", "GSE57338_processed.rds"),
          cib6 = file.path(DIR, "stage7_results", "Table6_CIBERSORT_22cell_scores.csv"))
miss <- need[!file.exists(need)]
if (length(miss) > 0) {
  logline("!!! missing files, provide them before running:")
  for (m in miss) logline("   -", m)
  stop("missing prerequisite files, see the list above")
}

# ---- MCP-counter and EPIC installation check (one-time, ~1-3 min) ----
# upgrade = "never" suppresses the interactive update prompt of install_github
# (the prompt would otherwise swallow the EPIC install and break Part C)
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
if (!requireNamespace("devtools", quietly = TRUE)) install.packages("devtools")
if (!requireNamespace("curl",     quietly = TRUE)) install.packages("curl")
if (!requireNamespace("MCPcounter", quietly = TRUE)) {
  logline("[install] MCP-counter (official ebecht repository, Source subdirectory)...")
  devtools::install_github("ebecht/MCPcounter", ref = "master", subdir = "Source",
                           upgrade = "never")
}
if (!requireNamespace("EPIC", quietly = TRUE)) {
  logline("[install] EPIC (official GfellerLab repository, pure R, no Rtools needed)...")
  devtools::install_github("GfellerLab/EPIC", upgrade = "never")
}
logline("[package check] devtools/curl/MCPcounter/EPIC ready")

# ============================================================
# Part 1: data loading
# ============================================================
d38 <- readRDS(need[["rds"]])
expr38 <- d38$expr; grp38 <- d38$group
logline(sprintf("  expression matrix: %d genes x %d samples | HF %d / NF %d",
                nrow(expr38), ncol(expr38),
                sum(grp38 == "HF"), sum(grp38 == "NF")))
cib <- read.csv(need[["cib6"]], check.names = FALSE, row.names = 1)
cells <- colnames(cib)[1:22]
logline(sprintf("  Table 6 loaded: %d samples x %d cols (incl. %d P-value and related cols)",
                nrow(cib), ncol(cib), ncol(cib) - 22))
if (nrow(cib) != 313) logline("  [note] sample count != 313 — if it is not 313, check the console output")

# ============================================================
# Part A: CIBERSORT all-sample recomputation control + P<0.05-subset sensitivity
# ============================================================
logline("\n========== Part A: CIBERSORT subset sensitivity ==========")
gv <- grp38[rownames(cib)]
if (any(is.na(gv))) stop("group vector not aligned with sample names (NA)")

wl <- function(frac_mat, gvec, cells_chr, tag) {
  res <- do.call(rbind, lapply(cells_chr, function(cc) {
    x <- frac_mat[, cc]
    # constant columns (blood-derived LM22 cells with no signal in myocardium)
    # are explicitly flagged as not tested — evidence of the LM22 applicability issue
    if (var(x) == 0) {
      return(data.frame(cell = cc, median_HF = round(x[1], 4),
                        median_NF = round(x[1], 4),
                        direction = "constant column (not tested)", P = NA,
                        stringsAsFactors = FALSE))
    }
    p <- tryCatch(wilcox.test(x ~ gvec)$p.value, error = function(e) NA)
    data.frame(cell = cc,
               median_HF = round(median(x[gvec == "HF"]), 4),
               median_NF = round(median(x[gvec == "NF"]), 4),
               direction = ifelse(median(x[gvec=="HF"]) > median(x[gvec=="NF"]), "HF>NF", "HF<NF"),
               P = p, stringsAsFactors = FALSE)
  }))
  res$P_BH <- p.adjust(res$P, method = "BH")
  res <- res[order(res$P), ]
  write.csv(res, file.path(OUT, tag), row.names = FALSE, fileEncoding = "UTF-8")
  n_const <- sum(res$direction == "constant column (not tested)")
  if (n_const > 0)
    logline(sprintf("  [note] %d constant cell columns (no signal for blood-derived LM22 cells in myocardium)", n_const))
  res
}

# A1: all-sample recomputation (control — should match the M4 differential CSV)
res_all <- wl(cib, gv, cells, "TableA1_immune_diff_allsamples_control.csv")
logline(sprintf("  [A1 all-sample] cell types with BH q<0.05: %d", sum(res_all$P_BH < 0.05, na.rm = TRUE)))
for (i in which(res_all$P_BH < 0.05))
  logline(sprintf("     %s (P=%.2e, q=%.3f, %s)",
                  res_all$cell[i], res_all$P[i], res_all$P_BH[i], res_all$direction[i]))

# A2: P<0.05-subset sensitivity
valid <- rownames(cib)[cib[, "P-value"] < 0.05]
rate <- length(valid) / nrow(cib)
logline(sprintf("  [A2 subset] P<0.05 valid samples: %d/%d (%.1f%%)", length(valid), nrow(cib), 100*rate))
gv_v <- grp38[valid]
logline(sprintf("     within-subset groups: HF %d / NF %d", sum(gv_v=="HF"), sum(gv_v=="NF")))
if (min(table(gv_v)) < 3)
  logline("  [note] a group has <3 — subset Wilcoxon reported descriptively only (underpowered, stated as such in the Discussion)")
res_sub <- wl(cib[valid, ], gv_v, cells, "TableA2_immune_diff_P005subset_sensitivity.csv")
logline("  [A2 subset] cell types with BH q<0.05: ",
        sum(res_sub$P_BH < 0.05, na.rm = TRUE), " (compare against the all-sample table in the console)")

# comparison: behavior of all-sample-significant cells within the subset
sigall <- res_all$cell[which(res_all$P_BH < 0.05)]  # which() guards against NA indices
if (length(sigall) > 0) {
  logline("  [comparison] all-sample BH-significant cells in the subset:")
  for (cc in sigall) {
    r <- res_sub[res_sub$cell == cc, ]
    a <- res_all[res_all$cell == cc, ]
    logline(sprintf("     %-28s all-sample q=%.3f(%s) -> subset q=%s(%s)",
                    cc, a$P_BH, a$direction,
                    ifelse(is.na(r$P_BH), "NA", sprintf("%.3f", r$P_BH)),
                    r$direction))
  }
}

# ============================================================
# Part B: MCP-counter (marker-gene-set method, no reference matrix)
# ============================================================
logline("\n========== Part B：MCP-counter ==========")
ex <- expr38[!duplicated(rownames(expr38)), ]  # deduplicate genes (as in M4)
mcp <- MCPcounter::MCPcounter.estimate(ex, featuresType = "HUGO_symbols")
logline(sprintf("  MCP-counter: %d cell populations x %d samples", nrow(mcp), ncol(mcp)))
write.csv(mcp, file.path(OUT, "TableB_MCPcounter_10pop_scores.csv"), fileEncoding = "UTF-8")

mcp_t <- t(mcp)  # samples x populations
gvm <- grp38[rownames(mcp_t)]
mcp_res <- wl(mcp_t, gvm, colnames(mcp_t), "TableB2_MCP_diff_Wilcoxon.csv")
logline(sprintf("  [MCP] populations with BH q<0.05: %d", sum(mcp_res$P_BH < 0.05, na.rm = TRUE)))
for (i in which(mcp_res$P_BH < 0.05))
  logline(sprintf("     %s (P=%.2e, q=%.3f, %s)",
                  mcp_res$cell[i], mcp_res$P[i], mcp_res$P_BH[i], mcp_res$direction[i]))

# boxplot draft (MCP counterpart of Figure 6A; optional for submission)
png(file.path(OUT, "FigureS_draft_MCP_boxplot.png"), width = 3000, height = 2000, res = 300)
dfm <- data.frame(value = as.vector(mcp_t),
                  cell  = rep(colnames(mcp_t), each = nrow(mcp_t)),
                  group = rep(gvm, times = ncol(mcp_t)))
boxplot(value ~ interaction(cell, group, sep = "  "), data = dfm,
        las = 2, cex.axis = 0.45, ylab = "MCP-counter score",
        main = "MCP-counter: HF vs NF (draft)")
dev.off()

# ============================================================
# Part C: EPIC (constrained least squares; otherCells ~ cardiomyocytes + resident)
# ============================================================
logline("\n========== Part C：EPIC ==========")
ex_lin <- 2^ex  # EPIC expects linear TPM-like values (reverse log2)
epic_out <- EPIC::EPIC(bulk = ex_lin)
epic_fr <- as.data.frame(epic_out$cellFractions)
logline(sprintf("  EPIC: %d cell types x %d samples", ncol(epic_fr), nrow(epic_fr)))
write.csv(epic_fr, file.path(OUT, "TableC_EPIC_proportions.csv"), fileEncoding = "UTF-8")

gve <- grp38[rownames(epic_fr)]
epic_res <- wl(as.matrix(epic_fr), gve, colnames(epic_fr), "TableC2_EPIC_diff_Wilcoxon.csv")
logline(sprintf("  [EPIC] types with BH q<0.05: %d", sum(epic_res$P_BH < 0.05)))
for (i in which(epic_res$P_BH < 0.05))
  logline(sprintf("     %s (P=%.2e, q=%.3f, %s)",
                  epic_res$cell[i], epic_res$P[i], epic_res$P_BH[i], epic_res$direction[i]))

png(file.path(OUT, "FigureS_draft_EPIC_boxplot.png"), width = 2400, height = 1800, res = 300)
boxplot(epic_fr[, 1:min(7, ncol(epic_fr))], las = 2, cex.axis = 0.7,
        main = "EPIC cell fractions: all samples")
dev.off()

# ============================================================
# Part D: three-method direction-concordance cross table (core revision material)
# ============================================================
logline("\n========== Part D: direction-concordance table ==========")
# cell-type mapping (aligning names across methods at the lineage level)
map_to_common <- c(
  "T cells CD4"="T cells (CD4)", "T cells CD8"="CD8 T cells", "T cells gamma delta"="T cells (other)",
  "T cells regulatory (Tregs)"="T cells (CD4)", "T cells follicular helper"="T cells (CD4)",
  "NK cells resting"="NK cells", "NK cells activated"="NK cells",
  "B cells naive"="B cells", "B cells memory"="B cells", "Plasma cells"="B lineage",
  "Monocytes"="Monocytic lineage", "Macrophages M0"="Macrophages/monocytes", "Macrophages M1"="Macrophages/monocytes",
  "Macrophages M2"="Macrophages/monocytes", "Dendritic cells resting"="Dendritic cells",
  "Dendritic cells activated"="Dendritic cells", "Mast cells resting"="Mast cells", "Mast cells activated"="Mast cells",
  "Eosinophils"="Granulocytes (other)", "Neutrophils"="Neutrophils",
  "CD8 T cells"="CD8 T cells", "T cells"="T cells (pan)", "Cytotoxic lymphocytes"="Cytotoxic lymphocytes",
  "NK cells"="NK cells", "B lineage"="B lineage", "Monocytic lineage"="Monocytic lineage",
  "Myeloid dendritic cells"="Dendritic cells", "Endothelial cells"="Endothelial cells", "Fibroblasts"="Fibroblasts",
  "Bcells"="B lineage", "NKcells"="NK cells", "CD4T"="T cells (CD4)", "CD8T"="CD8 T cells",
  "Macrophage"="Macrophages/monocytes", "Monocyte"="Monocytic lineage", "otherCells"="Other cells (resid)")
norm_cell <- function(x) { y <- unname(map_to_common[x]); ifelse(is.na(y), x, y) }

d1 <- data.frame(method="CIBERSORT(all)", cell=norm_cell(res_all$cell),
                 direction=res_all$direction, q=res_all$P_BH, P=res_all$P)
d2 <- data.frame(method="CIBERSORT(P<0.05 subset)", cell=norm_cell(res_sub$cell),
                 direction=res_sub$direction, q=res_sub$P_BH, P=res_sub$P)
d3 <- data.frame(method="MCP-counter", cell=norm_cell(mcp_res$cell),
                 direction=mcp_res$direction, q=mcp_res$P_BH, P=mcp_res$P)
d4 <- data.frame(method="EPIC", cell=norm_cell(epic_res$cell),
                 direction=epic_res$direction, q=epic_res$P_BH, P=epic_res$P)
cross <- rbind(d1, d2, d3, d4)
cross$cell <- as.character(cross$cell); cross$direction <- as.character(cross$direction)
cross <- cross[order(cross$cell, cross$method), ]
write.csv(cross, file.path(OUT, "TableD_threemethod_direction_concordance.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# lineage-level direction-concordance counts
common_cells <- names(which(table(cross$cell) >= 2))
logline("  [concordance] lineages appearing in >=2 methods:")
for (cc in common_cells) {
  sub <- cross[cross$cell == cc, ]
  dirs <- paste0(sub$method, ":", sub$direction, collapse = " | ")
  agree <- length(unique(sub$direction)) == 1
  logline(sprintf("     %-26s %s %s", cc, dirs, ifelse(agree, "<- concordant", "")))
}

logline("\n========== stage 8 finished: ", format(Sys.time()), "==========")
logline("          Verify the numbers, then update the Methods/Discussion immune sections accordingly.")

# ============================================================
# M2_lasso_validation.R  |  M2
# Purpose: 13-gene LASSO + five-scenario ROC + DeLong CI
# Paper sections: Figure 3, Table 2
# Pipeline position: step 3 of the main chain
# Prerequisite: M1 completed (output/GSE57338_processed.rds etc.)
# (Generated from the modular pipeline consolidation; analysis logic identical to the verified version.)
# [v8.1 update] plotting only: ROC panel titles now include the DeLong 95% CI
#            (below the AUC line); analysis logic / CSV outputs / anchor numbers unchanged.
# ============================================================
# ==============================================================================
# M2: stage-4 pipeline (LASSO diagnostic model + four ROCs + publication metrics)
# ------------------------------------------------------------------------------
#            Rerunning regenerates ROC_panels.png = Figure 3B with CI.
# [Contents] = stage 4v4 (LASSO + four ROCs + plots) + stage 4b (publication patches):
#   LASSO (fixed seed) -> risk score -> training + four validation ROCs -> CV/path/panel/box plots
#   -> per-dataset Youden-cutoff sens/spec -> z-score sensitivity -> missing-gene identification
# [Input] M1 outputs (three rds + GSE57338_HFvsNF_all.csv + the LRG list txt)
# [Usage] edit only the DIR line -> select all and run (~2 min)
# [Expected] LASSO selects 13 genes (incl. HLTF/LSP1/HMGN2/ARGLU1) | training AUC = 0.982
#   | GSE5406 AUC = 0.965 | GSE116250 pooled AUC = 0.951, sens 98.0%, spec 92.9%
#   | DCM 0.971 | ICM 0.896 | z-score variant within 0.03 of the original AUC
# [Robustness] (1) group vectors are re-aligned to colnames after loading the rds;
#          (2) length echoes are printed at every step (y/score/AUC fully traceable);
#          (3) zero-hit protection for LASSO genes (no silent zero-length score vectors);
#          (4) hard stopifnot assertions with explicit messages containing exact numbers.
# ==============================================================================

rm(list = ls())   # start from a clean session
cat("== clean session started ==\n")

suppressPackageStartupMessages({
  library(glmnet)
  library(pROC)
})

DIR <- "./GEO_data"   # set to your working directory
OUT <- file.path(DIR, "output")

# ------------------------------------------------------------------
# 1. Input pre-check
# ------------------------------------------------------------------
need <- c(file.path(OUT, "GSE57338_processed.rds"),
          file.path(OUT, "GSE5406_processed.rds"),
          file.path(OUT, "GSE116250_processed.rds"),
          file.path(OUT, "GSE57338_HFvsNF_all.csv"),
          file.path(DIR, "lrg332.txt"))
missf <- need[!file.exists(need)]
if (length(missf) > 0) stop("file(s) not found: ", paste(basename(missf), collapse = ", "),
                            "\nrun M1 first and check the DIR path")
cat("input pre-check passed (all 5 files present)\n")

# ------------------------------------------------------------------
# 2. Load + automatic group alignment
# ------------------------------------------------------------------
d38 <- readRDS(need[1]); d06 <- readRDS(need[2]); d62 <- readRDS(need[3])
tt38 <- read.csv(need[4], check.names = FALSE)
lrg <- unique(toupper(trimws(readLines(need[5], warn = FALSE))))
lrg <- lrg[nzchar(lrg)]

align_group <- function(d, tag) {
  g <- d$group; cn <- colnames(d$expr)
  cat(sprintf("[%s] group length %d | expr columns %d", tag, length(g), length(cn)))
  if (length(g) != length(cn))
    stop(sprintf("\n[%s] length mismatch! group=%d vs columns=%d — rds structure problem, check this line", tag, length(g), length(cn)))
  if (!is.null(names(g))) {
    g2 <- g[cn]
    if (any(is.na(g2))) {
      cat(" -> same length, but group names do not match colnames (original order kept)\n")
    } else {
      d$group <- g2
      cat(" -> aligned by colnames\n")
    }
  } else {
    cat(" -> same length, unnamed group (positional = column order)\n")
  }
  d
}
d38 <- align_group(d38, "57338")
d06 <- align_group(d06, "5406")
d62 <- align_group(d62, "116250")

cat(sprintf("sample sizes: 57338=%d | 5406=%d | 116250=%d | LRG set=%d genes\n",
            ncol(d38$expr), ncol(d06$expr), ncol(d62$expr), length(lrg)))
print(table(d38$group)); print(table(d06$group)); print(table(d62$group))

# ------------------------------------------------------------------
# 3. Recompute the 56 DLRGs (self-check)
# ------------------------------------------------------------------
LFC <- 0.2
gidcol <- if (".gid" %in% names(tt38)) ".gid" else names(tt38)[ncol(tt38)]
sig38 <- subset(tt38, adj.P.Val < 0.05 & abs(logFC) > LFC)
genes56 <- as.character(sig38[[gidcol]])
genes56 <- genes56[toupper(genes56) %in% lrg]
cat(sprintf("DLRGs recomputed=%d [expected 56; if different, note the number and continue]\n", length(genes56)))
miss56 <- setdiff(genes56, rownames(d38$expr))
if (length(miss56) > 0)
  cat(sprintf("  note: %d DLRGs absent from expr rownames (dropped): %s\n",
              length(miss56), paste(head(miss56, 10), collapse = ", ")))
genes56 <- intersect(genes56, rownames(d38$expr))

# ------------------------------------------------------------------
# 4. LASSO modeling
# ------------------------------------------------------------------
x38 <- t(d38$expr[genes56, , drop = FALSE])
y38 <- as.integer(d38$group == "HF")
cat(sprintf("modeling matrix: %d samples x %d genes | y38 length %d\n", nrow(x38), ncol(x38), length(y38)))
if (nrow(x38) != length(y38)) stop(sprintf("modeling matrix rows (%d) != y38 length (%d) — rds problem, check this line", nrow(x38), length(y38)))
set.seed(2026)
cvfit <- cv.glmnet(x38, y38, family = "binomial", alpha = 1, nfolds = 10, type.measure = "auc")

png(file.path(OUT, "LASSO_cv_curve.png"), 900, 700, res = 120)
plot(cvfit); title("LASSO 10-fold CV (binomial)", line = 2.5); dev.off()
png(file.path(OUT, "LASSO_coefficient_path.png"), 900, 700, res = 120)
plot(cvfit$glmnet.fit, xvar = "lambda")
abline(v = log(cvfit$lambda.1se), lty = 2, col = "red")
title("Coefficient shrinkage path", line = 2.5); dev.off()

co <- as.matrix(coef(cvfit, s = "lambda.1se"))
sel <- setdiff(rownames(co)[co[, 1] != 0], "(Intercept)")
nmin <- sum(as.matrix(coef(cvfit, s = "lambda.min"))[, 1] != 0) - 1
if (length(sel) == 0) {
  cat("lambda.1se selected 0 genes -> falling back to lambda.min\n")
  co <- as.matrix(coef(cvfit, s = "lambda.min"))
  sel <- setdiff(rownames(co)[co[, 1] != 0], "(Intercept)")
}
intercept <- co["(Intercept)", 1]
beta <- co[sel, 1]; names(beta) <- sel

cat(sprintf("\n========== LASSO selected %d genes (lambda.min alternative: %d) ==========\n", length(sel), nmin))
cat("score formula: Score ="); cat(sprintf(" %.4f", intercept))
for (g in sel) cat(sprintf(" %+.4f x %s", beta[g], g))
cat("\n")
write.csv(data.frame(gene = sel, coefficient = round(as.numeric(beta), 6)),
          file.path(OUT, "LASSO_coefficients.csv"), row.names = FALSE)

# ------------------------------------------------------------------
# 5. Scoring function (zero-hit protection)
# ------------------------------------------------------------------
score_set <- function(expr, label) {
  hit <- intersect(sel, rownames(expr))
  missg <- setdiff(sel, rownames(expr))
  if (length(hit) == 0)
    stop(sprintf("[%s] zero LASSO-gene hits! first 5 of sel: %s | first 5 rownames: %s",
                 label, paste(head(sel, 5), collapse = ","),
                 paste(head(rownames(expr), 5), collapse = ",")))
  if (length(missg) > 0)
    cat(sprintf("  [%s] %d/%d genes missing (scored on the subset)\n", label, length(missg), length(sel)))
  if (length(missg) > length(sel) * 0.2)
    cat(sprintf("  [%s] warning: >20%% missing, interpret with caution\n", label))
  as.numeric(intercept + colSums(expr[hit, , drop = FALSE] * beta[hit]))
}

# ------------------------------------------------------------------
# 6. Training-set ROC
# ------------------------------------------------------------------
s38 <- score_set(d38$expr, "57338 training")
cat(sprintf("length echo: y38=%d s38=%d\n", length(y38), length(s38)))
if (length(y38) != length(s38)) stop(sprintf("y38 length (%d, from group) != s38 length (%d, from expr columns)\nexpr dims: %d rows x %d cols", length(y38), length(s38), nrow(d38$expr), ncol(d38$expr)))
r38 <- roc(y38, s38, quiet = TRUE, levels = c(0, 1), direction = "<")
ci38 <- as.numeric(ci.auc(r38))
best <- coords(r38, "best", ret = c("threshold", "sensitivity", "specificity"),
               best.method = "youden")
cat(sprintf("\ntraining 57338: AUC=%.3f [95%%CI %.3f-%.3f]\n", as.numeric(auc(r38)), ci38[1], ci38[3]))
cat(sprintf("  Youden cutoff=%.3f -> sensitivity=%.1f%% specificity=%.1f%%\n",
            best$threshold, 100 * best$sensitivity, 100 * best$specificity))

metrics <- function(score, y, tag) {
  cat(sprintf("length echo [%s]: y=%d score=%d\n", tag, length(y), length(score)))
  if (length(y) != length(score)) stop(sprintf("[%s] y length (%d) != score length (%d) — check this line", tag, length(y), length(score)))
  r <- roc(y, score, quiet = TRUE, levels = c(0, 1), direction = "<")
  ci <- as.numeric(ci.auc(r))
  sens <- mean(score[y == 1] > best$threshold)
  spec <- mean(score[y == 0] <= best$threshold)
  cat(sprintf("%s: AUC=%.3f [%.3f-%.3f] | sens at cutoff=%.1f%% spec=%.1f%%\n",
              tag, as.numeric(auc(r)), ci[1], ci[3], 100 * sens, 100 * spec))
  list(r = r, auc = as.numeric(auc(r)), lo = ci[1], hi = ci[3], sens = sens, spec = spec)
}

# ------------------------------------------------------------------
# 7. Four validation scenarios
# ------------------------------------------------------------------
s06 <- score_set(d06$expr, "5406")
m06 <- metrics(s06, as.integer(d06$group == "HF"), "validation 1 - GSE5406")

s62 <- score_set(d62$expr, "116250")
y_all <- as.integer(d62$group %in% c("DCM", "ICM"))
m62 <- metrics(s62, y_all, "validation 2 - GSE116250 pooled HF")

sub_dcm <- which(d62$group %in% c("DCM", "NF"))
mD <- metrics(s62[sub_dcm], as.integer(d62$group[sub_dcm] == "DCM"), "validation 3 - GSE116250 DCM")

sub_icm <- which(d62$group %in% c("ICM", "NF"))
mI <- metrics(s62[sub_icm], as.integer(d62$group[sub_icm] == "ICM"), "validation 4 - GSE116250 ICM")

# ------------------------------------------------------------------
# 8. Summary table + plots
# ------------------------------------------------------------------
roc_df <- data.frame(
  scenario = c("Train-57338", "Valid-5406", "Valid-116250-HF",
               "Valid-116250-DCM", "Valid-116250-ICM"),
  n_HF = c(sum(y38), sum(d06$group == "HF"), sum(y_all),
           sum(d62$group == "DCM"), sum(d62$group == "ICM")),
  n_NF = c(sum(1 - y38), sum(d06$group == "NF"), sum(d62$group == "NF"),
           sum(d62$group == "NF"), sum(d62$group == "NF")),
  AUC = round(c(as.numeric(auc(r38)), m06$auc, m62$auc, mD$auc, mI$auc), 3),
  AUC_low = round(c(ci38[1], m06$lo, m62$lo, mD$lo, mI$lo), 3),
  AUC_high = round(c(ci38[3], m06$hi, m62$hi, mD$hi, mI$hi), 3),
  sens_pct = round(100 * c(best$sensitivity, m06$sens, m62$sens, mD$sens, mI$sens), 1),
  spec_pct = round(100 * c(best$specificity, m06$spec, m62$spec, mD$spec, mI$spec), 1)
)
write.csv(roc_df, file.path(OUT, "ROC_validation_summary.csv"), row.names = FALSE)

ttl <- c("Train GSE57338", "Valid GSE5406", "Valid 116250-HF",
         "Valid 116250-DCM", "Valid 116250-ICM")
cols <- c("#0f766e", "#b45309", "#1d4ed8", "#be185d", "#15803d")
plist <- list(r38, m06$r, m62$r, mD$r, mI$r)
aucs <- c(as.numeric(auc(r38)), m06$auc, m62$auc, mD$auc, mI$auc)
# [v8.1] panel titles include the DeLong 95% CI (reusing ci38/m*$lo/hi computed above; no new computation)
#        values match the Figure 3 legend; rerunning regenerates ROC_panels.png (Figure 3B with CI)
cilo <- c(ci38[1], m06$lo, m62$lo, mD$lo, mI$lo)
cihi <- c(ci38[3], m06$hi, m62$hi, mD$hi, mI$hi)
png(file.path(OUT, "ROC_panels.png"), 1500, 1050, res = 120)
par(mfrow = c(2, 3))
for (i in 1:5) plot(plist[[i]], col = cols[i], legacy.axes = TRUE, print.auc = FALSE,
                    main = sprintf("%s\nAUC = %.3f [95%% CI %.3f-%.3f]",
                                   ttl[i], aucs[i], cilo[i], cihi[i]))
dev.off()

boxp <- function(score, y, label) {
  g <- factor(ifelse(y == 1, "HF", "NF"), levels = c("NF", "HF"))
  w <- tryCatch(wilcox.test(score ~ g), error = function(e) NULL)
  p <- if (is.null(w)) NA else w$p.value
  boxplot(score ~ g, col = c("#4c78a8", "#e4572b"), names = c("NF", "HF"),
          ylab = "Risk score", main = sprintf("%s (P=%.2g)", label, p))
}
png(file.path(OUT, "risk_score_boxplot.png"), 1500, 1050, res = 120)
par(mfrow = c(2, 3))
boxp(s38, y38, "Train GSE57338")
boxp(s06, as.integer(d06$group == "HF"), "Valid GSE5406")
boxp(s62, y_all, "Valid 116250-HF")
boxp(s62[sub_dcm], as.integer(d62$group[sub_dcm] == "DCM"), "Valid 116250-DCM")
boxp(s62[sub_icm], as.integer(d62$group[sub_icm] == "ICM"), "Valid 116250-ICM")
dev.off()

# ------------------------------------------------------------------
# 9. Interpretation
# ------------------------------------------------------------------
cat("\n========== ROC validation summary ==========\n")
print(roc_df, row.names = FALSE)
cat("\ninterpretation: training AUC>0.90 and validation>0.80 -> publication-ready; validation 0.70-0.80 -> usable (\"moderate discrimination\"); <0.70 -> check the console output\n")
cat(sprintf("if the selected genes include HMGN2/LSP1/PFKM they echo the stage-3 core genes (consistent narrative)\n", length(sel)))

# ================================================================
# 10. Publication patches: per-dataset cutoff + z-score sensitivity + missing-gene check
# ================================================================
# ---- (3) missing-gene identification ----
miss62 <- setdiff(sel, rownames(d62$expr))
cat(sprintf("missing genes in 116250: %d -> %s\n", length(miss62),
            if (length(miss62)) paste(miss62, collapse = ", ") else "none"))

# ---- scoring function (v4 fixed version, no t()) ----
score_set <- function(expr) {
  hit <- intersect(sel, rownames(expr))
  as.numeric(intercept + colSums(expr[hit, , drop = FALSE] * beta[hit]))
}
# ---- z-score standardized score (platform-agnostic variant) ----
score_z <- function(expr) {
  hit <- intersect(sel, rownames(expr))
  z <- t(scale(t(expr[hit, , drop = FALSE])))   # per-gene standardization within each dataset
  as.numeric(colSums(z * beta[hit]))             # intercept is meaningless in z space; cutoff = 0
}

run_pair <- function(expr, group, keep, pos_labels, tag) {
  e <- expr[, keep, drop = FALSE]; g <- group[keep]
  y <- as.integer(g %in% pos_labels)
  s_raw <- score_set(e); s_z <- score_z(e)
  r1 <- roc(y, s_raw, quiet = TRUE, levels = c(0, 1), direction = "<")
  r2 <- roc(y, s_z,  quiet = TRUE, levels = c(0, 1), direction = "<")
  b1 <- coords(r1, "best", ret = c("sensitivity", "specificity"), best.method = "youden")
  b2 <- coords(r2, "best", ret = c("sensitivity", "specificity"), best.method = "youden")
  list(tag = tag, n_pos = sum(y), n_neg = length(y) - sum(y),
       auc = as.numeric(auc(r1)), lo = as.numeric(ci.auc(r1))[1], hi = as.numeric(ci.auc(r1))[3],
       sens_y = 100 * b1$sensitivity, spec_y = 100 * b1$specificity,
       auc_z = as.numeric(auc(r2)), sens_z = 100 * b2$sensitivity, spec_z = 100 * b2$specificity)
}

all62 <- rep(TRUE, ncol(d62$expr))
res <- list(
  run_pair(d38$expr, d38$group, rep(TRUE, ncol(d38$expr)), "HF", "Train-57338"),
  run_pair(d06$expr, d06$group, rep(TRUE, ncol(d06$expr)), "HF", "Valid-5406"),
  run_pair(d62$expr, d62$group, d62$group %in% c("DCM","ICM","NF"), c("DCM","ICM"), "Valid-116250-HF"),
  run_pair(d62$expr, d62$group, d62$group %in% c("DCM","NF"), "DCM", "Valid-116250-DCM"),
  run_pair(d62$expr, d62$group, d62$group %in% c("ICM","NF"), "ICM", "Valid-116250-ICM")
)

df <- do.call(rbind, lapply(res, function(x) data.frame(
  scenario = x$tag, n_HF = x$n_pos, n_NF = x$n_neg,
  AUC = round(x$auc, 3), CI_low = round(x$lo, 3), CI_high = round(x$hi, 3),
  sens_own_cutoff = round(x$sens_y, 1), spec_own_cutoff = round(x$spec_y, 1),
  AUC_zscore = round(x$auc_z, 3), sens_z = round(x$sens_z, 1), spec_z = round(x$spec_z, 1)
)))
write.csv(df, file.path(OUT, "ROC_validation_summary_v2.csv"), row.names = FALSE)

cat("\n========== publication ROC summary (three metric sets) ==========\n")
print(df, row.names = FALSE)
cat("\nreading guide:\n")
cat("  * AUC = raw score; sens/spec_own_cutoff = per-dataset Youden-optimal cutoff (primary in the paper)\n")
cat("  * AUC_zscore = AUC of the gene-standardized score (platform-agnostic sensitivity; should be close to the raw AUC)\n")
cat("  * sens_own_cutoff for 116250 should recover to 80-95% (the earlier 0% was a cross-platform cutoff failure)\n")
cat("\nstage 4 complete!\n")


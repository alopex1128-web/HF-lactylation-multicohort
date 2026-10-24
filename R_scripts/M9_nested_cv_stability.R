# ============================================================
# M9_nested_cv_stability.R
# Purpose: nested-CV outer-loop AUC + 200x bootstrap + CI summary
# Paper section: Results, LASSO paragraph
# Pipeline position: step 10 of the main chain
# Prerequisite: M2 completed (needs the M1 expression matrix + ROC_validation_summary_v2.csv)
# (Generated from the modular pipeline consolidation; analysis logic identical
#  to the verified version.)
# ============================================================
# ============================================================
# M9: nested CV + bootstrap stability + DeLong CI summary
# ============================================================
# [Background] The training-set AUC of 0.982 is optimistically biased because
#   feature selection and evaluation used the same data (Varma & Simon 2006).
#   This module adds three quantities:
#   (1) nested-CV outer-loop AUC (optimism-corrected performance estimate);
#   (2) bootstrap gene-selection stability (re-selection frequency of each of
#   the 13 genes);
#   (3) 95% CIs of the five-scenario AUCs (DeLong CIs already computed in M2;
#   summarized here).
# [Logic] (1) rebuild the M2 modelling matrix (56 DLRGs x 313 samples,
#   seed = 2026); (2) nested CV: 10 outer folds, inner cv.glmnet re-selecting
#   lambda.1se within each fold, pooled AUC from out-of-fold predictions;
#   (3) bootstrap: 200 resamples with replacement, re-running cv.glmnet gene
#   selection each time and recording per-gene selection frequency;
#   (4) read ROC_validation_summary_v2.csv and print the five-scenario CIs.
# [Expected] nested-CV outer-loop AUC 0.93-0.98 (a value below 0.982 is the
#   normal correction; 0.982 is the internal optimistic estimate — both are
#   reported in the paper); stability: HLTF/ARGLU1/LSP1/HMGN2 hubs expected
#   >=80%; weak genes may drop to 30-60% — report as observed.
# [Input] output/GSE57338_processed.rds + GSE57338_HFvsNF_all.csv
#   + the LRG gene list (data/lrg332.txt) + output/ROC_validation_summary_v2.csv (from M2)
# [Runtime] ~3-6 min (the 200 cv.glmnet fits dominate)
# [Reference] Varma S, Simon R. BMC Bioinformatics 2006;7:91 (PMID 16504092)
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# Set DIR to your working directory (forward slashes on Windows)
DIR <- "./GEO_data"
OUT <- file.path(DIR, "output")

LOG <- file.path(OUT, "M9_log.txt")
logline <- function(...) {
  msg <- paste0(..., "")
  cat(msg, "\n"); cat(msg, "\n", file = LOG, append = TRUE)
}
if (file.exists(LOG)) file.remove(LOG)
logline("========== M9: nested CV + bootstrap | start: ", format(Sys.time()), "==========")

for (p in c("glmnet", "pROC")) {
  if (!requireNamespace(p, quietly = TRUE)) install.packages(p)
}
suppressPackageStartupMessages({ library(glmnet); library(pROC) })

# ---- 1. Data loading (reproduces the M2 logic) ----
need <- c(rds = file.path(OUT, "GSE57338_processed.rds"),
          tt  = file.path(OUT, "GSE57338_HFvsNF_all.csv"))
miss <- need[!file.exists(need)]
if (length(miss) > 0) stop(paste("missing prerequisite file(s):", paste(miss, collapse="\n   ")))
d38 <- readRDS(need[["rds"]])
tt38 <- read.csv(need[["tt"]], check.names = FALSE)

lrg_file <- file.path(DIR, "lrg332.txt")
if (!file.exists(lrg_file)) lrg_file <- list.files(DIR, pattern = "lrg332\\.txt$", recursive = TRUE, full.names = TRUE)[1]
if (is.na(lrg_file)) stop("lrg332.txt (LRG gene list) not found")
lrg <- toupper(trimws(readLines(lrg_file, warn = FALSE)))
lrg <- lrg[lrg != ""]

LFC <- 0.2
gidcol <- if (".gid" %in% names(tt38)) ".gid" else names(tt38)[ncol(tt38)]
sig38 <- subset(tt38, adj.P.Val < 0.05 & abs(logFC) > LFC)
genes56 <- as.character(sig38[[gidcol]])
genes56 <- genes56[toupper(genes56) %in% lrg]
genes56 <- intersect(genes56, rownames(d38$expr))
logline(sprintf("[data] 56-DLRG recomputation = %d genes | samples = %d", length(genes56), ncol(d38$expr)))
if (length(genes56) != 56) logline("  note: DLRG count != 56; continuing")

x38 <- t(d38$expr[genes56, , drop = FALSE])
y38 <- as.integer(d38$group == "HF")

# ---- 2. Nested cross-validation (10 outer folds) ----
logline("[Part A] nested CV (10 outer folds x inner cv.glmnet)...")
set.seed(2026)
folds <- sample(rep(1:10, length.out = nrow(x38)))
pred_oof <- rep(NA_real_, nrow(x38))
for (k in 1:10) {
  te <- folds == k
  tr <- !te
  cvfit <- cv.glmnet(x38[tr, ], y38[tr], family = "binomial", alpha = 1,
                     nfolds = 10, type.measure = "auc")
  lam <- cvfit$lambda.1se
  if (length(coef(cvfit, s = "lambda.1se")@x) <= 1) lam <- cvfit$lambda.min
  p <- as.numeric(predict(cvfit, x38[te, , drop = FALSE], s = lam, type = "response"))
  pred_oof[te] <- p
}
roc_oof <- roc(y38, pred_oof, quiet = TRUE)
ci_oof <- as.numeric(ci.auc(roc_oof))
logline(sprintf("[nested CV] pooled out-of-fold AUC=%.3f [95%%CI %.3f-%.3f]", as.numeric(auc(roc_oof)), ci_oof[1], ci_oof[3]))
logline("  [reference] internal training AUC=0.982 (optimistic: same data for selection + evaluation)")

# ---- 3. Bootstrap stability (200 replicates) ----
logline("[Part B] bootstrap x200 (gene re-selection each time; several minutes)...")
set.seed(2026)
B <- 200
sel_count <- setNames(rep(0L, length(genes56)), genes56)
for (b in 1:B) {
  idx <- sample(nrow(x38), replace = TRUE)
  cvb <- cv.glmnet(x38[idx, ], y38[idx], family = "binomial", alpha = 1,
                   nfolds = 10, type.measure = "auc")
  co <- as.matrix(coef(cvb, s = "lambda.1se"))
  selb <- setdiff(rownames(co)[co[, 1] != 0], "(Intercept)")
  if (length(selb) == 0) {
    co <- as.matrix(coef(cvb, s = "lambda.min"))
    selb <- setdiff(rownames(co)[co[, 1] != 0], "(Intercept)")
  }
  sel_count[selb] <- sel_count[selb] + 1L
  if (b %% 40 == 0) logline(sprintf("   progress %d/200", b))
}
stab <- sort(sel_count / B, decreasing = TRUE)
stab_df <- data.frame(gene = names(stab), freq = round(as.numeric(stab), 3))
write.csv(stab_df, file.path(OUT, "M9_bootstrap_stability.csv"), row.names = FALSE)
logline("[stability top 10]")
for (i in 1:10) logline(sprintf("   %-8s %.0f%%", stab_df$gene[i], stab_df$freq[i] * 100))
n80 <- sum(stab >= 0.8)
logline(sprintf("[summary] genes with >=80%% stability: %d / 13-gene baseline", n80))

# ---- 4. DeLong CI summary (already available from M2) ----
logline("[Part C] five-scenario DeLong CIs (from M2 ROC_validation_summary_v2.csv)")
v2 <- file.path(OUT, "ROC_validation_summary_v2.csv")
if (file.exists(v2)) {
  rv <- read.csv(v2, check.names = FALSE)
  print(rv)
  for (i in seq_len(nrow(rv))) {
    logline(sprintf("   %s: AUC=%.3f [95%%CI %.3f-%.3f]",
                    rv$scenario[i], rv$AUC[i], rv$AUC_low[i], rv$AUC_high[i]))
  }
} else {
  logline("  !! ROC_validation_summary_v2.csv not found — rerun M2 first")
}

# ---- 5. Self-verification ----
logline("========== self-verification ==========")
logline(sprintf("nestedCV_AUC=%.3f [%.3f-%.3f]", as.numeric(auc(roc_oof)), ci_oof[1], ci_oof[3]))
logline(sprintf("stability_ge80 = %d genes", n80))
logline(paste0("top_stability = ", paste(sprintf("%s(%.0f%%)", stab_df$gene[1:5], stab_df$freq[1:5]*100), collapse=", ")))
logline("========== M9 done ==========")

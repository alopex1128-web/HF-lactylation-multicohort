# ============================================================
# M8_meff_correction.R
# Purpose: effective number of tests correction for inter-gene correlation
# Paper section: Table 3 note
# Pipeline position: step 9 of the main chain
# Prerequisite: M6 completed (P values are embedded; only the M1 expression matrix is needed)
# (Generated from the modular pipeline consolidation; analysis logic identical
#  to the verified version.)
# ============================================================
# ============================================================
# M8: Meff correction for the 26 MR-tested genes (inter-gene-correlation-
# adjusted multiple testing)
# ============================================================
# [Background] The 26 MR-tested genes are co-expressed and correlated, so a
#   26-test BH/Bonferroni correction is conservative. This module computes
#   the effective number of tests (Meff) via an eigenvalue approach and
#   recalculates adjusted q values.
# [Logic] (1) extract expression of the 26 testable genes from
#   GSE57338_processed.rds; (2) Pearson correlation matrix -> eigendecomposition;
#   (3) two Meff estimators: Nyholt 2004 [1+(M-1)(1-Var(lambda)/M)];
#   Li-Ji 2005 [sum(floor(lambda)+I(frac>=0.5))] — Li-Ji is the primary report;
#   (4) the 26 P values are the authoritative values of Table 3 (verified in M6,
#   embedded here); (5) BH-m=Meff q values + Bonferroni-Meff threshold.
# [Expected] Meff(Li-Ji) 18-24 (strong co-expression among 26 genes compresses
#   the effective count); the minimum adjusted q remains >0.05 (the global null
#   is unchanged; only the q numbers shrink). If any q<0.05 appears, that is a
#   genuine result.
# [Input] output/GSE57338_processed.rds (from M1)
# [Runtime] ~30 s
# [References] Li J, Ji L. Heredity 2005;95:221-227 (PMID 16077740);
#   Nyholt DR. Am J Hum Genet 2004;74:765-769 (PMID 14997420)
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# Set DIR to your working directory (forward slashes on Windows)
DIR <- "./GEO_data"
OUT <- file.path(DIR, "output")

LOG <- file.path(OUT, "M8_meff_log.txt")
logline <- function(...) {
  msg <- paste0(..., "")
  cat(msg, "\n"); cat(msg, "\n", file = LOG, append = TRUE)
}
if (file.exists(LOG)) file.remove(LOG)
logline("========== M8: Meff correction | start: ", format(Sys.time()), "==========")

if (!requireNamespace("pROC", quietly = TRUE)) install.packages("pROC")
suppressPackageStartupMessages(library(pROC))

# ---- 1. The 26 testable genes + authoritative Table 3 P values ----
genes26 <- c("SOD1","H2AZ1","ACAT2","LDHB","TP53","HLTF","MBD2","RB1","PFKM",
             "GATA2","HMGN1","TPR","WARS2","STAT2","NSUN2","ARGLU1","PPP1CC",
             "LSP1","LRPPRC","ZNFX1","DDX17","HMGN2","JAK2","EIF4G1","JMJD1C","STAT4")
P26 <- c(0.012,0.014,0.017,0.025,0.055,0.059,0.079,0.137,0.162,0.203,
         0.250,0.463,0.507,0.566,0.587,0.738,0.775,0.819,0.884,0.949,
         0.959,0.970,0.991,0.025,0.059,0.591)
stopifnot(length(genes26) == 26, length(P26) == 26)
logline("[input] 26 genes + P values embedded (Table 3 values, verified in M6)")

# ---- 2. Expression matrix ----
rds_path <- file.path(OUT, "GSE57338_processed.rds")
if (!file.exists(rds_path)) stop("GSE57338_processed.rds missing — run M1 first")
d38 <- readRDS(rds_path)
expr <- d38$expr
miss <- setdiff(genes26, rownames(expr))
if (length(miss) > 0) stop(sprintf("%d gene(s) not in expr rownames: %s", length(miss), paste(miss, collapse=",")))
X <- expr[genes26, , drop = FALSE]
logline(sprintf("[data] expression matrix %d genes x %d samples", nrow(X), ncol(X)))

# ---- 3. Meff computation ----
R <- cor(t(X), method = "pearson")
eig <- eigen(R, symmetric = TRUE)$values
M <- length(genes26)

# Nyholt 2004
meff_nyholt <- 1 + (M - 1) * (1 - var(eig) / M)
# Li-Ji 2005 (primary report)
fl <- floor(eig); frac <- eig - fl
meff_liji <- sum(fl + as.integer(frac >= 0.5))
# reference upper bound (all eigenvalues kept)
logline(sprintf("[Meff] Nyholt=%.1f | Li-Ji=%.1f (26 nominal tests)", meff_nyholt, meff_liji))

# ---- 4. Recalculated adjustment ----
Meff <- round(meff_liji)
q_bh26 <- p.adjust(P26, method = "BH")          # reproduces Table 3
q_bhMeff <- p.adjust(P26, method = "BH") * 26 / Meff  # BH denominator replaced by Meff
bonf_Meff <- pmin(1, P26 * Meff)
res <- data.frame(gene = genes26, P = P26,
                  q_BH26 = round(q_bh26, 3),
                  q_BH_Meff = round(q_bhMeff, 3),
                  P_Bonf_Meff = round(bonf_Meff, 3))
write.csv(res, file.path(OUT, "M8_meff_results.csv"), row.names = FALSE)

sig <- res[res$q_BH_Meff < 0.05, ]
logline(sprintf("[result] Meff(Li-Ji)=%d | min q_BH_Meff=%.3f | genes with q<0.05: %d",
                Meff, min(res$q_BH_Meff), nrow(sig)))
if (nrow(sig) > 0) {
  logline("[note] gene(s) below the Meff-adjusted threshold — genuine result:")
  for (i in seq_len(nrow(sig))) logline(sprintf("   %s P=%.3f q_Meff=%.3f", sig$gene[i], sig$P[i], sig$q_BH_Meff[i]))
} else {
  logline("[conclusion] global null unchanged — no significant gene after correction; the five suggestive hits only have smaller q")
}

# ---- 5. Self-verification ----
logline("========== self-verification ==========")
logline(sprintf("Meff_Nyholt=%.1f | Meff_LiJi=%d", meff_nyholt, Meff))
logline(sprintf("min_q_BH_Meff=%.3f | sig_count=%d", min(res$q_BH_Meff), nrow(sig)))
logline("========== M8 done ==========")

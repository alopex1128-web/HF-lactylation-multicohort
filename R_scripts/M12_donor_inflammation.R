# M12 · Quantifying donor-inflammation stratification of immune lineage scores
# Added in v9.2 (response to pre-submission revision R3).
# Output = Supplementary Table S9 + underlying data for Supplementary Figure S1
# Inputs:
#   1) expr_sym.csv           — M1 preprocessed output (rows = symbols, cols = 313 samples, GSE57338)
#   2) mcp_lineage_scores.csv — MCP-counter 1.2.0 lineage scores from M4/M5 (rows = 10 lineages, cols = 313 samples)
# Design: the 136 NF donors are split by the median of a composite inflammatory
#         score (mean of within-donor z-scores of IL1B/TNF/NFKB1) into high- vs
#         low-inflammation groups (68 each); (i) high vs low inflammation;
#         (ii) conservative sensitivity of HF vs low-inflammation donors.
# Expected values (from the reference run; see manuscript Suppl. Table S9):
#   Neutrophils high-vs-low q=0.005 (myeloid elevation); all lymphoid lineages q>0.14
#   HF vs low-inflammation: T cells/Cyto q<1e-4, Endothelial q=0.049,
#   Fibroblasts q=0.009 (robust); Monocytic q=0.218 and Neutrophils q=0.580 no longer significant
suppressMessages({library(data.table)})

expr <- fread('expr_sym.csv')
mcp  <- fread('mcp_lineage_scores.csv')
mat  <- as.matrix(expr[, -1]); rownames(mat) <- expr[[1]]
mcpM <- as.matrix(mcp[, -1]); rownames(mcpM) <- mcp[[1]]

# Group assignment: same source as M5_immune_sensitivity.R (group vector from the
# M1 preprocessing output, 177 HF / 136 NF). If your expr_sym.csv column names do
# not carry group labels, provide group.csv (two columns: sample,group) instead.
if (file.exists('group.csv')) {
  gp <- read.csv('group.csv', stringsAsFactors = FALSE)
  hf_cols <- colnames(mat) %in% gp$sample[gp$group == 'HF']
} else {
  hf_cols <- grepl('^HF|^GSM.*HF', colnames(mat))
}
nf_cols  <- !hf_cols
stopifnot(sum(hf_cols) == 177, sum(nf_cols) == 136)

# ---- Donor composite inflammatory score (mean of within-donor z-scores) ----
zscore <- function(x) (x - mean(x)) / sd(x)
infl <- rep(NA_real_, sum(nf_cols))
nf_names <- colnames(mat)[nf_cols]
for (g in c('IL1B','TNF','NFKB1')) {
  v <- mat[g, nf_names]
  infl <- infl + zscore(v)        # accumulate within-donor standardized values, then divide by 3
}
infl <- infl / 3
names(infl) <- nf_names
hi <- names(infl)[infl > median(infl)]
lo <- names(infl)[infl <= median(infl)]
stopifnot(length(hi) == 68, length(lo) == 68)

# ---- (i) High- vs low-inflammation donors (10 lineages, Wilcoxon + BH) ----
p_hi_lo <- apply(mcpM, 1, function(r) wilcox.test(r[hi], r[lo])$p.value)
q_hi_lo <- p.adjust(p_hi_lo, 'BH')

# ---- (ii) HF vs low-inflammation donors (conservative sensitivity) ----
hf_names <- colnames(mat)[hf_cols]
p_hf_lo  <- apply(mcpM, 1, function(r) wilcox.test(r[hf_names], r[lo])$p.value)
q_hf_lo  <- p.adjust(p_hf_lo, 'BH')

out <- data.frame(lineage = rownames(mcpM),
                  median_lo = apply(mcpM[, lo], 1, median),
                  median_hi = apply(mcpM[, hi], 1, median),
                  median_HF = apply(mcpM[, hf_names], 1, median),
                  P_hi_vs_lo = p_hi_lo, q_hi_vs_lo = q_hi_lo,
                  P_HF_vs_lo = p_hf_lo,  q_HF_vs_lo  = q_hf_lo)
write.csv(out, 'S9_donor_inflammation.csv', row.names = FALSE)

# ---- Self-verification ----
cat('Neutrophils q(hi vs lo):', sprintf('%.4f', q_hi_lo['Neutrophils']), '(expected 0.005)\n')
cat('T cells q(HF vs lo):', sprintf('%.2e', q_hf_lo['T cells']), '(expected <1e-4)\n')
cat('Endothelial q(HF vs lo):', sprintf('%.4f', q_hf_lo['Endothelial cells']), '(expected 0.0487)\n')
cat('Fibroblasts q(HF vs lo):', sprintf('%.4f', q_hf_lo['Fibroblasts']), '(expected 0.0090)\n')
cat('Monocytic q(HF vs lo):', sprintf('%.3f', q_hf_lo['Monocytic lineage']), '(expected 0.218 — no longer significant)\n')
cat('Neutrophils q(HF vs lo):', sprintf('%.3f', q_hf_lo['Neutrophils']), '(expected 0.580 — no longer significant)\n')
cat('Lymphoid lineages q(hi vs lo) all >0.14:', all(q_hi_lo[c('T cells','CD8 T cells','Cytotoxic lymphocytes','B lineage','NK cells')] > 0.14), '\n')

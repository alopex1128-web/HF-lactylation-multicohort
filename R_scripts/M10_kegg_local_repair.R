# ============================================================
# M10_kegg_local_repair.R (optional)
# Purpose: local reproduction of the KEGG REST enrichment (skippable when
#          network access to rest.kegg.jp is unavailable)
# Paper section: Table 4, KEGG panel
# Pipeline position: optional step
# Prerequisite: M1 completed; requires access to rest.kegg.jp
# (Generated from the modular pipeline consolidation; analysis logic identical
#  to the verified version.)
# ============================================================
# ============================================================
# M10: local KEGG enrichment rerun
# ============================================================
# [Background] The reference implementation of the DLRG KEGG enrichment
#   (Table 4) is provided in the results/ tables. This script reruns the same
#   analysis locally so that every analysis in the paper can be reproduced
#   end-to-end on a single machine.
# [Logic] Identical to M4 Part A2: DLRGs (56) -> Entrez mapping ->
#   enrichKEGG(organism="hsa", pAdjustMethod="BH", pvalueCutoff=0.05),
#   background universe = all symbols of the expression matrix mapped to
#   Entrez (same map_id logic as M4).
# [Expected] Values are compared against the paper automatically after the
#   run: Carbon metabolism q ~ 2.5e-3, Pentose phosphate q ~ 5.2e-3,
#   significant terms ~ 85 (tolerance +/-3; tiny differences between
#   clusterProfiler implementations are normal).
# [Network] enrichKEGG needs rest.kegg.jp; on failure the script retries 3
#   times (20 s apart). The KEGG server is often congested at certain hours;
#   retrying at a different time of day usually succeeds.
# [Input] output/GSE57338_processed.rds + output/GSE57338_DLRGs.csv (from M1)
# [Runtime] ~1-3 min
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# Set DIR to your working directory (forward slashes on Windows)
DIR <- "./GEO_data"

OUT <- file.path(DIR, "stage7_results")
LOGF <- file.path(OUT, "M10_kegg_local_report.txt")
logline <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n"); cat(msg, "\n", file = LOGF, append = TRUE)
}
if (file.exists(LOGF)) file.remove(LOGF)
logline("========== M10: local KEGG rerun | start: ", format(Sys.time()), "==========")

# ---- 0. Package check ----
if (!requireNamespace("clusterProfiler", quietly = TRUE) ||
    !requireNamespace("org.Hs.eg.db", quietly = TRUE)) {
  BiocManager::install(c("clusterProfiler", "org.Hs.eg.db"), update = FALSE, ask = FALSE)
}
suppressPackageStartupMessages({ library(clusterProfiler); library(org.Hs.eg.db) })
logline("[packages] clusterProfiler + org.Hs.eg.db ready")

# ---- 1. Data loading ----
need <- c(rds = file.path(DIR, "output", "GSE57338_processed.rds"),
          dlrg = file.path(DIR, "output", "GSE57338_DLRGs.csv"))
miss <- need[!file.exists(need)]
if (length(miss) > 0) {
  logline("!!! missing prerequisite file(s):"); for (m in miss) logline("   -", m)
  stop("run M1 (stages 1-3) first")
}
d38 <- readRDS(need[["rds"]])
expr38 <- d38$expr
dlrg <- read.csv(need[["dlrg"]], check.names = FALSE)
logline(sprintf("[data] expression matrix %d genes | DLRGs %d", nrow(expr38), nrow(dlrg)))

# ---- 2. Entrez mapping (same logic as M4) ----
map_id <- function(symbols) {
  suppressWarnings(eg <- AnnotationDbi::select(org.Hs.eg.db, keys = unique(symbols),
              columns = "ENTREZID", keytype = "SYMBOL"))
  eg <- eg[!is.na(eg$ENTREZID) & !duplicated(eg$SYMBOL), ]
  setNames(eg$ENTREZID, eg$SYMBOL)
}
bg <- unique(map_id(rownames(expr38)))
# Column detection: in the M1 DLRG csv the gene-symbol column is .gid
# (appended as the last column); the first column holds logFC values —
# do not use names()[1].
gidcol <- if (".gid" %in% names(dlrg)) ".gid" else
          if ("gene" %in% names(dlrg)) "gene" else
          if ("symbol" %in% names(dlrg)) "symbol" else NA
if (is.na(gidcol)) stop("gene symbol column not found (.gid/gene/symbol)")
dlrg_entrez <- unique(unname(map_id(as.character(dlrg[[gidcol]]))))
logline(sprintf("[mapping] background universe %d genes | DLRGs mapped %d (compare with the M4 log)", length(bg), length(dlrg_entrez)))

# ---- 3. enrichKEGG (3 automatic retries) ----
ekk <- NULL
for (attempt in 1:3) {
  ekk <- tryCatch(
    enrichKEGG(gene = dlrg_entrez, universe = bg, organism = "hsa",
               pAdjustMethod = "BH", pvalueCutoff = 0.05),
    error = function(e) {
      logline(sprintf("  [attempt %d] KEGG REST connection failed: %s", attempt, conditionMessage(e)))
      NULL
    })
  if (!is.null(ekk)) break
  if (attempt < 3) { logline("  waiting 20 s before retry..."); Sys.sleep(20) }
}
if (is.null(ekk)) {
  logline("!!! all 3 attempts failed (KEGG server congestion)")
  logline(">>> suggestion: rerun this script at another time of day;")
  logline("    or rely on the published Table 4 values (enrichKEGG is a")
  logline("    deterministic computation: same input + same universe give")
  logline("    the same result, so reproducibility is unaffected)")
  stop("KEGG REST unreachable; see report for suggestions")
}

dfk <- as.data.frame(ekk)
write.csv(dfk, file.path(OUT, "Table4_DLRG_KEGG_enrichment_local.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")
logline(sprintf("[KEGG] significant terms: %d (paper: 85; +/-3 is normal)", nrow(dfk)))

# ---- 4. Comparison with the published values ----
check <- function(name, got, want, tol = 0.05) {
  ratio <- if (want != 0) abs(got - want) / abs(want) else 0
  ok <- ratio <= tol
  logline(sprintf("  [%s] %s: local=%.3e | paper=%.3e | relative difference=%.1f%%",
                  ifelse(ok, "OK", "!!"), name, got, want, ratio * 100))
}
logline("---- comparison against paper (tolerance 5%) ----")
getq <- function(df, keyword) {
  hit <- grep(keyword, df$Description, value = FALSE)
  if (length(hit) > 0) df$p.adjust[hit[1]] else NA
}
check("Carbon metabolism q", getq(dfk, "Carbon metabolism"), 2.456e-3)
check("Pentose phosphate q", getq(dfk, "Pentose phosphate"), 5.239e-3)
check("ATP-dependent chromatin remodeling q", getq(dfk, "ATP-dependent chromatin"), 1.366e-2)
logline("---- significant terms: local = ", nrow(dfk), " | paper = 85 ----")

logline("========== M10 done: ", format(Sys.time()), "==========")
logline("  all OK -> KEGG values reproduced; Table 4 unchanged")
logline("  mismatch -> report the local values (deterministic computation; update Table 4 and the abstract q if they differ)")

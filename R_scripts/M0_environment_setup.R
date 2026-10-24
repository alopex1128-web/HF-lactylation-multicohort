# ============================================================
# M0_environment_setup.R  |  Module 0
# Purpose: install all R package dependencies for the pipeline
# Paper section: —
# Pipeline position: step 1 of the main chain
# Usage: run once on a fresh R installation (already-installed
#        packages are skipped automatically)
# ============================================================
# ==============================================================================
# Installs every package required by stages 1-3 plus the prerequisites of
# stage 4, from a clean R environment.
# ------------------------------------------------------------------------------
# [Prerequisites] R (>= 4.x) installed; no packages required beforehand.
# [Mirrors] CRAN and Bioconductor both point to the Tsinghua mirror — replace
#           with a mirror close to you (or the defaults) if preferred.
# [Runtime] ~20-40 min for a full install (the ~100 MB annotation packages
#           dominate); <1 min when everything is already installed.
# [Usage] Open in RStudio -> select all -> run (no paths to edit).
# [Next] Open M1_preprocessing_limma.R -> set DIR -> run.
# ==============================================================================

# ------------------------------------------------------------------
# 1. Mirrors (CRAN + Bioconductor)
# ------------------------------------------------------------------
options(repos = c(CRAN = "https://mirrors.tuna.tsinghua.edu.cn/CRAN/"))
options(BioC_mirror = "https://mirrors.tuna.tsinghua.edu.cn/bioconductor")
cat("Mirrors set: CRAN Tsinghua / Bioconductor Tsinghua\n")

# ------------------------------------------------------------------
# 2. Package list
# ------------------------------------------------------------------
cran_pkgs <- c(
  "BiocManager",   # bridge for Bioconductor installs
  "data.table",    # fast reading of large files (125 MB GSE57338 matrix)
  "ggplot2",       # volcano / PCA plotting engine
  "ggrepel",       # non-overlapping gene labels on volcano plots
  "pheatmap",      # core-DEG heatmaps
  "VennDiagram",   # three-dataset Venn intersection
  "glmnet",        # [stage 4 prerequisite] LASSO feature selection
  "pROC"           # [stage 4 prerequisite] ROC curves / AUC
)
bioc_pkgs <- c(
  "limma",                         # differential expression core
  "sva",                           # surrogate-variable batch correction
  "AnnotationDbi",                 # probe annotation engine
  "hugene11sttranscriptcluster.db",# GSE57338 annotation (Human Gene 1.1 ST, not 2.0)
  "hgu133a.db"                     # GSE5406 annotation (HG-U133A)
)
cat(sprintf("Planned install: %d CRAN + %d Bioconductor packages\n\n", length(cran_pkgs), length(bioc_pkgs)))

# ------------------------------------------------------------------
# 3. CRAN installs (only missing packages)
# ------------------------------------------------------------------
cat("========== CRAN packages ==========\n")
fail_cran <- c()
for (p in cran_pkgs) {
  if (requireNamespace(p, quietly = TRUE)) {
    cat(sprintf("  [skip] %-32s installed (v%s)\n", p, as.character(packageVersion(p))))
    next
  }
  cat(sprintf("  [install] %s ...\n", p))
  ok <- tryCatch({ install.packages(p); requireNamespace(p, quietly = TRUE) },
                 error = function(e) FALSE)
  if (ok) cat(sprintf("  [ok] %s v%s\n", p, as.character(packageVersion(p)))) else {
    cat(sprintf("  [FAIL] %s\n", p)); fail_cran <- c(fail_cran, p)
  }
}

# ------------------------------------------------------------------
# 4. Bioconductor installs (only missing packages)
# ------------------------------------------------------------------
cat("\n========== Bioconductor packages ==========\n")
if (!requireNamespace("BiocManager", quietly = TRUE)) {
  cat("  BiocManager missing, installing first...\n"); install.packages("BiocManager")
}
fail_bioc <- c()
todo <- bioc_pkgs[!sapply(bioc_pkgs, requireNamespace, quietly = TRUE)]
if (length(todo) == 0) {
  cat("  all installed, skipping\n")
} else {
  cat(sprintf("  to install: %s\n", paste(todo, collapse = ", ")))
  cat("  (annotation packages are ~100 MB; this may take 5-15 minutes)\n")
  ok <- tryCatch({ BiocManager::install(todo, ask = FALSE, update = FALSE); TRUE },
                 error = function(e) { cat("  install error:", conditionMessage(e), "\n"); FALSE })
  for (p in bioc_pkgs) {
    if (requireNamespace(p, quietly = TRUE)) {
      cat(sprintf("  [ready] %-32s v%s\n", p, as.character(packageVersion(p))))
    } else {
      cat(sprintf("  [FAIL] %s\n", p)); fail_bioc <- c(fail_bioc, p)
    }
  }
}

# ------------------------------------------------------------------
# 5. Annotation spot checks (actually call select() once, not just check load)
# ------------------------------------------------------------------
cat("\n========== Annotation package spot checks ==========\n")
tryCatch({
  suppressPackageStartupMessages(library(hugene11sttranscriptcluster.db))
  m1 <- AnnotationDbi::select(hugene11sttranscriptcluster.db,
        keys = "7892501", columns = "SYMBOL", keytype = "PROBEID")
  cat(sprintf("  hugene11st probe 7892501 -> %s  (expected TSPAN6)\n", m1$SYMBOL[1]))
}, error = function(e) cat("  hugene11st spot check failed:", conditionMessage(e), "\n"))
tryCatch({
  suppressPackageStartupMessages(library(hgu133a.db))
  m2 <- AnnotationDbi::select(hgu133a.db,
        keys = "1007_s_at", columns = "SYMBOL", keytype = "PROBEID")
  cat(sprintf("  hgu133a    probe 1007_s_at -> %s  (expected DDR1)\n", m2$SYMBOL[1]))
}, error = function(e) cat("  hgu133a spot check failed:", conditionMessage(e), "\n"))

# ------------------------------------------------------------------
# 6. Summary
# ------------------------------------------------------------------
cat("\n==================================================\n")
if (length(c(fail_cran, fail_bioc)) == 0) {
  cat("  All packages ready. Next:\n")
  cat("    open M1_preprocessing_limma.R\n")
  cat("    -> set the DIR variable at the top of each script (forward slashes)\n")
  cat("    -> select all and run (~15-25 min)\n")
} else {
  cat("  Some packages failed to install:", paste(c(fail_cran, fail_bioc), collapse = ", "), "\n")
  cat("  Review the console output above for error details\n")
}
cat("==================================================\n")

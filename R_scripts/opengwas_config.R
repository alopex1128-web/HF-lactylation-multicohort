# ============================================================
# OpenGWAS API configuration (MR modules only)
# Usage: add  source("opengwas_config.R")  at the top of any MR script.
# Updated 2026-08-23: compatible with ieugwasr 2.x (set_jwt_token is
# deprecated; the new version reads the OPENGWAS_JWT environment
# variable), while 1.x still uses set_jwt_token.
# ============================================================

# Your JWT token. Obtain one by registering at https://api.opengwas.io
# (free) and paste it below; tokens expire periodically and can be
# refreshed by logging in again.
OPENGWAS_TOKEN <- "PASTE_YOUR_NEW_TOKEN_HERE"

# Token setup (dual compatibility: ieugwasr 2.x via environment variable, 1.x via set_jwt_token)
Sys.setenv(OPENGWAS_JWT = OPENGWAS_TOKEN)
if ("set_jwt_token" %in% getNamespaceExports("ieugwasr")) {
  suppressWarnings(ieugwasr::set_jwt_token(OPENGWAS_TOKEN))
  cat("[OpenGWAS config] ieugwasr 1.x path: set_jwt_token() applied\n")
} else {
  cat("[OpenGWAS config] ieugwasr 2.x path: OPENGWAS_JWT environment variable set\n")
}

# ============================================================
# Key datasets for the MR layer (verified 2026-08-21 / re-verified
# per gene 2026-08-23)
# ============================================================
# Outcome — HERMES heart-failure GWAS:
#   ID: ebi-a-GCST009541
#   Shah S et al. 2020 Nat Commun (PMID: 31919418)
#   47,309 cases / 930,014 controls
#
# Exposure — eQTLGen cis-eQTL (blood tissue, 31,684 samples):
#   Per-gene dataset format: eqtl-a-ENSG00000xxxxx
#   Ensembl IDs for the 22 candidate genes verified individually.
#   Note three mapping pitfalls (handled in M6):
#     1. LSP1 != ENSG00000269099 (that is the pseudogene LSP1P1); the true ID is ENSG00000130592
#     2. BRD4  != ENSG00000143294 (that is PRCC); the true ID is ENSG00000141867
#     3. eqtl-a- datasets contain trans-eQTL signals; apply a cis filter
#        of gene region +/- 1 Mb (built into M6 as cis_filter)
#
# Quick test (standalone):
#   library(ieugwasr); source("opengwas_config.R")
#   associations(variants = "rs12143842", id = "ebi-a-GCST009541")
#   Returning association statistics for that SNP means the configuration works
#   (ieugwasr 1.x uses: associations(rsid = "rs12143842", id = "ebi-a-GCST009541"))
# ============================================================

cat("[OpenGWAS config done] token in use (replace when expired)\n")

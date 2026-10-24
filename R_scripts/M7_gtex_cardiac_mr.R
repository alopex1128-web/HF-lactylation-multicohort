# ============================================================
# M7_gtex_cardiac_mr.R  |  M7
# Purpose: three-gene cardiac cis-eQTL sensitivity MR (requires the R3 data bundle)
# Paper section: Suppl. Table S7
# Pipeline position: step 8 of the main chain
# Prerequisite: the R3_local data bundle (eQTL/HERMES/1000G packaged)
# (Generated from the modular pipeline consolidation; analysis logic identical to the verified version.)
# ============================================================
# ============================================================
# R3 local reproduction: GTEx cardiac-eQTL MR sensitivity analysis (runs end-to-end locally)
# ============================================================
# [Background] the R3 sensitivity analysis (GTEx cardiac eQTL x HERMES).
#   This script reimplements the same logic in R so the analysis can be reproduced
#   end-to-end on a local R installation; the local numbers are authoritative.
#
# [Workflow (maps one-to-one to the R3 paragraph of the Methods)]
#   1) extract SNP-type instruments at P<5e-8 for three genes from the GTEx left-ventricle eQTL
#   2) bridge rsID -> GRCh37 coordinates (dbSNP mapping table)
#   3) match the 1000G EUR panel by coordinate + alleles
#   4) PLINK v1.90 clump (r2<0.1, 2500 kb, p1=5e-8)
#   5) match the HERMES GWAS by coordinate + allele set
#   6) harmonize (allele alignment + palindromic frequency check)
#   7) Wald ratio (1 IV) / fixed-effect IVW (>=2 IVs) + BH
#
# [Prerequisites]
#   (1) place the R3_local folder (1000G/plink/HERMES data) in the working directory
#   (2) download PLINK v1.90 (https://www.cog-genomics.com/plink/1.9/
#       -> Windows 64-bit -> extract plink.exe -> put it in R3_local)
#   (3) the folder should already contain: hits_heart_lv.tsv / dbsnp_map.tsv /
#       hermes_window.tsv / blood_reference.csv /
#      plink_WARS2.bed|.bim|.fam / plink_LRPPRC.* / plink_NSUN2.*
#
# [Expected values] checked automatically step by step; all OK = successful reproduction
# [Runtime] ~2-4 min (read.delim on large files dominates)
# [Reproducibility] fully deterministic (no random numbers); PLINK is cross-platform consistent
# ============================================================

rm(list = ls())
options(stringsAsFactors = FALSE)

# Set DIR to your local R3_local folder (forward slashes)
DIR <- "./R3_local"
setwd(DIR)

LOGF <- file.path(DIR, "R3_local_report.txt")
logline <- function(...) {
  msg <- paste0(...)
  cat(msg, "\n"); cat(msg, "\n", file = LOGF, append = TRUE)
}
if (file.exists(LOGF)) file.remove(LOGF)
logline("========== R3 local reproduction: GTEx cardiac-eQTL MR, started: ", format(Sys.time()), "==========")

check <- function(label, got, want) {
  ok <- isTRUE(all.equal(got, want)) || identical(as.character(got), as.character(want))
  mark <- ifelse(ok, "OK", "!! mismatch")
  logline(sprintf("  [%s] %s: local=%s | expected=%s", mark, label,
                  paste(got, collapse=","), paste(want, collapse=",")))
  invisible(ok)
}

GENES <- c(WARS2 = "ENSG00000116874", LRPPRC = "ENSG00000138095", NSUN2 = "ENSG00000037474")
GENE_CHR <- c(WARS2 = "1", LRPPRC = "2", NSUN2 = "5")

# ---- 0. PLINK check ----
PLINK <- NULL
if (file.exists(file.path(DIR, "plink.exe"))) {
  PLINK <- file.path(DIR, "plink.exe")
} else if (file.exists(file.path(DIR, "plink"))) {
  PLINK <- file.path(DIR, "plink")
}
if (is.null(PLINK)) {
  logline("  [error] plink.exe not found! Download from https://www.cog-genomics.com/plink/1.9/")
  logline("        (Windows 64-bit), extract plink.exe into this folder, and rerun.")
  stop("plink.exe missing")
}
logline("  [0] PLINK ready:", PLINK)

# ---- 1. extract P<5e-8 SNP instruments for the three genes ----
logline("---- 1/7 extracting GTEx left-ventricle eQTL instruments (P<5e-8, type=SNP) ----")
logline("  reading hits_heart_lv.tsv (236,063 rows, ~30-60 s)...")
hits <- read.delim("hits_heart_lv.tsv", colClasses = c(
  variant="NULL", r2="NULL", pvalue="numeric", molecular_trait_object_id="NULL",
  molecular_trait_id="NULL", maf="NULL", gene_id="character", median_tpm="NULL",
  beta="numeric", se="numeric", an="integer", ac="integer",
  chromosome="character", position="integer", ref="character", alt="character",
  type="character", rsid="character"))
logline("  eQTL rows loaded:", nrow(hits), "rows x", ncol(hits), "cols")

hits <- hits[hits$gene_id %in% GENES, ]
instr <- hits[hits$pvalue < 5e-8 & hits$type == "SNP", ]
n_by_gene <- table(factor(sapply(instr$gene_id, function(g) names(GENES)[GENES==g]), levels=names(GENES)))
logline("  SNP instrument counts: WARS2=", n_by_gene["WARS2"], " LRPPRC=", n_by_gene["LRPPRC"],
        " NSUN2=", n_by_gene["NSUN2"])
check("instrument count WARS2", as.integer(n_by_gene["WARS2"]), 252)
check("instrument count LRPPRC", as.integer(n_by_gene["LRPPRC"]), 219)
check("instrument count NSUN2", as.integer(n_by_gene["NSUN2"]), 78)
instr$gene <- sapply(instr$gene_id, function(g) names(GENES)[GENES==g])

# ---- 2. rsID -> GRCh37 coordinate bridging (dbSNP mapping merge) ----
logline("---- 2/7 rsID -> GRCh37 bridging (dbSNP b151) ----")
dmap <- read.delim("dbsnp_map.tsv")
logline("  dbSNP mapping table:", nrow(dmap), "rows (263 of 549 rsIDs mapped)")
check("dbSNP mapped rows", nrow(dmap), 263)

m <- merge(instr, dmap, by.x = "rsid", by.y = "rsid", suffixes = c("", ".d"))
# check: GTEx alleles consistent with dbSNP + chromosome consistent
m$ok <- (m$ref == m$dbSNP_ref) & (m$alt == m$dbSNP_alts) & (m$chr37 == GENE_CHR[m$gene])
m <- m[m$ok, ]
logline("  usable instruments after bridging:", nrow(m))
n_m <- table(factor(m$gene, levels=names(GENES)))

# ---- 3. match the 1000G EUR panel by coordinate + alleles ----
logline("---- 3/7 1000G (EUR) matching by GRCh37 coordinate + alleles ----")
match_bim <- function(gene) {
  bim <- read.table(file.path(DIR, paste0("plink_", gene, ".bim")), sep="\t",
                    colClasses=c("character","character","NULL","integer","character","character"))
  names(bim) <- c("chr","id","bp","a1","a2")
  # index both bim allele columns
  keys <- rbind(data.frame(bp=bim$bp, a=bim$a1, id=bim$id),
                data.frame(bp=bim$bp, a=bim$a2, id=bim$id))
  kidx <- paste(keys$bp, keys$a, sep="_")
  gi <- m[m$gene == gene, ]
  gi$matched_bim <- NA
  k1 <- paste(gi$pos37, gi$alt, sep="_")
  k2 <- paste(gi$pos37, gi$ref, sep="_")
  for (i in seq_len(nrow(gi))) {
    hit <- which(kidx == k1[i] | kidx == k2[i])
    if (length(hit) > 0) gi$matched_bim[i] <- keys$id[hit[1]]
  }
  gi
}
matched <- do.call(rbind, lapply(names(GENES), match_bim))
for (g in names(GENES)) {
  n_g <- sum(matched$gene == g & !is.na(matched$matched_bim))
  logline(sprintf("  %s: 1000G coordinate matches %d/%d", g, n_g, sum(matched$gene==g)))
}
check("WARS2 match count", sum(matched$gene=="WARS2" & !is.na(matched$matched_bim)), 142)
check("LRPPRC match count", sum(matched$gene=="LRPPRC" & !is.na(matched$matched_bim)), 94)
check("NSUN2 match count", sum(matched$gene=="NSUN2" & !is.na(matched$matched_bim)), 25)
matched <- matched[!is.na(matched$matched_bim), ]

# ---- 4. PLINK clump ----
logline("---- 4/7 PLINK clump（r2<0.1, 2500kb, p1=5e-8）----")
clumped <- list()
for (g in names(GENES)) {
  gi <- matched[matched$gene == g, ]
  fin <- file.path(DIR, paste0("clumpin_", g, ".tsv"))
  write.table(data.frame(SNP=gi$matched_bim, P=sprintf("%.3e", gi$pvalue)),
              fin, sep="\t", row.names=FALSE, quote=FALSE)
  outf <- file.path(DIR, paste0("clump_", g))
  logline("  ", g, ": running plink clump...")
  invisible(system2(PLINK, args = c("--bfile", shQuote(file.path(DIR, paste0("plink_", g))),
    "--clump", shQuote(fin), "--clump-p1", "5e-8", "--clump-r2", "0.1",
    "--clump-kb", "2500", "--clump-snp-field", "SNP", "--clump-field", "P",
    "--out", shQuote(outf), "--silent"), wait = TRUE))
  clfile <- paste0(outf, ".clumped")
  if (!file.exists(clfile)) {
    logline("  [warning]", g, ": clump output missing (no independent IV for this gene)")
    clumped[[g]] <- gi[0, ]
    next
  }
  cl <- readLines(clfile)
  cl <- cl[!startsWith(trimws(cl), "#") & nzchar(trimws(cl)) & !grepl("^CHR\\s+F\\s+SNP", trimws(cl))]
  # .clumped columns: CHR F SNP BP P ... (SNP is column 3)
  snps <- sapply(cl, function(l) strsplit(trimws(l), "\\s+")[[1]][3])
  logline("  ", g, ": independent IVs after clump:", length(snps), ":", paste(snps, collapse=", "))
  clumped[[g]] <- gi[gi$matched_bim %in% snps, ]
}
check("WARS2 independent IVs", nrow(clumped[["WARS2"]]), 4)
check("LRPPRC independent IVs", nrow(clumped[["LRPPRC"]]), 2)
check("NSUN2 independent IVs", nrow(clumped[["NSUN2"]]), 1)
ivs <- do.call(rbind, clumped)

# ---- 5. match HERMES by coordinate + allele set ----
logline("---- 5/7 HERMES GWAS matching (GRCh37 coordinate + allele set) ----")
her <- read.delim("hermes_window.tsv", colClasses = c(
  SNP="character", CHR="character", BP="integer", A1="character", A2="character",
  freq="numeric", b="numeric", se="numeric", p="numeric", N="NULL"))
her$A1 <- toupper(her$A1); her$A2 <- toupper(her$A2)
her$key <- paste(her$BP, apply(cbind(her$A1, her$A2), 1, function(x) paste(sort(x), collapse="_")), sep="|")
ivs$refU <- toupper(ivs$ref); ivs$altU <- toupper(ivs$alt)
ivs$key <- paste(ivs$pos37, apply(cbind(ivs$refU, ivs$altU), 1, function(x) paste(sort(x), collapse="_")), sep="|")
hidx <- her[!duplicated(her$key), ]
ivs$hermes_row <- hidx[match(ivs$key, hidx$key), ]
n_hit <- sum(!is.na(ivs$hermes_row$SNP))
logline("  HERMES matches:", n_hit, "/", nrow(ivs), "IV loci hit")
check("HERMES hit count", n_hit, 7)

# ---- 6. harmonize (alignment + palindrome check) ----
logline("---- 6/7 allele alignment + palindromic frequency check ----")
comp <- function(a) chartr("ACGT", "TGCA", a)
ivs$beta_y <- NA; ivs$se_y <- NA; ivs$p_y <- NA; ivs$flip <- NA
ivs$palindromic <- FALSE; ivs$F_stat <- NA
for (i in seq_len(nrow(ivs))) {
  h <- ivs$hermes_row[i, ]
  if (is.na(h$SNP)) next
  e_alt <- ivs$altU[i]; e_ref <- ivs$refU[i]
  a1 <- h$A1; a2 <- h$A2
  beta <- NA
  if ((a1==e_alt && a2==e_ref) || (a1==comp(e_alt) && a2==comp(e_ref))) beta <- h$b
  else if ((a1==e_ref && a2==e_alt) || (a1==comp(e_ref) && a2==comp(e_alt))) beta <- -h$b
  if (is.na(beta)) next
  pali <- (e_alt=="A" && e_ref=="T") || (e_alt=="T" && e_ref=="A") ||
          (e_alt=="C" && e_ref=="G") || (e_alt=="G" && e_ref=="C")
  alt_freq <- ivs$ac[i] / ivs$an[i]
  if (pali) {
    d_direct <- abs(h$freq - alt_freq)
    d_flip <- abs(h$freq - (1 - alt_freq))
    if (min(d_direct, d_flip) > 0.2) {
      logline("    [drop] palindrome", ivs$rsid[i], ": frequency difference >0.2, cannot orient")
      next
    }
    if (beta == h$b && d_flip < d_direct) beta <- -h$b
    else if (beta == -h$b && d_direct < d_flip) beta <- h$b
  }
  ivs$beta_y[i] <- beta
  ivs$se_y[i] <- h$se
  ivs$p_y[i] <- h$p
  ivs$flip[i] <- (beta != h$b)
  ivs$palindromic[i] <- pali
  ivs$F_stat[i] <- (ivs$beta[i] / ivs$se[i])^2
}
ivs$used <- !is.na(ivs$beta_y) & !ivs$palindromic
excl <- ivs[!ivs$used, ]
logline("  palindromes dropped:", nrow(excl), "(", paste(excl$rsid, collapse=", "), ")")
check("palindrome drop count", nrow(excl), 2)

# ---- 7. Wald/IVW + BH ----
logline("---- 7/7 Wald ratio / fixed-effect IVW estimation ----")
res <- list()
for (g in names(GENES)) {
  keep <- ivs[ivs$gene == g & ivs$used, ]
  if (nrow(keep) == 0) {
    res[[g]] <- list(status="not estimable", nIV=0); next
  } else if (nrow(keep) == 1) {
    u <- keep[1, ]
    beta <- u$beta_y / u$beta
    se <- abs(u$se_y / u$beta)
    res[[g]] <- list(status="Wald", nIV=1, beta=beta, se=se)
  } else {
    num <- sum(keep$beta_y * keep$beta / keep$se_y^2)
    den <- sum(keep$beta^2 / keep$se_y^2)
    res[[g]] <- list(status="IVW (FE)", nIV=nrow(keep), beta=num/den, se=sqrt(1/den))
  }
}
blood <- read.csv("blood_reference.csv")
for (g in names(GENES)) {
  r <- res[[g]]; b <- blood[blood$Gene == g, ]
  if (!is.null(r$beta)) {
    z <- r$beta / r$se
    p <- 2 * pnorm(-abs(z))
    r$p <- p; res[[g]] <- r
    logline(sprintf("  %s: %s beta=%.6f P=%.4f (%d IV) | blood-eQTL beta=%s P=%s",
                    g, r$status, r$beta, p, r$nIV, b$beta_blood, b$P_blood))
  } else {
    logline("  ", g, ": not estimable | blood-eQTL beta=", b$beta_blood, " P=", b$P_blood)
  }
}
# BH across the 3 testable genes
testable <- names(res)[sapply(res, function(r) !is.null(r$p))]
if (length(testable) > 0) {
  ps <- sapply(testable, function(g) res[[g]]$p)
  q <- p.adjust(ps, method="BH")
  for (g in testable) { res[[g]]$q_BH <- unname(q[g]) }
  logline("  BH adjustment:", paste(sprintf("%s q=%.3f", testable, unname(q)), collapse="  "))
}
check("WARS2 beta", round(res[["WARS2"]]$beta, 6), -0.003045)
check("LRPPRC beta", round(res[["LRPPRC"]]$beta, 6), -0.020276)
check("NSUN2 beta", round(res[["NSUN2"]]$beta, 6), 0.008021)
check("WARS2 P", round(res[["WARS2"]]$p, 4), 0.8609)
check("LRPPRC P", round(res[["LRPPRC"]]$p, 4), 0.7683)
check("NSUN2 P", round(res[["NSUN2"]]$p, 4), 0.8118)

# ---- output ----
out <- data.frame()
for (g in names(GENES)) {
  keep <- ivs[ivs$gene == g, ]
  for (i in seq_len(nrow(keep))) {
    u <- keep[i, ]
    out <- rbind(out, data.frame(
      gene=g, method=res[[g]]$status, nIV=res[[g]]$nIV, rsID=u$rsid, pos37=u$pos37,
      beta_eQTL_alt=round(u$beta,4), se_eQTL=round(u$se,4),
      F_stat=ifelse(is.na(u$beta),NA,round((u$beta/u$se)^2,1)),
      eqtl_p=sprintf("%.2e", u$pvalue),
      beta_HF=ifelse(is.na(u$beta_y),"no HERMES match",round(u$beta_y,4)),
      se_HF=ifelse(is.na(u$se_y),"",round(u$se_y,4)),
      hermes_p=ifelse(is.na(u$p_y),"",round(u$p_y,3)),
      flip=u$flip, palindromic=u$palindromic,
      used_in_estimate=u$used))
  }
}
write.table(out, "R3_local_details.tsv", sep="\t", row.names=FALSE, quote=FALSE)
summ <- data.frame(
  Gene=testable,
  heart_method=sapply(testable, function(g) res[[g]]$status),
  heart_nIV=sapply(testable, function(g) res[[g]]$nIV),
  heart_beta=round(sapply(testable, function(g) res[[g]]$beta),4),
  heart_P=round(sapply(testable, function(g) res[[g]]$p),3),
  heart_qBH=round(sapply(testable, function(g) res[[g]]$q_BH),3),
  blood_beta=blood$beta_blood[match(testable, blood$Gene)],
  blood_P=blood$P_blood[match(testable, blood$Gene)])
write.csv(summ, "R3_local_summary.csv", row.names=FALSE)
logline("output: R3_local_details.tsv + R3_local_summary.csv")
logline("========== R3 local reproduction finished: ", format(Sys.time()), "==========")
logline("          Verify the anchors, then use the local numbers as authoritative.")

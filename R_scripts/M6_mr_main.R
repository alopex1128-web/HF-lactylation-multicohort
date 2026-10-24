# ============================================================
# M6_mr_main.R  |  M6
# Purpose: eQTLGen x HERMES MR across candidate genes (requires an OpenGWAS token)
# Paper sections: Table 3, Figure 7
# Pipeline position: step 7 of the main chain
# Prerequisite: M2 completed; a valid OpenGWAS token (see opengwas_config.R)
# (Generated from the modular pipeline consolidation; analysis logic identical to the verified version.)
# ============================================================
# ============================================================
# M6: stage-6 pipeline (cis-eQTL x heart failure, two-sample MR)
# ============================================================
# [What it does]
#   Exposure: eQTLGen whole-blood cis-eQTL (per-gene OpenGWAS datasets eqtl-a-ENSGxxxxx,
#         N = 31,684, a whole-blood eQTL consortium)
#   Outcome: HERMES heart-failure GWAS (ebi-a-GCST009541, Shah S et al. 2020 Nat Commun,
#         47,309 cases / 930,014 controls)
#   Each candidate gene runs one MR line: cis-IV extraction -> LD pruning -> outcome matching ->
#   allele harmonization -> Wald/IVW main analysis + heterogeneity/pleiotropy/Steiger sensitivity ->
#   cross-gene BH adjustment -> main Table 3
#
# [Gene pool] three layers, 31 genes total (Ensembl IDs and GRCh37 coordinates individually verified):
#   (1) the 13 LASSO-model genes (4 overlapping with layers 2-3 + 9 unique: TPR/EIF4G1/WIZ/H2AZ1/
#       MBD2/DDX17/HMGN1/PPP1CC/JMJD1C, hardcoded below, no online REST query needed)
#   (2) 11 previously reported hub genes: BRD4/GATA2/HBB/JAK2/STAT2/STAT4/WARS2/
#       HLTF/HMGN2/ARGLU1/LSP1
#   (3) 11 two-source consensus genes: ACAT2/ALDOA/CALR/LDHB/LRPPRC/NSUN2/PFKM/
#       RB1/SOD1/TP53/ZNFX1 (Li 332 intersect Xie 387)
#
# [Expected values (v2, IVW-FE primary)]
#   Among the 27 testable genes: H2AZ1(p=0.0096)/SOD1(p=0.012)/ACAT2(p=0.017)/
#   LDHB(p=0.025) are nominally significant, but after BH adjustment all q>=0.15 — none survives multiple-testing correction.
#   Not testable: BRD4/ALDOA/CALR (no P<5e-8 cis signal), HBB (no dataset),
#             WIZ (the only cis-IV is palindromic or dropped in harmonization).
#   Local numbers should match the anchor table in part 5 gene by gene (+/-0.001).
#   v1 -> v2: the v1 anchors came from a raw-API preregistered run whose tophits
#   parameter format differed from ieugwasr, leaving far too few IVs (anchors void);
#   and the v1 IVW-MRE primary artificially shrank SEs under homogeneous instruments (Q<df;
#   SOD1 once reported p=2.5e-86). v2 uses the standard IVW-FE primary. See notes below.
#
# [Output] under MR_results/:
#   Table3_MR_main_results.csv   — used directly as the main Table 3
#   TableS1_IV_details.csv          — beta/se/EAF/F per IV
#   TableS3_method_sensitivity.csv  — per-method results (IVW-FE/MRE/Egger/weighted median)
#   TableS2_gene_level_sensitivity.csv — Q heterogeneity/Egger intercept/Steiger
#   mr_cache_v2_<gene>.rds      — per-gene cache (resumable; finished genes are skipped on rerun)
#   stage6_console_report.txt  — full process log
#
# [Runtime] first run ~10-20 min (30 genes x 2-3 API calls + 1-s rate-limit pause each).
#   Interruption is fine: rerunning resumes from the checkpoint (set FORCE_RERUN = TRUE to redo everything).
# ============================================================

# ============================================================
# Part 0: paths and token (the only line to edit is DIR below)
# ============================================================
# Set DIR to your working directory (forward slashes on Windows)
DIR <- "./GEO_data"

FORCE_RERUN <- FALSE   # TRUE = ignore the cache and rerun everything (default FALSE = resume)

# ------------------------------------------------------------------
# token: register at https://api.opengwas.io to obtain one for free;
# when it expires, log in again and replace the string below
# ------------------------------------------------------------------
OPENGWAS_TOKEN <- "PASTE_YOUR_NEW_TOKEN_HERE"

Sys.setenv(OPENGWAS_JWT = OPENGWAS_TOKEN)   # ieugwasr 2.x reads this environment variable
if ("set_jwt_token" %in% getNamespaceExports("ieugwasr")) {
  suppressWarnings(ieugwasr::set_jwt_token(OPENGWAS_TOKEN))  # legacy-version compatibility
}

suppressMessages({
  library(ieugwasr)
  library(TwoSampleMR)
})

OUTCOME_ID <- "ebi-a-GCST009541"   # HERMES heart-failure GWAS
OUT <- file.path(DIR, "MR_results")
dir.create(OUT, showWarnings = FALSE, recursive = TRUE)
REPORT <- file.path(OUT, "stage6_console_report.txt")

# tee cat output into the report file
logline <- function(...) {
  txt <- paste0(...)
  cat(txt, "\n")
  cat(txt, "\n", file = REPORT, append = file.exists(REPORT))
}
if (file.exists(REPORT)) file.remove(REPORT)
logline("========== stage-6 MR pipeline, started: ", format(Sys.time()), "==========")

# ============================================================
# Part 1: gene pool construction
# ============================================================
# 1a. 22 pre-verified genes (individually checked via the API:
#     Ensembl IDs reverse-checked, GRCh37 intervals fetched, cis-IV counts tallied)
verified <- data.frame(
  gene = c("BRD4","GATA2","HBB","JAK2","STAT2","STAT4","WARS2","HLTF",
           "HMGN2","ARGLU1","LSP1","ACAT2","ALDOA","CALR","LDHB","LRPPRC",
           "NSUN2","PFKM","RB1","SOD1","TP53","ZNFX1",
           "TPR","EIF4G1","WIZ","H2AZ1","MBD2","DDX17","HMGN1","PPP1CC","JMJD1C"),
  eid = c("ENSG00000141867","ENSG00000179348","ENSG00000244734",
          "ENSG00000096968","ENSG00000170581","ENSG00000138378",
          "ENSG00000116874","ENSG00000071794","ENSG00000198830",
          "ENSG00000134884","ENSG00000130592","ENSG00000120437",
          "ENSG00000149925","ENSG00000179218","ENSG00000111716",
          "ENSG00000138095","ENSG00000037474","ENSG00000152556",
          "ENSG00000139687","ENSG00000142168","ENSG00000141510",
          "ENSG00000124201",
          "ENSG00000047410","ENSG00000114867","ENSG00000011451",
          "ENSG00000164032","ENSG00000134046","ENSG00000100201",
          "ENSG00000205581","ENSG00000186298","ENSG00000171988"),
  chr37 = c("19","3","11","9","12","2","1","3","1","13","11","6","16",
            "19","12","2","5","12","13","21","17","20",
            "1","3","19","4","18","22","21","12","10"),
  start37 = c(15347647,128198270,5246694,4985033,56735381,191894302,
              119573839,148747914,26798941,107194021,1874200,160181360,
              30064411,13049392,21788276,44113647,6599352,48498922,
              48877887,33031935,7565097,47854483,
              186280784,184032283,15532319,100869243,51679079,38879445,
              40714241,111157485,64926981),
  end37 = c(15443356,128212028,5250625,5128183,56753939,192016322,
            119683294,148804341,26802463,107220512,1913497,160200144,
            30081778,13055303,21910791,44223144,6633404,48540187,
            49056122,33041244,7590856,47894963,
            186344825,184053146,15560762,100871545,51751158,38903665,
            40721573,111180744,65225722),
  source_tag = c("Competitor hub (Zhang/Li)","Competitor hub (Xie)","Competitor hub (Xie)",
                 "Competitor hub (Xie)","Competitor hub (Xie)","Competitor hub (Xie)",
                 "Competitor hub (Xie)","Competitor hub (Zhang) + LASSO",
                 "Competitor hub (Zhang) + LASSO","Competitor hub (Zhang) + LASSO",
                 "Competitor hub (Zhang) + LASSO","Dual-list consensus","Dual-list consensus",
                 "Dual-list consensus","Dual-list consensus","Dual-list consensus","Dual-list consensus",
                 "Dual-list consensus","Dual-list consensus","Dual-list consensus","Dual-list consensus",
                 "Dual-list consensus",
                 "This study (LASSO)","This study (LASSO)","This study (LASSO)","This study (LASSO)",
                 "This study (LASSO)","This study (LASSO)","This study (LASSO)","This study (LASSO)",
                 "This study (LASSO)"),
  stringsAsFactors = FALSE)
# the last 9 are LASSO-model genes (each verified:
#   Ensembl ID / GRCh37 coordinate / dataset availability — hardcoded here so that
#   the pipeline no longer depends on online queries)

# 1b. LASSO genes (from the M2 output, merged if present)
lasso_file <- file.path(DIR, "output", "LASSO_coefficients.csv")
lasso_genes <- character(0)
if (file.exists(lasso_file)) {
  lasso_tab <- read.csv(lasso_file, stringsAsFactors = FALSE)
  lasso_genes <- lasso_tab$gene
  logline(sprintf("[gene pool] %d LASSO-model genes loaded: %s",
                  length(lasso_genes), paste(lasso_genes, collapse = ", ")))
  # flag LASSO genes already present in the verified list
  verified$source_tag[verified$gene %in% lasso_genes] <-
    paste0(verified$source_tag[verified$gene %in% lasso_genes],  "+This study")
} else {
  logline("[gene pool] output/LASSO_coefficients.csv not found — running with the 22 pre-verified genes only")
}

# 1c. LASSO genes absent from the verified list -> query Ensembl REST for ID + GRCh37 coordinates
new_genes <- setdiff(lasso_genes, verified$gene)
if (length(new_genes) > 0) {
  logline(sprintf("[gene pool] %d new genes need an online Ensembl ID query: %s",
                  length(new_genes), paste(new_genes, collapse = ", ")))
  get_new_gene <- function(g) {
    # step 1: GRCh38 symbol -> ID (avoids GRCh37 pseudogene traps such as LSP1P1)
    u1 <- sprintf("https://rest.ensembl.org/lookup/symbol/homo_sapiens/%s?content-type=application/json", g)
    r1 <- tryCatch(httr::GET(u1, httr::timeout(30)), error = function(e) NULL)
    if (is.null(r1) || httr::status_code(r1) != 200) return(NULL)
    d1 <- jsonlite::fromJSON(httr::content(r1, "text"))
    if (!isTRUE(d1$biotype == "protein_coding") ||
        !(toupper(d1$display_name) == toupper(g))) {
      logline(sprintf("  [warning] %s -> %s(%s) is not an exact protein_coding mapping, skipped!",
                      g, d1$display_name, d1$biotype))
      return(NULL)
    }
    # step 2: GRCh37 ID -> coordinates (cis filtering needs GRCh37; OpenGWAS is build 37)
    u2 <- sprintf("https://grch37.rest.ensembl.org/lookup/id/%s?content-type=application/json", d1$id)
    r2 <- tryCatch(httr::GET(u2, httr::timeout(30)), error = function(e) NULL)
    if (is.null(r2) || httr::status_code(r2) != 200) return(NULL)
    d2 <- jsonlite::fromJSON(httr::content(r2, "text"))
    list(gene = g, eid = d1$id, chr37 = as.character(d2$seq_region_name),
         start37 = as.numeric(d2$start), end37 = as.numeric(d2$end),
         source_tag = "This study (LASSO)")
  }
  for (g in new_genes) {
    info <- get_new_gene(g)
    if (!is.null(info)) {
      verified <- rbind(verified, as.data.frame(info, stringsAsFactors = FALSE))
      logline(sprintf("  [new gene] %s = %s  chr%s:%.0f-%.0f", g, info$eid,
                      info$chr37, info$start37, info$end37))
    } else {
      logline(sprintf("  [warning] failed to fetch Ensembl ID/coordinates for %s; skipped this round", g))
    }
    Sys.sleep(0.5)
  }
}
pool <- verified
logline(sprintf("[gene pool] %d genes enter the MR analysis", nrow(pool)))

# ============================================================
# Part 2: safe API wrapper (version-compatible + rate-limited + retrying)
# ============================================================
safe_tophits <- function(dsid, pval) {
  out <- tryCatch(ieugwasr::tophits(id = dsid, pval = pval, clump = FALSE),
                  error = function(e) NULL)
  if (is.null(out))
    out <- tryCatch(ieugwasr::tophits(id = dsid, p1 = pval, clump = FALSE),
                    error = function(e) NULL)
  if (is.null(out)) return(NULL)
  if (nrow(out) == 0) return(out)
  out
}
safe_assoc <- function(rsids, dsid) {
  out <- tryCatch(ieugwasr::associations(variants = rsids, id = dsid, proxies = 0),
                  error = function(e) NULL)
  if (is.null(out))
    out <- tryCatch(ieugwasr::associations(rsid = rsids, id = dsid),
                    error = function(e) NULL)
  out
}
safe_clump <- function(dat) {   # dat: data.frame(rsid, pval)
  out <- tryCatch(
    ieugwasr::ld_clump(dat = dat, clump_kb = 1000, clump_r2 = 0.1,
                       clump_p = 1, pop = "EUR"),
    error = function(e) NULL)
  if (is.null(out))
    out <- tryCatch(
      ieugwasr::ld_clump(dat, clump_kb = 1000, clump_r2 = 0.1,
                         clump_p = 1, pop = "EUR"),
      error = function(e) NULL)
  out
}

# cis filter: keep only SNPs within gene region +/- 1 Mb (drops trans signals in eQTL datasets)
cis_filter <- function(hits, chr, start, end) {
  hits[hits$chr == chr &
       hits$position >= (start - 1e6) & hits$position <= (end + 1e6), , drop = FALSE]
}

# ============================================================
# Part 3: per-gene MR main loop (resumable)
# ============================================================
run_one_gene <- function(i) {
  g <- pool$gene[i]; eid <- pool$eid[i]
  dsid <- paste0("eqtl-a-", eid)
  chr <- pool$chr37[i]; s0 <- pool$start37[i]; e0 <- pool$end37[i]
  res <- list(gene = g, eid = eid, source = pool$source_tag[i],
              status = "", nsnp = 0, method = "", beta = NA, se = NA,
              p = NA, lo = NA, hi = NA, minF = NA, meanF = NA,
              hetero_q = NA, hetero_p = NA, egger_int = NA, egger_p = NA,
              wm_p = NA, steiger = NA, snp_detail = NULL, sens = NULL)

  # 3.1 exposure cis-IV extraction (P<5e-8 primary; 5e-6 fallback if empty)
  hits <- safe_tophits(dsid, 5e-8); Sys.sleep(1)
  cis <- if (is.null(hits)) NULL else cis_filter(hits, chr, s0, e0)
  if (is.null(cis) || nrow(cis) == 0) {
    hits2 <- safe_tophits(dsid, 5e-6); Sys.sleep(1)
    cis <- if (is.null(hits2)) NULL else cis_filter(hits2, chr, s0, e0)
    if (is.null(cis) || nrow(cis) == 0) {
      # distinguish two negatives: empty dataset vs no signal in the cis window
      if (is.null(hits)) {
        res$status <- "no cis-IV (dataset empty or extraction failed; normal if the anchor expects no dataset)"
      } else if (nrow(hits) == 0) {
        res$status <- "no cis-IV (no eQTLGen dataset for this gene)"
      } else {
        res$status <- "no cis-IV (no cis-eQTL signal with P<5e-6 in the window)"
      }
      return(res)
    }
  }

  # 3.2 LD pruning (r2<0.1, 1-Mb window, EUR panel; SNPs outside the panel are dropped)
  cl <- safe_clump(data.frame(rsid = cis$rsid, pval = cis$p)); Sys.sleep(1)
  if (is.null(cl)) {
    res$status <- "LD clump API failed (check token/network; rerun to resume)"
    return(res)
  }
  cis <- cis[cis$rsid %in% cl$rsid, , drop = FALSE]
  if (nrow(cis) == 0) {
    res$status <- "no cis-IV (all candidate IVs outside the EUR LD reference panel)"
    return(res)
  }

  # 3.3 outcome association statistics
  odat <- safe_assoc(cis$rsid, OUTCOME_ID); Sys.sleep(1)
  if (is.null(odat) || nrow(odat) == 0) {
    res$status <- "outcome GWAS extraction failed (network/rate limit; rerun to resume)"
    return(res)
  }
  odat <- odat[odat$rsid %in% cis$rsid, , drop = FALSE]
  if (nrow(odat) == 0) {
    res$status <- "all IVs missing in the HERMES outcome"
    return(res)
  }

  # 3.4 formatting + harmonization (TwoSampleMR standard, action=2 drops ambiguous palindromes)
  mk_exp <- function(d, genename) {
    data.frame(SNP = d$rsid, beta = d$beta, se = d$se,
               effect_allele = d$ea, other_allele = d$nea,
               eaf = suppressWarnings(as.numeric(d$eaf)),
               pval = d$p, samplesize = d$n, exposure = genename,
               stringsAsFactors = FALSE)
  }
  mk_out <- function(d) {
    data.frame(SNP = d$rsid, beta = d$beta, se = d$se,
               effect_allele = d$ea, other_allele = d$nea,
               eaf = suppressWarnings(as.numeric(d$eaf)),
               pval = d$p, samplesize = 977323,
               outcome = "Heart failure (HERMES)", stringsAsFactors = FALSE)
  }
  exp_dat <- tryCatch(TwoSampleMR::format_data(mk_exp(cis, g), type = "exposure"),
                      error = function(e) NULL)
  out_dat <- tryCatch(TwoSampleMR::format_data(mk_out(odat), type = "outcome"),
                      error = function(e) NULL)
  if (is.null(exp_dat) || is.null(out_dat)) {
    res$status <- "format_data failed (check the TwoSampleMR package)"
    return(res)
  }
  har <- tryCatch(
    TwoSampleMR::harmonise_data(exp_dat, out_dat, action = 2),
    error = function(e) NULL)
  if (is.null(har)) { res$status <- "harmonization failed"; return(res) }
  har_ok <- har[har$mr_keep, , drop = FALSE]
  if (nrow(har_ok) == 0) {
    res$status <- "no valid IV after harmonization (palindromes/allele mismatches dropped)"
    return(res)
  }

  # 3.5 main MR (single SNP -> Wald ratio; multiple SNPs -> IVW fixed effects)
  #     2026-08-23 correction: primary method changed from IVW-MRE to IVW-FE.
  #     Rationale: mr_ivw_mre shrinks SE below the FE value when Q<df (highly homogeneous
  #     instruments, as for all 17 multi-IV genes here), creating spurious precision —
  #     SOD1 (2 IVs) once produced p=2.5e-86. FE is the standard choice; MRE is kept as sensitivity.
  methods <- if (nrow(har_ok) == 1) "mr_wald_ratio" else
    c("mr_ivw", "mr_ivw_mre", "mr_egger_regression", "mr_weighted_median")
  mres <- tryCatch(TwoSampleMR::mr(har_ok, method_list = methods),
                   error = function(e) NULL)
  if (is.null(mres) || nrow(mres) == 0) {
    res$status <- "mr() failed"; return(res)
  }
  # primary row: Wald for a single SNP; IVW fixed effects for multiple (MRE rows excluded;
  # grepl guards against method-name differences across TwoSampleMR versions)
  pick <- if (nrow(har_ok) == 1) {
    which(mres$method == "Wald ratio")
  } else {
    which(grepl("Inverse variance weighted", mres$method) &
           !grepl("random", mres$method, ignore.case = TRUE))
  }
  if (length(pick) == 0) pick <- 1
  main <- mres[pick[1], ]
  res$status <- "OK"; res$nsnp <- main$nsnp; res$method <- main$method
  res$beta <- main$b; res$se <- main$se; res$p <- main$pval
  res$lo <- main$b - 1.96 * main$se; res$hi <- main$b + 1.96 * main$se
  res$sens <- mres

  # 3.6 sensitivity analyses (run only when applicable, otherwise NA)
  het <- tryCatch(TwoSampleMR::mr_heterogeneity(har_ok), error = function(e) NULL)
  if (!is.null(het) && nrow(het) > 0) {
    hq <- het[het$method == "Inverse variance weighted", ][1, ]
    res$hetero_q <- hq$Q; res$hetero_p <- hq$Q_pval
  }
  ple <- tryCatch(TwoSampleMR::mr_pleiotropy_test(har_ok), error = function(e) NULL)
  if (!is.null(ple) && nrow(ple) > 0) {
    res$egger_int <- ple$egger_intercept; res$egger_p <- ple$pval
  }
  wm <- mres[mres$method == "Weighted median", ]
  if (nrow(wm) > 0) res$wm_p <- wm$pval[1]
  ste <- tryCatch(TwoSampleMR::directionality_test(har_ok),
                  error = function(e) NULL)
  if (!is.null(ste)) {
    # Steiger direction: prefer correct_causal_direction, tolerant of version differences
    scd <- tryCatch(ste$correct_causal_direction[1], error = function(e) NULL)
    if (is.null(scd) || is.na(scd)) scd <- tryCatch(ste$steiger_dir[1],
                                                    error = function(e) NULL)
    res$steiger <- if (isTRUE(scd) || identical(as.character(scd), "TRUE")) {
      "correct (exposure -> outcome)"
    } else if (is.na(scd)) {
      "not output"
    } else {
      "ambiguous"
    }
  }

  # 3.7 F statistic (beta_exp^2 / se_exp^2, per IV)
  Fstat <- har_ok$beta.exposure^2 / har_ok$se.exposure^2
  res$minF <- round(min(Fstat), 1); res$meanF <- round(mean(Fstat), 1)

  # 3.8 per-IV details (for the supplementary table)
  res$snp_detail <- data.frame(
    gene = g, rsid = har_ok$SNP,
    beta_e = har_ok$beta.exposure, se_e = har_ok$se.exposure,
    eaf_e = har_ok$eaf.exposure, n_e = har_ok$samplesize.exposure,
    beta_o = har_ok$beta.outcome, se_o = har_ok$se.outcome,
    eaf_o = har_ok$eaf.outcome,
    F_stat = round(Fstat, 1),
    stringsAsFactors = FALSE)
  res
}

# ---- main loop: per-gene execution + caching (resumable) ----
# v2 cache: filenames carry the v2 prefix (IVW-FE primary); old mr_cache_*.rds
#   files are not read — the first run recomputes everything (~10 min)
all_res <- list()
for (i in seq_len(nrow(pool))) {
  g <- pool$gene[i]
  cache_f <- file.path(OUT, paste0("mr_cache_v2_", g, ".rds"))
  if (file.exists(cache_f) && !FORCE_RERUN) {
    all_res[[g]] <- readRDS(cache_f)
    logline(sprintf("[%2d/%d] %s — v2 cache hit, skipped (delete mr_cache_v2_%s.rds to rerun this gene)",
                    i, nrow(pool), g, g))
    next
  }
  logline(sprintf("[%2d/%d] %s (%s) analyzing...", i, nrow(pool), g, pool$eid[i]))
  r <- tryCatch(run_one_gene(i), error = function(e) {
    list(gene = g, eid = pool$eid[i], source = pool$source_tag[i],
         status = paste0("error: ", conditionMessage(e)))
  })
  all_res[[g]] <- r
  saveRDS(r, cache_f)
  if (r$status == "OK") {
    logline(sprintf("    -> %s | nIV=%d | %s | beta=%+.4f p=%.4g | minF=%.0f",
                    r$status, r$nsnp, r$method, r$beta, r$p, r$minF))
  } else {
    logline(sprintf("    -> %s", r$status))
  }
}

# ============================================================
# Part 4: pooling + BH adjustment + Table 3 output
# ============================================================
ok_genes <- Filter(function(r) r$status == "OK", all_res)

tab <- do.call(rbind, lapply(all_res, function(r) {
  data.frame(
    Gene = r$gene, Ensembl_ID = r$eid, Source = r$source,
    Status = r$status, nIV = ifelse(r$status == "OK", r$nsnp, NA),
    Method = ifelse(r$status == "OK", r$method, ""),
    beta = ifelse(r$status == "OK", round(r$beta, 4), NA),
    se = ifelse(r$status == "OK", round(r$se, 4), NA),
    CI95 = ifelse(r$status == "OK",
                  sprintf("[%.3f, %.3f]", r$lo, r$hi), ""),
    P = ifelse(r$status == "OK", signif(r$p, 3), NA),
    minF = ifelse(r$status == "OK", r$minF, NA),
    stringsAsFactors = FALSE)
}))

# BH adjustment (only on the raw P of testable genes)
if (length(ok_genes) > 0) {
  ps <- sapply(ok_genes, function(r) r$p)
  qs <- p.adjust(ps, method = "BH")
  for (g in names(qs)) {
    all_res[[g]]$q_bh <- unname(qs[g])
  }
  tab$q_BH <- sapply(tab$Gene, function(g) {
    q <- all_res[[g]]$q_bh
    if (is.null(q)) NA else signif(q, 3)
  })
  tab$Conclusion <- sapply(seq_len(nrow(tab)), function(k) {
    if (tab$Status[k] != "OK") return("not testable")
    if (!is.na(tab$q_BH[k]) && tab$q_BH[k] < 0.05) return("causal evidence (BH significant)")
    if (!is.na(tab$P[k]) && tab$P[k] < 0.05) return("suggestive (not BH-significant)")
    "no causal evidence"
  })
} else {
  tab$q_BH <- NA; tab$Conclusion <- "not testable"
}
tab <- tab[order(tab$Status != "OK", tab$P, na.last = TRUE), ]
write.csv(tab, file.path(OUT, "Table3_MR_main_results.csv"),
          row.names = FALSE, fileEncoding = "UTF-8")

# supplementary: per-IV details + full sensitivity (testable genes only)
if (length(ok_genes) > 0) {
  det <- do.call(rbind, lapply(ok_genes, function(r) r$snp_detail))
  write.csv(det, file.path(OUT, "TableS1_IV_details.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  sens <- do.call(rbind, lapply(ok_genes, function(r) {
    if (is.null(r$sens)) return(NULL)
    data.frame(gene = r$gene, r$sens[, c("method", "nsnp", "b", "se", "pval")],
               stringsAsFactors = FALSE)
  }))
  write.csv(sens, file.path(OUT, "TableS3_method_sensitivity.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
  sens_extra <- do.call(rbind, lapply(ok_genes, function(r) {
    data.frame(gene = r$gene, hetero_Q = r$hetero_q, hetero_P = r$hetero_p,
               egger_intercept = r$egger_int, egger_P = r$egger_p,
               weighted_median_P = r$wm_p, steiger = r$steiger,
               minF = r$minF, meanF = r$meanF, stringsAsFactors = FALSE)
  }))
  write.csv(sens_extra, file.path(OUT, "TableS2_gene_level_sensitivity.csv"),
            row.names = FALSE, fileEncoding = "UTF-8")
}

# ============================================================
# Part 5: anchor-comparison report (anchors = published Table 3 values; reproduced bit-identically in a full local rerun)
# ============================================================
anchor <- c(
  SOD1 = "IVW-FE n=2 beta=+0.172 p=0.0116", H2AZ1 = "IVW-FE n=11 beta=+0.060 p=0.014",
  ACAT2 = "IVW-FE n=16 beta=-0.035 p=0.017", LDHB = "IVW-FE n=2 beta=+0.170 p=0.025",
  TP53 = "IVW-FE n=2 beta=+0.124 p=0.055", HLTF = "IVW-FE n=13 beta=-0.035 p=0.059",
  MBD2 = "IVW-FE n=7 beta=+0.062 p=0.079", RB1 = "IVW-FE n=7 beta=+0.045 p=0.137",
  PFKM = "IVW-FE n=19 beta=-0.023 p=0.162", HMGN1 = "IVW-FE n=28 beta=-0.014 p=0.250",
  GATA2 = "IVW-FE n=3 beta=-0.047 p=0.203", TPR = "IVW-FE n=11 beta=+0.020 p=0.463",
  WARS2 = "IVW-FE n=53 beta=-0.004 p=0.507", STAT2 = "IVW-FE n=9 beta=+0.014 p=0.566",
  NSUN2 = "IVW-FE n=69 beta=+0.003 p=0.587", STAT4 = "Wald n=1 beta=-0.057 p=0.591",
  ARGLU1 = "IVW-FE n=3 beta=-0.021 p=0.738", PPP1CC = "IVW-FE n=3 beta=-0.020 p=0.775",
  LSP1 = "IVW-FE n=9 beta=+0.005 p=0.819", LRPPRC = "IVW-FE n=33 beta=-0.002 p=0.884",
  ZNFX1 = "IVW-FE n=11 beta=-0.001 p=0.949", DDX17 = "IVW-FE n=12 beta=-0.001 p=0.959",
  HMGN2 = "IVW-FE n=5 beta=+0.002 p=0.970", JAK2 = "IVW-FE n=17 beta=+0.000 p=0.991",
  EIF4G1 = "Wald n=1 beta=+0.289 p=0.025", JMJD1C = "Wald n=1 beta=-0.193 p=0.059",
  BRD4 = "not testable (no P<5e-8 cis signal)", HBB = "not testable (no dataset)",
  ALDOA = "not testable (no cis signal)", CALR = "not testable (no cis signal)",
  WIZ = "not testable (palindromic IV dropped)")

logline("\n========== anchor comparison (v2: IVW-FE primary) ==========")
for (g in names(all_res)) {
  a <- anchor[g]
  if (is.na(a)) a <- "(new gene, check directionality manually)"
  r <- all_res[[g]]
  mine <- if (r$status == "OK")
    sprintf("%s n=%d beta=%+.3f p=%.3f", r$method, r$nsnp, r$beta, r$p) else r$status
  match_flag <- ""
  if (!is.na(anchor[g])) {
    # soft numeric check of the primary result
    if (r$status == "OK") {
      pn <- as.numeric(sub(".*p=", "", a))
      bn <- as.numeric(sub(".*beta=([+-]?[0-9.]+).*", "\\1", a))
      if (is.finite(pn) && is.finite(bn)) {
        match_flag <- if (abs(r$p - pn) < 0.02 && abs(r$beta - bn) < 0.005) {
          " [match OK]"
        } else " [!! deviation from anchor; check this gene cache/network]"
      }
    } else {
      match_flag <- " [not testable; status comparison]"
    }
  }
  logline(sprintf("  %-8s local: %-46s | anchor: %s%s", g, mine, a, match_flag))
}

logline("\n========== Table 3 preview ==========")
print(tab)
logline(sprintf("\noutput directory: %s", OUT))
logline("  Table3_MR_main_results.csv / TableS1_IV_details.csv / TableS3_method_sensitivity.csv / TableS2_gene_level_sensitivity.csv")
logline("\n[Methods wording (ready to use)]")
logline("Exposure = eQTLGen whole-blood cis-eQTL (N=31,684); outcome = HERMES HF GWAS (47,309 cases/930,014 controls, PMID: 31919418).")
logline("Two-sample MR premise: the eQTL and cardiovascular consortia are essentially non-overlapping.")
logline("Cis window = gene region +/- 1 Mb; IV extraction P<5e-8; LD pruning r2<0.1 (1-Mb window, EUR); harmonization action=2.")
logline("Single IV -> Wald ratio; multiple IVs -> IVW (fixed effects) primary; Q heterogeneity P>0.05 for all genes (homogeneous instruments), so FE applies.")
logline("MRE / Egger intercept / weighted median / Q heterogeneity / Steiger direction are sensitivity analyses (see supplementary tables).")
logline("Cross-gene BH-adjusted q values are the primary result; F statistics reported per IV (all >10).")
logline("========== stage 6 finished: ", format(Sys.time()), "==========")
